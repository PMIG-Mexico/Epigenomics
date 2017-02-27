__author__ = 'jorge'

import csv, re, argparse, os, copy, time, sqlite3
from multiprocessing import JoinableQueue, cpu_count, Process, Queue, active_children

parser = argparse.ArgumentParser()
parser.add_argument('-g', action='store', dest='gene_name',
                    help='The name of the gene to be studied')
parser.add_argument('-G', action='store', dest='gene_info',
                    help='Path of the Gencode (gene annotation) gene_info file')
parser.add_argument('-f', action='store', dest='gtf',
                    help='Path of the Gencode (exons, introns, etc.) gtf file')
parser.add_argument('-p', action='store', dest='prom_ann',
                    help='Path of the FANTOM promoter annotation file')
parser.add_argument('-a', action='store', dest='association',
                    help='Path of the FANTOM enhancer-promoter association and correlation file')
parser.add_argument('-e', action='store', dest='eg_name',
                    help='Path of the EG reference file with the tissues IDs')
parser.add_argument('-x', action='store', dest='expression_dir',
                    help='Path of the directory with the Roadmap Epigenomics expression files')
parser.add_argument('-t', action='store', dest='hist_dir',
                    help='Path of the directory with the Roadmap Epigenomics '
                         'histone tagAlign files in directories named by tissue')
parser.add_argument('-m', action='store', dest='meth_dir',
                    help='Path of the directory with the Roadmap Epigenomics '
                         'methylation bedgraph files in directories named by tissue')
parser.add_argument('-D', default=float('inf'), type=int, action='store', dest='dist',
                    help='Distance value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes the promoters with a SMALLER distance than the one given (default: no filter)')
parser.add_argument('-C', default=0, type=float, action='store', dest='corr',
                    help='Correlation value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a HIGHER correlation than the one given (default: no filter)')
parser.add_argument('-P', default=float('inf'), type=float, action='store', dest='pval',
                    help='P-value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a SMALLER p-value than the one given (default: no filter)')
parser.add_argument('-F', default=float('inf'), type=float, action='store', dest='fdr',
                    help='FDR value to filter promoters from the '
                         'enhancer-promoter association file. '
                         'Takes promoters with a SMALLER FDR than the one given (default: no filter)')
parser.add_argument('-l', default=None, action='store', dest='permitted_file',
                    help='File with the list of permitted genes, otherwise it will search for it')
parser.add_argument('-d', default='RFEpigenomics2.db', action='store', dest='db_name',
                    help='Path to the database to be used')
pargs = parser.parse_args()


def filteredDatabase(gene_info, prom_ann, association, dist, corr, pval, fdr):
    # Readers for the files
    gene_reader = csv.reader(open(gene_info, 'r'), delimiter='\t')
    prom_reader = csv.reader(open(prom_ann, 'r'), delimiter='\t')
    asso_reader = csv.reader(open(association, 'r'), delimiter='\t')

    # Clean the first lines (no information there)
    for next_line in range(1, 9):
        prom_reader.next()
        if next_line == 7:  # Just one line here
            asso_reader.next()

    # Get a dictionary of all promoters associated to a gene
    promoters = {}  # Keys are genes and values are promoter coordinates
    for pr in prom_reader:
        # In case there are more than 1 associated gene
        split = pr[1].split(',')
        coord = pr[0].split(',')[0]

        for s in split:
            found = re.search(r'p(\d)@(.*)', s)

            if bool(found):
                # Add the entries to the promoters dictionary
                if found.group(2) not in promoters.keys():
                    promoters[found.group(2)] = [coord]

                else:
                    promoters[found.group(2)].append(coord)

    # Get a dictionary of the genes in the Gencode file
    gencode_genes = {} # Keys are genes and values are Encode IDs
    for gr in gene_reader:
        if gr[6] != '' and gr[6] != 'NA':
            gencode_genes[gr[6]] = gr[0]

    # Get a list of the promoters and enhancers in the FANTOM association file with the corresponding filters
    promoter_enhancer = {}  # Promoter coordinates as keys and enhancer coordinates as values
    for ar in asso_reader:
        if float(ar[2]) < dist and float(ar[3]) > corr and float(ar[4]) < pval and float(ar[5]) < fdr:
            prom_coord = ar[1].split(',')[0]

            if prom_coord not in promoter_enhancer.keys():
                promoter_enhancer[prom_coord] = [ar[0]]

            else:
                promoter_enhancer[prom_coord].append(ar[0])

    # Make a deep copy of the promoters dictionary to iterate through it
    copy_promoters = copy.deepcopy(promoters)
    # Purge the promoters dictionary
    for key in copy_promoters.keys():
        # Remove the promoters not found in the associated promoters
        for coord_value in copy_promoters[key]:
            if coord_value not in promoter_enhancer.keys():
                promoters[key].remove(coord_value)

        # If the gene has no associated promoters delete its entry
        # Remove the promoter-associated genes not found in the Gencode file
        if len(promoters[key]) == 0 or key not in gencode_genes.keys():
            promoters.pop(key)

    # Create the txt tab delimited file with the list of genes
    filename = 'IMPermitted-D{}C{}P{}F{}.txt'.format(dist, corr, pval, fdr)
    with open(filename, 'a') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for gene in promoters.keys():
            # reduce() joins the list of lists of associated enhancers for each promoters associated to a gene
            enhancers = reduce(lambda x, y: x + y, [promoter_enhancer[prom] for prom in promoters[gene]])
            # Columns in the file will be Encode ID - gene name - promoters - enhancers
            writer.writerow([gencode_genes[gene]] + [gene] + [','.join(promoters[gene])] + [','.join(enhancers)])


class Producer(object):

    def __init__(self, gene, permitted, hist_queue, meth_queue):
        self.gtf_file = pargs.gtf
        self.gene = gene
        self.permitted = permitted
        self.hist_queue = hist_queue
        self.meth_queue = meth_queue
        self.poison_queue = Queue()

        # Start the processes
        pr_en_process = Process(target=self.promoters_enhancers)
        pr_en_process.start()
        ex_in_process = Process(target=self.exons_introns)
        ex_in_process.start()
        poison_process = Process(target=self.poison)
        poison_process.start()

    def promoters_enhancers(self):
        # Search for the gene name in the permitted genes file
        print 'Opening {} file'.format(self.permitted)
        with open(self.permitted, 'r') as perm_file:
            permitted_reader = csv.reader(perm_file, delimiter='\t')
            for line in permitted_reader:
                if line[1] == self.gene:
                    # Split the promoters and enhancers
                    promoters = line[2].split(',')
                    enhancers = line[3].split(',')

                    # Put every promoter in the queues
                    for pr in promoters:
                        print 'Sending promoters to the queues'
                        self.hist_queue.put({'promoter': pr})
                        self.meth_queue.put({'promoter': pr})

                    # Put every enhancer in the queues
                    for en in enhancers:
                        print 'Sending enhancers to the queues'
                        # Change the coordinate format from chr#:#-# to chr#:#..#
                        self.hist_queue.put({'enhancer': en.replace('-', '..')})
                        self.meth_queue.put({'enhancer': en.replace('-', '..')})

                    break
        # Put a message to the poison queue to see it has finished
        self.poison_queue.put('promoters_enhancers')

    def exons_introns(self):
        # Search for the gene name in the gtf file
        print 'Opening {} file'.format(self.gtf_file)
        with open(self.gtf_file) as handle:
            gtf_reader = csv.reader(handle, delimiter='\t')
            for line in gtf_reader:
                # The gene name is in the 8th column of each line 4th entry of the semicolon separated string
                gtf_info = line[8].split(';')
                gene_name = re.search(r' gene_name \"(.*)\"', gtf_info[4]).group(1)

                # If there is a match put the gene's elements in the queues with the corresponding label
                if gene_name == self.gene:
                    # For introns
                    if line[2] == 'intron':
                        print 'Sending introns to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'intron': coordinates})
                        self.meth_queue.put({'intron': coordinates})

                    # For exons
                    elif line[2] == 'exon' and re.search(r' transcript_type \"(.*)\"',
                                                         gtf_info[5]).group(1) != 'retained_intron':
                        print 'Sending exons to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'exon': coordinates})
                        self.meth_queue.put({'exon': coordinates})

                    # For retained introns
                    elif line[2] == 'exon' and re.search(r'transcript_type \"(.*)\"',
                                                         gtf_info[5]).group(1) == 'retained_intron':
                        print 'Sending retained introns to the queues'
                        coordinates = '{}:{}..{}'.format(line[0], line[3], line[4])
                        self.hist_queue.put({'retained_intron': coordinates})
                        self.meth_queue.put({'retained_intron': coordinates})

        # Put a message to the poison queue to signal it has finished
        self.poison_queue.put('exons_introns')

    def poison(self):
        producers = []
        while True:
            # Get the messages from the poison queue
            producers.append(self.poison_queue.get())

            # Once both producers have finished put the poison pills
            if ('promoters_enhancers' and 'exons_introns') in producers:
                for pill in xrange(cpu_count()):
                    self.hist_queue.put(None)
                    self.meth_queue.put(None)

                break


class Consumer(object):
    def __init__(self, hist_queue, meth_queue, res_queue):
        self.hist_dir = pargs.hist_dir
        self.meth_dir = pargs.meth_dir
        self.hist_queue = hist_queue
        self.meth_queue = meth_queue
        self.results_queue = res_queue

        # Start the processes
        hist_process = Process(target=self.hist_worker)
        hist_process.start()
        meth_process = Process(target=self.meth_worker)
        meth_process.start()

    def hist_worker(self):
        while True:
            print 'Is hist_queue empty? {}'.format(self.hist_queue.empty())
            # Get a job and break in case it is a poison pill
            job = self.hist_queue.get()
            print 'Starting with job {}'.format(job)
            if job is None:
                self.hist_queue.task_done()
                break

            # There should only be one key in the dictionary from the queue
            element_type = job.keys()[0]
            coords = job[element_type]
            element_info = re.search(r'(.*):(\d*)\.\.(\d*)', coords)

            # Extract the element information
            chromosome = element_info.group(1)
            element_start = float(element_info.group(2))
            element_end = float(element_info.group(3))

            # The folders in the histones directory should have the tissues names/ids
            hist_walker = os.walk(self.hist_dir)

            # Walk through the directories
            for h_dir_path, h_dirs, h_files in hist_walker:

                print 'Walking through {}'.format(self.hist_dir)
                # Search for tagAlign files in each folder
                for track in h_files:
                    found_tagAlign = re.search(r'(.*)-(.*)\.tagAlign', track)

                    # If found, read the file
                    if bool(found_tagAlign):
                        print 'Found file {}'.format(track)
                        file_path = h_dir_path + '/' + track
                        hist_modification = found_tagAlign.group(2)

                        with open(file_path) as handle:
                            tagAlign_reader = csv.reader(handle, delimiter='\t')

                            # It will count the amount of modifications found in the element
                            modification_counter = 0

                            for tr in tagAlign_reader:
                                # Check if the chromosome is the right one
                                if tr[0] == chromosome:
                                    hist_start = float(tr[1])
                                    hist_end = float(tr[2])

                                    # If the modification is inside or in one of the borders of the element count it
                                    if (element_start <= hist_start < hist_end <= element_end or
                                        hist_start < element_start <= hist_end or
                                            hist_start <= element_end < hist_end):

                                        # Add to the modification counter
                                        print 'Adding to the {} counter'.format(hist_modification)
                                        modification_counter += 1

                            # Make a list with the results
                            print 'Got results for {}'.format(hist_modification)
                            results = [element_type, coords, hist_modification, modification_counter]

                            # Put the results in the queue as a dictionary
                            tissue = h_dir_path.split('/')[-1]
                            self.results_queue.put({tissue: results})
            print 'Task done with job {}'.format(job)
            self.hist_queue.task_done()

    def meth_worker(self):
        while True:
            print 'Is meth queue empty? {}'.format(self.meth_queue.empty())
            job = self.meth_queue.get()
            print 'Starting with job {}'.format(job)
            if job is None:
                self.meth_queue.task_done()
                break

            # There should only be one key in every dictionary from the queue
            element_type = job.keys()[0]
            coords = job[element_type]
            element_info = re.search(r'(.*):(\d*)\.\.(\d*)', coords)

            # Extract the element information
            chromosome = element_info.group(1)
            element_start = float(element_info.group(2))
            element_end = float(element_info.group(3))

            # The folders in the methylation directory should have the tissues names/ids
            meth_walker = os.walk(self.meth_dir)

            # Walk through the directories
            for m_dir_path, m_dirs, m_files in meth_walker:

                print 'Walking through {}'.format(self.meth_dir)
                # Look for a bedgraph file in each folder
                for track in m_files:
                    found_bedgraph = re.search(r'(.*)\.bedgraph', track)

                    if bool(found_bedgraph):
                        file_path = m_dir_path + '/' + track
                        print 'Found file {}'.format(track)

                        # Read the file
                        with open(file_path) as handle:
                            bedgraph_reader = csv.reader(handle, delimiter='\t')

                            # List to get the sum of the methylation scores
                            meth_score_sum = []

                            for br in bedgraph_reader:
                                # Check that the chromosome is the right one
                                if br[0] == chromosome:
                                    meth_start = float(br[1])
                                    meth_end = float(br[2])
                                    meth_score = float(br[3])

                                    # If the methylated region is inside the element or in any border of the element
                                    if (element_start <= meth_start < meth_end <= element_end or
                                        meth_start < element_start <= meth_end or
                                            meth_start <= element_end < meth_end):

                                        # Add the score to the list
                                        print 'Calculating sums for {}'.format(job)
                                        meth_score_sum.append(meth_score)

                            # Calculate the sum of the scores
                            meth_score_total = sum(meth_score_sum)

                            # Make a list with the results
                            print 'Got results for {}'.format(job)
                            results = [element_type, coords, 'Methylation', meth_score_total]

                            # Put the results in the queue as a dictionary
                            tissue = m_dir_path.split('/')[-1]  # The containing folder has the name of the tissue
                            self.results_queue.put({tissue: results})

                        # It will only use the first bedgraph file found
                        break
            print 'Ending task with job {}'.format(job)
            self.meth_queue.task_done()


def expression(gene, permitted, gene_info, eg_name, expr_dir, tissue_dir, res_queue):
    # Match the gene name with the Ensembl ID
    with open(permitted) as perm_handle:
        perm_reader = csv.reader(perm_handle, delimiter='\t')

        for per in perm_reader:
            # In case the gene is found get its ID
            if per[1] == gene:
                print 'Getting Ensembl ID for the gene'
                gene_id = per[0]
                # Break to avoid changing the ID
                break

            # This is just to know if the name given is not valid
            else:
                gene_id = 'NA'

    print 'The gene ID is {}'.format(gene_id)

    # Get the coordinates for the gene
    with open(gene_info) as gi_handle:
        gi_reader = csv.reader(gi_handle, delimiter='\t')

        for gir in gi_reader:
            if gir[0] == gene_id:
                print 'Getting coordinates of the gene'
                coord = 'chr{}:{}..{}'.format(gir[1], gir[2], gir[3])
                print 'Coordinates are {}'.format(coord)
                break

    # Get the list of the tissues being used from the histones or methylation directory
    tissue_walker = os.walk(tissue_dir)
    tissues_names = tissue_walker.next()[1]

    # Match the corresponding tissue names with their IDs
    tissues = {} # Keys are IDs and values are names
    with open(eg_name) as tissue_handle:
        eg_reader = csv.reader(tissue_handle, delimiter='\t')

        # Iterate through the file to find the IDs that match the names
        for er in eg_reader:
            if er[1] in tissues_names:
                print 'Tissue {} matched'.format(er[1])
                # Save the matches in a dictionary
                tissues[er[0]] = er[1]

    print tissues

    # Find the expression level for the gene in each tissue
    for tissue_id in tissues.keys():
        # Set expression score as None
        expression_score = None

        # Walk through the expression directory
        expr_walker = os.walk(expr_dir)
        for dir_path, dirs, files in expr_walker:
            for track in files:
                found_pc = re.search(r'(.*)\.pc', track)
                found_nc = re.search(r'(.*)\.nc', track)
                found_rb = re.search(r'(.*)\.rb', track)

                if bool(found_pc):
                    print 'Found the protein-coding file'
                    file_path = dir_path + track
                    with open(file_path) as handle:
                        pc_reader = csv.reader(handle, delimiter='\t')

                        # Get the column number of the tissue
                        column = pc_reader.next().index(tissue_id)
                        print 'Protein-coding: tissue {}, index {}'.format(tissue_id, column)

                        for pcr in pc_reader:
                            # Find the row of the corresponding gene
                            if pcr[0] == gene_id:
                                print 'Gene found in protein-coding'
                                expression_score = float(pcr[column])
                                expression_result = ['gene', coord, 'Expression', expression_score]
                                res_queue.put({tissues[tissue_id]: expression_result})
                                break
                    # Break only in case of finding the gene in the file
                    if expression_score is not None:
                        break

                elif bool(found_nc):
                    print 'Found the non-coding file'
                    file_path = dir_path + track
                    with open(file_path) as handle:
                        nc_reader = csv.reader(handle, delimiter='\t')

                        # Get the columns number of the tissue
                        column = nc_reader.next().index(tissue_id)
                        print 'Non-coding: tissue {}, index {}'.format(tissue_id, column)

                        for ncr in nc_reader:
                            # Find the row of the corresponding gene
                            if ncr[0] == gene_id:
                                print 'Gene found in non-coding'
                                expression_score = float(ncr[column])
                                expression_result = ['gene', coord, 'Expression', expression_score]
                                res_queue.put({tissues[tissue_id]: expression_result})
                                break
                    # Break only in case of finding the gene in the file
                    if expression_score is not None:
                        break

                elif bool(found_rb):
                    print 'Found the ribosomal file'
                    file_path = dir_path + track
                    with open(file_path) as handle:
                        rb_reader = csv.reader(handle, delimiter='\t')

                        # Get the columns number of the tissue
                        column = rb_reader.next().index(tissue_id)
                        print 'Ribosomal: tissue {}, index {}'.format(tissue_id, column)

                        for rbr in rb_reader:
                            # Find the row of the corresponding gene
                            if rbr[0] == gene_id:
                                print 'Gene found in ribosomal'
                                expression_score = float(rbr[column])
                                expression_result = ['gene', coord, 'Expression', expression_score]
                                res_queue.put({tissues[tissue_id]: expression_result})
                                break
                    # Break only in case of finding the gene in the file
                    if expression_score is not None:
                        break

def relief(rel_queue):
    # Make the dictionaries for tissue and chromosomes
    tissue_ids = {'Small_Intestine': 'SMI', 'Pancreas': 'PAN', 'Adult_Liver': 'ADL',
                  'Spleen': 'SPL', 'Thymus': 'THY', 'Lung': 'LUN', 'Esophagus': 'ESO'}
    chromosome_num = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10,
                      '11': 11, '12': 12, '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19,
                      '20': 20, '21': 21, '22': 22, 'X': 23, 'Y': 24, 'M': 25}
    # Schema of the database
    schema = '(Gene TEXT, Tissue TEXT, Feature TEXT, Site TEXT, Value REAL, Chr INT, Start INT, End INT)'
    with sqlite3.connect(pargs.db_name) as con:
        cursor = con.cursor()
        # Creates the table in case it doesn't exist
        table_name = 'D{}C{}P{}F{}'.format(pargs.dist, pargs.corr, pargs.pval, pargs.fdr)
        cursor.execute("CREATE TABLE IF NOT EXISTS " + table_name + schema)

        while True:
            # If the queue is not empty get a job
            if not rel_queue.empty():
                print 'Relief queue is not empty'
                rel_job = rel_queue.get()
                # Append it to the database
                output_job = re.search(r"\{'(.*)': \['(.*)', 'chr(.*):(\d*)\.\.(\d*)', '(.*)', (.*)\]\}", str(rel_job))
                cursor.execute("INSERT INTO " + table_name + " VALUES(?, ?, ?, ?, ?, ?, ?, ?)",
                               (pargs.gene_name, tissue_ids[output_job.group(1)], output_job.group(6),
                                output_job.group(2), float(output_job.group(7)), chromosome_num[output_job.group(3)],
                                int(output_job.group(4)), int(output_job.group(5))))

            # If the queue is empty check for active children; 0 children means its not going to get any more jobs
            elif len(active_children()) == 0:
                print 'Relief queue is empty and there are no more children'
                break

if __name__ == '__main__':

    start_time = time.time()

    # Create the required queues
    histone_queue = JoinableQueue(maxsize=100)
    methylation_queue = JoinableQueue(maxsize=100)
    results_queue = Queue()

    # Search for the permitted genes file if it is None
    permitted_file = pargs.permitted_file
    if permitted_file is None:
        print 'Searching for the permitted file'
        permitted_walker = os.walk('./')
        for p_dir_path, p_dirs, p_files in permitted_walker:
            for file_name in p_files:
                found_permitted = re.search(r'IMPermitted-D(.*)C(.*)P(.*)F(.*)\.txt', file_name)

                if bool(found_permitted):
                    permitted_file = file_name
                    break

        # If the permitted file is not found, make one and use it
        if permitted_file is None:
            print 'File not found, creating one'
            filteredDatabase(pargs.gene_info, pargs.prom_ann, pargs.association,
                             pargs.dist, pargs.corr, pargs.pval, pargs.fdr)

            permitted_file = 'IMPermitted-D{}C{}P{}F{}.txt'.format(pargs.dist, pargs.corr, pargs.pval, pargs.fdr)

    # The amount of real processes will be cpu_count() * 2 (see the Consumer class)
    consumer_processes = cpu_count()

    # Start the producers
    Producer(pargs.gene_name, permitted_file, histone_queue, methylation_queue)
    print 'Creating {} Consumers'.format(consumer_processes * 2)
    # Start the consumers
    for i in xrange(consumer_processes):
        Consumer(histone_queue, methylation_queue, results_queue)

    print active_children()

    # Start a process to get expression levels
    print 'Getting expression scores'
    expression_process = Process(target=expression, args=(pargs.gene_name, permitted_file, pargs.gene_info,
                                                          pargs.eg_name, pargs.expression_dir, pargs.hist_dir,
                                                          results_queue,))
    expression_process.start()

    # Relief the results queue to avoid overloading
    print 'Dropping the results'
    relief(results_queue)

    # Wait for the consumers to stop working
    histone_queue.join()
    methylation_queue.join()
    print 'Queues joined'

    print (time.time() - start_time) / 60.0, 'minutes'