"""
Script computing nucleotide diversity over a range of window sizes
Input file is an Unzipped VCF, containing exclusively individuals to be considered for analysis
"""
import time
import sys

#Class to hold output of grab_chunk, etc.
class big_chunk():
    def __init__(self, snp_pos_lis, haplo_lis):
        self.snp_pos_lis = snp_pos_lis
        self.haplo_lis = haplo_lis

#Function to pull list of transcribed sites in an interval from a gene table
#Use inside of grab_chunk
def get_tx_win(path, start, stop, chr):
    tx_start = 0
    tx_end = 0
    pos_lis = list()
    with open(path) as gtable:
        for line in gtable:
            line = line.rstrip().split('\t')
            if line[2] == chr:
                #Starts within window and ends within window
                if int(line[4]) >= start and int(line[5]) <= stop:
                    tx_start = int(line[4])
                    tx_end = int(line[5])
                #Starts before window and ends within window
                elif int(line[4]) < start and int(line[5]) >= start and int(line[5]) <= stop:
                    tx_start = start
                    tx_end = int(line[5])
                #Starts before window and ends outside of window
                elif int(line[4]) < start and int(line[5]) > stop:
                    tx_start  = start
                    tx_end = stop
                #Starts within window and ends outside of window
                elif int(line[4]) >= start and int(line[4]) <= stop and int(line[5]) > stop:
                    tx_start = int(line[4])
                    tx_end = stop
                #Falls completely outside of window
                #Greater than(break from loop and return list) -- this requires a sorted gene table!
                elif int(line[4]) >= stop:
                    break
                else:
                    continue
                add_pos_lis = list(range(tx_start, tx_end + 1))
                pos_lis.extend(add_pos_lis)
    #cast list into set to remove duplicate entries
    pos_set = set(pos_lis)
    #return as sorted list
    pos_lis = list(pos_set)
    pos_lis.sort()
    return(pos_lis)

#Function to grab 5Mb from VCF
def grab_chunk(path, start, stop, chr, exclude_tx = False, g_table_path = ''):
    #if excluding tx sites, start by retrieving list of transcribed sites
    if exclude_tx is True:
        tx_pos_lis = get_tx_win(path = g_table_path, start = start, stop = stop, chr = chr)
    elif exclude_tx is False:
        tx_pos_lis = []
    with open(path) as inf:
        for line in inf:
            if line.startswith('##'):
                continue
            #Count number of individuals in VCF
            elif line.startswith('#CHROM'):
                line = line.rstrip().split('\t')
                #Get ind number by subtracting number of non-genotype columns
                n_ind = len(line) - 9
                #initialize haplotype list
                haplo_lis = []
                for i in range(n_ind * 2):
                    haplo_lis.append('')
                #initialize snp position list
                snp_pos_lis = []
            else:
                line = line.rstrip().split('\t')
                #check chromosome match
                if line[0] == chr and int(line[1]) >= start and int(line[1]) <= stop:
                    #check whether site is transcribed
                    #if transcribed, insert 'X' into haplotype sequence
                    if int(line[1]) in tx_pos_lis:
                        for i in range(len(haplo_lis)):
                            haplo_lis[i] += 'X'
                        #record snp position for these as well
                        snp_pos_lis.append(int(line[1]))
                    else:
                        geno_lis = line[9:]
                        #tidy up genotype list
                        #also split on '/' before appending to haplotype list
                        haplo_lis_append = []
                        for i in geno_lis:
                            haplo_lis_append.extend(i.split(':')[0].split('/'))
                        #append snp genotype info to haplotype list
                        for i in range(len(haplo_lis)):
                            haplo_lis[i] += haplo_lis_append[i]
                        #also record snp position
                        snp_pos_lis.append(int(line[1]))
                #break if on proper chromosome and position is greater than window stop position
                elif line[0] == chr and int(line[1]) > stop:
                    break
                else:
                    continue
    snp_pos_haplo_lis = big_chunk(snp_pos_lis, haplo_lis)
    return(snp_pos_haplo_lis)

#Missing data counter (single site)
def count_miss(haplo_lis, pos):
    #count missing data
    miss_data = 0
    for i in haplo_lis:
        if i[pos] == '.':
            miss_data += 1
    return(miss_data)

#Missing data threshold/removal (single site)
def rm_miss(max_miss_prop, haplo_lis):
    rm_pos_lis = []
    num_chroms = len(haplo_lis)
    #all haplotype strings in haplo_lis are same length, just use first to get site count
    num_sites = len(haplo_lis[0])
    for i in range(num_sites):
        miss_data = count_miss(haplo_lis = haplo_lis, pos = i)
        prop_miss = miss_data/num_chroms
        if prop_miss > max_miss_prop:
            rm_pos_lis.append(i)
        else:
            continue
    #edit haplo_lis to flag sites above missing data threshold
    haplo_lis_new = []
    for i in range(num_chroms):
        haplo_lis_new.append('')
    for i in range(num_sites):
        if i in rm_pos_lis:
            for j in range(num_chroms):
                haplo_lis_new[j] += 'M'
        else:
            for j in range(num_chroms):
                haplo_lis_new[j] += haplo_lis[j][i]
    return(haplo_lis_new)

#Function counting number of sites with insufficient data
def count_miss_rm(haplo_lis, haplo):
    miss_num = haplo_lis[haplo].count('M')
    return(miss_num)

#Function to count number of pairwise differences for a single site
#returns a tuple (diff_count, num_pairwise_comparisons)
def site_pair_diff(snp_lis, exclude_tx_miss = True):
    diff_count = 0
    num_pairwise_comparisons = 0
    for i in range(len(snp_lis)):
        #exclude sites above missing data threshold or transcribed sites. Also exclude individuals with missing data
        if exclude_tx_miss is True and snp_lis[i] == 'X' or snp_lis[i] == 'M' or snp_lis[i] == '.':
            continue
        else:
            #to avoid double-counting pairs, only use item in list greater than first
            for j in range(i + 1, len(snp_lis)):
                #exclude sites above missing data threshold or transcribed sites
                if exclude_tx_miss is True and snp_lis[j] != 'X' and snp_lis[i] != 'M' and snp_lis[i] != '.':
                    num_pairwise_comparisons += 1
                    if snp_lis[i] != snp_lis[j]:
                        diff_count += 1
    return((diff_count, num_pairwise_comparisons))

#Function iterating single-site pairwise diff counter over a chunk of sites
#Returns a list of tuples: [(diff_count, num_pairwise_comparisons), (diff_count, num_pairwise_comparisons), ... ]
#Entires correspond to values for individual SNPs
def sum_pair_diff(haplo_lis, exclude_tx_miss = True):
    diff_counts = []
    for i in range(len(haplo_lis[0])):
        snp_lis = []
        for j in haplo_lis:
            snp_lis.append(j[i])
        diff_counts.append(site_pair_diff(snp_lis = snp_lis, exclude_tx_miss = exclude_tx_miss))
    return(diff_counts)

#Function zipping diff count lis and snp pos lis into dictionary
def get_pos_diff_dict(snp_pos_lis, diff_counts):
    return(dict(zip(snp_pos_lis, diff_counts)))

#Function finding denominator for nucleotide diversity computation (number of sites considered in computation)
#Get number of sites in a window
#Excludes tx and missing data sites
def get_denom(haplo_lis, g_table_path, start, stop, chr):
    #look through haplo_lis, count number of transcribed sites and sites beyond missing data threshold
    #using get_tx_win, find number of TX sites
    #inefficient, should call get_tx_win outside of this function (redundant calls of get_tx_win otherwise)
    #counts all transcribed sites, including invariant transcribed sites. Difference between window size and this value is the proper denominator for nucletide diversity in non-transcribed sites.
    #include exception for empty haplo_lis
    win_size = stop - start + 1
    tx_lis = get_tx_win(path = g_table_path, start = start, stop = stop, chr = chr)
    try:
        miss_count = haplo_lis[0].count('M')
    except IndexError:
        miss_count = 0
    denom = win_size - len(tx_lis) - miss_count
    return(denom)

#function to subset a haplotype object, given start and stop positions
def subset_haplo_obj(haplo_obj, start, stop):
    #convert genomic start position to positoin in haplo_obj.haplo_lis
    #get minimum and maximum positions in haplo_obj.snp_pos_lis that are within defined interval boundaries
    start_counter = 0
    #index sometimes out of range for below line
    #if the defined window tries to start past the snp position list, will throw an error.
    #need to return zeros for when this happens
    try:
        while haplo_obj.snp_pos_lis[start_counter] < start:
            start_counter += 1
    except IndexError: #error when loop tries to start at a position past the furthest snp, return zeros
        sub_haplo_obj = big_chunk(snp_pos_lis = [], haplo_lis = [])
        return(sub_haplo_obj)
        exit()
    stop_counter = start_counter
    try:
        while haplo_obj.snp_pos_lis[stop_counter] < stop:
            stop_counter += 1 #one postition past the final snp position w/in defined interval
    except IndexError:
        stop_counter -= 1
    sub_haplo_lis = []
    for i in haplo_obj.haplo_lis:
        sub_haplo_lis.append(i[start_counter:stop_counter])
    sub_snp_pos_lis = haplo_obj.snp_pos_lis[start_counter:stop_counter]
    sub_haplo_obj = big_chunk(snp_pos_lis = sub_snp_pos_lis, haplo_lis = sub_haplo_lis)
    return(sub_haplo_obj)

#Function computing nucleotide diversity for a single window
#check division step
#include functionality to list several window sizes
#in result dictionary, include number transcirbed sites and number of SNPs considered for analysis
def get_nuc_divers(path, start, stop, chr, g_table_path = '', mult_wins = False, win_size_lis = []):
    start_time = time.time()
    #get snp chunk from VCF
    haplo_obj = grab_chunk(path = path, start = start, stop = stop, chr = chr, exclude_tx = True, g_table_path = g_table_path)
    #don't include code to remove missing missing data
    diff_counts = sum_pair_diff(haplo_lis = haplo_obj.haplo_lis)
    #Return a single value if not computing for multiple window sizes
    if mult_wins is False:
        n_nonTX = get_denom(haplo_lis = haplo_obj.haplo_lis, g_table_path = g_table_path, start = start, stop = stop, chr = chr)
        #divide mean diff count by number of nonTX sites
        sum = 0
        for i in diff_counts:
            #number of differences divided by number of comparisons (avg per-site number)
            #some entries will have zero comparisons, depending on upstream filtering
            if i[1] > 0:
                sum += i[0] / i[1]
        try:
            nuc_divers = sum / n_nonTX
        except ZeroDivisionError:
            nuc_divers = 'NA'
    #Return a dictionary if computing for multiple window sizes
    elif mult_wins is True:
        #get dictionary of diff counts and snp position
        #position is key, diff count is value
        pos_diff_dict = get_pos_diff_dict(snp_pos_lis = haplo_obj.snp_pos_lis, diff_counts = diff_counts)
        #initialize result dictionary
        nuc_divers = {}
        for win_size in win_size_lis:
            #initialize list of nucleotide diversity results for current window size
            nuc_divers_lis = []
            #compute number of sub-windows
            win_num = int((stop - start + 1) / win_size)
            #compute nucleotide diversity for each sub-window
            for subwin in range(win_num):
                subwin_start = start + (win_size * (subwin))
                subwin_stop = subwin_start + win_size - 1
                #subset haplo_lis
                sub_haplo_obj = subset_haplo_obj(haplo_obj = haplo_obj, start = subwin_start, stop = subwin_stop)
                n_nonTX = get_denom(haplo_lis = sub_haplo_obj.haplo_lis, g_table_path = g_table_path, start = subwin_start, stop = subwin_stop, chr = chr)
                #subset diff_counts
                sub_pos_diff_dict = {}
                for k,v in pos_diff_dict.items():
                    if k >= subwin_start and k <= subwin_stop:
                        sub_pos_diff_dict[k] = v
                #divide mean diff count by number of nonTX sites
                sum = 0
                snp_counter = 0
                for i in sub_pos_diff_dict.values():
                    if i[1] > 0:
                        sum += i[0] / i[1]
                        if i[0] > 0:
                            snp_counter += 1
                try:
                    #append tubple to list: (substart, stubstop, n_nontx, n_snps, nuc_divers)
                    nuc_divers_lis.append((subwin_start, subwin_stop, n_nonTX, snp_counter, (sum / n_nonTX)))
                except ZeroDivisionError: #throws ZeroDivisionError when all sites in a window are transcribed
                    nuc_divers_lis.append((subwin_start, subwin_stop, n_nonTX, snp_counter, 'NA'))
            nuc_divers[win_size] = nuc_divers_lis
    return(nuc_divers)

#function to initialize outfiles
#follows format of PREnucdivers_vcfPi_winsize.txt where 'PRE' is replaced by file prefix
def init_outfiles(win_size_lis, pre = ''):
    for size in win_size_lis:
        with open(pre + 'nucdivers_vcfPi_' + str(size) + '.txt', 'w') as outf:
            outf.writelines('chr\tstart\tstop\tn_nonTX_sites\tn_snp_nonTX\tnuc_divers\n')
            continue

#function to loop over window sizes, writing info in nuc_divers output dictionary to outfiles
#def write_res(n):
def write_res(res_dict, chr, pre = ''):
    #dictionary keys are window window sizes
    win_size_lis = list(res_dict.keys())
    for size in win_size_lis:
        with open(pre + 'nucdivers_vcfPi_' + str(size) + '.txt', 'a') as outf:
            for res in res_dict[size]:
                outf.writelines(chr + '\t' + str(res[0]) + '\t' + str(res[1]) + '\t' + str(res[2]) + '\t' + str(res[3]) + '\t' + str(res[4]) + '\n')

#function to get list of window positions for a given window size
#default start at pos 3000000
def define_wins(win_size, path, chr, start = 3000000):
    pos_lis = []
    with open(path) as vcf:
        for line in vcf:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                if line[0] == chr:
                    pos_lis.append(int(line[1]))
    last_pos = pos_lis[len(pos_lis) - 1]
    start_lis = list(range(start, last_pos, win_size))
    stop_lis = list(range(start + win_size - 1, last_pos, win_size))
    stop_lis.append(last_pos)
    res_lis = []
    for pos in range(len(start_lis)):
        res_lis.append((start_lis[pos], stop_lis[pos]))
    return(res_lis)

#function to loop nuc_divers over genome, writing continuously to outfiles (loops over entire VCF)
#assumes that a vcf contains data for single chromosome
def loop_nuc_divers(path, chr, g_table_path = '', mult_wins = False, win_size_lis = [], outfile_prefix = '', major_win = 1000000):
    #initialize init_outfiles
    print('Initiating Outfiles')
    init_outfiles(win_size_lis = win_size_lis, pre = outfile_prefix)
    #scrape VCF to get list of major window positions
    print('Finding Major Windows')
    major_win_pos = define_wins(win_size = major_win, path = path, chr = chr)
    #looping over major windows, run get_nuc_divers
    for win in major_win_pos:
        print('Window: ' + str(win))
        res_dict = get_nuc_divers(path = path, start = win[0], stop = win[1], chr = chr, g_table_path = g_table_path, mult_wins = mult_wins, win_size_lis = win_size_lis)
        #update output outfiles
        write_res(res_dict = res_dict, chr = chr, pre = outfile_prefix)


loop_nuc_divers(path = sys.argv[1], chr = sys.argv[2], g_table_path = sys.argv[3], mult_wins = True, win_size_lis = [2500, 5000, 50000, 100000, 500000, 1000000], outfile_prefix = sys.argv[4])
