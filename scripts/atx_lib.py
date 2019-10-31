#Compute divergence from chained & netted alignments e.g.: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsRn6/
"""
rm_ref_gaps(), get_divergence(), scrape_vcf() and rm_segsites() required for ENTRY class
ENTRY class required for ATX class
"""

"""
Function stripping reference gaps from an ENTRY_OBJ. Returns a new ENTRY_OBJ.
"""

def rm_ref_gaps(entry_obj):
    seq_1 = ''
    seq_2 = ''
    for i in range(len(entry_obj.seq_1)):
        if entry_obj.seq_1[i] == '-':
            continue
        else:
            seq_1 += entry_obj.seq_1[i]
            seq_2 += entry_obj.seq_2[i]
    new_obj = entry(seq_1 = seq_1, seq_2 = seq_2, ref_start = entry_obj.ref_start, ref_end = entry_obj.ref_end, ref_len = entry_obj.ref_len, align_start = entry_obj.align_start, align_end = entry_obj.align_end, align_len = entry_obj.align_len)
    return(new_obj)
    
"""
Outputs list of positions that are different between alignments in dictionary: {genomic_position:(ref_bp, align_bp)}
Assumes gaps filtered from ref. If not, will filter gaps from ref.
"""

def get_divergence(entry_obj):
    diff_dict = {}
    if '-' in entry_obj.seq_1:
        entry_obj = rm_ref_gaps(entry_obj)
    for i in range(entry_obj.ref_len):
        if entry_obj.seq_1[i].upper() != entry_obj.seq_2[i].upper():
            if entry_obj.seq_2[i] != '-':
                diff_dict[entry_obj.ref_start + i] = (entry_obj.seq_1[i], entry_obj.seq_2[i])
    return(diff_dict)
    
"""
Obtains list of variant positions from VCF
Assumes that VCF contains only variant sites, no indels, chromosome 
mouse -- sites polymorphic in any three populations
"""  

def scrape_vcf(VCF_path, chr):
    var_pos_ls = []
    with open(VCF_path) as vcf:
        for line in vcf:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                if line[0] == chr:
                    pos = int(line[1])
                    var_pos_ls.append(pos)
    return(var_pos_ls)

"""
Given dictionary of differing site positions, remove sites that are polymorphic
SET of variant positions obtained via scrape_vcf
"""

#THIS VER SLOWER THAN CASTING INTO SETS (by ~3X)
# def rm_segsites(diff_dict, var_pos_ls):
#     for pos in list(diff_dict.keys()):
#         if pos in var_pos_ls:
#             diff_dict.pop(pos)
#     return(diff_dict)

#Note, requires coercion of variant list to set before calling 
def rm_segsites(diff_dict, var_pos_set):
    diff_pos_set = set(list(diff_dict.keys()))
    not_var_diff_pos = diff_pos_set - var_pos_set
    var_diff_pos = diff_pos_set - not_var_diff_pos
    for pos in var_diff_pos:
        diff_dict.pop(pos)
    return(diff_dict)


"""
Counts number of non-indel sites after removal of reference gaps via rm_ref_gaps
This is the number that the count of divergent sites will be divided by to obtain divergent proportion
"""

def get_site_counts(entry_obj):
    site_count = 0
    if '_' in entry_obj.seq_1:
        entry_obj = rm_ref_gaps(entry_obj)
    for i in range(entry_obj.ref_len):
        if entry_obj.seq_2[i] != '_':
            site_count += 1
    return(site_count)
            
"""
Gets divergence proportion by dividing the number of divergent sites (length of diff_dict output by get_divergence / count output by get_site_counts)
Assumes that reference gaps and polymorphic sites have been removed from diff dict
"""

def get_div_prop(diff_dict, site_count):
    n_div = len(diff_dict)
    div_prop = n_div / site_count
    return(div_prop)

"""
Class storing information for a single .atx entry. Retuned by ATX_OBJ.get_entry().
for var_pos_ls, submit a list representing the set of variant positions for all populations of interest
"""

class entry(object):
    def __init__(self, seq_1, seq_2, ref_start, ref_end, ref_len, align_start, align_end, align_len):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_len = ref_len
        self.align_start = align_start
        self.align_end = align_end
        self.align_len = align_len
    def get_diffs(self, var_pos_ls):
        diff_dict = get_divergence(self)
        diff_dict = rm_segsites(diff_dict, var_pos_ls)
        self.diff_dict = diff_dict
        site_count = get_site_counts(self)
        self.site_count = site_count
        div_prop = get_div_prop(diff_dict, site_count)
        self.div_prop = div_prop         

"""
Class storing information from .atx file (format detailed at: http://genome.ucsc.edu/goldenPath/help/axt.html).
"""

class atx(object):
    def __init__(self, path):
        with open(path, 'r') as inf:
            entry_count = 0
            seq_len = 0
            entry_seq_len_ls = []
            #dictionary of entries -- keys: entry number. Values: tuple containing entry header followed by aligned sequences (3 elements in this tuple)
            entry_dict = {}
            for line in inf:
                if line.startswith('#'):
                    continue
                elif len(line.rstrip().split(' ')) > 1:
                    entry_count += 1
                    seq = inf.readline().rstrip()
                    seq_len += len(seq)
                    entry_seq_len_ls.append(len(seq))
                    seq_2 = inf.readline().rstrip()
                    entry_dict[int(line.rstrip().split(' ')[0])] = (line.rstrip(), seq, seq_2)
            self.entrycount = entry_count
            self.seqlen = seq_len
            self.seqlen_ls = entry_seq_len_ls
            self.entries = entry_dict
    
    """
    .get_entry() returns ENTRY_OBJ for a single ATX_OBJ entry
    """
            
    def get_entry(self, entry_number):
        seq_1 = self.entries[entry_number][1]
        seq_2 = self.entries[entry_number][2]
        header = self.entries[entry_number][0].split(' ')
        ref_start = int(header[2])
        ref_end = int(header[3])
        align_start = int(header[5])
        align_end = int(header[6])
        ref_len = ref_end - ref_start + 1
        align_len = align_end - align_start + 1
        return(entry(seq_1 = seq_1, seq_2 = seq_2, ref_start = ref_start, ref_end = ref_end, ref_len = ref_len, align_start = align_start, align_end = align_end, align_len = align_len))
    
    """
    .rm_entry() removes an entry from ATX_OBJ
    """
    
    def rm_entry(self, entry_number):
        self.entries.pop(entry_number)

"""
Function subsetting an ENTRY_OBJ
entry_obj -- an ENTRY_OBJ obtained via ATX_OBJ.get_entry()
    NOTE: expects an entry object w/ reference gaps removed via rm_ref_gaps()
    if reference sequence contains gaps, rm_ref_gaps will be run
stpos -- start position for the reference(gaps in reference alignment not counted)
endpos -- end position for the reference
"""

def subset_entry(entry_obj, stpos, endpos):
    if '-' in entry_obj.seq_1:
        entry_obj = rm_ref_gaps(entry_obj)
    correction = entry_obj.ref_start #by subtracting this correction from stpos and endpos, can obtain postion indexes for subsetting reference and alignment sequences
    slice_start = stpos - correction #distance of subset start from original genomic starting position
    slice_end = endpos - correction #distance of subset end from original genomic starting position
    seq_1 = entry_obj.seq_1[slice_start:slice_end + 1] #add one to subset over closed interval stpos:endpos
    seq_2 = entry_obj.seq_2[slice_start:slice_end + 1]
    ref_start = stpos
    ref_end = endpos
    align_start = entry_obj.align_start + slice_start
    align_end = entry_obj.align_start + slice_end
    ref_len = endpos - stpos + 1
    align_len = align_end - align_start + 1
    new_obj = entry(seq_1 = seq_1, seq_2 = seq_2, ref_start = ref_start, ref_end = ref_end, ref_len = ref_len, align_start = align_start, align_end = align_end, align_len = align_len)
    return(new_obj) 

"""
Pool entries containing data along a specified interval
Outputs a dictionary of entry objects (for now)
remove = False
if remove == True, remove entries completely falling w/in specified interval
NOTE: only removes entry from entries dictionary. Doesn't modify other ATX_OBJ attributes. Should be used with caution. 
"""

def pool_entries(atx_obj, int_start, int_end, remove = False):
    entries_in_interval = {}
    rm_ls = []
    key_ls = list(atx_obj.entries.keys()) #need to do this b/c dictionary size changes during loop
    for i in key_ls:
        curr_ent = rm_ref_gaps(atx_obj.get_entry(i))
        #if an entry begins after interval start
        if curr_ent.ref_start > int_start and curr_ent.ref_start < int_end:
            #if no entry truncation is necessary
            if curr_ent.ref_end < int_end:
                entries_in_interval[i] = curr_ent
                if remove == True:
                    atx_obj.rm_entry(i)
                continue
            #if entry extends beyond interval end (truncation required)
            elif curr_ent.ref_end > int_end:
                curr_ent_trunc = subset_entry(curr_ent, curr_ent.ref_start, int_end)
                entries_in_interval[i] = curr_ent_trunc
                continue
        #if an entry begins before interval start, but ends after interval start      
        elif curr_ent.ref_start < int_start and curr_ent.ref_end > int_start:
            curr_ent_trunc = subset_entry(curr_ent, int_start, curr_ent.ref_end)
            entries_in_interval[i] = curr_ent_trunc
            if remove == True:
                atx_obj.rm_entry(i)
            continue
        #if entry begins after interval end, break out of loop
        elif curr_ent.ref_start > int_end:
            break
    return(entries_in_interval)
    
"""
Function to compute divergence proportion for pooled entries
pool_dict -- output from pool_entries
"""
                
def compute_pooled_div(pool_dict, var_pos_ls):
    diff_count = 0
    site_count = 0
    for v in pool_dict.values():
        v.get_diffs(var_pos_ls = var_pos_ls)
        diff_count += len(v.diff_dict)
        site_count += v.site_count
    try:
        pool_div = diff_count / site_count
    except ZeroDivisionError:
        pool_div = 0
    return(diff_count, site_count, pool_div)
        











