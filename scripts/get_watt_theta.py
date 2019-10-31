import sys

#1st argument -- input filename
#2nd argument -- number of column containing snp number (0-based)
#3rd argument -- output file prefix
#4th argument -- number of individuals
#Takes file contaning newline-delimited snp numbers
#Returns input file with appended watterson's theta (unstandardizes) values
#Assumes diploidy and space-delimited input file
#Changed to tab-delimited
"""
give: sample size (number of individuals, not number of chromosomes
assumes diploidy
"""
def get_watt_theta(n_segsites, n_inds):
	n_chrom = n_inds * 2
	a = 0
	for i in range(n_chrom - 1):
		i += 1
		a += 1/i
	watt_theta = n_segsites / a
	return(watt_theta)

col = int(sys.argv[2])

with open(sys.argv[1]) as inf, open(sys.argv[3] + '.txt', 'w') as outf:
	outf.writelines(inf.readline().rstrip() + '\twatt_theta\n')
	for line in inf:
		line = line.rstrip()
		n_snps = int(line.split('\t')[col])
		watt_t = get_watt_theta(n_segsites = n_snps, n_inds = int(sys.argv[4]))
		outf.writelines(line + '\t' + str(watt_t) + '\n')
