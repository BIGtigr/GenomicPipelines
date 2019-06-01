##########################################################################################
#About this script: 
#	Version : 20180330/001
#	Author  : Dahu Zou
#	Description : this script is used to find olfactory receptor genes in genomes
##########################################################################################
import sys,os,glob
#import NGStool
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--path', help=('absolute path to genome sequence file'))
parser.add_option('--query', help=('amino acid sequence'))
parser.add_option('--cpu', type='int', default=30,help='numbers of threads')
usage = u'''`Usage:python %prog [options...] [path] [query] [cpu]'''
(options, args) = parser.parse_args()


numOpts = len(sys.argv)
if numOpts < 2:
    parser.print_help()
    sys.exit()

###########################################################################
CODE = {
        "standard":{
                                'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A',                                 # Alanine
                                'TGC' : 'C', 'TGT' : 'C',                                                           # Cysteine
                                'GAC' : 'D', 'GAT' : 'D',                                                           # Aspartic Acid
                                'GAA' : 'E', 'GAG' : 'E',                                                           # Glutamic Acid
                                'TTC' : 'F', 'TTT' : 'F',                                                           # Phenylalanine
                                'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G',                                 # Glycine
                                'CAC' : 'H', 'CAT' : 'H',                                                           # Histidine
                                'ATA' : 'I', 'ATC' : 'I', 'ATT' : 'I',                                              # Isoleucine
                                'AAA' : 'K', 'AAG' : 'K',                                                           # Lysine
                                'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L', 'TTA' : 'L', 'TTG' : 'L',       # Leucine
                                'ATG' : 'M',                                                                        # Methionine
                                'AAC' : 'N', 'AAT' : 'N',                                                           # Asparagine
                                'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',                                 # Proline
                                'CAA' : 'Q', 'CAG' : 'Q',                                                           # Glutamine
                                'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R', 'AGA' : 'R', 'AGG' : 'R',       # Arginine
                                'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S', 'AGC' : 'S', 'AGT' : 'S',       # Serine
                                'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',                                 # Threonine
                                'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V',                                 # Valine
                                'TGG' : 'W',                                                                        # Tryptophan
                                'TAC' : 'Y', 'TAT' : 'Y',                                                           # Tyrosine
                                'TAA' : 'X', 'TAG' : 'X', 'TGA' : 'X'                                               # Stop
                                }
        ## more translate table could be added here in future
}

#ATCGatcgRYMKrymk
#TAGCtagcYRKMyrkm
complement_table = {'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','R':'Y','Y':'R','M':'K','K':'M','r':'y','y':'r','m':'k','k':'m','W':'W','N':'N','n':'n'}
stop_codon_p = ['TAA','TAG','TGA']
start_codon_p = ['ATG']
stop_codon_n = ['AAT','GAT','AGT']
start_codon_n = ['GTA']
##############################################################################
# read seq length of fasta file into dict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fas2length(fasfile):
	seqlen = {}
	i = 0
	with open(fasfile) as fh:
		for line in fh.readlines():
			i+=1
			#print (i)
			line = line.strip()
			if line.startswith('>'):
				seq_title = line.replace('>','')
				#print (seq_title)
				seqlen[seq_title] = 0
			else:
				seqlen[seq_title] += len(line)
	return seqlen

#############################################################################
# read fasta file into dict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fasta2dict(fFile):
        (defLine, seq) = ('', '')
        d = {}
        for l in fFile:
		l = l.strip()
                if '>' == l[0]: # a new def line
                        if defLine: # record previous def line and seq
                                d[defLine[1:].split()[0]] = (defLine, seq)
                        (defLine, seq) = (l, '')
                else: # more of the sequence
                        seq += l
        if defLine: d[defLine[1:].split()[0]] = (defLine, seq)
        return d

##############################################################################
# parsing blast result
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def blast_parser(tab,seq_len):
	dicta = {}
	with open(tab) as TAB:
		for line in TAB:
			lines = line.split()
			pos = {}
			start = int(lines[8])
			end = int(lines[9])
			extend_length = 3*seq_len[lines[0]] - 3*(int(lines[7]) - int(lines[6]) + 1)+30
			if lines[1] in dicta:
				if lines[8] < lines[9]:
					dicta[lines[1]].append([start,end,extend_length])
				else:
					dicta[lines[1]].append([start,end,extend_length])
			else:
				dicta[lines[1]] = []
				if lines[8] < lines[9]:
					dicta[lines[1]].append([start,end,extend_length])
				else:
					dicta[lines[1]].append([start,end,extend_length])
	return dicta

##############################################################################
# reverse
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def reverse(seq):
	reverse_out = seq[::-1]
	return(reverse_out)

##############################################################################
# complement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def complement(seq):
	complement_out = ''.join([complement_table[i] for i in seq])
	return complement_out

##############################################################################
#  extract seq seed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extractor(seq,start,end):
	son_seq = seq[start:end]
	return son_seq


##############################################################################
# sort 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def sort(inlist):
	out = sorted(inlist,cmp = lambda x,y: cmp(x[2],y[2]))
	return out

##############################################################################
# are two locus belongs to one gene ?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def isone(locus1,locus2):
	if abs(locus1[2] - locus2[2]) < 1500:
		return True
	else:
		return False
##############################################################################
#  merge two locus
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def merge(locus1,locus2):
	if abs(locus1[0]-locus1[1]) > abs(locus2[0]-locus2[1]):
		return locus1
	else:
		return locus2
#############################################################################
# remove duplicate locus of on a scaffold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter(inlist):
	list0 = []
	list1 = sort(inlist)
	temp = list1[0]
	for i in range(1,len(list1)):
		if isone(temp,list1[i]):
			temp = merge(temp,list1[i])
		else:
			list0.append(temp)
			temp = list1[i]
	list0.append(temp)
	return list0
##############################################################################
# find codon that was given leftward
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def leftcodon(seq,point,codon,length):
	point = int(point)
	n_site = point
	end = point - length
	if end < 0:
		end = 0
	for i in range(point-1,end,-3):
		#print (i),
		NNN = seq[i-3:i]
		#print ('\t'),
		#print (NNN)
		if NNN.upper() in codon:
			n_site = i-2
	return n_site

##############################################################################
# go right to search codon which was given
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rightcodon(seq,point,codon,length):
	n_site = point
	start = int(point)
	end = point + length
	if point + length > len(seq):
		end = len(seq)
	for i in range(point,end,3):
		NNN = seq[i:i+3]
		if NNN in codon:
			n_site = i+3
			#print ('\t\t%s\t%s' %(n_site,NNN)),
	return n_site


############################################################################
# extend for start or stop codon
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extend(seq,hit):
	new_hit = [0,0,0]
	seq_p = seq
	seq_n = complement(seq_p)
	start = hit[0]
	end = hit[1]
	extend_len = hit[2]
	if int(start) < int(end):
		hit[0] = leftcodon(seq_p,start+3,start_codon_p,extend_len)
		hit[1] = rightcodon(seq_p,end,stop_codon_p,extend_len)
	else:
		hit[0] = leftcodon(seq_n,end,stop_codon_n,extend_len)
		hit[1] = rightcodon(seq_n,start-3,start_codon_n,extend_len)
	new_hit[2] = max(hit[0],hit[1])-abs(dicta[scaffold][i][1] - dicta[scaffold][i][0])/2
	return new_hit
	
	
#############################################################################
# main body
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################

genomes = glob.glob('%s/*fasta' %options.path)
seqlen = fas2length(options.query)

for genome in genomes:
	if not genome.endswith('fasta'):
		continue
	FH = open(genome)
	name2seq = fasta2dict(FH)
	FH.close()
	genome_name = genome.split('/')[-1]
	os.mkdir(genome_name)
	os.chdir(genome_name)
	mkdbCMD = 'makeblastdb -in %s -dbtype nucl -title tblastdb -parse_seqids -out ./genome -logfile makeblastdb.log' %genome
	tblastnCMD = 'tblastn -num_threads %s -evalue 1e-10 -db genome -outfmt 6 -query %s -out tblastn.f6' %(options.cpu,options.query)
	os.system(mkdbCMD)
	os.system(tblastnCMD)
	
	dicta = blast_parser('./tblastn.f6',seqlen)
	for scaffold in dicta:
		seq_p = name2seq[scaffold][1]
		i = 0
		for locus in dicta[scaffold]:
			dicta[scaffold][i] = extend(seq_p,locus)
			i += 1

	out = open('result.cds','w')
	for scaffold in dicta:
		seq_p = name2seq[scaffold][1]
		list_filtered = filter(dicta[scaffold])
		for locus in list_filtered:
			if locus[0] < locus[1]:
				or_seq = extractor(seq_p,locus[0]-1,locus[1])
				out.write('>%s|%s|%s' %(scaffold,locus[0],locus[1]))
				out.write(or_seq)
			if locus[0] > locus[1]:
				or_seq = extractor(seq_p,locus[1]-1,locus[0])
				or_seq = reverse(complement(or_seq))
				out.write('>%s|%s|%s' %(scaffold,locus[0],locus[1]))
				out.write(or_seq)
	out.close()
	os.chdir('../')

#########################################################################
# the end 
#########################################################################
