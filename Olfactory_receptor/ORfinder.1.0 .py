##########################################################################################
#About this script: 
#	Version : 20180330/001
#	Author  : Dahu Zou
#	Description : this script is used to find olfactory receptor genes in genomes
##########################################################################################
import sys,os,glob,re
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
# parafly blast
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def paraBlast(query, database = 'genome', evalue = 1e-10, cpu = 30):
	os.mkdir('temp')
	i = 0
	temp = open('temp/query0.fas','w')
	cmd = open('blast.sh','w')
	with open(query) as QUERY:
		for line in QUERY:
			line = line.strip()
			if line.startswith('>'):
				i += 1
				temp.close()
				title = line.lstrip('>')
				temp = open('temp/query%s.fas' %i,'w')
				temp.write('>%s\n' %title)
				cmd.write('tblastn -num_threads 1 -evalue %s -db %s -outfmt 6 -query ./temp/query%s.fas -out ./temp/tblastn.%s\n' %(evalue,database,i,i))
			else:
				temp.write('%s\n' %line)
	cmd.close()
	temp.close()
	os.system('ParaFly -c blast.sh -CPU %s' %cpu)
	os.system('cat ./temp/tblastn* > tblastn.f6')

##############################################################################
# parsing blast result
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def blast_parser(tab,seq_len,add = 90):
	dicta = {}
	with open(tab) as TAB:
		for line in TAB:
			lines = line.split()
			pos = {}
			start = int(lines[8])
			end = int(lines[9])
			coordinate = max(start, end)-abs(start - end)/2
			extend_length1 = int(lines[6]) * 3 + add
			extend_length2 = (seq_len[lines[0]] - int(lines[7])) * 3 + add
			#extend_length = 3*seq_len[lines[0]] - 3*(int(lines[7]) - int(lines[6]) + 1)+30
			if lines[1] in dicta:
				if lines[8] < lines[9]:
					dicta[lines[1]].append([start,end,coordinate,extend_length1,extend_length2])
				else:
					dicta[lines[1]].append([start,end,coordinate,extend_length1,extend_length2])
			else:
				dicta[lines[1]] = []
				if lines[8] < lines[9]:
					dicta[lines[1]].append([start,end,coordinate,extend_length1,extend_length2])
				else:
					dicta[lines[1]].append([start,end,coordinate,extend_length1,extend_length2])
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
# is the same direction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def is_same_direction(locus1,locus2):
	val = (locus1[0] - locus1[1])*(locus2[0] - locus2[1])
	if val > 0:
		return True
	else:
		return False
	
##############################################################################
# are two locus belongs to one gene ?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def isone(locus1,locus2,dis = 900):
	if is_same_direction(locus1,locus2):
		if abs(locus1[2] - locus2[2]) < dis:
			return True
		else:
			return False
	else:
		return False

##############################################################################
# are two locus overlaped ? (another method)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def isoverlap(locus1,locus2,dis=0):
	if is_same_direction(locus1,locus2):
		if max(locus1[0],locus1[1]) > min(locus2[0],locus2[1]):
			return True
		else:
			return False
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
# remove duplicate locus of on a scaffold 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter(inlist,dis = 900):
	list0 = []
	list1 = sort(inlist)
	temp = list1[0]
	for i in range(1,len(list1)):
		if isone(temp,list1[i],dis):
			temp = merge(temp,list1[i])
		else:
			list0.append(temp)
			temp = list1[i]
	list0.append(temp)
	return list0

#############################################################################
#remove duplicate locus of on a scaffold 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def filter1(inlist,dis = 900):
        list0 = []
        list1 = sort(inlist)
        temp = list1[0]
        for i in range(1,len(list1)):
                if isoverlap(temp,list1[i],dis):
                        temp = merge(temp,list1[i])
                else:
                        list0.append(temp)
                        temp = list1[i]
        list0.append(temp)
        return list0
##############################################################################
# find codon that was given leftward
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def leftcodon(seq,point,codon,length,direction):
	if direction == '+':
		point = int(point)-3
	else:
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
def rightcodon(seq,point,codon,length,direction):
	if direction == '+':
		point = int(point)
	else:
		point = int(point)+3
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
	new_hit = [0,0,0,0,0]
	seq_p = seq
	seq_n = complement(seq_p)
	start = hit[0]
	end = hit[1]
	extend_len1 = hit[3]
	extend_len2 = hit[4]
	if int(start) < int(end):
		direction = '+'
		new_hit[0] = leftcodon(seq_p,start+3,start_codon_p,extend_len1,direction)
		new_hit[1] = rightcodon(seq_p,end,stop_codon_p,extend_len2,direction)
	else:
		direction = '-'
		new_hit[1] = leftcodon(seq_n,end,stop_codon_n,extend_len2,direction)
		new_hit[0] = rightcodon(seq_n,start-3,start_codon_n,extend_len1,direction)
	new_hit[2] = max(hit[0],hit[1])-abs(dicta[scaffold][i][1] - dicta[scaffold][i][0])/2
	return new_hit

#############################################################################	
# translate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def translate(seq):
	#aa = ''.join([CODE['standard'][seq[x:x+3].upper()] for x in range(0,len(seq),3)])
	a = []
	for x in range(0,len(seq),3):
		if seq[x:x+3].upper() in CODE['standard']:
			a.append(CODE['standard'][seq[x:x+3].upper()])
		else:
			a.append('N')
		aa = ''.join(a)
	return aa


#############################################################################
# is psudo or not
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ispsudo(seq):
	aa = translate(seq)
	if re.search('X',aa[1:-1]):
		return True
	else:
		return False

#############################################################################
#tag: intact, patial, or psudo
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tag(seq):
	tag = ''
	aa = translate(seq)
	if ispsudo(seq):
		tag = 'psudo'
	elif aa[0] is 'M' and aa[-1] is 'X':
		tag = 'intact'
	else:
		tag = 'partial'
	return tag

#############################################################################
# main body
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################

genomes = glob.glob('%s/*fasta' %options.path)
seqlen = fas2length(options.query)

dictc = {}
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
	#tblastnCMD = 'tblastn -num_threads %s -evalue 1e-10 -db genome -outfmt 6 -query %s -out tblastn.f6' %(options.cpu,options.query)
	os.system(mkdbCMD)
	#os.system(tblastnCMD)
	paraBlast(options.query)
	

	dicta = blast_parser('./tblastn.f6',seqlen)
	out = open('result.cds1','w')
	dictb = {}
	for scaffold in dicta:
		dna = name2seq[scaffold][1]
		dictb[scaffold] = filter1(dicta[scaffold],150)
		i = 0
		j = 0
		for locus in dictb[scaffold]:
			dictb[scaffold][i] = extend(dna,locus) #the first three element of locus were updated, and the last two change to zero
			print (locus),
			print (dictb[scaffold][i])
			i += 1
			j += 1

		list_filtered = filter(dictb[scaffold])
		for locus in list_filtered:
			if locus[0] < locus[1]:
				cds = extractor(dna,locus[0]-1,locus[1])
				cds_tag = tag(cds)
				out.write('>%s|%s|%s|%s\n' %(scaffold,locus[0],locus[1],cds_tag))
				out.write('%s\n' %cds)
			if locus[0] > locus[1]:
				cds = extractor(dna,locus[1]-1,locus[0])
				cds = reverse(complement(cds))
				cds_tag = tag(cds)
				out.write('>%s|%s|%s|%s\n' %(scaffold,locus[0],locus[1],cds_tag))
				out.write('%s\n' %cds)
	out.close()
	os.chdir('../')
		
#########################################################################
# the end 
#########################################################################
