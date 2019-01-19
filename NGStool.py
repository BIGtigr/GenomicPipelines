#!/usr/bin/python

import glob
import sys
import os
import re
from optparse import OptionParser
import argparse
from Bio import SeqIO
import random


parser = OptionParser()

#parser = argparse.ArgumentParser(description="Remote tool to convert fasta to phy for paml")
#parser.add_argument("-i",dest="i", required=True, help="Provide the input file")
#parser.add_argument("-o",dest="o", help="Output path")
#parser.add_argument("-s",dest="stop",default=False, help="weather or not to remove stop codons")#not used

#arg = parser.parse_args()
#pathIn=arg.i
#pathOut=arg.o
#rmStop=arg.stop

#parser.add_option("-v", dest="vcf_file",
#                  help="VCF File to Annotate", metavar="VCF")
#parser.add_option("-o", dest="vcf_output",
                  #help="VCF File to Annotate", metavar="VCF")
#(options, args) = parser.parse_args()

#vcffile = options.vcf_file
prank_dir="/home/zdh/Software/prank/bin"
muscle_dir="/home/zdh/Software/Muscle"
ParaFly_dir = "/home/zdh/Software/parafly-r2013-01-21/bin"
genblastg_dir = "/home/zdh/Software/genBlast_v138_linux_x86_64"
blast_dir = "/home/zdh/Software/ncbi-blast-2.6.0+-src/c++/bin"
trimal_dir = "/home/zdh/Software/trimAl/source"
gblock_dir = "/home/zdh/Software/Gblocks_0.91b"

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

class Fasta:
	def __inity__(self):
		pass
	
	def fas2hash(self,fasfile):
		name2seq = {}
		with open(fasfile) as fh:
			for line in fh.readlines():
				line = line.translate(None,'\n')
				if line.startswith('>'):
					seq_title = line.replace('>','')
					name2seq[seq_title] = ''
				else:
					name2seq[seq_title] += line
			return name2seq
	def fasta2dict(fFile):
		(defLine, seq) = ('', '')
		d = {}
		for l in fFile:
			if '>' == l[0]: # a new def line
				if defLine: # record previous def line and seq
					d[defLine[1:].split()[0]] = (defLine, seq)
				(defLine, seq) = (l, '')
			else: # more of the sequence
            			seq += l
		if defLine: d[defLine[1:].split()[0]] = (defLine, seq)
		return d 
		
	def align(self, fasfile, software = 'prank'):
		print ("Aligning reads...")
		basename = (fasfile.split('.'))[0]
		if software is 'prank':
			cmd = "%s/prank -d=%s -o=%s -codon +F -termgap -f=fasta" %(prank_dir, fasfile, basename)
			os.system(cmd)
		
		if software is 'muscle':
			cmd = "%s/muscle3.8.31_i86linux64 -in %s -out %s.muscle -quiet " %(muscle_dir, fasfile, basename)
			os.system(cmd)

	def trimal(self, fasfile):
		basename = (fasfile.split('.'))[0]
		print ("Run trimal...")
		cmd = "%s/trimal -in %s -out %s.trimal -automated1" %(trimal_dir, fasfile, basename)
		os.system(cmd)
		
	def gblock(self, fasfile):
		basename = (fasfile.split('.'))[0]
		print ("Run gblock...")
		cmd = "%s/Gblocks %s -e=.fa -b5=n -t=c -b0=4" %(trimal_dir, fasfile)
		os.system(cmd)

	def convert_format(self, fasfile,software='prank'):
		basename = (fasfile.split('.'))[0]
                if software is 'prank':
                        cmd = "%s/prank -d=%s -o=%s -f=paml -convert" %(prank_dir, fasfile, basename)
                        os.system(cmd)

	def getSeq(self,fasfile,seqid):
		name2seq = self.fas2hash(fasfile)
		return name2seq[seqid]

	def cds2aa(self,fasfile):
		name2seq = self.fas2hash(fasfile)
		for title in name2seq:
			seq = name2seq[title].upper()
			i = 0
			sequence = ''
			while i+2 < len(seq):
				NNN = seq[i:i+3]
				aa = CODE['standard'][NNN]
				sequence = sequence + aa
				i = i+3
			print ('>%s\n%s' %(title,sequence))
	
	def triplize(self,fasfile):
		name2seq = self.fas2hash(fasfile)
		for title in name2seq:
			aex = {}
			seq = name2seq[title].upper()
			for i in range(3):
				j = i
				position = 0
				aex_num = 0
				sequence = []
				while i+2 < len(seq):
					NNN = seq[i:i+3]
					sequence.append(NNN)
					if NNN in CODE['standard']:
						aa = CODE['standard'][NNN]
						if aa is 'X':
							aex_num = aex_num + 1
							position = i
					i = i + 3
				#if aex_num == 0 or aex_num == 1:
				#if aex_num == 0:
				#if aex_num == 1:
				triple_seq = ''.join(sequence)
				print (">%s|%s|%d|%d|%d\n%s" %(title,aex_num,position,len(seq),j,triple_seq))

	def triplize1(self,fasfile):
		fasfile_core = (fasfile.split('/'))[-1]
		print (fasfile_core)
		name2seq = self.fas2hash(fasfile)
		#fasfile_core = (fasfile.split('\/'))[-1]
		#print (fasfile_core)
		out = open(fasfile_core,'w')
		for title in name2seq:
			seq = name2seq[title].upper()
			for i in range(3):
				new_seq = seq[i:]
				#print (">%s|%s\n%s\n" %(title,i,new_seq))
				out.write(">%s|%s\n%s\n" %(title,i,new_seq))
		out.close()
	
	def extract_seqs(self,fasfile,*titles):
		fasfile_core = (fasfile.split('/'))[-1]
		name2seq = self.fas2hash(fasfile)
		out = open(fasfile_core,'w')
		for i in titles:
			out.write('>%s\n%s\n' %(i,name2seq[i]))
		out.close()

	def filter(self,fasfile,min_length):
		name2seq = self.fas2hash(fasfile)
		outfile = "%s.filtered" %fasfile
		outfh = open(outfile,'w')
		for title in name2seq:
			if len(name2seq[title]) > int(min_length):
				outfh.write('>%s\n%s\n' %(title,name2seq[title]))
		outfh.close()
	
	def codeml(self, phy, tree, model):
		phy_file = phy.split('/')[-1]
		models = {'Model01ratio':[0, 0, 10, 0, 0.4],
			'Model1Neutral':[0, 1, 10, 0, 0.4],
        		'Model2Selection':[ 0, 2, 10, 0, 0.4],
        		'Model3Discrtk2':[0, 3, 2, 0, 0.4],
        		'Model3Discrtk3':[0, 3, 3, 0, 0.4],
        		'Model7beta':[0, 7, 10, 0, 0.4],
        		'Model8beta':[0, 8, 10, 0, 0.4],
        		'ModelAbranchSite':[2, 2, 3, 0, 0.4],
        		'NullModelAbranchSite':[ 2, 2, 3, 1, 1],
        		'ModelBbranchSite':[2, 3, 3, 0, 0.4],
        		'CladeModelC':[3, 2, 3, 0, 0.4],
        		'NullCladeModelC':[3, 2, 3, 1, 1],
        		'Model2RatioBranch':[2, 0, 10, 0, 0.4],
        		'ModelFreeRatio':[1, 0, 10, 0, 0.4],
			}
		os.makedirs('%s/%s' %(phy_file,model))
		os.chdir('%s/%s' %(phy_file,model))
		
		ctl = open('codeml.ctl','w')
		ctl.write('''
    seqfile  = %s
    treefile = %s
    outfile  = mlc.out
       noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1  * 0: concise; 1: detailed, 2: too much
     runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
     seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
      aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
  aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                    * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
       model = %s * models for codons:
                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                    * models for AAs or codon-translated AAs:
                    * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                    * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
     NSsites = %s * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
       icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
   fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
       kappa = 3  * initial or fixed kappa
   fix_omega = %s  * 1: omega or omega_1 fixed, 0: estimate
       omega = %s * initial or fixed omega, for codons or codon-based AAs
   fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
       alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
      Malpha = 0  * different alphas for genes
       ncatG = %s * # of categories in dG of NSsites models
       clock = 0  * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
       getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .5e-6
''' %(phy, tree, models[model][0], models[model][1], models[model][3], models[model][4], models[model][2]))
		ctl.close()
		os.system("codeml codeml.ctl")

	def kaks(self, sp1, sp2):
		pass

	def addSuffix(self, fasfile):
		basename = (fasfile.split('.'))[0]
		cmd = "perl -pi -e 's/>/>%s|/' %s" %(basename,fasfile)
		os.system(cmd)

	def NcbiLongestIso(self,fasta_file):
		length = {}
		dict_a = {}
		raw_title = {}
		with open (fasta_file) as fh:
			for line in fh:
				line.strip()
				if line.startswith('>'):
					lines = re.findall("\(\S*\)",line)
					acession_num = line.split()[0]
					gene = lines[-1]
					raw_title[gene] = line
					if length.has_key(gene):
						length[gene][acession_num] = 0
					else:
						length[gene] = {}
						length[gene][acession_num] = 0
					dict_a[acession_num] = ''
				else:
					dict_a[acession_num] = dict_a[acession_num] + line
					length[gene][acession_num] = length[gene][acession_num] + len(line)
		for g in length:
			#print (length[g])
			max_key = max(length[g], key=length[g].get)
			print ('%s|%s|%s\n%s' %(max_key,g,raw_title[g],dict_a[max_key])),
	def NcbiLongestIso1(self,fasta_file):
		length = {}
		dict_a = {}
		name2seq = self.fas2hash(fasta_file)
		for key in name2seq:
			if key.endswith('mRNA'):
				bracket = re.findall("\(\S*\)",key)
				acession_num = key.split()[0]
				gene = bracket[-1]
				#print (gene)
				if length.has_key(gene):
					length[gene][acession_num] = len(name2seq[key])
				else:
					length[gene] = {}
					length[gene][acession_num] = len(name2seq[key])
				dict_a[acession_num] = name2seq[key]
		for g in length:
			max_key = max(length[g], key=length[g].get)
			print ('>%s|%s\n%s' %(max_key,g,dict_a[max_key]))

	def geneAnotate(self,query,target):
		pass

	def blast(self, query, target, cpu = 30):
		parameters = "-evalue 1e-5 -outfmt 6 -max_hsps 1 -max_target_seqs 5"
		target_basename = (target.split('.'))[0]
		query_basename = (query.split('.'))[0]
		if re.match(r'[^ATGC]',(open(target).read())[1],re.I):
			mkdb_cmd = "makeblastdb -in %s -dbtype prot -out %s" %(target, target_basename)
			os.system(mkdb_cmd)
			if re.match(r'[^ATGC]',(open(query).read())[1],re.I):
				blast_cmd = "blastp -num_threads %d -query %s -db %s %s -out %s_%s.blastp" %(cpu,query,target_basename,parameters,query_basename,target_basename)
				os.system(blast_cmd)
			else:
				blast_cmd = "blastx -num_threads %d -query %s -db %s %s -out %s_%s.blastx" %(cpu,query,target_basename,parameters,query_basename,target_basename)
				os.system(blast_cmd)
		else:
			mkdb_cmd = "makeblastdb -in %s -dbtype nucl -out %s" %(target, target_basename)
			os.system(mkdb_cmd)
			if re.match(r'[^ATGC]',(open(query).read())[1],re.I):
				blast_cmd = "tblastn -num_threads %d -query %s -db %s %s -out %s_%s.tblastn" %(cpu,query,target_basename,parameters,query_basename,target_basename)
				os.system(blast_cmd)
			else:
				blast_cmd = "blastn -num_threads %d -query %s -db %s %s -out %s_%s.blastn" %(cpu,query,target_basename,parameters,query_basename,target_basename)
				os.system(blast_cmd)
	def remoteblastp(self,fasfile,out='remoteblastp.out'):
		from Bio.Blast import NCBIWWW
		from Bio.Blast import NCBIXML
		my_perc_ident = 'none'
		my_blast_program = 'blastp'
		my_evalue_treshold = 0.00001
		my_hitlist_size = 10
		
		filein = open(fasfile, "r")
		fileout = open(out, 'w')
		for record in SeqIO.parse(filein, format="fasta"):
			result = NCBIWWW.qblast(my_blast_program, "nr", record.format("fasta"), hitlist_size = my_hitlist_size, expect = my_evalue_treshold, perc_ident = my_perc_ident)
			fileout.write(result.read())
		filein.close()
		fileout.close()
			
		

	def split(self, fasfile, outdir = 'Split'):
		os.system('mkdir %s' %outdir)
		basename = (fasfile.split('.'))[0]
		for record in SeqIO.parse(fasfile,"fasta"):
			ofile = open (outdir + '/' + record.id + '.fa', 'w')
			SeqIO.write(record,ofile,"fasta")
	def complement(self, seq):
		complement = {'A':'T','C':'G','G':'C','T':'A','a':'t','t':'a','c':'g','g':'c','R':'Y','Y':'R','M':'K','K':'M','r':'y','y':'r','m':'k','k':'m'}
		bases = list(seq)
		
		for i in range(len(bases)):
			bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
		return ''.join(bases)
	def reverse_complement(self, seq):
		return self.complement(seq[::-1])
	
	def randseq(self,n, l, gc):
		"""
    		Prints all the random sequences be-they generated from adhoc or template
    		arguments. The resultant sequences are sent to standard-out.
    		"""
		if l <= 0:
			raise IOError('Positive sequence length (-l) required [error].')
		if n <= 0:
			raise IOError('Positive number of sequence (-n) required [error].')
		if gc > 100 or gc < 0:
			raise IOError('GC percentage (-gc) must be between 0 .. 100 [error].')

		num_seqs, seq_len, gc_perc = n, l, gc / 100.0
		seqs = []
		for _ in range(num_seqs):
			# begin by making an AT repeat-sequence of the user-desired length
			seq_list = list('AT' * seq_len)[:seq_len]
			num_gc_reqd = int(len(seq_list) * gc_perc)  # number of GCs required
			# create list of unique indices
			gc_positions = list(range(0, len(seq_list)))
			random.shuffle(gc_positions)  # jumble their positions and add G or C
			gc_positions = gc_positions[: num_gc_reqd]
			for position in gc_positions:
				g_or_c = random.choice(['G', 'C'])
				seq_list[position] = g_or_c  # insert either a G or C
			seq_str = ''.join(seq_list)
			seqs.append(seq_str)  # save as FASTA
		for i, seq in enumerate(seqs):
			# shuffle bases so that if a sequence with 0 GC% are not only AT dimers
			seq = list(seq)
			random.shuffle(seq)
			print('>sequence_' + str(i + 1) + '\n' + ''.join(seq))
	def get_n_perc(self,seq):
		"""
		Derive the N % of a sequence.
		@param seq: DNA sequence.
		@return: percentage representing the N percentage.
		"""
		n_count = float(str(seq).upper().count('N'))
		return n_count / len(seq) * 100

	def slidingwindow(self,f, w, o, n):
		"""
		Runs the logic of the sliding window algorithm with overlapping segments.
		@param f: User-provided input FASTA file.
		@param w: Sliding window size.
		@param o: Sliding window overlap size.
		@param n: Hard-mask threshold; ignore sequences exceeding this value.
		"""
		entries = SeqIO.parse(f, 'fasta')
		
		for entry in entries:
			seq = str(entry.seq)
			d = entry.description  # sequence descriptor
			chunk1 = seq[0: w]  # the first chunk has no overlaps
			start, end = 0, w
			if self.get_n_perc(chunk1) < n:
				print('>' + d + '|' + str(start) + '|w|' + str(w) + '|o|' + str(o))
				print(chunk1)
			while True:
				start = end - o
				end = start + w
				win = seq[start: end]
				if start > len(seq):
					break
				if start != len(seq) and self.get_n_perc(win) < n:
					print('>' + d + '|' + str(start) + '|w|' + str(w) + '|o|' + str(o))
				print(win)

	def swissprot(self,fas):
		blast_cmd = "blastp -num_threads 40 -query %s -db %s -outfmt 5 -evalue 1e-5 -out swissprot.out" %(fas,swissdb)
		os.system(blast_cmd)
	
class MCMCtree(Fasta):
	def four_fold_site(self,paml_path,sp_num):
		four_fold_nuc = glob.glob('%s/*/*/4fold.nuc' %paml_path)
		four_fold = {}
		for nuc_file in four_fold_nuc:
			seq = {}
			with open(nuc_file) as nuc:
				for line in nuc:
					line = line.strip()
					line = line.replace("\s*$","")
					if re.match('[A-Z]',line):
						if re.findall(r'[^ATGCN]',line):
							seq_title = line
						else:
							if seq_title in seq:
								seq[seq_title] = seq[seq_title] + line
							else:
								seq[seq_title] = line
				if len(seq) == int(sp_num):
					for seq_title in seq:
						if seq_title in four_fold:
							four_fold[seq_title] = four_fold[seq_title] + seq[seq_title]
						else:
							four_fold[seq_title] = seq[seq_title]

		out_4fold = open('4fold_concatenated.fas','w')
		for i in four_fold:
			out_4fold.write('>%s\n%s\n' %(i,four_fold[i]))
		out_4fold.close()

		phycmd='clustalw -INFILE=4fold_concatenated.fas -CONVERT -TYPE=DNA -OUTPUT=PHYLIP -OUTFILE=4fold_concatenated.phylip'
		ModelTest='java -jar /home/wangkai/home/tools/jmodeltest2/jModelTest-2.1.4/jModelTest.jar -d 4fold_concatenated.phylip -o 4fold_modeltest.out -s 11 -f -i -g 4 -w -AIC -AICc -BIC'
		#fasta = Fasta()
		#fasta.convert_format('4fold_concatenated.fas')
		self.convert_format('4fold_concatenated.fas')
		os.system(phycmd)
		os.system(ModelTest)
	def four_fold_site1(self,path):
		four_fold = {}
		four_fold_nuc = glob.glob('4fold.nuc','%s/Model01ratio/' %path)
			

	def mcmctree(self,phy = '4fold_concatenated.phy',tree = 'mcmctree.input',out = 'mcmctree.out'):
		ctl = open('mcmctree.ctl','w')
                ctl.write('''
	seed = -1
     seqfile = %s
    treefile = %s
     outfile = %s
       ndata = 1
     seqtype = 0 * 0: nucleotodes; 1: codons; 2: AAs
     usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
       clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
     RootAge = < 3.304  * safe constraint on root age, used if no fossil for root.

       model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
		    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
 
       alpha = 3.6060    * alpha for gamma rates at sites
       ncatG = 5    * No. categories in discrete gamma

   cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

     BDparas = 1 1 0    * birth, death, sampling
 kappa_gamma = 6 2      * gamma prior for kappa
 alpha_gamma = 1 1      * gamma prior for alpha

 rgene_gamma = 2 2     * gamma prior for overall rates for genes
sigma2_gamma = 1 10    * gamma prior for sigma^2     (for clock=2 or 3)

    finetune = 0: 0.09 0.11 0.28 0.19 0.3  * times, rates, mixing, paras, RateParas

       print = 1
      burnin = 2000
    sampfreq = 5
     nsample = 10000''' %(phy,tree,out))
		ctl.close()
		os.system('mcmctree mcmctree.ctl')
		
			
class Bam:
	pass

class Fastq:
	pass
class VCF:
	def ReadsMapping(self,genome,read1,read2):#fit for one samplle back_mapping to genome
		index_cmd = 'bwa index -a bwtsw %s -p genome' %genome
		map_cmd = 'bwa mem -t 40 genome %s %s | samtools view -b - | samtools rmdup -S - rmdup.bam' %(read1,read2)
		sort_cmd = 'samtools sort -o rmdup_sort.bam rmdup.bam'
		list_cmd = 'ls *sort.bam > sort_bam.list'
		mpileup_cmd = 'samtools mpileup -q 20 -Q 20 -C 50 -t DP,SP -m 2 -F 0.002 -u -g -f %s -b sort_bam.list | bcftools call -vcO v -o raw.vcf -' %genome
	
		#os.system(index_cmd) #genome index
		#os.system(map_cmd) #read mapping, remove PCR duplicates
		os.system(sort_cmd) #sort bam file
		os.system(list_cmd) #creat bam list for mpileup
		os.system(mpileup_cmd) #parsing reads aligment,snp/indel calling
	def MappingStat(self,genome,vcffile):
		pass		
		
		
class Gff:
	def statistics(self,gff):
		gene_num = 0.0
		exon_num = 0.0
		total_exon_length = 0
		intron_num = 0.0
		intron_length = 0
		total_intron_length = 0
		total_transcript_length = 0
		with open(gff) as gf:
			for line in gf:
				if line.startswith('//') :
					continue
				lines = line.split('\t')
				if len(lines) < 5:
					continue
				if lines[2] == "match" or lines[2] == "gene":
					gene_num = gene_num + 1
					transcript_length = abs(int(lines[3]) - int(lines[4])) + 1
					total_transcript_length = total_transcript_length + transcript_length
				if lines[2] == "cds" or lines[2] == "CDS":
					exon_num = exon_num + 1
					exon_length = abs(int(lines[3]) - int(lines[4])) + 1
					total_exon_length = total_exon_length + exon_length
				if lines[2] == "intron":
					intron_num = intron_num + 1
					intron_length = abs(int(lines[3]) - int(lines[4])) + 1
					total_intron_length = total_intron_length + intron_length
		if total_intron_length == 0:
			total_intron_length = total_transcript_length - total_exon_length
			intron_num = exon_num - gene_num		

		total_intron_length = total_intron_length + intron_length
		average_transcript_length = total_transcript_length / gene_num
		average_cds_length = total_exon_length / gene_num
		average_exon_length = total_exon_length / exon_num
		average_intron_length = total_intron_length / intron_num
		average_exon_per_gene = exon_num / gene_num
		
		print("gene_num:%s" %gene_num)
		print("average_transcript_length:%s" %average_transcript_length)
		print("average_cds_length:%s" %average_cds_length)
		print("average_exon_per_gene:%s" %average_exon_per_gene)
		print("average_exon_length:%s" %average_exon_length)
		print("average_intron_length:%s" %average_intron_length)


class Statistics:
	pass

class Draw:
	pass

class Machine_learn:
	pass
