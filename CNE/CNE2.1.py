import sys,os,glob
from optparse import OptionParser


parser = OptionParser()
parser.add_option('--gm', help=('Absolute path to genomes'))
parser.add_option('--ref', help=('reference genome'))
parser.add_option('--gff', help=('gff file for ref'))
#parser.add_option('--cpu', type='int', default=30,help='numbers of threads')

usage = u'''`Usage:python %prog [options...] [genomes] [ref]'''
(options, args) = parser.parse_args()

numOpts = len(sys.argv)
if numOpts < 3:
    parser.print_help()
    sys.exit()

genomes = os.listdir(options.gm)
ref = options.ref
GFF = options.gff

g, s = ref.split('_')
ABBR_ref = g[0:3].upper() + s[0:2].upper()

#here are something you need to specify
##tree
tree = '(((((((GYPHI,(HALLE,BUTHE)),AQUCH),GYPBA),CATAU),(PICPU,(TAEGU,PSEHU))),NIPNI),GALGA)'
tree1 = tree.replace(',',' ')

## species in different group
TARGETSPECIES = 'GYPHI;GYPBA;CATAU'
OUTGROUP = 'BUTHE;HALLE;AQUCH'
CONSERVE = 'GALGA;PICPU;TAEGU;PSEHU;NIPNI'

## path to several important tools
bin_path = '/home/zdh/Projects/Scavenge/68.CNE/01.genome_alignment/bin'
multi_bin = '/home/zdh/Projects/Scavenge/68.CNE/01.genome_alignment/02.CNE/01.multiple/bin'
msaview_bin = '/home/zdh/Tools/05.48birds/Highly_conserved_elements.bak/bin'


################################################################################
def maf_read(mafs):
	startONscaf = {}
        dict1 = {}
        dict2 = {}
        num = 0
        pwd = os.getcwd()
        if not os.path.exists('split_maf'):
                os.mkdir('split_maf')
        for maf in mafs:
                scaf = maf.split('/')[-2]
                dict1[scaf] = {}
                dict2[scaf] = {}
                startONscaf[scaf] = []
                with open(maf) as MAF:
                        for line in MAF:
                                line = line.strip()
                                if line.startswith('a score'):
                                        num += 1
                                        maf_tem = open('split_maf/%s.maf' %num,'w')
                                        maf_tem.write('%s\n' %line)
					flag = 0
                                elif line.startswith('s '):
                                        maf_tem.write('%s\n' %line)
                                        if flag == 0 :
						flag = 1
                                                lines = line.split()
                                                dict1[scaf][lines[2]] = '%s/split_maf/%s.maf' %(pwd,num)
                                                #maf_tem.write('%s\n' %line)

                                                length = len(lines[6])
                                                start = int(lines[2]) - 1
                                                startONscaf[scaf].append(lines[2])
                                                for i in range(length):
                                                        pos = i + 1
                                                        if lines[6][i] != '-':
                                                                start = start + 1
                                                                dict2[scaf][start] = pos
        return dict1,dict2,startONscaf
#######################################################################################
def maf_read0(maf):
	dict_ = {}
	species = []
	num = 0
	with open(maf) as MAF:
		for line in MAF:
			line = line.strip()
				if line.startswith('a score'):
					dict_[num] = {}
					flag = 0
					num += 1
				elif line.startswith('s'):
					lines = line.split()
					if flag == 0:
						dict_[num]['ref'] = line
					dict_[num][lines[1]] = lines[-1]
					if lines[1] not in species:
						species.append(lines[1])
	return dict_,species
					

#######################################################################################
def maf_read2(mafs):
	startONscaf = {}
	dict1 = {}
	dict2 = {}
	num = 0
	pwd = os.getcwd()
	species = {}
	if not os.path.exists('split_maf'):
		os.mkdir('split_maf')
	for maf in mafs:
		temp,species = maf_read0(maf)
		for i in temp:
			seq_len = len(temp[i]['ref'].split()[-1])
			for sp in species:
				if sp not in temp[i]:
					temp[i][sp] = '-'*seq_len
		current = {}
		for i in temp:
			info = temp[i]['ref']
			start  = info[2]
			length = info[3]
			strand = info[4]
			scaf_len = info[5]

			if not current:
				for sp in species:
					current['start'] = start

		
#############################################################################################
def bed_read(beds):
        dict_ = {}
        num = 0
        for bed in beds:
                scaf_ = bed.split('/')[-1]
                scaf = scaf_.split('.')[0]
                if scaf not in dict_:
                        dict_[scaf] = {}
                with open(bed) as BED:
                        for line in BED:
                                num += 1
                                dict_[scaf][num] = []
                                lines = line.split()
                                dict_[scaf][num].append(int(lines[1]))
                                dict_[scaf][num].append(int(lines[2])-1)
        return dict_

##############################################################################################
def FasRead(faspath):
        dict_ = {}
	list_ = []
        CNEfass = os.listdir(faspath)
        for fas in CNEfass:
                dict_[fas] = {}
                pfas = os.path.join(faspath,fas)
                with open(pfas) as FAS:
                        for line in FAS:
                                line = line.strip()
                                if line.startswith('>'):
                                        title = line.split()[-1]
                                        if title not in list_:
                                                list_.append(title)
                                        dict_[fas][title] = ''
                                else:
                                        dict_[fas][title] += line
	return dict_,list_

##############################################################################################
print('Starting pairwise alignment......')
#pairwise genome alignment
if not os.path.exists('00.pairwise'):
	os.mkdir('00.pairwise')
os.chdir('00.pairwise')

pairwise_path = os.getcwd()

ref_path = os.path.join(sys.argv[1],ref)
for genome in genomes:
	if not genome.endswith('fasta'):
		continue
	if genome != ref:
		if not os.path.exists(genome):
			os.mkdir(genome)
		os.chdir(genome)
		genome_file = os.path.join(sys.argv[1],genome)
		lastz_cmd = 'perl %s/lastz_CNM.pl --hspthresh 2200 --inner 2000 --ydrop 3400 --gappedthresh 10000 --scores HoxD55 --chain --linearGap loose --cuts 50 --cpu 50 %s %s' %(bin_path,ref_path,genome_file)
		#os.system(lastz_cmd)
		os.chdir('../')
os.chdir('../')

if not os.path.exists('01.multiple'):
	os.mkdir('01.multiple')
os.chdir('01.multiple')

# multi genome alignment
mafs = glob.glob('%s/*.fasta/7.maf/all.maf' %pairwise_path)
dicta = {}
for maf in mafs:
	genus_species = maf.split('/')[-3]
	genus, species = genus_species.split('_')
	abbr = genus[0:3] + species[0:2]
	ABBR = abbr.upper()
	with open(maf) as MAF:
		for line in MAF:
			if line.startswith('##'):
				info = line
			if line.startswith('a score'):
				a = line
			if line.startswith('s target'):
				line = line.lstrip('s target')
				scaf = line.split()[0]
				line_ = 's %s.'%ABBR_ref + line
				a += line_
				if scaf not in dicta:
					dicta[scaf] = {}
			if line.startswith('s query'):
				line = line.strip('s query')
				line_ = 's %s.%s' %(ABBR,line)
				a += line_
				if ABBR not in dicta[scaf]:
					dicta[scaf][ABBR] = []
				dicta[scaf][ABBR].append(a)
#print alignment by scaffold
for scaf in dicta:
	if not os.path.exists('%s.maf' %scaf):
        	os.mkdir('%s.maf' %scaf)
        LIST = open('%s.list' %scaf, 'w')
        os.chdir('%s.maf' %scaf)
        pwd = os.getcwd()
        for sp in dicta[scaf]:
                LIST.write('%s\t%s\t%s/%s.%s.%s.maf\n' %(ABBR_ref,sp,pwd,scaf,ABBR_ref,sp))
                out = open('%s.%s.%s.maf' %(scaf,ABBR_ref,sp), 'w')
                out.write('%s' %info)
                for align in dicta[scaf][sp]:
                        out.write('%s\n' %align)
                out.close()
        LIST.close()
        os.chdir('../')
# multiple alignment
print('Starting multiple alignment......')
multi_path = os.getcwd()
for list_ in os.listdir(multi_path):
	if list_.endswith('list'):
		scaf = list_.split('.')[0]
		if not os.path.exists(scaf):
			os.mkdir(scaf)
		generate_run_multiz = 'python %s/run_multiz.py --pair_align %s --tree "%s" --out %s' %(multi_bin,list_,tree1,scaf)
		os.system(generate_run_multiz)
		os.chdir(scaf)
		#os.system('sh run_multiz.sh')
		splitCMD = 'msa_split %s_ref.maf --in-format MAF  --windows 1000000,0 --out-root %s.msa_split --out-format SS --min-informative 1000 --between-blocks 5000' %(ABBR_ref,scaf)
		os.system(splitCMD)
		##obtain the 4-fold generate site and estimate the nonconserve model for each chromosome
		FD_sitesCMD = '/home/zdh/Software/phast/bin0/msa_view %s_ref.maf --in-format MAF --4d --features %s > %s.codon.ss' %(ABBR_ref,GFF,scaf)
		FD2faCMD = 'msa_view %s.codon.ss  --in-format SS --out-format FASTA  --tuple-size 1 > %s.sites.fa' %(scaf,scaf)
		phyloFitCMD = 'phyloFit --tree "%s" --msa-format FASTA --out-root %s.nonconserved-4d %s.sites.fa' %(tree,scaf,scaf)
		os.system(FD_sitesCMD)
		os.system(FD2faCMD)
		os.system(phyloFitCMD)
		os.chdir('../')
os.system('ls ./*/*nonconserved-4d.mod  > all_mod.lista')
mods = ''
with open('all_mod.lista') as LIST:
	for line in LIST:
		mods += (line[0:-1]+',')
	mods = mods.rstrip(',')
nonconsCMD = 'phyloBoot --read-mods %s --output-average noncons.mod' %mods
os.system(nonconsCMD)
split_mafs = glob.glob('./*/*msa_split*.ss')
ConsCMD = open('cons.sh','w')
for split_maf in split_mafs:
	conscmd = 'phastCons --target-coverage 0.25 --expected-length 20 --estimate-rho %s --no-post-probs %s noncons.mod' %(split_maf,split_maf)
	ConsCMD.write('%s\n' %conscmd)
ConsCMD.close()
os.system('ParaFly -c cons.sh -CPU 30')

os.system('ls ./*/*.ss.cons.mod > cons.txt')
mods = ''
with open('cons.txt') as LIST:
	for line in LIST:
		mods += (line[0:-1]+',')
	mods = mods.rstrip(',')
consBoot = 'phyloBoot --read-mods %s --output-average cons.mod' %mods
os.system(consBoot)
if not os.path.exists('ELEMENTS'):
	os.mkdir('ELEMENTS')
if not os.path.exists('SCORES'):
	os.mkdir('SCORES')
for i in glob.glob('./*/*msa_split*ss'):
	i_core = i.split('/')[-1]
	phastConsCMD = 'phastCons --most-conserved ELEMENTS/%s.bed --score %s cons.mod,noncons.mod > SCORES/%s.wig' %(i_core,i,i_core)
	os.system(phastConsCMD)

#exclude coding regions from bed file
if not os.path.exists('Elements_1'):
	os.mkdir('Elements_1')
beds = os.listdir('./ELEMENTS')
for bed in beds:
	if not bed.endswith('bed'):
		continue
	bedCMD = 'bedtools subtract -a ./ELEMENTS/%s -b %s > Elements_1/%s' %(bed,GFF,bed)
	os.system(bedCMD)


# Go to ../ to extract CNE sequence
mafs = glob.glob('./*/*ref.maf')
beds = glob.glob('./Elements_1/*bed')

mafpath,site2pos,startONscaf = maf_read(mafs)
bedinfo = bed_read(beds)

if not os.path.exists('CNE_raw'):
	os.mkdir('CNE_raw')
for scaf in bedinfo:
	starts = startONscaf[scaf]
	for num in bedinfo[scaf]:
		start_ = bedinfo[scaf][num][0]
		end_ = bedinfo[scaf][num][1]
		start = site2pos[scaf][start_]
		end = site2pos[scaf][end_]
		if int(start_) > int(starts[-1]):
			maf_path = mafpath[scaf][starts[-1]]
		else:
			for i in range(1,len(starts)):
				if int(start_) >= int(starts[i-1]) and int(start_) < int(starts[i]):
					maf_path = mafpath[scaf][starts[i-1]]
		extractCNECMD = 'msa_view %s --in-format MAF --start %s --end %s --refidx 0 > CNE_raw/%s_%s_%s.fas' %(maf_path,start,end,scaf,start_,end_)
		os.system(extractCNECMD)

os.chdir('../')

# this is the part to prepare for PhyloAcc input files: .bed/.fas/ID.txt
os.mkdir('02.PhyloAcc')
os.chdir('02.PhyloAcc')

para_file = '''PHYTREE_FILE ./input.mod
SEG_FILE ./input.bed
ALIGN_FILE ./input.fas
ID_FILE ./ID.txt
RESULT_FOLDER ./result_tmp/
PREFIX result
BURNIN 500
MCMC 1000
CHAIN 1
TARGETSPECIES %s
OUTGROUP %s
CONSERVE %s
GAPCHAR -
NUM_THREAD 4
VERBOSE 1
'''%(TARGETSPECIES,OUTGROUP,CONSERVE)

para = open('param.txt','w')
para.write(para_file)
para.close()

CNEpath = '../01.multiple/CNE_raw'
CNEseq,Species = FasRead(CNEpath)

num = 0
bed = []
start = 0
end = 0
confas = {}
for fas in CNEseq:
	if ABBR_ref not in CNEseq[fas]:
		continue
	for sp in Species:
		if sp == ABBR_ref:
			lenth = len(CNEseq[fas][ABBR_ref])
			end += lenth
			bed.append('%s\t%s\t%s\t%s\t%s' %(num,start,end,fas,lenth))
			start = end
			num += 1

		if sp in CNEseq[fas]:
			if sp not in confas:
				confas[sp] = ''
			confas[sp] += CNEseq[fas][sp]
		else:
			if sp not in confas:
				confas[sp] = ''
			confas[sp] += '-'*lenth

F = open('input.fas','w')
for sp in confas:
	F.write('>%s\n%s\n' %(sp,confas[sp]))
F.close()

F = open('input.bed','w')
ID = open('ID.txt','w')
for b in bed:
	id_ = b.split()[0]
	F.write('%s\n' %b)
	ID.write('%s\n' %id_)
F.close()
ID.close()

os.system('tree_doctor --name-ancestors ../01.multiple/noncons.mod > input.mod')

os.mkdir('result_tmp')
PhyloAcc = 'PhyloAcc param.txt > PhyloAcc.log'
os.system(PhyloAcc)


############################# The end #######################################

