#!/usr/bin/perl
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;

$usage = <<USAGE;
About this script:
    Version : 20180625/001
    Trimal was added in this version
Usage:
    perl paml.pl -i=seq_path -t=tree_file -m=mode
    3 parameters should be given
	parameters mode can be set as follows:
		Model01ratio
        	Model1Neutral
       		Model2Selection
        	Model3Discrtk2
	        Model3Discrtk3
        	Model7beta
        	Model8beta
        	ModelAbranchSite
        	NullModelAbranchSite
        	ModelBbranchSite
        	CladeModelC
        	NullCladeModelC
        	Model2RatioBranch
        	ModelFreeRatio

For example:
    perl $0 -i=/home/zdh/Software/Genome_tools/Paml_pip/fas -t=/home/zdh/Software/Genome_tools/Paml_pip/8sp.tre -m=Model01ratio -m=ModelFreeRatio -m=Model2RatioBranch

USAGE

if (@ARGV < 3){die $usage}

GetOptions('in|i=s' => \$path,
           'tree|t=s' => \$tree,
           'model|m=s' => \@models,
	   'cpu|cpu=s' => \$cpu,
	  );

$cpu ||= 30;
@suffixlist=qw(.txt .fa .fasta .fas .nuc.fas.fa .fas.fa .trimal.fa);

print ">>> running prank...\n";
mkdir "prank", 0755 or warn "cannot make $orthologid directory: $!";
chdir "prank";

##生成比对命令文件
@ortholog = `ls $path`;
open (PRANK_CMD,">prank.sh");
foreach (@ortholog) {
	chomp;
	$fas = $_;
	($name_a, $na, $suffix)=fileparse($_, @suffixlist);
	$prefix = s/\.\w+//;
	#print PRANK_CMD "prank \-d\=$path\/$fas \-o\=$name_a \-translate \-F \-termgap \-f\=fasta \n";
	print PRANK_CMD "prank \-d\=$path\/$fas \-o\=$name_a \-codon \+F \-termgap \-f\=fasta \n";
}
close PRANK_CMD;

`chmod 777 prank.sh`;
unless (-e "prank.ok"){
	system ("ParaFly -c prank.sh -CPU $cpu");
	open OUT, ">", "prank.ok" or die $!; close OUT;
}

$prank_out_path = `pwd`;
chomp $prank_out_path;
opendir PRANK_OUT,$prank_out_path;
@prank_out = (readdir PRANK_OUT);

print ">>> running trimal...\n";
chdir "../";
mkdir "trimal", 0755 or warn "cannot make $orthologid directory: $!";
chdir "trimal";
open (TRIMAL_CMD,">trimal.sh");
foreach $file(@prank_out) {
        next if $file !~ /best\.fas$/ && $file !~ /best\.nuc\.fas$/;
        chomp;
        print TRIMAL_CMD "trimal -in $prank_out_path/$file -out $file\.trimal -automated1\n";
}
close TRIMAL_CMD;
`ParaFly -c trimal.sh -CPU $cpu`;

$trimal_out_path = `pwd`;
chomp $trimal_out_path;
opendir TRIMAL_OUT,$trimal_out_path;
@trimal_out = (readdir TRIMAL_OUT);

open TIM,">new_title.sh";
foreach (@trimal_out){
        next if $_ !~ /trimal$/;
        $cmd = "perl -pi -e 's/\\s\\d+\\sbp//g' $trimal_out_path/$_";
        print TIM "$cmd\n";
}
close TIM;
`ParaFly -c new_title.sh -CPU $cpu`;

print ">>> running Gblocks...\n";
#gblock
chdir "../";
mkdir "gblock", 0755 or warn "cannot make $orthologid directory: $!";
chdir "gblock";
open (GBLOCK_CMD,">gblock.sh");
foreach (@trimal_out) {
	next if $_ !~ /trimal$/;
	chomp;
	print GBLOCK_CMD "Gblocks $trimal_out_path\/$_ \-e\=\.fa \-b5\=N \-t\=c \-b0\=4\n";
	`Gblocks $trimal_out_path\/$_ \-e\=\.fa \-b5\=N \-t\=c \-b0\=4\n`;
}
close GBLOCK_CMD;
#system ("ParaFly -c gblock.sh -CPU $cpu");

print ">>> convert to paml format...\n";
#格式转换
chdir "../";
mkdir "convert_paml", 0755 or warn "cannot make $orthologid directory: $!";
chdir "convert_paml";
open (CONVERT_CMD,">format_paml.sh");
opendir PRANK, $trimal_out_path;
foreach (readdir PRANK) {
	chomp;
	if (/\.trimal\.fa$/) {
		($name_b, $na, $suffix)=fileparse($_, @suffixlist);
		print CONVERT_CMD "prank -d=$trimal_out_path/$_ \-o=$name_b \-f\=paml \-convert\n";
	}
}
close CONVERT_CMD;
closedir PRANK;

system ("ParaFly -c format_paml.sh -CPU $cpu");
$convert_paml_path = `pwd`;
chomp $convert_paml_path;
chdir "../";

print ">>> running paml...\n";
mkdir "run_paml", 0755 or warn "cannot make $orthologid directory: $!";
chdir "run_paml";
$paml_path = `pwd`;
chomp $paml_path;
	
opendir CONVERT_PAML,$convert_paml_path;
$pm = new Parallel::ForkManager($cpu);
foreach $ortholog_phy (readdir CONVERT_PAML) {
	if ($ortholog_phy =~ /phy$/) {
		$count++;
		chomp $ortholog_phy;
		$fasta_file = "$convert_paml_path\/$ortholog_phy";
		mkdir "$ortholog_phy", 0755 or warn "cannot make $orthologid directory: $!";
		chdir "$ortholog_phy";
		run_codeml(@models);
		chdir "../";
	}	
}
$pm->wait_all_children;

##########################################################################
#                             subroutine
##########################################################################
sub run_codeml {
	@models_in = @_; 

	%models = (
	'Model01ratio' 		=> [ 0, 0, 10, 0, 0.4], 
	'Model1Neutral' 	=> [ 0, 1, 10, 0, 0.4],
	'Model2Selection' 	=> [ 0, 2, 10, 0, 0.4],
	'Model3Discrtk2' 	=> [ 0, 3, 2, 0, 0.4],
	'Model3Discrtk3' 	=> [ 0, 3, 3, 0, 0.4],
	'Model7beta' 		=> [ 0, 7, 10, 0, 0.4],
	'Model8beta' 		=> [ 0, 8, 10, 0, 0.4],
	'ModelAbranchSite' 	=> [ 2, 2, 3, 0, 0.4],
	'NullModelAbranchSite' 	=> [ 2, 2, 3, 1, 1],
	'ModelBbranchSite' 	=> [ 2, 3, 3, 0, 0.4],
	'CladeModelC'		=> [ 3, 2, 3, 0, 0.4],
	'NullCladeModelC'	=> [ 3, 2, 3, 1, 1],
	'Model2RatioBranch'	=> [ 2, 0, 10, 0, 0.4],
	'ModelFreeRatio' 	=> [ 1, 0, 10, 0, 0.4],
	);

	foreach $model (@models_in){
		mkdir "$model", 0755 or warn "cannot make $model directory: $!";
                chdir "$model" or die "cannot chdir to $model: $!";
                $current_dir = `pwd`;
                chomp $current_dir;

                open (CTL,">$model.ctl")||die;
  		print CTL "#be carefully!
    seqfile  = $fasta_file
    treefile = $tree
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
       model = $models{$model}[0] * models for codons:
                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
       	            * models for AAs or codon-translated AAs:
               	    * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                    * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
     NSsites = $models{$model}[1] * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
       	            * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
               	    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
       icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
   fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
       kappa = 3  * initial or fixed kappa
   fix_omega = $models{$model}[3]  * 1: omega or omega_1 fixed, 0: estimate
       omega = $models{$model}[4] * initial or fixed omega, for codons or codon-based AAs
   fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
       alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
      Malpha = 0  * different alphas for genes
       ncatG = $models{$model}[2] * # of categories in dG of NSsites models
       clock = 0  * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
       getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .5e-6\n";
			
			close CTL;
			$pm->start and chdir "../" and next;
			system ("codeml $model.ctl");
			$pm->finish;
	}
}

###########################################################################################
