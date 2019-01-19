#!/usr/bin/perl
use Getopt::Long;
use FileHandle;
use Storable qw(store retrieve);
use List::MoreUtils qw(mesh);
use List::Util qw(max min);
use utf8;


#input: a, masked genome;b, cds_pep sequencs of query species
#ouput: gene model -> homology.gff
#method: #1. genblastg
	 #2. extract candidate gene regions from genome
	 #3. blastx toward to respective pep_cds
	 #4. find best hit of each candidate region
	 #5. genewise between best hit cds_pep and corresponding candidate gene region
	 #6. integrate genewise results of different kind of species
	 #7. remove defective gene models
	 #8. ouptut final gene models

$usage = <<USAGE;
About this script:
    Version : 20181108/001
    Author  : bigtigerzou@163.com
Usage:
    perl <$0> -g=genome -p=cds_path -t=10
    3 parameters should be given
For example:
    perl $0 -g=/home/zdh/Software/Genome_tools/HBGP.1.02/genome/vulture.genome -p=/home/zdh/Software/Genome_tools/HBGP.1.02/test_pep -t=10 evalue=1e-10

USAGE

if (@ARGV != 3){die $usage}

GetOptions('genome|g=s' => \$genome,
           'path|p=s' => \$path,
	'threads|t=s' => \$threads,
	'evalue|e=f' => \$evalue
          );
$threads ||= 30;
$evalue ||= 1e-10;

print ">>> reading genome into hash...\n";
open (FAS, "$genome")||die;
while (<FAS>){
        chomp;
	if (/^>(\S+)/) {$title = $1;}
	else           {$sequence{$title} .= $_; $length{$title} = $length += length $_;}
}
close FAS;

store \%sequence,"genome_seq";
#%sequence = {};

#####################################################
print ">>> running genblastg...\n";

opendir PATH,$path;
foreach $fasta_file (readdir PATH){
	$i = 0;
	$cds_title = '';
	$cds_seq = '';
	@commonds = ();
        chomp $fasta_file;
        if($fasta_file =~ /fasta$/ or $fasta_file =~ /fas$/ or $fasta_file =~ /fa$/){
                mkdir "$fasta_file",0755 or warn "cannot make $orthologid directory: $!";
                `cp alignscore.txt blastall $fasta_file`;
                chdir "$fasta_file";
                `ln -s $genome genome.fasta`;
		`formatdb -i genome.fasta -p F`;
		`makeblastdb -in $path/$fasta_file -dbtype prot -out protein`;
                open CMD,">$fasta_file.sh";

		open FASTA,"$path/$fasta_file";
		while ($line = <FASTA>){
			chomp $line;
			if ($line =~ s/^>//){
				$cds_title = $line;
			}
			else {
				$pep{$cds_title} .= $line;
				$pep1{$cds_title} .= $line;
			}
		}
		foreach $pep_name (keys %pep1){
			$i++;
			open LINE,">temp_$pep_name\_$i";
			print CMD "genblast_v138_linux_x86_64 -p genblastg -q temp_$pep_name\_$i -t genome.fasta -e $evalue -g T -gff -o temp_$pep_name\_$i\n";
			if (length($pep{$pep_name}) > 10000){
				$len_of_pep = length($pep{$pep_name});
				$half_len = $len_of_pep/2;
				$seq1 = substr($pep{$pep_name},0,$half_len+500);
				$seq2 = substr($pep{$pep_name},$half_len-500,$len_of_pep+1);
				print LINE ">$pep_name.a\n";
				print LINE "$seq1\n";
				print LINE ">$pep_name.b\n";
				print LINE "$seq2\n";
				close LINE;
				
			}
			else{
				print LINE ">$pep_name\n";
				print LINE "$pep{$pep_name}\n";
				close LINE;
			}
		}
		#store \%pep,"cds_pep";
		%pep1 = ();
                close CMD;
                close FASTA;
	
		system("nohup ParaFly -c $fasta_file.sh -CPU $threads");
                chdir "../";
        }
}
closedir PATH;


##################################################################################################
print ">>> analyzing results of genblastg...\n";
$geneblastg_out_path = `pwd`;
chomp $geneblastg_out_path;
opendir PATH,$geneblastg_out_path;

foreach $file (readdir PATH) {
	$flag = 0;
	chomp $file;
	next if $file =~ /^\.+/ or -f $file;
	next if $file !~ /fasta$/ && $file !~ /fas$/ && $file !~ /fa$/;
        chdir "$geneblastg_out_path/$file";
	print "-----------------$file----------------------\n";
        $pwd = `pwd`;
        chomp $pwd;
	opendir GFF_DIR,$pwd;
	@gffs = readdir GFF_DIR;
	open TOTAL_GFF,">$file.gff";
	foreach $gff (@gffs){
		if ($gff =~ /gff$/){
			open GFF,"$pwd/$gff";
			local $/ = undef;
			while (<GFF>){
				print TOTAL_GFF "$_";
			}
			close GFF;
		}
	}
	$/ = "\n";
	closedir GFF_DIR;
	close TOTAL_GFF;
	foreach $gff (@gffs){
		`rm $gff` if $gff=~ /genome/ or $gff =~ /temp/
	}
	print ">>> dealing with integrated gff...\n";
        open GFF,"$pwd/$file.gff";
	open CANDIDATE_SEQ,"> $file.fas";
        while (<GFF>) {
                next if /^\#/;
                #GFF TO TAB;
                if ($_ =~ /transcript/) {
                        $_=~s/\s+$//;
                        @line=split(/\t/,$_);
                        $sca=$line[0];
                        $location = "$line[3]\t$line[4]\t$line[6]";
                        $sca{$sca} = "love";
                        push @{$sca},$location;
                }
        }
	close GFF;

        @key = keys %sca;
	%sca = {};
	@scaffords = sort {$a <=> $b} @key;
	foreach $scafford (@scaffords) {
                $flag++;
                $scafford_length = $length{$scafford};
                $scafford_seq = $sequence{$scafford};
                @location = ();
                @raw_location = ();
                @location_merged = ();
                foreach $location (@{$scafford}) {
                        push @raw_location,$location;
                }
		foreach $location (@raw_location){
			($a,$b,$d) = split /\t/,$location;
			$start = $a - 1000;
			$end = $b + 1000;
			$start = 1 if $start <= 0;
			$end = $scafford_length if $end >= $scafford_length;
			$gene_location = $start."\t".$end;
			push @location,$gene_location;
		}

                @location_merged = merge_group(@location);
                foreach (@location_merged){
                        ($start,$end) = split /\t/,$_;
                        $sub_seq = substr ($scafford_seq, ($start - 1), ($end - $start + 1));
                        $gene_title = "$scafford\_$start\_$end";
                        print CANDIDATE_SEQ ">$gene_title\n$sub_seq\n";
                        $gene_mother{$gene_title} = $sub_seq;
                }
		



        }
	close CANDIDATE_SEQ;
	print ">>> running blastx...\n";

	`blastx -num_threads $threads -db $geneblastg_out_path/$file/protein -evalue $evalue -outfmt 6 -max_hsps 1 -max_target_seqs 4 -query $pwd/$file.fas -out $file\_blastx.out`;
	
	print ">>> resolving blastx out...\n";
	open BLASTX_OUT,"$pwd/$file\_blastx.out" or die;
	while ($blastx = <BLASTX_OUT>){
		@blastx = split /\t/,$blastx;
		if ($blastx[0] ne $temp_q){
			if ($blastx[6] < $blastx[7]){
				$best_pos_pro{$blastx[0]} = $blastx[1];
			}
			if ($blastx[6] > $blastx[7]){
				$best_neg_pro{$blastx[0]} = $blastx[1];

			}
			$temp_q = $blastx[0];
			$temp_t = $blastx[1];
		}
#		if ($blastx[0] eq $temp_q && $blastx[1] ne $temp_t){
#			print "---\n";
#			if (!exists $best_neg_pro{$blastx[0]} && $blastx[6] > $blastx[7]){
#				$best_neg_pro{$blastx[0]} = $blastx[1];
#			}
#			if (exists $best_pos_pro{$blastx[0]} && !exists $sub_pos_pro{$blastx[0]} && $blastx[6] < $blastx[7]){
#                               $sub_pos_pro{$blastx[0]} = $blastx[1];
#                               next;
#                       }
#			if (!exists $best_pos_pro{$blastx[0]}&& $blastx[6] < $blastx[7]){
#				$best_pos_pro{$blastx[0]} = $blastx[1];
#			}
#			if (exists $best_neg_pro{$blastx[0]} && !exists $sub_neg_pro{$blastx[0]} && $blastx[6] > $blastx[7]){
#                              $sub_neg_pro{$blastx[0]} = $blastx[1];
#				next;
#			}
#		}		
	}
	
	print ">>> running genewise...\n";
	mkdir "genewise";
	chdir "genewise";
	$genewise_path = `pwd`;
	chomp $genewise_path;
	open GENEWISE_CMD,">genewise.sh";
	foreach $nuc(keys %best_neg_pro){
		mkdir "$nuc\_n";
		chdir "$nuc\_n";
		$prot = $best_neg_pro{$nuc};
		open NUC,">$nuc.fas";
		print NUC ">$nuc\n$gene_mother{$nuc}";
		close NUC;
		open PROT,">$prot.fas";
		print PROT ">$prot\n$pep{$prot}";
		close PROT;
		$pwd = `pwd`;
		chomp $pwd;
		print GENEWISE_CMD "genewise $pwd/$prot.fas $pwd/$nuc.fas -gff -trev -pseudo -nosplice_gtag -silent > $pwd/genewise.out\n";
		chdir "../";
	}
	foreach $nuc (keys %best_pos_pro){
		mkdir "$nuc\_p",0755 or warn "cannot make $orthologid directory: $!";
		chdir "$nuc\_p";
		$prot = $best_pos_pro{$nuc};
		open NUCL,">$nuc.fas";
		print NUCL ">$nuc\n$gene_mother{$nuc}";
		close NUCL;
		open PROT,">$prot.fas";
		print PROT ">$prot\n$pep{$prot}";
		close PROT;
		$pwd = `pwd`;
		chomp $pwd;
		print GENEWISE_CMD "genewise $pwd/$prot.fas $pwd/$nuc.fas -gff -tfor -pseudo -nosplice_gtag -silent > $pwd/genewise.out\n";
		chdir "../";
	}
	foreach $nuc (keys %sub_pos_pro){
                mkdir "$nuc\_sp",0755 or warn "cannot make $orthologid directory: $!";
                chdir "$nuc\_sp";
                $prot = $sub_pos_pro{$nuc};
                open NUCL,">$nuc.fas";
                print NUCL ">$nuc\n$gene_mother{$nuc}";
                close NUCL;
                open PROT,">$prot.fas";
                print PROT ">$prot\n$pep{$prot}";
                close PROT;
                $pwd = `pwd`;
                chomp $pwd;
                print GENEWISE_CMD "genewise $pwd/$prot.fas $pwd/$nuc.fas -gff -tfor -pseudo -nosplice_gtag -silent > $pwd/genewise.out\n";
                chdir "../";
       }
	foreach $nuc (keys %sub_neg_pro){
                mkdir "$nuc\_sn",0755 or warn "cannot make $orthologid directory: $!";
                chdir "$nuc\_sn";
                $prot = $sub_neg_pro{$nuc};
                open NUCL,">$nuc.fas";
                print NUCL ">$nuc\n$gene_mother{$nuc}";
                close NUCL;
                open PROT,">$prot.fas";
                print PROT ">$prot\n$pep{$prot}";
                close PROT;#
                $pwd = `pwd`;
                chomp $pwd;
                print GENEWISE_CMD "genewise $pwd/$prot.fas $pwd/$nuc.fas -gff -trev -pseudo -nosplice_gtag -silent > $pwd/genewise.out\n";
                chdir "../";
        }
	close GENEWISE_CMD;
	system("nohup ParaFly -c genewise.sh -CPU $threads");
	%best_neg_pro = ();
	%best_pos_pro = ();
	%sub_pos_pro = ();
	%sub_neg_pro = ();
	#join genewise result of different segment of genome
	opendir GENEWISE,$genewise_path;
	open WISE_OUT,">sub_genewise.gff3";
	foreach $file (readdir GENEWISE){
		open WISE_temp,"$file/genewise.out";
		local $/ = undef;
		while (<WISE_temp>){
			print WISE_OUT "$_";
		}
		close WISE_temp;		

	}
	close WISE_OUT;
	closedir GENEWISE;
	$/ = "\n";
	chdir "../";
}
chdir "../";
$pwd = `pwd`;
chomp $pwd;
#join genewise result of diferent kinds of species
`cat ./*/genewise/sub_genewise.gff3 > total_genewise.gff3`;

open INPUT,"$pwd/total_genewise.gff3";
while ($line = <INPUT>){
        next if $line =~ /\/\//;
        if ($line =~ /match/){
                $wise{$tag} = $wise_out;
                $wise_out = '';
                @line = split /\t/,$line;
                $contig = (split /\_/,$line[0])[0];
                push @raw_contigs,$contig;
                $score = $line[5];
                if ($line[6] =~ /\+/){
                        push @{$contig.'_forward'},$line[0];
                        push @{$line[0].'_forward'},$score;
                }
                if ($line[6] =~ /\-/){
                        push @{$contig.'_reverse'},$line[0];
                        push @{$line[0].'_reverse'},$score;
                }
                $tag = "$line[0]=$line[6]=$line[5]";
                $wise_out .= $line;
        }
        else {
                $wise_out .= $line;
        }
}
$wise{$tag} = $wise_out;
close INPUT;

@contig{@raw_contigs} = ();
@contigs = keys %contig;

#join overlaped segments
foreach $contig (@contigs){
        $orientation = 'forward';
        merge($orientation,@{$contig.'_forward'});
        $orientation = 'reverse';
        merge($orientation,@{$contig.'_reverse'});
}

open OUT,">raw_homology.gff3";
foreach $contig(@contigs){
        foreach $segment (@{$contig.'_forward'}){
                next if @{$segment.'_forward'} == ();
                $max_score = max (@{$segment.'_forward'});
                $tag = "$segment=+=$max_score";
                if (exists $wise{$tag}){
                        push @result,$wise{$tag};
                }
        }
        foreach $segment (@{$contig.'_reverse'}){
                next if @{$segment.'_reverse'} == ();
                $max_score = max (@{$segment.'_reverse'});
                $tag = "$segment=-=$max_score";
                if (exists $wise{$tag}){
                        push @result,$wise{$tag};
                }
        }
}
@hash{@result} = ();
@last = keys %hash;
foreach $last (@last){
        print OUT "$last//\n";
}
close OUT;
open RAW,"raw_homology.gff3";
open OUT,">homology.gff3";
while ($line = <RAW>){
	print OUT $line and next if $line =~ /\/\//;
	@lines=split(/\t/,$line);
	if($lines[0]=~/^(.+)\_(\d+)\_\d+$/){
		$scaf=$1;
		$ns=$2;
		$lines[0]=$scaf;
		$lines[3]=$lines[3]+$ns-1;
		$lines[4]=$lines[4]+$ns-1;
		$_=join("\t",@lines);
	}
	print OUT "$_";
}
close RAW;
close OUT;


############################# functions ################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub min_max {
        @input_a = @_;
        @input_sort_b = sort {$a <=> $b} @input_a;
        return ($input_sort_b[0],$input_sort_b[-1]);
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub merge_two {
        ($a,$b) = @_;
        ($a_start,$a_end) = split /\t/,$a;
        ($b_start,$b_end) = split /\t/,$b;
        ($c,$d) = min_max($a_start,$a_end,$b_start,$b_end);
        return $c."\t".$d;
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub merge_group {
        @input = @_;
        $count = @input;
	@out =();
	if ($count == 1){
		@out = @input;
	}
	if ($count > 1){
	@input_sort = sort {(split /\t/,$a)[0] <=> (split /\t/,$b)[0]} @input;
        for($i=0;$i<$count;$i++){
                if ($i == 0){
                        $temp = $input_sort[$i];
                }
                else {
                        $current = $input_sort[$i];
                        ($temp_start,$temp_end)=split /\t/,$temp;
                        ($current_start,$current_end)=split /\t/,$current;
                        if ($current_start <= $temp_end){
                                $temp = merge_two($temp,$current);
                                if ($i == $count -1){push @out,$temp;}
                        }
                        else {
                                push @out,$temp;
                                $temp = $current;
				if ($i == $count -1){push @out,$temp;}
                        }
                }
        }
	}
        return @out;
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub merge {
        ($orientation,@input) = @_;
        $temper = '';
        $current = '';
        @memery = ();
        $max_current_score = '';
        @input_sort = sort {(split /\_/,$a)[1] <=> (split /\_/,$b)[1]} @input;
        for($i = 0;$i < scalar @input_sort;$i++){
                if ($i == 0){
                        $temper = $input_sort[$i];
                        $temper_start = (split /\_/,$temper)[1];
                        $temper_end = (split /\_/,$temper)[2];
                }
                else {
                        $current = $input_sort[$i];
                        $current_start = (split /\_/,$current)[1];
                        $current_end = (split /\_/,$current)[2];
                        if ($current_start <= $temper_end){
                                $temper_start = min($temper_start,$temper_end,$current_start,$current_end);
                                $temper_end = max($temper_start,$temper_end,$current_start,$current_end);
                                $max_current_score = max (@{"$current\_$orientation"});
                                push @{"$temper\_$orientation"},$max_current_score;
                                push @memery,$current;
                        }
                        if ($current_start > $temper_end){
                                foreach $member (@memery){
                                        @{"$member\_$orientation"} = @{"$temper\_$orientation"};
                                }
                                @memery = ();
                                $temper = $input_sort[$i];
                                $temper_start = (split /\_/,$temper)[1];
                                $temper_end = (split /\_/,$temper)[2];
                        }
                }
        }
        foreach $member (@memery){
                @{"$member\_$orientation"} = @{"$temper\_$orientation"};
        }
}                      
