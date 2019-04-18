#!/usr/bin/perl
use Getopt::Long;
use Statistics::R;

GetOptions (\%opt,"table:s","help");

$R=Statistics::R->new();
open IN, "$ARGV[0]" or die "$!";
$R->startR;
#<IN>;
while($line = <IN>){
        chomp $line;
        next if ($line eq "");
        @line = split /\t/,$line;
        push @p,$line[-1];
        push @obs,$line;
}

$temp=join(",",@p);
#print "$temp\n";

$R->send("parray <- c($temp)");
$R->send("print (parray)");
$a=$R->read;
#print "P:$a\n";
$R->send("q <- p.adjust(parray, method='fdr')");
$R->send("print (q)");
$q=$R->read;
#print "Q:$q\n";
$R->stopR();

@fdr=split(" ",$q);
for($j=0;$j<@fdr;$j++){
        if ($fdr[$j]=~/\[|\]/){
                splice (@fdr,$j,1);
        }
}
for($j=0;$j<@fdr;$j++){
        $fdr = sprintf("%.8f",$fdr[$j]);
        $term = "$obs[$j]\t$fdr[$j]\t$fdr";
        #print "$obs[$j]\t$fdr[$j]\n";
        push @result,$term;
}
@result_sort = sort {(split /\t/,$a)[7] <=> (split /\t/,$b)[7]} @result;
$pwd = `pwd`;
chomp $pwd;
#`mkdir selected_ortholog`;
for($j=0;$j<@result_sort;$j++){
        print "$result_sort[$j]\n";
        if ((split /\t/,$result_sort[$j])[7] < 0.05){
                $selected_gene++;
                $ortholog_name = (split /\t/,$result_sort[$j])[0];
                $ortholog_name =~ s/\.phy$//;
                #`cp $pwd/prank/$ortholog_name.fas $pwd/selected_ortholog`;
                #`cp $pwd/prank/$ortholog_name.fas.fa $pwd/selected_ortholog`;
        }
}
