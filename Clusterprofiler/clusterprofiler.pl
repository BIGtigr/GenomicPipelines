#!/usr/bin/perl
use Getopt::Long;
use Statistics::R;
use Bio::TreeIO;
use Bio::Tree::TreeI;
#use Bio::Tree::TreeFunctionsIi;

$usage = <<USAGE;
About this script:
    Version : 20210301/001
    Mainly used for GO enrichment analysis
Usage:
    perl clusterprofiler.pl -a=forground -b=background
    pAdjustMethod: one of "holm", "hochberg", "hommel", "bonferroni", "BH","BY", "fdr", "none"
    2 parameters should be given

For example:
    perl $0 -a=111 -b=222

USAGE

if (@ARGV < 2){die $usage}

GetOptions('for|a=s' => \$file1,
           'back|b=s' => \$file2,
          );

`cp $file1 ./aaa.txt`;
`cp $file2 ./bbb.txt`;

$Rbin = '/home/zdh/Software/R-3.5.3/bin/R';
$R=Statistics::R->new(r_bin => $Rbin);
$R->startR;
$R->run(q'library(clusterProfiler)');
$R->run(q'library(org.Hs.eg.db)');

#$R->set('file1', "$file1");
#$R->set('file2', "$file2");
$R->send("degenes <- read.table('aaa.txt',head=FALSE)");
$R->send('genelist <- degenes$V1');
$R->send('genelist1 <- genelist[!duplicated(genelist)]');
$R->send("genelist2 <- as.vector(genelist1)");
$R->send("x1 <- c(genelist2)");
$R->send('eg1 <- bitr(x1, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")');
$R->send('genelist3 <- eg1$ENTREZID');

$R->send("background <- read.table('bbb.txt',head=FALSE)");
$R->send('backgene <- background$V1');
$R->send('backgene1 <- backgene[!duplicated(backgene)]');
$R->send("backgene2 <- as.vector(backgene1)");
$R->send("x2 <- c(backgene2)");
$R->send('eg2 <- bitr(x2, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")');
$R->send('backgene3 <- eg2$ENTREZID');

print"Conducting GO enrichment analysis\n";
$R->send("go <- enrichGO(genelist3, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',universe = backgene3,pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')");
$R->send('gos <- simplify(go, cutoff = 0.7, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)');

print"Writing results into cvs file\n";
$R->send('write.csv(summary(gos),"GO-enrich.csv",row.names =F)');

print"Drawing pictures\n";
$R->send('barplot(gos,showCategory=20,drop=T)');
$R->send('dotplot(gos,showCategory=50)');

$R->stopR();
`rm aaa.txt`;
`rm bbb.txt`;

