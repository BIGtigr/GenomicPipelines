$usage = <<USAGE;
About this script:
    Version : 20181206/001
    Author  : bigtigerzou@163.com
    Dependencies : orthMCL, BLAST+, MYSQL, MCL.
    This script is used for lazyly running orthomcl. To use it, DIGS-tool(https://github.com/giffordlabcvr/DIGS-tool)
    should be properly installed first. it comprise a lots of useful modules, part of which were used here.
    Before initialise your analysis, please set the following 4 environment variables:\$DIGS_HOME,\$DIGS_MYSQL_USER and \$DIGS_MYSQL_PASSWORD.
    Following is what i "echoed" in my ~/.bashrc file:
		export DIGS_MYSQL_USER=tigerzou
		export DIGS_MYSQL_PASSWORD=123456
		export DIGS_HOME=/home/zdh/Software/DIGS-tool
    Any problem please contact bigtigerzou\@163.com and you will not be answered.
Usage:
    perl $0 -p=<path to pep files>

USAGE
# Include the PERL module library for DIGS
use DBI;
use lib ($ENV{DIGS_HOME}) . '/modules/';
use Getopt::Long;
use Base::Console;
use Base::FileIO;

# Third party program interface modules
use Interface::MySQLtable;   # Interface to BLAST

if (@ARGV < 1){die $usage}
$help = undef;
GetOptions ('path|p=s'     => \$pep_path,
	    'help'           => \$help,
) or die $usage;
if ($help) {die $usage}
$pep_path =~ s/\/$//;
system('mkdir compliantFasta');
opendir PEP,$pep_path;
foreach $fas (readdir PEP){
	next if $fas !~ /aa\.longest$/;
	$sp = (split(/\./,$fas))[0];
	open Afas,"$pep_path/$fas";
	open Bfas,">compliantFasta/$sp.fasta";
	while (<Afas>){
		next if $_ =~ /^$/ or $_ =~ /^>$/;
		if ($_ =~ s/^>//){
			print Bfas ">$sp|$_";
		}
		else{print Bfas "$_";}
	}
	close Afas;
	close Bfas;
}
# filter sequence
`orthomclFilterFasta compliantFasta/ 30 20`;
all2allBLASTP('goodProteins.fasta');
`orthomclBlastParser orthomcl.out compliantFasta > similarSequences.txt`;
$db_name = 'orthomcl'.time;
# MySQL database connection parameters
%mysql = (
        mysql_username         => $ENV{DIGS_MYSQL_USER},
        mysql_password         => $ENV{DIGS_MYSQL_PASSWORD},
        db_name                => $db_name,   # Obtained from control file or console
        mysql_server           => 'localhost'   # Obtained from control file or console
);

create_screening_db($db_name);

$orthomcl_config = "# this config assumes a mysql database named 'orthomcl'.  adjust according
# to your situation.
dbVendor=mysql
dbConnectString=dbi:mysql:$db_name
dbLogin=$mysql{'mysql_username'}
dbPassword=$mysql{'mysql_password'}
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
dbConnectString=dbi:mysql:$db_name:mysql_local_infile=1
";
open CONFIG,">orthomcl.config";
print CONFIG "$orthomcl_config";
close CONFIG;

`orthomclInstallSchema orthomcl.config`;
`orthomclLoadBlast orthomcl.config similarSequences.txt > orthomclLoadBlast.log`;
`orthomclPairs orthomcl.config orthomcl_pairs.log cleanup=no`;
`orthomclDumpPairsFiles orthomcl.config`;
`mcl mclInput --abc -I 1.5 -o mclOutput > mcl.log`;
`orthomclMclToGroups OCG 1 < mclOutput > groups.txt`;

#parsing clustering algorithm output and generate file for cafe(comes from chenlianfu)
open IN, "groups.txt" or die $!;
while (<IN>) {
	if (/^([^:]+?):/) {
		$group_ID = $1; push @group_ID, $group_ID;
		@genes = m/(\w+)\|/g;
		foreach (@genes) {
			$taxon{$_} = 1;
			$gene_family_size{$group_ID}{$_} ++;
		}
	
	}
}
close IN;
open OUT,">orthomcl.cafe.tab";
@taxon = sort { $a cmp $b } keys %taxon;
print OUT "Description\tID";
foreach (@taxon) {
	print OUT "\t$_";
}
print "\n";
foreach my $group_ID (@group_ID) {
	$out = '';
	$out .= "NA\t$group_ID";
	foreach (@taxon) {
		if (exists $gene_family_size{$group_ID}{$_}) {
			$num += $gene_family_size{$group_ID}{$_};
			$out .= "\t$gene_family_size{$group_ID}{$_}";
			$number ++;
		}
		else{
			$out .= "\t0";
		}
	}
	$out .= "\n";
	print OUT $out if $num >= (@taxon / 2) && $number >= 3;
}

############################################################################
# Subroutines
############################################################################

#####################################################################################
# Subroutine: create_screening_db
# Description: creat mysql datebase
#####################################################################################
sub create_screening_db {

        ($db_name) = @_;

        # Get connection variables from self
        $server   = $mysql{'mysql_server'};
        $username = $mysql{'mysql_username'};
        $password = $mysql{'mysql_password'};

        # CREATE THE DB
        print "\n\t ### Creating '$db_name' screening database\n\n";
           $drh = DBI->install_driver("mysql");
           $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
           $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
        unless ($dbh) {
                die "\n\t # Couldn't create screening database '$db_name'\n\n";
        }
}

#####################################################################################
# Subroutine: all2allBLASTP
# Description: all-to-all blastp
#####################################################################################
sub all2allBLASTP {
	$goodProteins = @_;
	`makeblastdb -in goodProteins.fasta -dbtype prot -title orthomcl -parse_seqids -out orthomcl -logfile orthomcl.log > makeblastdb.log`;
	`blastp -db orthomcl -query goodProteins.fasta -seg yes -out orthomcl.out -num_threads 30 -evalue 1e-7 -outfmt 6 > blastp.log`;
}

############################################################################
# End of file
############################################################################
