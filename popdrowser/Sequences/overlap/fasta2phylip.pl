use Bio::AlignIO;
use Getopt::Long;
use warnings;

my $infile = "";
my $outfile = "";

GetOptions('in=s' => \$infile,  'out=s' => \$outfile);
die("Usage = fasta2phylip.pl -in <fasta file>  -out <phylip>  \n") if(!$infile || !$outfile);

$in  = Bio::AlignIO->new(-file => "$infile" ,
                         -format => 'fasta');
$out = Bio::AlignIO->new(-file => ">$outfile",
                         -format => 'phylip');
while ( my $aln = $in->next_aln() ) { $out->write_aln($aln); }

