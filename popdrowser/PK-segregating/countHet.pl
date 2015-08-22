	#!/usr/bin/perl

#Neutrality finder.
#
#Input: Multi-aligned Fasta (outgrup first), RAW Recoded sequence
#Output: GFF3

### LIBS
use Getopt::Long;
use Tie::File;


# INPUT Fasta file with command line options
my $aln = "";
my $chr = "";

GetOptions('aln=s' => \$aln, 
	   'chr=s' => \$chr, 
	  );
#die("Usage = ??? \n") if(!$aln || !$seq);




#test Tie::File
tie my @fileArray, 'Tie::File', $aln or die "$!\n";

my %seqHash   = {};
my $fasta_id  = "";
my $fasta_seq = "";

for (@fileArray) {
	#chomp;	
	if ($_=~/^>(\S+)/){
		$fasta_id           = $1;
	}else{
		$seqHash{$fasta_id} = $_;
	}
}

#print header
print "chr\tline\thet\n";

foreach my $key (sort (keys %seqHash)){
	my $countHet = ($seqHash{$key} =~ tr/(M|R|W|S|Y|K)//);
	print $chr."\t".$key."\t".$countHet."\n";
}


exit 0;











