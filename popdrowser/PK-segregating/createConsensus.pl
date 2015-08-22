#!/usr/bin/perl

#Neutrality finder.
#
#Input: Multi-aligned Fasta (outgrup first), RAW Recoded sequence
#Output: GFF3

### LIBS
use Bio::SeqIO;
use Getopt::Long;
use Math::BigInt;


# INPUT Fasta file with command line options
my $aln = "";
my $winsize = 1000;  #default
my $seq = "";
my $chr = "";
my $lab = "";
GetOptions('aln=s' => \$aln, 'seq=s' => \$seq, 'winsize=i' => \$winsize , 'chr=s' => \$chr, 'lab=s' => \$lab);
die("Usage = neutral2Ps.pl -aln <MAF file with outgrup first> -seq <RAW Recoded seq> [-winsize <int> (default=1000) -chr <chromosome> -lab <label>] > output\n") if(!$aln || !$seq);

#$seq = "test/test_seq.txt";
##$winsize = int(length($seq)/4);
#$aln = "test/test_aln.fasta";

#open recoded sequence
open SEQ, $seq or die $!;
my $recoded_seq = "";
	while (<SEQ>) {
		next if /^(\s|\t)+$/;  # skip blank lines
		chomp;              # remove trailing newline characters
		$recoded_seq = $_;
	}
close SEQ;

##put each sequence in a new array 
(!&readfasta($aln,\@seqArray)) && do{die "$!\n";};


#Sliding windows loop variables
my $len = length($seqArray[0]);
die ("Window size must be smaller than total seq. lenght(Total lenght:$len bp)\n") if ($len < $winsize);


my $whole_outgroup = ">".$chr.".dsim\n";
my $whole_consensus = ">".$chr.".dgrp_consensus\n";

########

#Sliding window Loop for Pa
SWINDOW: for(my $i = 0; $i < $len; $i += $winsize){
	

	my $sw_outgroup = "";
	$sw_outgroup = substr($seqArray[0], $i, $winsize );
	
	my $sw_preconsensus = "";
	for my $cns (1..$#seqArray){
		$sw_preconsensus .= substr($seqArray[$cns], $i, $winsize );
		$sw_preconsensus .= "\n";
	}


	my ($sw_consensus) = doConsensus($sw_preconsensus);
	
	$whole_outgroup .= $sw_outgroup;
	$whole_consensus .= $sw_consensus;
	
} #-> END WINDOW loop

$whole_outgroup .= "\n";
$whole_consensus .= "\n";

print $whole_outgroup.$whole_consensus;

exit 0;


######## ____________FUNCTIONS____________ ########
###################################################

sub readfasta($$){
	my ($file,$seqs)=@_;
	(!open(FASTA,$file)) && do{return 0;};	
	my $id="";
	my $seq="";
	while(<FASTA>){
		chomp;	
		if ($_=~/^>(\S+)/){
			push(@$seqs,$seq) if ($id ne "");
			$seq="";
			$id=$1;
		}else{
			$seq.=$_;
		}
	}
	push(@$seqs,$seq) if ($id ne "");
	close(FASTA);
	return 1;
}


sub doConsensus {
	use consensus;
	use vars qw(@ISA $_NO_CONSENSUS_SYMBOL);
	use UnivAln;
	@ISA                  = qw(Bio::UnivAln);
	$_NO_CONSENSUS_SYMBOL = $Bio::UnivAln::_NO_CONSENSUS_SYMBOL;
	
	my ($raw)=@_;

	my $self = Bio::UnivAln->new( -seqs => $raw );
	my $consensus = greedy( $self, 0.5, 1 );
	
	return $consensus;
}
