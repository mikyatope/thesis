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


#combinatorial number
my $x = Math::BigInt->new(2);
$x->bnok(2);

$n = 2;
$nbin2 = $n*($n-1)/2;


########

#Sliding window Loop for Pa
SWINDOW: for(my $i = 0; $i < $len; $i += $winsize){
	
	#get 4fold positions (array_ref)	
	#$sw_recoded = substr($recoded_seq, $i, $winsize );	
	my $positions= getPos($i, $recoded_seq, '4', $winsize);
	
	my %diff = ();

	my $sw_outgroup = "";
	$sw_outgroup = substr($seqArray[0], $i, $winsize );
	my $sw_preconsensus = "";
	for my $cns (1..$#seqArray){
		$sw_preconsensus .= substr($seqArray[$cns], $i, $winsize );
		$sw_preconsensus .= "\n";
	}


	my ($sw_consensus) = doConsensus($sw_preconsensus);
	
	my @consArray = ();	
	push (@consArray, $sw_outgroup, $sw_consensus); 
#	$sw_recoded = substr($recoded_seq, $i, $winsize );
#	print "\n\n".$sw_recoded."\n".$consArray[0]."\n".$consArray[1]."\n";

	#comparisons count to 0
	for my $j (0 .. $#consArray){
		for my $k (0 .. $#consArray){
			next if $j <= $k;
			my $hname = $j."-".$k;             #print $hname."\n";
			$diff{$hname}=0;
		}	
	}

	#Count diferences for each 4-fold position
	foreach my $pos (@$positions){

		#count missing and heterozigous sites in position
		my $col ="";		
		for (0 .. $#consArray){
			my $nuc=uc(substr($consArray[$_],$pos,1));
			$col .= $nuc;
		}	
		my $d = ($col =~ tr/(N|M|R|W|S|Y|K)//);
		next if $d > 8;                        #print $d." ";
				
		#pairwise comparision
		for my $j (0 .. $#consArray){
			for my $k (0 .. $#consArray){
				next if $j <= $k;
				my $hname = $j."-".$k;
				my $nucj=uc(substr($consArray[$j],$pos,1));
				my $nuck=uc(substr($consArray[$k],$pos,1));
				my $pair = $nucj.$nuck;
				my $Ns = ($pair =~ tr/(N|M|R|W|S|Y|K)//);
				next if $Ns > 0;                        #print "$pair $Ns \n";
				$diff{$hname}++ if $nucj ne $nuck;
			}	
		}
	}
	
	my $sum = 0;
	for my $key ( keys %diff ) {
        my $value = $diff{$key};
        #print "$key => $value\n";
        $sum += $value;
    }
	my $m = $#$positions+1;
	my $ks = 0;
	if ($m == 0) {
		$ks = -0.001;
	} else {
		$ks = ($sum/$nbin2)/$m;
		$ks = sprintf("%.6f",$ks);
	}
	
	print "# $ks = ($sum/$nbin2)/$m \n";
	
	my $ks_JC=JC_correction($m,$ks);
	$ks_JC = sprintf("%.6f",$ks_JC);
	
	print "# $ks_JC = ($sum/$nbin2)/$m \n";
	
	my $end = $i+$winsize;
	$end = $len if $end > $len;
	
    print $chr."\t".$lab."\tremark\t".$i."\t".$end."\t".$ks_JC."\t.\t.\tName=".$lab."\n";

} #-> END WINDOW loop


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

sub getPos {
	my $i           = shift;
	my $recoded_seq = shift;
	my $regex       = shift;
	my $winsize     = shift;

	my $slide_recoded = substr( $recoded_seq, $i, $winsize );
	my @positions = ();
	while ( $slide_recoded =~ m/($regex)/gi ) {
		my $pos = length($`);                                  #lenght of string before the match, equivalent to relative position of the match starting with 0
		#$pos = $pos + $i if $i > 0;
		push( @positions, $pos );
	}
	return \@positions;
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

sub JC_correction {
	my ($numSites,$divergence)=@_;
	
	if ($numSites != '0') {
		#my $divergence2correct = $divergence/$numSites;

		#correccion JC
		if (((1-((4/3)*$divergence2correct)) > 0) ){
			#$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence2correct))); 
			$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence))) or $divergence_JC = 5 ; 
			#$divergence_JC = ($divergence_JC*$numSites);
		} else {
			$divergence_JC = ''; 
		}
	} else {
		$divergence_JC='';
	}
	
	return ($divergence_JC);
}