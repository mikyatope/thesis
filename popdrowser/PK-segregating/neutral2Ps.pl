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
my $xfold = "4"; #this can be a regexp
my $discard = 8;
GetOptions('aln=s' => \$aln, 'seq=s' => \$seq, 'winsize=i' => \$winsize , 'chr=s' => \$chr, 'lab=s' => \$lab, 'fold=s' => \$xfold, 'discard=i' => \$discard);
#die("Usage = neutral2Ps.pl -aln <MAF file with outgrup first> -seq <RAW Recoded seq> [-fold <int/regexp> -winsize <int> (default=1000) -chr <chromosome> -lab <label>] > output\n") if(!$aln || !$seq);

#$seq = "test/test_seq.txt";
#$winsize = int(length($seq)/4);
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
my $numArray = $#seqArray;   # number of aligned sequences;



########

#combinatorial number
my $x = Math::BigInt->new($#seqArray+1);
$x->bnok(2);

$n = $numArray+1;
$nbin2 = $n*($n-1)/2;

#Sliding window Loop for Pa
SWINDOW: for(my $i = 0; $i < $len; $i += $winsize){
	
	#get 4fold positions (array_ref)		
	my $positions= getPos($i, $recoded_seq, $xfold, $winsize);
	
	my %diff = ();
	
	#comparisons count to 0
	for my $j (0 .. $#seqArray){
		for my $k (0 .. $#seqArray){
			next if $j <= $k;
			my $hname = $j."-".$k;             #print $hname."\n";
			$diff{$hname}=0;
		}	
	}

	my $miss = 0;
	my $het  = 0;
	my $gaps = 0;
	my $analyzed = 0;
	#Count diferences for each 4-fold position
	POSITION: foreach my $pos (@$positions){
	
		#count missing and heterozigous sites in position
		my $col ="";		
		for (0 .. $#seqArray){
			my $nuc=uc(substr($seqArray[$_],$pos,1));
			#$gaps++ if $nuc eq "-" && $_ == 0;
			#next POSITION if $nuc eq "-" && $_ == 0;
			#$miss++ if ($nuc =~ tr/(N|M|R|W|S|Y|K)//) > 0;
			#next POSITION if ($nuc =~ tr/(N|M|R|W|S|Y|K)//) > 0;		
			$col .= $nuc;
		}	
		
		#discard sites with many N's and gaps from the analysis
		my $count_miss = ($col =~ tr/N//);
		my $count_het  = ($col =~ tr/(M|R|W|S|Y|K)//);
		my $count_gaps = ($col =~ tr/-//);
		$miss += $count_miss;
		$gaps += $count_gaps;
		$het  += $count_het;
		
		next POSITION if $count_gaps+$count_miss+$count_het > $discard;                        

		$analyzed++; #filters passed, count column as valid
				
		#pairwise comparision
		for my $j (0 .. $#seqArray){
			for my $k (0 .. $#seqArray){
				next if $j <= $k;
				my $hname = $j."-".$k;
				my $nucj=uc(substr($seqArray[$j],$pos,1));
				my $nuck=uc(substr($seqArray[$k],$pos,1));
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
    
	my $m_total = $#$positions+1;
	my $m = $m_total-$invalid ;
	my $ps = 0;
    
    
    
    if ($m == 0 || $analyzed < ($m_total/2) ) {
		print "\n# $ks = ($sum/$nbin2)/$m | m_total=$m_total miss=$miss gaps=$gaps \n";
		$ps = -0.01 ;
	} else {
		$ps = ($sum/$nbin2)/$m;
		$ps = sprintf("%.6f",$ps);
		print "\n# $ks = ($sum/$nbin2)/$m | m_total=$m_total miss=$miss gaps=$gaps \n";
	}
	
	print "# $ps = ($sum/$nbin2)/$m \n";
	
	my $end = $i+$winsize;
	$end = $len if $end > $len;
	
    print $chr."\t".$lab."\tremark\t".$i."\t".$end."\t".$ps."\t.\t.\tName=".$lab.";total_".$xfold."fold_sites=$m_total;analyzed=$analyzed;Ns=$miss;heterozigous_sites=$het;gaps=$gaps\n";
    

} #-> END WINDOW loop


exit 0;
###########################
######## FUNCTIONS ########
###########################

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
		$pos = $pos + $i if $i > 0;
		push( @positions, $pos );
	}
	return (\@positions);
}
