#!/usr/bin/perl

#Neutrality finder.
#
#Input: Multi-aligned Fasta (outgrup first), RAW Recoded sequence
#Output: GFF3

### LIBS
use Bio::SeqIO;
use Getopt::Long;
use Tie::File;


# INPUT Fasta file with command line options
my $aln = "";
my $winsize = 1000;  #default
my $seq = "";
my $chr = "";
my $lab = "";
my $xfold = "4"; #this can be a regexp
my $fixnum = 150;
my $discard = 8; ####### Functionality not tested! based on fixnum now.
GetOptions('aln=s' => \$aln, 'seq=s' => \$seq, 'winsize=i' => \$winsize , 'chr=s' => \$chr, 'lab=s' => \$lab, 'fold=s' => \$xfold, 'fixnum=i' => \$fixnum, , 'discard=i' => \$discard);
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
#(!&readfasta($aln,\@seqArray)) && do{die "$!\n";};

#test Tie::File
tie my @fileArray, 'Tie::File', $aln or die "$!\n";

my @seqArray = ();

my $fasta_id="";
my $fasta_seq="";
for (@fileArray) {
	#chomp;	
	if ($_=~/^>(\S+)/){
		push(@seqArray,$fasta_seq) if ($fasta_id ne "");
		$fasta_seq="";
		$fasta_id=$1;
	}else{
		$fasta_seq.=$_;
	}
}
push(@seqArray,$fasta_seq) if ($fasta_id ne "");


#Sliding windows loop variables
my $len = length($seqArray[0]);
die ("Window size must be smaller than total seq. lenght(Total lenght:$len bp)\n") if ($len < $winsize);
my $numArray = $#seqArray;   # number of aligned sequences;

########
#combinatorial number
my $n = $fixnum;
my $nbin2 = $n*($n-1)/2;

#number of discarded sites to consider
my $discard = $numArray-$fixnum if $fixnum <= $numArray or die("discarded($fixnum) larger than total($n)");

#Sliding window Loop for Pa
SWINDOW: for(my $i = 0; $i < $len; $i += $winsize){
	
	#get 4fold positions (array_ref)		
	my $positions= getPos($i, $recoded_seq, $xfold, $winsize);
	
	my %diff = ();
	
	#comparisons count to 0
	for my $j (0 .. $fixnum-1){
		for my $k (0 .. $fixnum-1){
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
	
		#get the column in the alignment for the current position
		my $col ="";		
		for (0 .. $numArray){
			my $nuc=uc(substr($seqArray[$_],$pos,1));
			$col .= $nuc;
		}	
		
		#count missing and heterozigous sites in position
		#discard sites with many N's and gaps from the analysis
		my $count_miss = ($col =~ tr/N//);
		my $count_het  = ($col =~ tr/(M|R|W|S|Y|K)//);
		my $count_gaps = ($col =~ tr/-//);
		$miss += $count_miss;
		$gaps += $count_gaps;
		$het  += $count_het;
		
		my $total_non_valid = $count_gaps + $count_miss + $count_het;
				
		next POSITION if $total_non_valid > $discard;                        
		
		$col =~ s/(N|M|R|W|S|Y|K|-)//g ;   #remove N's, gaps or heterozigous sites if any
		next POSITION if length($col) < $fixnum;
		my @col_split = split(//, $col);  #transform scalar into array
		shuffle(\@col_split);
		my @slice = @col_split[0..$fixnum-1]; #get the first $fixnum elements of the shuffled array
		
		#filters passed, count column as valid and proceed to analyze
		$analyzed++; 
				
		#pairwise comparision for each element in the array
		for my $j (0 .. $fixnum-1){
			for my $k (0 .. $fixnum-1){
				next if $j <= $k;
				my $hname = $j."-".$k;
				my $nucj=uc($slice[$j]);
				my $nuck=uc($slice[$k]);
				my $pair = $nucj.$nuck;
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
	my $m = $analyzed ;
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
	
	print $chr."\t".$lab."\tremark\t".$i."\t".$end."\t".$ps."\t.\t.\tName=".$lab.";total_".$xfold."fold_sites=$m_total;analyzed=$analyzed;Ns=$miss;heterozigous_sites=$het;gaps=$gaps;combinations=$nbin2\n";

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
		my $pos = length($`);                   #lenght of string before the match, equivalent to relative position of the match starting with 0
		$pos = $pos + $i if $i > 0;
		push( @positions, $pos );
	}
	return (\@positions);
}

# randomly permutate @array in place
sub shuffle
{
	#fisher yates shuffle
    my $array = shift;
    my $i = @$array;
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}
