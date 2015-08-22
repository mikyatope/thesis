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
my $jc = 0;
my $low = 0;
GetOptions('aln=s' => \$aln, 
	'seq=s' => \$seq, 
	'winsize=i' => \$winsize , 
	'chr=s' => \$chr, 
	'lab=s' => \$lab, 
	'fold=s' => \$xfold, 
	'fixnum=i' => \$fixnum, 
	'discard=i' => \$discard,
	'jc=i' => \$jc,
	'low' => \$low
	);
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

##put all sequences in a new array 
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

#$winsize = 10;
#$fixnum = 10;
#my @seqArray = ('ATGCATCGATGC',
#		'ATGCATCGATGC',
#		'ATGCATCGATGC',
#		'ATGCATCGATGC',
#		'ATGCATCGATGT',
#		'ATGCATCGATGT',
#		'ATGCATCGATGT',
#		'ATGCATCGATGT',
#		'ATGCATCGACGA',
#		'ATGCATCGACGA',
#		'ATGCATCGACGA',
#		'AAGCATCGACGA',
#		'AAGTATCGAAGA'
#		);


#Sliding windows loop variables
my $len = length($seqArray[0]);
die ("Window size must be smaller than total seq. lenght(Total lenght:$len bp)\n") if ($len < $winsize);
my $numArray = $#seqArray+1;   # number of aligned sequences ;


#number of discarded sites to consider
my $discard = $numArray-$fixnum if $fixnum <= $numArray or die("fixnum($fixnum) larger than total($numArray)");

#print header
print "chr\tlabel\tstart\tend\tm_total\tanalyzed\tsegregating\tlowfreq\tmiss\thet\tgaps\n";

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
	my $segregating = 0;
	my $lowfreq = 0;
	
	#Count diferences for each 4-fold position
	
        #@$positions = ( 1,2,3,4,5,6,7,8,9,10);

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
		
		$col =~ s/(N|M|R|W|S|Y|K|-)//g ;                #remove N's, gaps or heterozigous sites if any
		next POSITION if length($col) < $fixnum;
		my @col_split = split(//, $col);                #transform scalar into array
		shuffle(\@col_split);                           #randomize array
		my @slice = @col_split[0..$fixnum-1];           #get the first $fixnum elements of the shuffled array
		$col = join('',@slice);                         #array to string again
		
		#filters passed, count column as valid and proceed to analyze
		$analyzed++; 

		#Count nucleotides
		my @countNuc = ();
		$countNuc[0] = ($col =~ tr/A//);
		$countNuc[1] = ($col =~ tr/C//);
		$countNuc[2] = ($col =~ tr/T//);
		$countNuc[3] = ($col =~ tr/G//);
		
#print $col."\n";
#print "  A ".$countNuc[0] ;
#print "    C ".$countNuc[1] ;
#print "    T ".$countNuc[2] ;
#print "    G ".$countNuc[3]."\n";

		my $zeroes = 0;
		my $lowfreq_value = 7;
		my $countdiscard = 0;
		
		#check absent nucleotides in column
		for (0..3) {
			#check for low freq variants
			if ($countNuc[$_] < $lowfreq_value  && $countNuc[$_] != 0 && $low == 1) {
				$countNuc[$_] = 0;
				$countdiscard = 1;
			}
			$zeroes++ if $countNuc[$_] == 0;
		}
		$segregating++ if $zeroes < 3;    #count as segregating if 2 or more nucleotide counts are not 0
		$lowfreq +=  $countdiscard ;
		
#print "zeroes : $zeroes \n";
#print "discard : $countdiscard \n";

	}
	

	my $m_total = $#$positions+1;
	my $m = $analyzed ;
	my $ps = "";

	if ($m == 0 || $analyzed < ($m_total/2) ) {
		$ps = "NA" ;
	} else {
		$ps = $segregating ;
		$ps = sprintf("%.4f", JC_correction($m,$segregating) ) if $jc == 1;
	}
	
	
	my $end = $i+$winsize;
	$end = $len if $end > $len;

	print $chr."\t".$lab."\t".$i."\t".$end."\t".$m_total."\t".$analyzed."\t".$segregating."\t".$lowfreq."\t".$miss."\t".$het."\t".$gaps."\n";
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
sub shuffle {
	#fisher yates shuffle
	my $array = shift;
	my $i = @$array;
	while ( --$i ) {
		my $j = int rand( $i+1 );
		@$array[$i,$j] = @$array[$j,$i];
	}
}


sub JC_correction {
	my ($numSites,$divergence)=@_;
	if ($numSites != '0') {
		my $divergence2correct = $divergence/$numSites;
		#correccion JC
		if (((1-((4/3)*$divergence2correct)) > 0) ){
			$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence2correct))); 
			$divergence_JC = ($divergence_JC*$numSites);
		} else {
			$divergence_JC = "NA"; 
		}
	} else {
		$divergence_JC= "NA";
	}
	return ($divergence_JC);
}
