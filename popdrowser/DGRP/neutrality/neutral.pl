#!/usr/bin/perl

#Neutrality finder.
#
#Input: Multi-aligned Fasta (outgrup first), RAW Recoded sequence
#Output: GFF3

### LIBS
use Bio::SeqIO;
use Getopt::Long;


# INPUT Fasta file with command line options
my $aln = "";
my $winsize = 1000;  #default
my $seq = "";
my $chr = "";
my $lab = "";
GetOptions('aln=s' => \$aln, 'seq=s' => \$seq, 'winsize=i' => \$winsize , 'chr=s' => \$chr, 'lab=s' => \$lab);
die("Usage = neutral.pl -aln <MAF file with outgrup first> -seq <RAW Recoded seq> [-winsize <int> (default=1000) ] > output\n") 
if(!$aln || !$seq);

#create Bio::Seq objects for all the aligned sequences
my $input_aln = Bio::SeqIO->new( '-file' => $aln, '-format' => 'Fasta');

#open recoded sequence
open SEQ, $seq or die $!;
my $input_seq = "";
	while (<SEQ>) {
		next if /^(\s)*$/;  # skip blank lines
		chomp;              # remove trailing newline characters
		$input_seq = $_;
	}
close SEQ;

#put each sequence in a new array 
my @seqArray;
push(@seqArray, $input_seq );                                #Recoded sequence first (array[0])
	while (my $seq_input = $input_aln->next_seq() ) {    #Sequences in alignment (outgrup first , array[1])
		push(@seqArray, $seq_input->seq );
	}


#Sliding windows loop variables
my $len = length($seqArray[0]);
die ("Window size must be smaller than total seq. lenght(Total lenght:$len bp)\n") if ($len < $winsize);
my $numseq = scalar(@seqArray) -2;   # number of aligned sequences without recoded and outgrup sequence;

######
#Ps and Ds from the whole sequence
my $total_ps = 0;
my $total_ds = 0;

    COLUMN: for(my $i = 0; $i <= $len-1; $i++) {
				my @column = ();
				foreach my $seq (@seqArray) {
					my $pos = substr($seq,$i,1);
					push(@column,$pos);
				}
				
				if ($column[0] == 4 && $column[1] ne "-") {    #Count only 4 fold elements / discard positions with gaps in outgrup
				
					shift @column;                             #delete first (4fold number)
					my $outgrup = shift @column;               #get outgrup nucleotide and delete from array
					my $numseqs = scalar @column;              #number of polymorphic seq
					
					#count nucleotides only in polymorphics sequences.
					my %nucleotides = (
							'A' => scalar(grep {uc($_) eq "A"} @column),
							'T' => scalar(grep {uc($_) eq "T"} @column),
							'G' => scalar(grep {uc($_) eq "G"} @column),
							'C' => scalar(grep {uc($_) eq "C"} @column),
							'N' => scalar(grep {uc($_) eq "N"} @column),
					);
					

					next COLUMN if $nucleotides{N} == $numseq;     #Check if all positions are not N
					foreach my $key (%nucleotides){                #Loop for each nucleotide count
						next if $key eq "N";
						next if $nucleotides{$key} == 0;
						if($nucleotides{$key} == $numseqs ||
						   $nucleotides{$key}+$nucleotides{N} == $numseqs) {    #No polymorphism (All sequences with the same nuc. or only one nucl. with Ns)
							if ($key ne $outgrup){                              #check divergence
								$total_ds++;
								next COLUMN;
							}      
						} else {                                                #We have polymorphism
							if ($key eq $outgrup){                              #check divergence
								$total_ps++;                                    #if no divergence, count polymorphism 
							}
							next COLUMN;
						}
					}
				}
	} #-> END COLUMN loop
$total_ds = JC_correction($len,$total_ds);
#print "Ds: $total_ds - Ps: $total_ps \n\n";
#print "Length: $len\n\n";

########
#Sliding window Loop for pn and dn
	SWINDOW: for(my $i = 0; $i <= $len-1; $i = $i+$winsize) {
		
		my @working_seqs = ();
		my $pn = 0;
		my $dn = 0;

		#get slide for each sequence
		foreach my $seq (@seqArray) {
			my $window = substr($seq,$i,$winsize);
			my @splited = split("",$window);
			push(@working_seqs,[@splited]);
		}

		#Compare each ($p)position in the slide for each ($s)sequence (Columns)
		WCOLUMN: for $p ( 0 .. $#{$working_seqs[0]} ){
			
			my @column = ();
			for $s (0..$#working_seqs){
				push (@column,$working_seqs[$s][$p]);
			} 
			
			#Compare nucleotides
			if ($column[0] != 4 && $column[1] ne "-") {    #Count only NON-4fold elements / discard positions with gaps in outgrup
			
				shift @column;                             #delete array[0] (4fold number)
				my $outgrup = shift @column;               #get outgrup nucleotide (array[1]) and delete from array
				my $numseqs = scalar @column;              #number of polymorphic seq
				
				#count nucleotides (only polymorphic sequences remain in array).
				my %nucleotides = (
						'A' => scalar(grep {uc($_) eq "A"} @column),
						'T' => scalar(grep {uc($_) eq "T"} @column),
						'G' => scalar(grep {uc($_) eq "G"} @column),
						'C' => scalar(grep {uc($_) eq "C"} @column),
						'N' => scalar(grep {uc($_) eq "N"} @column),
				);
				
				#Nucleotide Comparation
				next WCOLUMN if $nucleotides{N} == $numseq;    #Check if all positions are not N
				foreach my $key (%nucleotides){                #Loop for each nucleotide count
					next if $key eq "N";
					next if $nucleotides{$key} == 0;
					if($nucleotides{$key} == $numseqs ||
					   $nucleotides{$key}+$nucleotides{N} == $numseqs) {    #No polymorphism (All sequences with the same nuc. or only one nucl. with Ns)
						if ($key ne $outgrup){                              #check divergence
							$dn++;
							next WCOLUMN;
						}      
					} else {                                                #We have polymorphism
						if ($key eq $outgrup){                              #check divergence
							$pn++;                                    #if no divergence, count polymorphism 
						}
						next WCOLUMN;
					}
				}
			}
		} #-> END Wcolumn
		
		#print "Dn: $dn - Pn: $pn \n";
		#Calculate Neutrality
		my $numsites = $#{$working_seqs[0]}+1;
		$dn = JC_correction($numsites,$dn);
		
		my ($chi_square, $p_value, $NI, $alpha) = CalculateNeutrality($total_ps, $pn, $total_ds, $dn);
		
		#output data
		#print "$ps, $ds, $dn, $pn \n";
		#$p_value < 0.05 ? print $NI."\n" : print "NULL\n";
		print $i."\t".($i+$winsize)."\t".$NI."\t".$p_value."\t".$total_ps."\t".$total_ds."\t".$dn."\t".$pn."\n";
		
	} #-> END WINDOW loop
	



###########################
######## FUNCTIONS ########
###########################

sub CalculateNeutrality {
	
	#input
	my ($ps, $pn, $ds, $dn) = @_;
	
	#output
	my ($chi_square, $p_value, $NI, $alpha)= "" ; 

	if ( (($dn == 0) and ($pn == 0)) or (($ds == 0) and ($ps == 0)) or 
		(($dn == 0) and ($ds == 0))  or (($pn == 0) and ($ps == 0)) ){
		$chi_square='NULL';
		$p_value='NULL'; 
	} else {
		my $intermediate=($dn*$ps)-($pn*$ds);
		$chi_square=(($dn+$ds+$pn+$ps)*($intermediate*$intermediate))/(($dn+$pn)*($ds+$ps)*($dn+$ds)*($pn+$ps));
		$p_value = chisqrprob(1,$chi_square);
	}
	if (($dn == 0) or ($ps == 0)) {
		$NI='NULL';
		$alpha='NULL';
	} else {
		$NI=($pn*$ds)/($ps*$dn);
		$alpha= 1-$NI;
	}
	
	return ($chi_square, $p_value, $NI, $alpha);

}

sub JC_correction {

	my ($numSites,$divergence)=@_;
	
	if ($numSites != '0') 
		{
		my $divergence2correct = $divergence/$numSites;

		#correccion JC
		if (((1-((4/3)*$divergence2correct)) > 0) )
			{
			$divergence_JC = (-(3/4)*log(1-((4/3)*$divergence2correct))); 
			$divergence_JC = ($divergence_JC*$numSites);
			}
		else
			{
			$divergence_JC = ''; 
			}
		}
	else {
		$divergence_JC='';
		}
	return ($divergence_JC);
}


## From statistics::distributions ##
####################################

sub chisqrprob {                                              # Upper probability   X^2(x^2,n)
	my ($n,$x) = @_;
	if (($n <= 0) || ((abs($n) - (abs(int($n)))) != 0)) {
		die "Invalid n: $n\n";                                # degree of freedom
	}
	return precision_string(_subchisqrprob($n, $x));
}

sub _subchisqrprob {
	use constant PI => 3.1415926536;

	my ($n,$x) = @_;
	my $p;

	if ($x <= 0) {
		$p = 1;
	} elsif ($n > 100) {
		$p = _subuprob((($x / $n) ** (1/3)
				- (1 - 2/9/$n)) / sqrt(2/9/$n));
	} elsif ($x > 400) {
		$p = 0;
	} else {   
		my ($a, $i, $i1);
		if (($n % 2) != 0) {
			$p = 2 * _subuprob(sqrt($x));
			$a = sqrt(2/PI) * exp(-$x/2) / sqrt($x);

			$i1 = 1;
		} else {
			$p = $a = exp(-$x/2);
			$i1 = 2;
		}

		for ($i = $i1; $i <= ($n-2); $i += 2) {
			$a *= $x / $i;
			$p += $a;
		}
	}
	return $p;
}

sub _subuprob {
	my ($x) = @_;
	my $p = 0; # if ($absx > 100)
	my $absx = abs($x);

	if ($absx < 1.9) {
		$p = (1 +
			$absx * (.049867347
			  + $absx * (.0211410061
			  	+ $absx * (.0032776263
				  + $absx * (.0000380036
					+ $absx * (.0000488906
					  + $absx * .000005383)))))) ** -16/2;
	} elsif ($absx <= 100) {
		for (my $i = 18; $i >= 1; $i--) {
			$p = $i / ($absx + $p);
		}
		$p = exp(-.5 * $absx * $absx) 
			/ sqrt(2 * PI) / ($absx + $p);
	}

	$p = 1 - $p if ($x<0);
	return $p;
}

sub precision {
	use constant SIGNIFICANT => 5; # number of significant digits to be returned
	my ($x) = @_;
	return abs int(log10(abs $x) - SIGNIFICANT);
}

sub precision_string {
	my ($x) = @_;
	if ($x) {
		return sprintf "%." . precision($x) . "f", $x;
	} else {
		return "0";
	}
}

sub log10 {
	my $n = shift;
	return log($n) / log(10);
}

