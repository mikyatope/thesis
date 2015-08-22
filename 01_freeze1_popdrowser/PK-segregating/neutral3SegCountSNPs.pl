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
my $seq = "";
my $chr = "";

GetOptions('aln=s' => \$aln, 
	   'seq=s' => \$seq, 
	   'chr=s' => \$chr, 
	  );
#die("Usage = ??? \n") if(!$aln || !$seq);



#open recoded sequence
open SEQ, $seq or die $!;
my $recoded_seq = "";
	while (<SEQ>) {
		next if /^(\s|\t)+$/;  # skip blank lines
		chomp;              # remove trailing newline characters
		$recoded_seq = $_;
	}
close SEQ;


#test Tie::File
tie my @fileArray, 'Tie::File', $aln or die "$!\n";

my @seqArray  = ();
my $fasta_id  = "";
my $fasta_seq = "";

for (@fileArray) {
	#chomp;	
	if ($_=~/^>(\S+)/){
		push(@seqArray,$fasta_seq) if ($fasta_id ne "");
		$fasta_seq  = "";
		$fasta_id   = $1;
	}else{
		$fasta_seq .= $_;
	}
}
push(@seqArray,$fasta_seq) if ($fasta_id ne "");


#Sliding windows loop variables
my $len      = length($seqArray[0]);
my $numArray = $#seqArray+1;   # number of aligned sequences ;


#print header
print "chr\tposition\tA\tC\tT\tG\tclass\tvalid\tnon_valid\tNs\tgaps\thet\n";

#Sliding window Loop for Pa
POSITION: for(my $i = 0; $i < $len; $i++){

        #get the column in the alignment for the current position
        my $col ="";		
        for (0 .. $numArray){
	        my $nuc=uc(substr($seqArray[$_],$i,1));
	        $col .= $nuc;
        }	
        
        #duplicate column to count missing sites later
        my $col2 = $col;                       
	
	#remove N's, gaps or heterozigous sites if any	
        $col =~ s/(N|M|R|W|S|Y|K|-)//g ;       
        my $valid =  length($col);

        #Count nucleotides
#        my @countNuc = ();
#        $countNuc[0] = ($col =~ tr/A//);
#        $countNuc[1] = ($col =~ tr/C//);
#        $countNuc[2] = ($col =~ tr/T//);
#        $countNuc[3] = ($col =~ tr/G//);

        my %countNuc = ();
        $countNuc{"A"} = ($col =~ tr/A//);
        $countNuc{"C"} = ($col =~ tr/C//);
        $countNuc{"T"} = ($col =~ tr/T//);
        $countNuc{"G"} = ($col =~ tr/G//);

        my $zeroes = 0;

        #check alleles in column. 'Zeroes' sum+1 for each absent allele
        for (keys %countNuc) {
	        $zeroes++ if $countNuc{$_} == 0;
        }

        # Zeroes is less than 3 if the position is an SNP (2 or more alleles)
        if ($zeroes < 3) {
                #Recoded seq value
                my $class = substr($recoded_seq, $i, 1);
                
                if (uc($class) eq "F" || uc($class) eq "L" || $class == 3 ) {
                        $class = check2fold($class, \%countNuc, $valid);
                }
                
                #count missing and heterozigous sites in position
                my $miss  = ($col2 =~ tr/N//);
                my $het   = ($col2 =~ tr/(M|R|W|S|Y|K)//);
                my $gaps  = ($col2 =~ tr/-//);

                my $non_valid = $gaps + $miss + $het;
                  
                print "$chr\t".($i+1)."\t$countNuc{'A'}\t$countNuc{'C'}\t$countNuc{'T'}\t$countNuc{'G'}\t$class\t$valid\t$non_valid\t$miss\t$gaps\t$het\n";
        }
               

	
	
} #-> END loop

exit 0;



###################
## SUBRUTINES 
###################


sub check2fold {
        my $class        = shift;
        my $countNuc_ref = shift;
        my $valid        = shift;
        
        $class = uc $class;          #uper and lower case not used yet
        my $returnClass = "";
        
        #sort keys (nucleotides) descending by value
        my %counts   = %$countNuc_ref;
        my @keys     = sort { $counts{$b} <=> $counts{$a} } (keys %counts);
        
        if ( $counts{$keys[2]} > 0 && $counts{$keys[2]} > ($valid*0.05) ) {
                $returnClass = "Oth";     #special class site if more than 2 alleles in high frequency
        } else {
                
                $a1 = $keys[0];           #allele 1
                $a2 = $keys[1];           #allele 2
                
                if ($class == 3) {

                        $returnClass = 4;
                        $returnClass = 0 if $a1 eq "G" || $a2 eq "G";
                }
                
                if ($class eq "F") {

                        $returnClass = 0;
                        $returnClass = 4 if $a1 eq "G" && $a2 eq "A";
                        $returnClass = 4 if $a1 eq "A" && $a2 eq "G";
                }
                
                if ($class eq "L") {

                        $returnClass = 0;
                        $returnClass = 4 if $a1 eq "G" && $a2 eq "A";
                        $returnClass = 4 if $a1 eq "A" && $a2 eq "G";
                        $returnClass = 4 if $a1 eq "T" && $a2 eq "C";
                        $returnClass = 4 if $a1 eq "C" && $a2 eq "T";
                }
        }
        

        return $returnClass;
}









