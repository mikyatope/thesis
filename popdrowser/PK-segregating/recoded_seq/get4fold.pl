#!/usr/bin/perl

# Get Nfold degenerate coding nucleotides for whole chromosomes.

# LIBS
use Bio::SeqIO;
use Getopt::Long;

# CODONS
my %degeneracies = (	
	"TTT" => [0,0,2],	"TCT" => [0,0,4],	"TAT" => [0,0,2],	"TGT" => [0,0,2],
	"TTC" => [0,0,2],	"TCC" => [0,0,4],	"TAC" => [0,0,2],	"TGC" => [0,0,2],
	"TTA" => [2,0,2],	"TCA" => [0,0,4],	"TAA" => [0,0,2],	"TGA" => [0,0,0],
	"TTG" => [2,0,2],	"TCG" => [0,0,4],	"TAG" => [0,0,2],	"TGG" => [0,0,0],
	
	"CTT" => [0,0,4],	"CCT" => [0,0,4],	"CAT" => [0,0,2],	"CGT" => [0,0,4],
	"CTC" => [0,0,4],	"CCC" => [0,0,4],	"CAC" => [0,0,2],	"CGC" => [0,0,4],
	"CTA" => [2,0,4],	"CCA" => [0,0,4],	"CAA" => [0,0,2],	"CGA" => [2,0,4],
	"CTG" => [2,0,4],	"CCG" => [0,0,4],	"CAG" => [0,0,2],	"CGG" => [2,0,4],

	"ATT" => [0,0,2],	"ACT" => [0,0,4],	"AAT" => [0,0,2],	"AGT" => [0,0,2],
	"ATC" => [0,0,2],	"ACC" => [0,0,4],	"AAC" => [0,0,2],	"AGC" => [0,0,2],
	"ATA" => [0,0,2],	"ACA" => [0,0,4],	"AAA" => [0,0,2],	"AGA" => [2,0,2],
	"ATG" => [0,0,0],	"ACG" => [0,0,4],	"AAG" => [0,0,2],	"AGG" => [2,0,2],

	"GTT" => [0,0,4],	"GCT" => [0,0,4],	"GAT" => [0,0,2],	"GGT" => [0,0,4],
	"GTC" => [0,0,4],	"GCC" => [0,0,4],	"GAC" => [0,0,2],	"GGC" => [0,0,4],
	"GTA" => [0,0,4],	"GCA" => [0,0,4],	"GAA" => [0,0,2],	"GGA" => [0,0,4],
	"GTG" => [0,0,4],	"GCG" => [0,0,4],	"GAG" => [0,0,2],	"GGG" => [0,0,4]
	);


## INPUT 
my $cdsfile = "";
my $seqfile = "";
my $outfile = "";
my $size = "";
GetOptions('cds=s' => \$gff_file, 'seq=s' => \$seqfile, 'out=s' => \$outfile);
die("Usage = get4fold.pl -cds <gff file>  -seq <fasta>  \n") if(!$gff_file || !$seqfile);

my $seqinput = Bio::SeqIO->new( '-file' => $seqfile, '-format' => 'Fasta');
my $seq = $seqinput->next_seq;


## CREATE $size-LENGHT STRING WITH N's
my $fourfold = "N" x $seq->length ;

## PARSE ANNOTATION INFO ##
open (GFF_FILE,$gff_file) or die $!;

#create array of mRNA_IDs and array-of-hushes of cds_info with sequences
my @mrna_id = ();
my @cds_all = ();
my @utrs_all = ();

while (my $gffline = <GFF_FILE>) {
	
	#skip blank lines in file
	if ($gffline !~ /^(\s)*$/){  
		                                     
		my @gff_fields = split (/\t/,$gffline);
		if ($gff_fields[2] eq "mRNA"){
			$gff_fields[8] =~ /^.*ID=(\w+);.*$/;
			push @mrna_id, $1;
		}elsif ($gff_fields[2] eq "CDS"){
			$gff_fields[8] =~ /^.*Parent=(.*).*$/;
			my $start_cds = $gff_fields[3];
			my $end_cds = $gff_fields[4];
			push @cds_all, { 
							'parent' => $1,
							'start'  =>	$start_cds,
							'end'    =>	$end_cds,
							'frame'  => $gff_fields[7],
							'strand' => $gff_fields[6],
							'seq'    => $seq->subseq($start_cds,$end_cds)
							};
		}elsif ($gff_fields[2] eq "five_prime_UTR" || $gff_fields[2] eq "three_prime_UTR"){
			$gff_fields[8] =~ /^.*Parent=(.*).*$/;
			my $start_utr = $gff_fields[3];
			my $end_utr = $gff_fields[4];
			push @utrs_all, { 
							'parent' => $1,
							'start'  =>	$start_utr,
							'end'    =>	$end_utr,
							};	
		}

	}
}
close GFF_FILE;


## RE-CODE SEQUENCE ##
#Stuff for each mRNA and its CDSs
foreach my $id (@mrna_id){

	#create sub-array of CDSs for current mRNA
	my $whole_cds = "";
	my @current_cds = ();
	my @parents = ();
	
	for my $i (0..$#cds_all){
		#split if there are multiple parents
		@parents = split /,/,$cds_all[$i]{parent};
		for my $n (0..$#parents){
			if ($parents[$n] eq $id){
				push(@current_cds,$cds_all[$i]);
			}
		}
	}

	#Do stuff only if complete cds has length multiple of 3 
	foreach  (@current_cds){
		$whole_cds .= $_->{seq};
	}
	
	if ((length($whole_cds)%3) == 0){
				
		#Translate UTRs	
		@parents = ();
		for my $i (0..$#utrs_all){
			
			#split if there are multiple parents
			@parents = split /,/,$utrs_all[$i]{parent};
			for my $j (0..$#parents){
				if ($parents[$j] eq $id){
					$utr_length = ($utrs_all[$i]{end} - $utrs_all[$i]{start})+1;
					my $get_fold = substr ($fourfold,$utrs_all[$i]{start}-1,$utr_length);
					my $us = "";
		
					#Only substitute N's
					foreach my $char (split //, $get_fold) {
						if ($char eq "N"){
							$us .= "U";
						} else {
							$us .= $char;
						}
					}
					substr ($fourfold,$utrs_all[$i]{start}-1,$utr_length,$us);		
				}
			}
		}		
	
		#"translate" each exon. Save nucleotides for the next exon when not ended in triplet.
    	my $orphan = "";
    	
    	for my $i (0..$#current_cds){
			
			######
			my $pos = 0;
			while($pos < (length($current_cds[$i]{seq}))){
				
				#if not "orphan" nucleotides from inexact codons, get full codon 
				if($orphan eq ""){
					my $codon = substr($current_cds[$i]{seq},$pos,3);
					
					#if exact triplet -> do normally
					if((length($codon)) == 3){
						
						#translate codon
						my $subs = codonTranslate($current_cds[$i]{strand},$codon);        
						
						#get equivalent codon from the coded new sequence
						my $fold = substr ($fourfold,$current_cds[$i]{start}+$pos-1,3);    
						
						#compare the 2 codons and return the most conservate positions
						my $comp = codonCompare($subs,$fold);
						
						#substitute
						substr ($fourfold,$current_cds[$i]{start}+$pos-1,3,$comp);
						$pos+=3;

					#if not exact triplet: save orphan nucleotides for the next exon.
					} elsif ((length($codon)) < 3) {
						$orphan = $codon;
						$pos+=length($orphan);
					}
					
				#if orphans from the last exon, add to the next nucleotides to make an exact triplet.	
				} elsif ($orphan ne "") {
					
					#translate codon
					my $next = substr($current_cds[$i]{seq},$pos,3-(length($orphan)));
					my $codon = $orphan.$next;
					my $subs = codonTranslate($current_cds[$i]{strand},$codon);
					
					#get equivalent codon from the new coded sequence
				    my $fold_orphan = substr ($fourfold,$current_cds[$i-1]{end}-(length($orphan)),length($orphan) );
				    my $fold_next = substr ($fourfold,$current_cds[$i]{start}-1,length($next) ) ;
				    my $fold = $fold_orphan.$fold_next;
				    
				    #compare the 2 codons
				    my $comp = codonCompare($subs,$fold);
				    
				    #Split new codon codification
				    #previous cds translation
					my $subs_orphan = substr ($comp,0,length($orphan));
					substr ($fourfold,$current_cds[$i-1]{end}-(length($orphan)),length($orphan),$subs_orphan) ;
					#current cds translation
					my $subs_next = substr ($comp,length($orphan),length($next));		
					substr ($fourfold,$current_cds[$i]{start}-1,length($next),$subs_next) ;
										
					$orphan = "";
					$pos+= length($next);
				}
			}
			######
			
		}
	}
	
}


print "$fourfold";


print "\n\n";




####### SUBROUTINES ########

#Translate codon to degeneration codes taking strand into acount
sub codonTranslate {
	my ($strand,$codon) = @_;
	if($strand eq "+"){
		$subs = join("",@{ $degeneracies{$codon} });
	} elsif($strand eq "-") {
		$codon =~ tr/ATGC/TACG/ ;
		$codon= reverse $codon;
		my $pre_subs = join("",@{ $degeneracies{$codon} });
		$subs = reverse $pre_subs;	
	}
	return $subs;
}

#Compare codon nucleotides to translate and keep the most conserved ones
sub codonCompare {
	my ($current,$fromfold) = @_;
	my $subs = "";

	if ($fromfold eq "NNN"){
		$subs = $current;
	} elsif ($current ne $fromfold) {
		my @split_current = split //, $current;
		my @split_fromfold = split //, $fromfold;

	    ####
		for (my $i=0; $i<=2; $i++){		
			if ($split_current[$i] ne $split_fromfold[$i] && $split_fromfold[$i] ne "N" && $split_fromfold[$i] ne "U"){
				if ( int($split_current[$i]) > int($split_fromfold[$i]) ) {
					$subs .= $split_fromfold[$i];
				} else {
					$subs .= $split_current[$i];
				}
			} else {
				$subs .= $split_current[$i];
			}
		} 
	    #####
	 
	} else {
		$subs = $current;
	}
	return $subs;
}	
		
		
