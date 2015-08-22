#!/usr/bin/perl

=head

axt2fasta.pl 
   by Miquel Ràmia , miquel.ramia@uab.cat
   Universitat Autònoma de Barcelona

This script parses and transforms the "hit" (aligned) sequence of multiple blocks UCSC's .axt alignment 
into a linear sequence, given the "query" sequence. The positions in the aligned sequence corresponding 
to gap positions in the query sequence are deleted, to make the alignment consistent with the query se-
quence coordinates.

Input:
    - Axt file
	- Query fasta sequence
	
Output:
	- Fasta sequence
	



=cut

# LIBS
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

# Input
my $in = "";
my $ref = "";
my $id = "";
GetOptions('in=s' => \$in, 'seq=s' => \$ref, 'id=s' => \$id);
die("Usage = -id [string] -in <axt file>  -seq <fasta> > stdout \n") if(!$in || !$ref);

#Get ref sequence
my $seqinput = Bio::SeqIO->new( '-file' => $ref, '-format' => 'Fasta');
my $seq = $seqinput->next_seq;
##Create an only-gaps string (same lenght as ref sequence)
my $fake_aln = "-" x $seq->length ;


use Bio::SearchIO;
my $parser = Bio::SearchIO->new(-format => 'axt',
							    -file   => $in);
while( my $result = $parser->next_result ) {
	while( my $hit = $result->next_hit ) {
		while( my $hsp = $hit->next_hsp ) {
			my $ref = $hsp->query_string;
			my $aln = $hsp->hit_string;
			my $new = cutGapCoords($ref,$aln);
			my $start = $hsp->start('query');
			my $length = length($new);
			#don't know why, but $hsp->start seems to add +1 to the position.
			substr ($fake_aln, $start-2, $length, $new);
		}
	}
}

print ">".$id."\n";
print uc $fake_aln."\n";

### SUBRUTINES ###

sub cutGapCoords {
	my($refseq,$alnseq) = @_;
	while ( $refseq =~ m/(-+)/g ) {
		my $end = pos($refseq);
		my $start = $end-length($1);
		my $length = length($1);
		substr($refseq, $start, $length, "");
		substr($alnseq, $start, $length, "");
	}
	return $alnseq;
}
