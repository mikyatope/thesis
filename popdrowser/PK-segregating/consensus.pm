#A subroutine to open the file- given filename, set filehandle
sub open_file {

	my ($filename) = @_;
	my $fh;

	unless ( open( $fh, $filename ) ) {
		print "Cannot open file $filename\n";
		exit;
	}
	return $fh;
}

#A subroutine to get next record
sub get_next_record {

	my ($fh) = @_;

	my ($offset);
	my ($fasta)                = '';
	my ($save_input_separator) = $/;

	$/ = "//\n";

	$fasta = <$fh>;

	$/ = $save_input_separator;

	return $fasta;
}

# Fasta to Raw format and get "greedy" consensus

sub fasta2raw {
	my ($fasta) = @_;
	my $raw;

	$raw = '';

	#while ($fasta =~ /^>([A-Z]+seq\d+)[^\n]*\n([^>]+)/gms)
	while ( $fasta =~ /^([A-Z]+tr\d+)[^\n]*\n([^FB]+)/gms )

	  #while ($fasta =~ /^(FBtr\d+_\d+_.*\t)[^\n]*\n([^FB]+)/gms)
	{
		$poly_seq = uc($2);
		$poly_seq =~ s/[^A-Z-]//g;
		$raw .= $poly_seq . "\n";
	}

	return ($raw);
}

sub greedy {
	my ( $self, $threshold, $mincol ) = @_;
	$threshold = 0.75 unless defined($threshold);
	$mincol    = 1    unless defined($mincol);
	bless $self;

	my ($row);
	my ( @row_start, @row_stop );
	my ( $start,     $end );

	# measure the bounds for each row

	foreach $row ( split( /\n/, $self->layout('raw') ) ) {

		#    ($start) = $row =~ /^(\.*)/g;
		#    ($end) = $row =~ /^(.*[^.])\.*$/;
		($start) = $row =~ /^(-*)/g;
		($end)   = $row =~ /^(.*[^-])-*$/;
		push @row_start, 1 + length $start;
		push @row_stop,  length $end;
	}

	my $consfunc = &Bio::UnivAln::_c_consensus_of_array($threshold);

	my @consensus;
	my $width  = $self->width;
	my $height = $self->height;
	my @inside;
	my $seq;
	foreach my $ncol ( 1 .. $width ) {
		@inside = ();    #which rows are inside their bounds
		foreach my $i ( 0 .. $height - 1 ) {
			push( @inside, $i + 1 )
			  if ( $ncol >= $row_start[$i] )
			  && ( $ncol <= $row_stop[$i] );
		}
		if ( $#inside + 1 < $mincol ) {
			push @consensus, '#';    #my marker for insufficent lines
		}
		else {
			$seq = $self->seqs( [@inside], [$ncol] );
			push @consensus, &$consfunc( [ split( /\n/, $seq ) ] );
		}
	}
	return wantarray ? @consensus : join( "", @consensus );
}

#A subroutine to chage the format from non-format to fasta-format
#Para sec consensus (modificado): $formateado = $formateado . "\n" . substr($sequence, $pos, $length). "\n";
#Para count (como est): $formateado = $formateado . substr($sequence, $pos, $length). "\n";
sub format_fasta {
	my ( $sequence, $length ) = @_;

	my $formateado = "";
	for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
		$formateado = $formateado . substr( $sequence, $pos, $length ) . "\n";
	}
	return $formateado;
}

sub get_file_data {
	my ($filename) = @_;
	use strict;
	use warnings;

	#Inizialize variables
	my @filedata = ();
	unless ( open( GET_FILE_DATA, $filename ) ) {
		print STDERR "Cannot open file \"$filename\"\n\n";
		exit;
	}
	@filedata = <GET_FILE_DATA>;
	close GET_FILE_DATA;
	return @filedata;
}

#A Subroutine to counts
sub count {
	($seq) = @_;
	my ($count) = 0;

	$count = ( $seq =~ tr/[AGTCNagtcn]// );
	return $count;
}

1
