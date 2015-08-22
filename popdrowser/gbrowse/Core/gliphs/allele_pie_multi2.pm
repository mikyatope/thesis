package Bio::Graphics::Glyph::allele_pie_multi2;

# $Id: allele_logos.pm,v 1.4 2004/06/08 22:04:15 avs Exp $
# Glyph for drawing a pie chart with s slices representing allele frequency.

use strict;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';
use Bio::Graphics::Glyph::generic;

# Give enough height to fit in the pie chart, min height is 20
sub height {
    my $self = shift;
    my $height = $self->SUPER::height;
    if(defined $self->option('stacked')){
	$height = 19 if $height < 19;
	my $freq = defined($self->option('freq')) ? $self->option('freq') : 'NO';
	my @pops = split /;/, $freq;
	return $height * scalar (@pops);
    }else{
	return $height > 19 ? $height : 19;
    }
}


# Need to make room for the allele pies if there is room
sub pad_right {
  my $self = shift;
  my $right = $self->SUPER::pad_right;
  my $height = $self->height;
  my $freq = defined ($self->option('freq')) ? $self->option('freq') : "NO";
  my @pops = split /;/, $freq;
  if(defined $self->option('stacked')){
      return $right > $height + 6 ? $right : $height + 6 if $self->label;
  }else{
      return $right > $height * scalar(@pops) + 6  ? $right: $height * scalar(@pops) + 6  if $self->label;
  }
}

sub majorcolor {
    my $self = shift;
    my $color = $self->option('majorcolor') || '#0033FF';
    $self->factory->translate_color($color);
}

sub minorcolor {
    my $self = shift;
    my $color = $self->option('minorcolor') || '#FF0000';
    $self->factory->translate_color($color);
}

sub get_description{
    my $self = shift;
    my $feature = shift;
    my $freq = defined ($self->option('freq')) ? $self->option('freq') : undef;
    if($freq && $freq =~ /:/){
	return join(' ', map(/(.*):/, map {split /;/} $freq));
    }else{
	return join '; ', eval {$feature->notes};
    }
}

sub pad_bottom{
    my $self = shift;
    if (defined $self->option('stacked')){
	return 0;
    }else{
	return $self->SUPER::pad_bottom;
    }
}

sub draw_description {
    my $self = shift;
    my ($gd,$left,$top,$partno,$total_parts) = @_;
    my $label = $self->description or return;
    my @items = split /\s+/, $label;
    my $x = $self->left + $left;
    $x = $self->panel->left + 1 if $x <= $self->panel->left;
    my $pie_size;
    #if(defined $self->option('stacked')){
	my $freq = defined($self->option('freq')) ? $self->option('freq') : 'NO';
        my @pops = split /;/, $freq;
	$pie_size = ($self->height)/scalar(@pops);
    #}else{
	#$pie_size = $self->height;
    #}
    for(my $i = 0; $i < @items; $i++){

    	if($self->option('stacked')){
	        $gd->string($self->font,
			($x + $pie_size + 7),
			((($pie_size +1 )* $i) - 1 + $pie_size/2 + $self->font->height/2 + $self->top + $top),
			$items[$i],
			$self->font2color);

		}else{
			my $d = 7 + ($pie_size * (0.5 +  $i)) - (length($items[$i]) * ($self->font->width/2));
			$gd->string($self->font,
				($x + $d),
				$self->bottom - $self->pad_bottom + $top,
				$items[$i],
				$self->font2color);
		}

    }
}

sub draw_component {
    my $self = shift;
    my $gd = shift;

    # find the center and vertices
    my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);


    my $feature = $self->feature; 
    my @alleles = $feature->attributes('alleles');
    #@alleles    = split /\//,$self->option('alleles') unless @alleles == 2;
    @alleles    = split /\//,@alleles[0];


    if (@alleles) {
	# If it is on the minus strand
	if (my $strand = $self->option('ref_strand') <0){
	    foreach (@alleles) {
		tr/ACTG/TGAC/ if $self->option('complement');
	    }
	}
	my $height = $self->height;
	my $maf =  defined ($self->option('maf')) ? $self->option('maf') : "NO";
	my $freq = defined ($self->option('freq')) ? $self->option('freq') : "NO";
	my @pop_freqs = split /;/, $freq;

# to prevent the pie chart from flipping when the MAF is 0.50;
	$freq = 0.499 if $freq == 0.5;
        # write the alleles
        # Changed for DGRP
    my $align = $self->option('align') || 1;     
    $gd->string(GD::Font->MediumBold, $x1 - $align, -2 + $y1, @alleles[0], $self->majorcolor);
    $gd->string(GD::Font->Small, $x1 -1, -2 + 10 + $y1, @alleles[1], $self->minorcolor);

	
	
	
	# draw the pie charts
	if($self->label){
	    for(my $i = 0  ; $i < scalar(@pop_freqs); $i++){
		my $freq;
		$freq = $pop_freqs[$i] unless $pop_freqs[$i] =~ /:/;
		if($pop_freqs[$i] =~ /:/){
		    my @s = split /:/, $pop_freqs[$i];
		    $freq = $s[1];
		}
		$freq = 0.499 if $freq == 0.5;
		$freq = 'NO' unless $freq =~ /^[0-9.]+$/;
		my $ph = $height -1 ;
		my $ph = (defined $self->option('stacked')) ? $height/scalar(@pop_freqs)  : $height;
		my $xd = (defined $self->option('stacked')) ? 0 : ($height  * $i);
		my $yd = (defined $self->option('stacked')) ? (($ph + 1)  * $i) : 0;
		
		$gd->arc((($x1 + $ph/2)+6 + $xd) ,($y1+$ph/2 + $yd),($ph - 1),($ph - 1), 
			     0, 360, $self->majorcolor) if $freq eq 'NO';
		$gd->filledArc((($x1 + $ph/2)+6 + $xd) ,($y1+$ph/2 + $yd),($ph - 1),($ph - 1),
			       270, (360 * (1- $freq)) + 270, $self->minorcolor) unless $freq eq 'NO';
		$gd->filledArc((($x1 + $ph/2)+6 + $xd),($y1+$ph/2 +$yd  ),($ph - 1),($ph - 1),
			       (360*(1-$freq))+270, 270 +360, $self->majorcolor) unless $freq eq 'NO';
	    }
	}
    }
}


1

;

__END__

=head1 NAME

Bio::Graphics::Glyph::allele_tower - The "allele_tower" glyph

=head1 SYNOPSIS

  See <Bio::Graphics::Panel> and <Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a letter for each allele found at a SNP position, one above the other (i.e. in a column). For example:
    A
    G

See also http://www.hapmap.org/cgi-perl/gbrowse/gbrowse 'genotyped SNPs' for an example.

The common options are available (except height which is calculated
based on the number of alleles).  In addition, if you give the glyph
the minor allele frequency (MAF) and indicate which is the minor
allele, the glyph will display these differences.


=head2 GETTING THE ALLELES

To specify the alleles, create an "Alleles" attribute for the feature.
There should be two such attributes.  For example, for a T/G
polymorphism, the GFF load file should look like:

 Chr3  .  SNP   12345 12345 . . . SNP ABC123; Alleles T ; Alleles G

Alternatively, you can pass an "alleles" callback to the appropriate
section of the config file.  This option should return the two alleles
in any way it wishes.

  alleles = sub {
	my $snp = shift;
	my $d   = $snp->attributes('AllelePair');
	return split m!/!,$d;
    }

=head2 OPTIONS

 . Glyph Colour
 . Different colour for alleles on the reverse strand
 . Print out the complement for alleles on the reverse strand
 . Major allele shown in bold
 . Horizontal histogram to show allele frequency

=head3 GLYPH COLOR

The glyph color can be configured to be different if the feature is on the plus or minus strand.  Use fgcolor to define the glyph color for the plus strand and bgcolor for the minus strand.  For example:

   fgcolor     = blue
   bgcolor     = red

For this option to work, you must also set ref_strand to return the strand of the feature:
   ref_strand        = sub {shift->strand}

=head3 REVERSE STRAND ALLELES

If the alleles on the negative strand need to be the complement of what is listed in the GFF files, (e.g. A/G becomes T/C), set the complement option to have value 1

complement   = 1

For this option to work, you must also set ref_strand to return the strand of the feature:

ref_strand        = sub {shift->strand}

=head3 MAJOR/MINOR ALLELE

Use the 'minor_allele' option to return the minor allele for the SNP.  If you use this option, the major allele will appear in bold type.

=head3 ALLELE FREQUENCY HISTOGRAMS

Use the 'maf' option to return the minor allele frequency for the SNP.  If you use this option, a horizontal histogram will be drawn next to the alleles, to indicate their relative frequencies. e.g.

 A______
 C__

Note: The 'label' option must be set to 1 (i.e. on) for this to work.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::cds>,
L<Bio::Graphics::Glyph::crossbox>,
L<Bio::Graphics::Glyph::diamond>,
L<Bio::Graphics::Glyph::dna>,
L<Bio::Graphics::Glyph::dot>,
L<Bio::Graphics::Glyph::ellipse>,
L<Bio::Graphics::Glyph::extending_arrow>,
L<Bio::Graphics::Glyph::generic>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::heterogeneous_segments>,
L<Bio::Graphics::Glyph::line>,
L<Bio::Graphics::Glyph::pinsertion>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::rndrect>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::ruler_arrow>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::translation>,
L<Bio::Graphics::Glyph::allele_tower>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Albert Vernon Smith E<lt>albert.smith@cshl.eduE<gt> in Lincoln Stein's lab E<lt>steinl@cshl.eduE<gt>.

Copyright (c) 2003 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
