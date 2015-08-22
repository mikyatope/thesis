package Bio::Graphics::Browser2::Plugin::PolySeqDump2;
use Bio::Graphics::Browser2::Plugin;
use CGI ':standard';
@ISA = 'Bio::Graphics::Browser2::Plugin';

use warnings;
use strict;

# called by gbrowse to return name of plugin for popup menu
sub name        { 'Polymorphic DGRP sequences' }

# called by gbrowse to return the descriptive verb for popup menu
sub verb        { 'Download' }

# called by gbrowse to return description of plugin
sub description { 'This plugin downloads the sequence for the 158 lines of the DGRP project along with the D.simulans outgrup' }

# called by gbrowse to return type of plugin
sub type        { 'dumper' }

# called by gbrowse to configure default settings for plugin
sub config_defaults {
    my $self = shift;
    return {foo => '',
            bar => ''}
}

# called by gbrowse to reconfigure plugin settings based on CGI parameters
sub reconfigure {
  my $self = shift;
  my $current = $self->configuration;
  #$current->{foo} = $self->config_param('foo');
  #$current->{bar} = $self->config_param('bar');
}

# called by gbrowse to create a <form> fragment for changing settings
sub configure_form {
  my $self    = shift;
  my $current = $self->configuration;
  #my $form = textfield(-name  => $self->config_name('foo'),
   #                    -value => $current->{foo})
    #         .
     #        textfield(-name  => $self->config_name('bar'),
      #                 -value => $current->{bar});
  return "";
}

sub mime_type {
  my $self = shift;
  my $config = $self->configuration;
  return ('application/octet-stream','sequence_region.fasta');
}

# called by gbrowse to annotate the DNA, returning features
sub dump {
	my $time_start = time;
	
	my $self     = shift;
	my $segment  = shift;
	my $config   = $self->configuration;
	my $seg_start = $segment->start;
	my $seg_stop  = $segment->stop;
	my $chr       = $segment->seq_id;
	my $file     = "/home/share/dgrp.freeze1/".$chr.".alloutgroups.nodupl.fasta";

    my $aln = readfasta($file,$seg_start,$seg_stop);
    print "###Region ".$chr.":".$seg_start."..".$seg_stop."\n";
    print "###Size: ".(($seg_stop-$seg_start)+1)."\n";
    print $aln;
    my $time_end = time - $time_start;
    print "\n\n### Time: $time_end s\n";
}

sub readfasta($$$){
	my ($file,$start,$end)=@_;
	my $aln;
	(!open(FASTA,$file)) && do{return 0;};	
	while(<FASTA>){
		chomp;	
		if ($_=~/^>(\S+)/){
			$aln .= $_."\n";
		}else{
			$aln .= substr($_, $start, ($end-$start)+1);
			$aln .= "\n";
		}
	}
	close(FASTA);
	return $aln;
}

