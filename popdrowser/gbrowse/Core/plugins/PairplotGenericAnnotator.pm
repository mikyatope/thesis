package Bio::Graphics::Browser2::Plugin::PairplotGenericAnnotator;
use Bio::Graphics::Browser2::Plugin;
use CGI ':standard';
@ISA = 'Bio::Graphics::Browser2::Plugin';

use warnings;

# called by gbrowse to return name of plugin for popup menu
sub name        { 'LD Plot [Testing feature]' }

# called by gbrowse to return the descriptive verb for popup menu
sub verb        { 'Display' }

# called by gbrowse to return description of plugin
sub description { 'This is an example plugin' }

# called by gbrowse to return type of plugin
sub type        { 'annotator' }

# called by gbrowse to configure default settings for plugin
sub config_defaults {
    my $self = shift;
    return {foo => $value1,
            bar => $value2}
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
  return $form;
}

# called by gbrowse to annotate the DNA, returning features
sub annotate {
   my $self     = shift;
   my $segment  = shift;
   my $config   = $self->configuration;
   
   my $seg_start = $segment->start;
   my $seg_stop  = $segment->stop;
   my $chr       = $segment->seq_id;
   
   #Feature list
   my $feature_list = $self->new_feature_list;
   my $feature_type = "";
   $feature_list->add_type($feature_type => {
                                             glyph =>'pairplot',
                                             angle => 45,
                                             bgcolor => 'blue',
                                             point => 1
                                            }); 


   # Create one big feature using the PairFeature
   # glyph (see end of synopsis for an implementation)
   my $block = PairFeature->new(-seq_id => $chr,
                                -start  => $seg_start,
                                -end    => $seg_stop,
                                -type   => $feature_type
                                );

   # It will contain a series of subfeatures.
   #Search database

   my $dbh = connect2db("vs_ld");

   my $markers = get_markers($dbh, $chr, $seg_start, $seg_stop);
#   my $values = get_pair_values($dbh, $chr, $seg_start, $seg_stop);

   while (my @row = $markers->fetchrow_array()){
       $block->add_SeqFeature(Bio::Graphics::Feature->new(
                                                          -start=>$row[0],
                                                          -end  =>$row[0],
                                                         ));
   }
    
#   while (my @row = $values->fetchrow_array()){ 
#       $block->{"D"} = {};
#       $block->{"D"}->{$row[0]}->{$row[1]} = $row[6];
#   }

   $feature_list->add_feature($block => $feature_type);
   $dbh->disconnect();
   
   return $feature_list;

}


#Database connection
sub connect2db {
	my ($database) = @_;
	my $dsn = "dbi:mysql:$database:158.109.215.162:3306";
	my $user = "nobody";
	my $pw = "";
	my $dbh = DBI->connect($dsn, $user, $pw, { RaiseError => 1, AutoCommit => 1}) 
				   || die "failed connection to database ($databse): $DBI::errstr";
	return $dbh;
}

#Get pairwise markers
sub get_markers {
	my ($dbh, $chr, $start, $end) = @_;
	
	my $query = "(SELECT DISTINCT pos1 FROM pairs WHERE pos1 >= ".$start." AND pos2 <= ".$end." AND chr = \"".$chr."\")
                 UNION
                 (SELECT DISTINCT pos2 FROM pairs WHERE pos1 >= ".$start." AND pos2 <= ".$end." AND chr = \"".$chr."\")";
	
	my $sth = $dbh->prepare($query);
	$sth->execute();  
	return $sth;
}





###########################################################
###########################################################
#
##PairFeature object
#package PairFeature;
#    #use base 'Bio::SeqFeature::Generic'; 
#    use base 'Bio::Graphics::Feature';
#    
#    sub pair_score {
#        my $self = shift;
#        my ($sf1,$sf2) = @_;
#        
#
#        return $score;
#    }



###########################################
###########################################


#PairFeature object
package PairFeature;
    #use base 'Bio::SeqFeature::Generic'; 
    use base 'Bio::Graphics::Feature';
    
    sub pair_score {
        my $self = shift;
        my ($sf1,$sf2) = @_;
        my $score = '';
        my $chr = $self->seq_id ;

       my $dbh_pair = connect2db_pairfeat("vs_ld");
       
       my $pos1 = $sf1->start;
       my $pos2 = $sf2->start;

       my $values = get_pair_values_pairfeat($dbh_pair, $chr , $pos1, $pos2);

	   while (my @row = $values->fetchrow_array()){ 

	       $score = $row[4] ; 
	       $score = $score*10 if $score <= 0.01;
	       $score = 1 if $row[4] > 0.01;
	   }

   		$dbh_pair->disconnect();

	return $score;
    }


#Database connection
sub connect2db_pairfeat {
	my ($database) = @_;
	my $dsn = "dbi:mysql:$database:158.109.215.162:3306";
	my $user = "nobody";
	my $pw = "";
	my $dbh = DBI->connect($dsn, $user, $pw, { RaiseError => 1, AutoCommit => 1}) 
				   || die "failed connection to database ($databse): $DBI::errstr";
	return $dbh;
}


#Get pairwise values
sub get_pair_values_pairfeat {
	my ($dbh, $chr, $pos1, $pos2) = @_;

	my $query = "SELECT d, dabs, dprime, dprimeabs, r2 ". 
			    "FROM pairs ".
			    "WHERE pos1 = ".$pos1.
			    " AND pos2 = ".$pos2.
		   	    " AND chr = \"".$chr."\"";

	my $sth = $dbh->prepare($query);
	$sth->execute();  
	return $sth;
}





