package Bio::Graphics::Browser2::Plugin::PairplotDGVBAnnotator;
# $Id: PairplotAnnotator.pm,v 1.16 2008/08/12 18:50:51 wei Exp $

use strict;
use Bio::Graphics::Browser2::Plugin;
use CGI::Pretty qw(:standard *table);
use Bio::DB::GFF;
use Data::Dumper;
use POSIX;
#use test;
use vars '$VERSION','@ISA', '$DEBUG';
$DEBUG = 0;
$VERSION = '0.21';

@ISA = qw(Bio::Graphics::Browser2::Plugin);

my @markerlist ;
my @COLORS = qw(red grey green blue orange cyan black 
		turquoise brown indigo wheat yellow emerald);
my @LIGHTCOLORS = qw(lightgrey ivory palegreen peachpuff antiquewhite white yellow);  

sub name { "LD Plot" }

sub description {
    p("Linkage Disequilibrium (LD) is the non random association between alleles in the population",
      "The LD plot plugin generates a pairwise plot of marker-to-marker LD values,on the current view.",
      " The ticks denote the markers with respect to their physical location on the chromosome segment.",
      " The boxes represent the marker pair relationship and it is plotted <i>between</i> the two markers.",
      " The color of the box (or the  intensity of the color) is based on the raw score for that marker pair.", 
      "Three LD properties can be viewed as separate tracks more  detailed explaination about the tracks and how to configure them are available <a href='/gbrowse_help.html#LD'> here </a>. " ).
    p("This plugin is based on the plugin written by Lalitha Krishnan &amp; Lincoln Stein in <a href=\"http://hapmap.ncbi.nlm.nih.gov/\" target=\"_blank\">HapMap</a>."); 
}

sub type { 'annotator' }

sub config_defaults {
    my $self = shift;
    return{dprime_plot_color => 'red',
	   md_color              => 'lightgrey',
	   rsq_plot_color        => 'blue',
	   rsq_md_color          => 'lightgrey',
	   box_size              => 'Proportionate',
	   seg_len               => '250Kb',
	   pop_code              => 1,
	   orient                => 'normal',
	   ld_property           => 'lod',
	   no_of_SNP             => '200',
	   match_MAF             => '0.05',
	   dprime_cutoff         => '0.3',
	   rsquare_cutoff        => '1.0',
	   Dprime                => 1,
	   Rsquare               => 0,
       RAL_orient            => 'invert'
	   };
}

sub reconfigure {

    my $self = shift;
    my $current_config = $self->configuration;
    foreach my $param ( $self->config_param() ) {
		$current_config->{$param} = $self->config_param($param);
    }

}

sub configure_form {

    my $self = shift;
    my $current_config = $self->configuration;

    my @rows;
    push @rows,TR(td(['&nbsp;','&nbsp;']));
    push @rows,
    TR( 
	td({-align=>'LEFT',-colspan=>'1'},b("Max. Segment Size&nbsp;&nbsp;&nbsp;<br>"),
	   popup_menu(-name=>$self->config_name("seg_len"),
		      -values=> [qw( 100Kb 200Kb 250Kb 300Kb 500Kb 750Kb 1Mb)],
		      -default=> $current_config->{"seg_len"})),

	td({-align=>'LEFT',-colspan=>'1'},b("Max. # gt'd SNPs&nbsp;&nbsp;&nbsp;<br>"),
	   popup_menu(-name=>$self->config_name("no_of_SNP"),
		      -values=> [qw(200 250 300 350 400 450 500 550 600 650 700 750 850 1000)],
		      -default=> $current_config->{"no_of_SNP"})),
	
	td({-align=>'LEFT',-colspan=>'1'},b("Box Size&nbsp;&nbsp;&nbsp;<br>"),
	   popup_menu(-name=>$self->config_name("box_size"),
		      -default=>$current_config->{"box_size"},
		      -values=>[qw(Proportionate Uniform )],))	
	);
    push @rows,TR(td(['&nbsp;','&nbsp;']));
    push @rows,
    TR(  td({-align=>'LEFT',-class=>"searchtitle"},b("LD Properties:")),
	 td({-align=>'LEFT',-rowspan=>'1'},
	    popup_menu(-name=>$self->config_name("ld_property"),
		       -values=> [qw(dprime rsquare lod)],
		       -default=> $current_config->{"ld_property"})),
	 td({-align=>'LEFT',-rowspan=>'1',-class=>"searchtitle"}," greater than<br>",
	    popup_menu(-name=>$self->config_name("dprime_cutoff"),
		       -values=> [qw(1.0 0.9 0.8 0.75 0.7 0.5 0.4 0.3 0.2 0.1 0.0)],
		       -default=> $current_config->{"dprime_cutoff"}
		       )),
	 td({-align=>'LEFT',-rowspan=>'1',-class=>"searchtitle"},"and less than<br>",
	    popup_menu(-name=>$self->config_name("rsquare_cutoff"),
		       -values=> [qw(1.0 0.9 0.8 0.75 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0)],
		       -default=> $current_config->{"rsquare_cutoff"}
		       )),
         TR(td(['&nbsp;','&nbsp;'])),
	 TR(
	    td({-align=>'LEFT',-rowspan=>'1',-class=>"searchtitle"},b("Color:")),
	    td({-align=>'LEFT'}, "Pairwise plot<br>",
	       popup_menu(-name=>$self->config_name("dprime_plot_color"),
			  -values=> \@COLORS,
			  -default=> $current_config->{"dprime_plot_color"}
			  )) 
	    )
	 );
    
    push @rows,TR(td(['&nbsp;','&nbsp;']));  
    
    push @rows,
    TR(
       td({-align=>'LEFT',-class=>"searchtitle"},b("Orientation:")),
       td({-align=>'LEFT'},
	  popup_menu(-name=>$self->config_name('RAL_orient'),
		     -values  =>[qw(normal invert)],
		     -default => $current_config->{"RAL_orient"}
		     ))
       );
    
    return table({#-width=>'100%',
		  -border=>0},
		 @rows); 

}

sub annotate {

    my $self = shift;
    my $segment = shift;
    my $config  = $self->configuration;
    my $browser_config = $self->browser_config;
    
##    my $dsource = $browser_config->source; 
    
    # initialising values for the top level segment
    my $ref        = $segment->ref;
    my $abs_start  = $segment->start;
    my $abs_stop   = $segment->stop ;

    my $ld_property= $config->{"ld_property"};
    my $Dprime = $config->{"Dprime"};
    my $Rsquare = $config->{"Rsquare"};
    my $dprime_cutoff =$config->{"dprime_cutoff"};
    my $rsq_cutoff = $config->{"rsquare_cutoff"}; 
    $dprime_cutoff ||= 0.3 ; # force default to 0 if not defined
    $rsq_cutoff ||= 1; # force default to 0 if not defined
    my $Haplo = 1; # has to set it up later;
    my $seg_length = $config->{"seg_len"};
    my $box= $config->{"box_size"};
    my $num_of_snp =$config->{"no_of_SNP"};
    my $last_char = substr($seg_length,-2,2,"");
    my $box_size =1;
    my ($start_rs,$stop_rs)="";
    my @pop_list=( "RAL" );
    $box_size =0 if ($box eq "Proportionate");
    $seg_length = $seg_length * 1000 if ($last_char eq "Kb");
    $seg_length = $seg_length * 1000000 if ($last_char eq "Mb");


    # creating the feature file   
    my $feature_list = Bio::Graphics::FeatureFile->new;
#
#    # This is the top level feature
#    my $block = PairFeature->new(-start=>$abs_start, -end=>$abs_stop);
#
#    foreach my $pop(@pop_list){
#		next unless($config->{$pop});
#		my $orient =$pop."_orient" ;  
#		$block->{$pop}= "N" if ($config->{$orient} eq "normal");
#		$block->{$pop}="I" if ($config->{$orient} eq "invert");
#    }
#    
#    $block->{ld_property} = $config->{"ld_property"};
#    $block->{plot_color} =$config->{"dprime_plot_color"};
#    $block->{dprime_cutoff} = $dprime_cutoff;
#    $block->{rsquare_cutoff} = $rsq_cutoff;    
#    $block->{num_of_snp} =$num_of_snp;
#    $block->{md_color} =$config->{"md_color"};
##    $block->{source} =$dsource;
#    my $coordinate={};  
#
#    # check for the length of the segment currently displayed. If more than seg length cut off do not plot the pair plot.
#    if ($abs_stop - $abs_start > $seg_length){
#		$feature_list = pretty_plot($block,$feature_list,'box',$box_size);
#		return $feature_list;
#    }
#
#    my @snps = $segment->features(-types=>['remark:snp']);
#    my $totalsnps = @snps;
    
#    print "\ntotalsnps".$totalsnps."\n";

  
#    return if ($totalsnps == 0);
#    # check number of snps to view against total snps in the segment. If more trim from both ends of the segment.
#    my $ctr = 0;
#    my $dropctr=1;
#    if ($totalsnps >= $num_of_snp){
#		my $trimsnps =($totalsnps - $num_of_snp);
#		$trimsnps = 1 if ($trimsnps == 0);
#		$dropctr =int(($totalsnps / $trimsnps)+0.5);
#		# if SNP cut off is less than half the total snp present do not plot the pairplot
#		if ($dropctr <= 1){
#		    $feature_list=    pretty_plot($block,$feature_list,'box',$box_size) ;
#		    return $feature_list;
#		}
#    }
#    
#    my %seen = ();
#    my @sortsnp =sort {$a->start <=> $b->start} @snps;
#    my @snp_subfeats;
#    $stop_rs = $sortsnp[$totalsnps-1]->name ;
#    $start_rs =$sortsnp[0]->name ;
#
#    for (my $snp=0;$snp<@sortsnp;$snp++) {
#	    #add snp as subfeature to big block feature after checking for maf cutoff and duplicate snps and snp cut off
#		if ($dropctr >=2 and $ctr <= $num_of_snp) {
#		    $ctr++;
#		    next if ($snp % $dropctr == 0);
#		} 
#	
#		next if $seen{$sortsnp[$snp]}++;
#        my $maf = find_maf($sortsnp[$snp]);
#        next if ($maf <= 0.05 ); 
#        my $name= $sortsnp[$snp]->name;
#        my $pos =$sortsnp[$snp]->stop;
#
#        # warn "Pairplot SNP:$name,MAF:$maf \n"; 
#        push @snp_subfeats, new Bio::Graphics::Feature(-start => $sortsnp[$snp]->start,
#					               -stop  => $sortsnp[$snp]->end,
#					               -type  => 'snp',
#					               -display_name => $sortsnp[$snp]->name);
#
#        $coordinate->{$name}=$pos;
#        $start_rs = $sortsnp[$snp]->name if ($snp == 0);
#  
#    } # End for
#
#    $block->add_segment(@snp_subfeats);
#    $block->{coordinate} = $coordinate;
#    my $ds = $block->source;
#    my ($dsuffix) = $ds =~ /hapmap_(.+)/;
#
#    my $no_of_markerpair = get_markerpair_values($block,$segment->ref, $abs_start, $abs_stop,$ld_property,$Rsquare,$Haplo);
#
#    $feature_list = pretty_plot($block,$feature_list,'pairplot',$box_size) ;
#
#    return  $feature_list  if ($no_of_markerpair > 0);


}


sub find_maf {

#    my $snp= shift;
#    my @acounts =$snp->attributes('acounts') ;
#    my $ac = $acounts[0];
#    $ac =~ /^(\w+):/; 
#    $ac=~ s/$ac://;
#    my @m = split ' ', $ac;
#    my $maf = $m[1] <$m[4]?$m[1]:$m[4];
#    return $maf ;
   
}


sub get_markerpair_values
{

#    my ($block,$ref,$start,$stop,$ld_property,$rsquare_fl,$haploview_fl) = @_ ;
#    my ($counter, $dp_median_score, $rsq_median_score)=0; 
#    my $tbl_name = lc($ref) . "_ldscore";
#    $ref =lc($ref);
#    my $dprime_cutoff = $block->{"dprime_cutoff"};
#    my $rsquare_cutoff = $block->{"rsquare_cutoff"};
#    my $coordinate =$block->{"coordinate"}; 
#    $dprime_cutoff = 0 if ($ld_property eq "lod");     
#    $rsquare_cutoff = 100 if ($ld_property eq "lod");      
#    my @seg = $block->segments;
#
#
#    my $dsource = $block->source;
#    my ($dbsuffix) = $dsource =~ /hapmap_(.+)/;
#    $dbsuffix ||= 'current'; #default to current
#
###############!!!!!!!!!!!!!!!!!!!!!!!!
#    my $dsn;
#    $dsn = "DBI:mysql:database=hapmap_ldplot_phase3;host=wildcat;user=krishnan;password=krishnan" if( $dsource eq "hapmap3_B36" );
#    $dsn = "DBI:mysql:database=hapmap_ldplot_phase3_r2;host=wildcat;user=krishnan;password=krishnan" if( $dsource eq "hapmap3r2_B36" );
#    $dsn = "DBI:mysql:database=hapmap_ldplot_r27;host=wildcat;user=krishnan;password=krishnan" if( $dsource eq "hapmap27_B36" );
#
#
#    my $dbh = DBI->connect($dsn) or die "can't connect to $dsn ";
#
##############!!!!!!!!!!!!!!!!!!!!!!!!!
#    my $query = "SELECT fstart,fstop,pop_code,marker1,marker2,dprime,$ld_property,rsquare FROM $tbl_name WHERE fstart >= $start and fstart <= $stop and $ld_property >= $dprime_cutoff and $ld_property <=$rsquare_cutoff  and  pop_code =?";
#
#
#    my $sth = $dbh->prepare($query) || die $dbh->errstr;
#
#    my @pop_list=("RAL");
#    
#    #pairwise _lookuphash will have the chromosomal reference for each snp pair(sub feature) Dprime lookup table is needed for haplo plotas well as dprime plot hence will always have the dataset ;
#    #loading subfeatures in the lookuphash
#    foreach my $pop(@pop_list) {
#        next unless ($block->{$pop});
#        my $pop_code =lc($pop) if ( $block->{$pop});
#        my $hashname="$pop-lookuphash";
#        my $lod_hashname="$pop-dprime-lookuphash"; 
#
#        $block->{$hashname} = {};
#        $block->{$lod_hashname} = {} if ($ld_property eq "lod");
#        $sth->bind_param(1,$pop_code);
#        my $start_time = (times)[0]; 
#        $sth->execute;
#        while (my ($pos1, $pos2,$p_code,$snp1,$snp2, $dprime_score,$score,$r2,$fbin) = $sth->fetchrow_array()) {
#            $counter++;
#            $dp_median_score +=  $dprime_score;
#            $rsq_median_score += $score;
#            $block->{$lod_hashname}->{$pos1}->{$pos2} = $dprime_score if ($ld_property eq"lod")  ;
#            $block->{$hashname}->{$pos1}->{$pos2} = $score  ;
#		}
#        my $end_time =(times)[0];
#    }
#
#    $sth->finish;
#    $dbh->disconnect;
#    return $counter;

}


sub get_b34 {

#    my $abs_start=shift;
#    my $abs_stop = shift;
#    my $ref = shift;
#    my $coordinate ;
#    my $db = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt', 
#                                -dsn => 'dbi:mysql:database=gbrowse_hapmap_current',
#                                -host=> 'osceola',
#                                -user => 'nobody',
#                                -password => 'nobody') or die "I do not have a database";
#
#    my $seg = $db->segment($ref,$abs_start=>$abs_stop);
#    my @snps = $seg->features(-types=>['snp:HapMap_gt','snp:hapmap_gt']);
#    foreach my $snp (@snps)
#    {
#        my $name=$snp->name;
#        my $stop =$snp->stop;
#        $coordinate->{$name}=$stop;
#    }
#
#    return $coordinate;

}


sub find_b34position {

#    my $rs = shift;
#    my $ref = shift;
#    my $dbsuffix = shift;
#    $ref = lc($ref);
#    my $db;
#    $db = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt',
#                             -dsn => 'dbi:mysql:database=gbrowse_hapmap20_B35_v3:osceola.cshl.edu',
#                             -user => 'nobody',
#                             -password => 'nobody'
#                           ) or die "I do not have a database" if ( $dbsuffix eq "B35");
#
#    $db = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt',
#                             -dsn => 'dbi:mysql:database=gbrowse_hapmap23_B36:kinsman.cshl.org',
#                             -user => 'nobody',
#                             -password => 'nobody'
#                           ) or die "I do not have a database" if ( $dbsuffix eq "B36");
#
#    $db = Bio::DB::GFF->new( -adaptor => 'dbi::mysqlopt',
#                             -dsn => 'dbi:mysql:database=gbrowse_hapmap_current',
#                             -host=> 'osceola', 
#                             -user => 'nobody',
#                             -password => 'nobody'
#                           ) or die "I do not have a database" if ($dbsuffix eq "B34");
#
#
#    my @features = $db->fetch_feature_by_name('SNP'=> $rs) or return 0;
#    my @true_feats;
#    map{$_->abs_start== $_->abs_stop and push(@true_feats,$_);} @features;
#    my $dbsnp_feat ;
#
#    my $start= $true_feats[0]->abs_start;
#
#    return $start;

}


sub pretty_plot {

#    my ($block,$feature_list,$glyph_type,$box_size)=@_ ;
#    my $ld_property =$block->{"ld_property"};
#    my $ppcolor =$block->{"plot_color"};
#    my $mdcolor = $block->{"md_color"};
#    my $median_score = 0;
#    my $num_of_snp=$block->{"num_of_snp"};
#
#    my @pop_list=("RAL");
#
#    foreach my $pop(@pop_list){
#        my $glyph_orient = "N";
#        next unless($block->{$pop});
#        my $f_type =$pop  if($block->{$pop}) ;
#        my $track_name = "$pop:$ld_property";
#        my $feature_type= "$pop-$ld_property";
#        $glyph_orient = $block->{$pop};
#
#        if ($glyph_type eq "box"){
#            $f_type .= ":too many features..choose 'Annotate LD Plot' under 'Reports & Analysis' and press 'Configure' to adjust the segment size";
#            $feature_list->add_type( $feature_type=>{ glyph => 'box',
#					              key => $f_type,
#					              bgcolor => $ppcolor,
#					              bump => "0",
#					              citation=>'Linkage Disequilibrium' } );
#        }
#        if ($glyph_type eq 'pairplot'){
#            $feature_list->add_type( $feature_type=>{ glyph =>'myreverseplot',
#										              angle => 45,
#					                                  point => 1,
#										              bump =>0,
#										              glyph_orient=>$glyph_orient,
#										              key =>$track_name,
#					                                  maxdepth => 1,
#					                                  numsnp=>$num_of_snp,
#										              bgcolor =>$ppcolor,
#					                                  mdcolor=>'lightgrey',
#										              f_type =>$f_type,
#										              ld_property=>$ld_property,
#										              box_size=>$box_size,
#                                                      citation=>'Linkage Disequilibrium',
#                                                      med_score=>$median_score } );
#        }
#        $feature_list->add_feature($block => $feature_type);
#    }
#    return $feature_list;

}



package PairFeature ;
use Data::Dumper;
use base 'Bio::Graphics::Feature';
sub pair_score {
    my $self = shift;
    my ($sf1,$sf2,$ld_fl) = @_;
    my $score = "0.00";
    #actual pairwise lookup
    my $lookup_table = $ld_fl . "-lookuphash";
    my $e1 = $sf1->start;
    my $e2 = $sf2->start;
    
    if ($self->{$lookup_table}->{$e1}->{$e2}){ 
        $score = $self->{$lookup_table}->{$e1}->{$e2}; 
    }else{ 
        $score = "99.99";
    }
    return $score ;
}







