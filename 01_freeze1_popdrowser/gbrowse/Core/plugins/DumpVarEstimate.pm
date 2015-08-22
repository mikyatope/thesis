 package Bio::Graphics::Browser2::Plugin::DumpVarEstimate;
 use Bio::Graphics::Browser2::Plugin;
 use CGI ':standard';
 @ISA = 'Bio::Graphics::Browser2::Plugin';

 # called by gbrowse to return name of plugin for popup menu
 sub name        { 'Variation Estimates Download' }

 # called by gbrowse to return the descriptive verb for popup menu
 sub verb        { 'On-the-fly' }

 # called by gbrowse to return description of plugin
 sub description { 
	p('ON-THE- FLY VARIATION ESTIMATES DOWNLOADER');
	p('This plugin allows to recompute nucleotide variation estimates and lists the results in text format in a new window.');
 }

 # called by gbrowse to return type of plugin
 sub type        { 'dumper' }

 # called by gbrowse to configure default settings for plugin
	 sub config_defaults {
		my $self = shift;
		return {
			usemuts          => 1,
			completedeletion => 0,
			fixnum           => 0,
			numnuc           => 150,
			windowtype       => 3,
			ldsinglets       => 1,
			outgroup         => "no_out",
		}
		
	 }

 # called by gbrowse to reconfigure plugin settings based on CGI parameters
	 sub reconfigure {
		my $self = shift;
		my $current_config = $self->configuration;

		$current_config->{ws} = $self->config_param('ws');
		$current_config->{jump} = $self->config_param('jump');
		$current_config->{statistic} = $self->config_param('statistic');
		$current_config->{windowtype} = $self->config_param('windowtype');
		
		$current_config->{outgroup} = $self->config_param('outgroup');
		#$current_config->{startpos} = $self->config_param('startpos');
		#$current_config->{endpos} = $self->config_param('endpos');
		
		$current_config->{usemuts} = $self->config_param('usemuts');
		$current_config->{completedeletion} = $self->config_param('completedeletion');
		$current_config->{fixnum} = $self->config_param('fixnum');
		$current_config->{numnuc} = $self->config_param('numnuc');
		
		$current_config->{ldsinglets} = $self->config_param('ldsinglets');

	 }

 # called by gbrowse to create a <form> fragment for changing settings
	 sub configure_form {
		my $self     = shift;
		my $current  = $self->configuration;
		my $defaults = $self->config_defaults;
	
		my @order= ('s','s_inter','eta','eta_e','pi','theta','tajd','fulidast','fulifast','k','fulid','fulif','faywuh','ldsites','d','dabs','dprime','dprimeabs','r2','h','hd','fufs');
		my %statistic_labels= ( 's'          => 'Segregating Sites - S',
								's_inter'    => 'Number of sites differing between ingroup and outgroup - S_inter',
								'eta'        => 'Number of mutations - Eta',
								'eta_e'      => 'Number of singletons - Eta_E',
								'pi'         => 'Nucleotide diversity - Pi',
								'theta'      => 'Watterson\'s nucleotide diversity per site - Theta',
								'tajd'       => 'Tajima\'s D',
								'fulidast'   => 'Fu & Li\'s D*',
								'fulifast'   => 'Fu & Li\'s F*',
								'k'          => 'Divergence per site - K',
								'fulid'      => 'Fu & Li\'s D',
								'fulif'      => 'Fu & Li\'s F',
								'faywuh'     => 'Fay & Wu\'s H',
								'ldsites'    => 'Number of polymorphic sites used for LD test - LD_sites',
								'd'          => 'D value',
								'dabs'       => 'Absolute D value - |D|',
								'dprime'     => 'D\' value',
								'dprimeabs'  => 'Absolute D\' value - |D\'|',
								'r2'         => 'R square',
								'h'          => 'Number of haplotypes - h',
								'hd'         => 'Haplotype diversity - Hd',
								'fufs'       => 'Fu\'s Fs'
								);

		my $form = h2("Options:");
		
		#Statistic and Windows Size
		$form .= "<fieldset><legend>Statistic & Windows size</legend>";
		$form   .= p(strong("Statistic:"));
		$form   .= popup_menu( -name     => $self->config_name('statistic'), 
							   -values   => \@order,
							   -labels   => \%statistic_labels,
							   -onChange => 'changeDefaults(this)',
							   -default  => $current->{statistic}
								);
								
		$form   .= p(strong("Sliding-window size:").
		             " <br />(Default: Unique window with current region size)(Max = Current region size or maximum of 200kb)");
		$form   .= textfield(-name  => $self->config_name('ws'),
							 -value => $current->{ws}
							)."bp   ".em(" &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Leave blank for default ");
		
		$form   .= p(strong("Jump size:") . 
		             " <br />(Only if Sliding-Window size < Current region) (Default/Max = Sliding-Window size)");
		$form   .= textfield(-name  => $self->config_name('jump'),
							 -value => $current->{jump}
							)."bp   ".em(" &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp Leave blank for default ");
		
		
		$form .= p("<br /><strong>Note: </strong>Maximum region range analyzed is 1Mb.");
		
		$form .= "</fieldset>";
		

		#Sequences
		$form .= "<fieldset><legend>Sequences</legend>";
		
		$form   .= p(strong("Outgroup:"));
		$form   .= p("To estimate '<strong>K</strong>', '<strong>S inter</strong>', '<strong>Fu & Li D</strong>', '<strong>Fu & Li F</strong>' and '<strong>Fay& Wu H</strong>' it is mandatory to define an outgroup species. PopDrowser will use <em>D. simulans</em> as default.");
		
		my %outgroup_labels= ( 'no_out' => 'No outgroup',
						        'dsim'   => 'D. simulans',
						        'dyak'   => 'D. yakuba'
						      );
		$form   .= popup_menu( -name     => $self->config_name('outgroup'), 
							   -values   => ["no_out","dsim","dyak"],
							   -labels   => \%outgroup_labels,
							   -default  => $current->{outgroup}
								);
		
=head1
		$form .= p(strong("Select sequences:"));
		$form .= checkbox_group( -name=>$self->config_name('sequences'),
								 -values=>['outgrup'],
								 -default=>['outgrup'],
								 -linebreak=>'false',
								 -disabled=>[''],
								 -labels=>\%sequence_labels);

							
		$form   .= p(strong("Modifly coordenates to analyze... :")."(default: current region)");
		$form .= submit(-name    => 'button_coord_vs',
				-value   => '[show/hide]', 
				-onClick => 'displaySeqCoord()'); 
				
		$form .= "<div id=\"vs_seqcoord\" style=\"display:none\">";
			$form   .= "<p></p>";
			$form   .= "<span>Start position </span>".textfield(-name  => $self->config_name('startpos'),
								 -value => $current->{startpos} );
			
			$form   .= p(strong(""));
			$form   .= "<span>End position  </span>".textfield(-name  => $self->config_name('endpos'),
								 -value => $current->{endpos} );
								
		$form   .= "</div>";
		
		$form .= p(strong("...or select a region visualy:"));
		$form .= submit(-name    => 'button_select_region',
						-value   => 'Select Region', 
						-onClick => 'scrollTo("track_Region Scale");'); 
=cut

						
		$form .= "</fieldset>";


		#LDSinglets
		$form .= "<fieldset><legend>LD options</legend>";
		$form .= p(strong("Singletons usage:"));
		my %ldsinglets_labels = ( '0' => 'Ignore singletons for LD analysis.',
		                      '1' => 'Use singletons for LD analysis.');
		$form .= radio_group( -name=>$self->config_name('ldsinglets'),
								 -values=>['0','1'],
								 -default=>$current->{ldsinglets},
								 -linebreak=>'true',
								 -labels=>\%ldsinglets_labels);
		$form .= "</fieldset>";
		
		
		
		#ADVANCED
		$form .= h2("Advanced Options:");
		$form .= "<fieldset><legend id=\"advanced_vs_title\">Advanced options</legend>";

=head1
		$form .= submit(-name    => 'button_advanced_vs',
						-value   => '[show/hide]', 
						-style   => 'float:right',
						-onClick => 'displayAdvancedVS()'); 
		$form .= "<div id=\"advanced_vs_options\" style=\"display:none\">";
=cut
		
		$form .= p(strong("Window Type:"));
		my %windowtype_labels = ( '0' => 'Number of sites refers to number of columns in the alignment.',
		                          '1' => 'Number of sites refers to number of net sites (excluding discarded sites).',
		                          '2' => 'Number of sites refers to the number of polymorphic sites.',
		                          '3' => 'Number of sites refers to number of positions in the reference sequence.');
		$form .= radio_group( -name=>$self->config_name('windowtype'),
								 -values=>['0','1','2','3'],
								 -default=>$current->{windowtype},
								 -linebreak=>'false',
								 -labels=>\%windowtype_labels);


		#$form .= $jscript;
		$form .= p(strong("Segregating Sites / Number of Mutations:"));
		my %usemuts_labels = ( '0' => 'Use total number of Segregating Sites',
		                       '1' => 'Use total number of Mutations');
		$form .= radio_group( -name=>$self->config_name('usemuts'),
								 -values=>['0','1'],
								 -default=>$current->{usemuts},
								 -linebreak=>'true',
								 -labels=>\%usemuts_labels);
								 
		$form .= p(strong("Gap/Missing/Ambiguous sites usage:"));
				
		my $numnuc_field = textfield(-name  => $self->config_name('numnuc'),
								   -value => $current->{numnuc}
								  );  
		my %fixnum_labels = ( '0' => 'at least a',
		                      '1' => 'a defined');
		my $fixnum_field = popup_menu( -name=>$self->config_name('fixnum'),
									   -values=>['0','1'],
                                       -labels=>\%fixnum_labels,
									   -default=>$current->{fixnum} ,
								       -linebreak=>'true');     				  
		my %completedeletion_labels = ( '0' => 'Use all sites with:',
		                                '1' => 'Complete Deletion (Exclude Gap/Missing/Ambiguous sites. Allways used in LD analysis)');                               
		$form .= radio_group( -name=>$self->config_name('completedeletion'),
								 -values=>['1','0'],
								 -default=>$current->{completedeletion},
								 -linebreak=>'true',
								 -labels=>\%completedeletion_labels).$fixnum_field.span(" number of ").$numnuc_field.span(" (A/T/G/C) sequences.");          # JS in /gbrowse/js/dgvb.js 
      
#$form .= "</div>";
		$form .= "</fieldset>";


		return $form;
	 }

 # called by gbrowse to annotate the DNA, returning features
	 sub dump {
		my $self    = shift;
		my $segment = shift;
		my $config = $self->configuration;
		
		my $time_start = time;

   #######Check form parameters:
		
	#Sequence start/end
   		my $chr         = $segment->seq_id;
		my $start       = $segment->start;      #odd, segment->start 'starts' at 0 instead of 1 ?
		my $end         = $segment->end;        #odd, segment->end adds +1 when start is not 1 ?
		my $region_size = ($end-$start);
		        
		#print "Gbrowse coord: ".$start."..".$end."\n";

	#Window size
		my $ws    = $config->{ws};
		my $jump  = $config->{jump};

		if ($ws >= $region_size || $ws == "" || $ws == 0) { 
			$ws = ($region_size+1); 
			$jump = $ws;
		} elsif ($jump > $ws || $jump == 0) {
			$jump = $ws;
		}
		
	#Max region and window size
		$ws = 200000 if $ws > 200000;
		$jump = $ws if $jump > 200000;
		$end = $start + 1000000 if $region_size > 1000000;
		

	#window type
		$windowtype       = $config->{windowtype};
		
	#Sites usage
		$usemuts          = $config->{usemuts};
		$completedeletion = $config->{completedeletion};
		$fixnum           = $config->{fixnum};
		$numnuc           = $config->{numnuc};
		
	#LDsinglets
		$ldsinglets       = $config->{ldsinglets};

#TEMP FIX!!!!!
#		#window type
#		$windowtype       = 3;
#		
#		#Sites usage
#		$usemuts          = 1;
#		$completedeletion = 0;
#		$fixnum           = 0;
#		$numnuc           = 150;
#		
#		#LDsinglets
#		$ldsinglets       = 1;


   #################### VariScan Options
		my (%options, %results, @cols)=();
		my @order=('s','s_inter','eta','eta_e','pi','theta','tajd','fulidast','fulifast','k','fulid','fulif','faywuh','ldsites','d','dabs','dprime','dprimeabs','r2','h','hd','fufs');

		die "Less than 22 stats in \@order" if (scalar(@order)!=22);

	#default VariScan values (Fixed = can't not be changed by parameters)
		$options{'RefSeq'}=2;
		$options{'RefSeq'}=1 if $config->{outgroup} eq "no_out";
		$options{'UseMuts'}=$usemuts;
		$options{'WindowType'}=$windowtype;
		$options{'Outgroup'}="first";
		$options{'Outgroup'}="none" if $config->{outgroup} eq "no_out";
		$options{'UseLDSinglets'}=$ldsinglets;
		$options{'WidthSW'}= $ws;
		$options{'JumpSW'}= $jump;
		$options{'SeqChoice'}="all";                               #fixed
		$options{'StartPos'}=$start; 
		$options{'BlockDataFile'}="none";                          #fixed
		$options{'CompleteDeletion'}=$completedeletion;            #fixed
		$options{'FixNum'}=$fixnum;                                #fixed, not used if CompleteDeletion=1
		$options{'NumNuc'}=$numnuc;                                #fixed, not used if CD=1
		$options{'IndivNames'}="a b c d";                          #not used
		$options{'RefPos'}=1;                                      #fixed, the StartPos and EndPos in the reference Seq (0 in the aln)
		$options{'EndPos'}=$end;                                   #if 0, the entire aln
		$options{'SlidingWindow'}=1;
		$options{'chr'}=$chr;
		$options{'PairLDFile'}="/tmp/LDpair.tsv";
		$options{'seed'}=5; #test pourposes

	#CHOOSE FILE PATH
		my $is_out = $config->{outgroup};
		
	#WITH-OUTGROUP FIXED PARAMs
		if ($config->{statistic} eq "k" || $config->{statistic} eq "s_inter" || $config->{statistic} eq "fulid" || $config->{statistic} eq "fulif" || $config->{statistic} eq "faywuh"){
			$is_out = "dsim" if $config->{outgroup} eq "no_out";
			$options{'Outgroup'}="first";
			$options{'RefSeq'}=2;
		}
		
	#LD FIXED PARAMS
		if ($config->{statistic} eq 'ldsites' || $config->{statistic} eq 'd' || $config->{statistic} eq 'dabs' || $config->{statistic} eq 'dprime' || $config->{statistic} eq 'dprimeabs' || $config->{statistic} eq 'r2' || $config->{statistic} eq 'h' || $config->{statistic} eq 'hd' || $config->{statistic} eq 'fufs') {
			$options{'CompleteDeletion'}=1 ;
			$is_out = "dsim" if $config->{outgroup} eq "no_out";
			$options{'Outgroup'}="first";
			$options{'RefSeq'}=2;
		}
		
		#$options{'exe'}="./home/dgrp/precomputed/variscan";
		#$options{'f'}="/home/share/dgrp.freeze1/".$chr.".all.".$is_out.".nodupl.fasta.phy";     #phylip aln

		my $resRef = chooseFile($start,$end,$chr,$is_out);
		my @res_arr = @$resRef;
		$options{'f'} = $res_arr[0];
		
	#convert coordenates to in-file coordenates
		my $alnStart = 0;
		my $alnEnd = 0;
		
		if ($res_arr[1] ne "NO"){
			$alnStart = $res_arr[1];
			$alnEnd = $res_arr[2];
			#print "File start-end: ".$alnStart."..".$alnEnd."\n";
			$options{'StartPos'} =  ($start - $alnStart) +1;
			$options{'EndPos'} =  ($end - $alnStart) +1;
			#$options{'EndPos'} =  ($end - $alnStart) + 1 if $options{'StartPos'} <= 1;
			#print "Modified coord: ".$options{'StartPos'}."..".$options{'EndPos'}."\n";
		}


	#Variscan Tests
		$options{'s'}="";
		$options{'s_inter'}="";
		$options{'eta'}="";
		$options{'eta_e'}="";
		$options{'pi'}="";
		$options{'theta'}="";
		$options{'tajd'}="";
		$options{'fulidast'}="";
		$options{'fulifast'}="";
		$options{'k'}="";
		$options{'fulid'}="";
		$options{'fulif'}="";
		$options{'faywuh'}="";
		$options{'ldsites'}="";
		$options{'d'}="";
		$options{'dabs'}="";
		$options{'dprime'}="";
		$options{'dprimeabs'}="";
		$options{'r2'}="";
		$options{'h'}="";
		$options{'hd'}="";
		$options{'fufs'}="";



	###### check selected statistics
		foreach my $stat (keys %options){
			if ($stat eq $config->{statistic}) { 
				$options{$stat} = $config->{statistic}; 
			}
		}
	
	###### create one-liner
		$options{'EndPos'} = $options{'EndPos'} -1 ;                   #fix different coordenate systems between VariScan and GBrowse
		my $variscanopt=&options2string(\%options, \@order, \@cols);
		$options{'EndPos'} = $options{'EndPos'} +1 ;                   #fix different coordenate systems
		
		
		#&Usage if (!$options{'f'} || scalar(@cols)==0);
		unlink ($options{'PairLDFile'});

		my $command="/home/share/variscan/variscan ".$options{"f"}.' '.$variscanopt;
#		print STDERR "Running $command\n";
#		print $command."\n\n";

		(!&runAndSave(\@cols, $command, \%results)) && do{die "Error computing variscan\n$command\n";};
		(!&printTEXT(\%results, \%options, $time_start, $alnStart, $is_out)) && do{die;};
				
	}





################################################
###############  VARISCAN  #####################
###############  FUNCTIONS #####################
################################################

#### Choose sequence file

	sub chooseFile {
		my ($start_cord, $end_cord, $chr, $outgr) = @_;
		my $region_size = $end_cord - $start_cord;
		my $folder = "/home/share/dgrp.freeze1/overlap/$chr/$outgr/";
		my @results = ();
		
		if ($region_size < 125000) {
			$folder .= "250Kb/";
		} elsif ($region_size < 250000) {
			$folder .= "500Kb/";
		} elsif ($region_size < 500000) {
			$folder .= "1Mb/";
		} elsif ($region_size < 1500000) {
			$folder .= "3Mb/";
		} elsif ($region_size < 3500000) {
			$folder .= "7Mb/";
		} else {
			$folder = "/home/share/dgrp.freeze1/".$chr.".all.".$outgr.".nodupl.fasta.phy";	
			push(@results, $folder);
			push(@results, "NO");
			push(@results, "NO");
		}
		
		#print $folder."\n";
		opendir(ODIR, $folder) or die("Cannot open directory $folder"); 
		my @dirfiles = readdir(ODIR);
		close ODIR;

		my $previous_1 = 0;
		foreach (@dirfiles) {
			next if $_ eq ".";
			next if $_ eq "..";
			$_ =~ m/.*\.(\d+)-(\d+)\.fa\.phy/;
		
			if ( $1 <= $start_cord && $2 >= $end_cord  ){ 
					my $path = $folder.$_;
					
					if ($1 > $previous_1) {
						push(@results, $path);
						push(@results, $1);
						push(@results, $2);
					}
					
					$previous_1 = $1;
			}
		}
		return \@results;
	}



#### Variscan options string
	sub options2string($$$){
		my ($options, $order, $cols)=@_;
		my $string="\"RunMode = ";
		
		#add statistic options
		foreach my $option (@$order){
		unless($$options{$option}){
			$$options{$option}=0;
		}else{
			$$options{$option}=1;
			push(@$cols,$option);
		}
		$string.="$$options{$option} ";
		}

		chop $string; $string.="|";

		#add the other variscan options
		foreach my $option (keys %$options){
		next if (grep(/$option/,@$order) || $option eq "inputfile" || $option eq "chr");
		$string.="$option = $$options{$option}\|";
		
		}
		chop $string; $string.="\"";;
		return $string;
	}

####EXECUTE (from modification)
	sub runAndSave($$$){
		my ($cols,$command, $results)=@_;

		(!open (COMMAND,"$command |")) && do{warn "$!\n";return 0;};
			while(<COMMAND>){
					chomp;
					next if ($_=~/^#/);
					$_=~s/^\s+//;	
					my @line=split(/\s+/,$_);
					for (my $i=8;$i<scalar(@line);$i++){
					push(@{$$results{$$cols[$i-8]}},{'startrefpos'=>$line[0],'endrefpos'=>$line[1], 'midrefpos'=>$line[2],'alnpos'=>$line[5],'val'=>$line[$i]});
				}
			#	print "$_\n";
			}
		close(COMMAND);
		if (!defined($$results{$$cols[0]})){return 0;}
		return 1;
	}

####PRINT
	sub printTEXT(){
		my ($results, $options, $time_start, $alnStart, $is_out)=@_;
		

		
		my %ldsinglets_labels       = ( '0' => 'Ignore singletons for LD analysis.',
					                    '1' => 'Use singletons for LD analysis.');
		my %windowtype_labels       = ( '0' => 'Number of sites refers to number of columns in the alignment.',
						                '1' => 'Number of sites refers to number of net sites (excluding discarded sites).',
						                '2' => 'Number of sites refers to the number of polymorphic sites.',
						                '3' => 'Number of sites refers to number of positions in the reference sequence.');
		my %usemuts_labels          = ( '0' => 'Use total number of Segregating Sites',
					                    '1' => 'Use total number of Mutations');
		my %fixnum_labels           = ( '0' => 'at least a',
						                '1' => 'a defined');
		my %completedeletion_labels = ( '0' => 'Use all sites with '.$fixnum_labels{$$options{'FixNum'}}.' number of '.$$options{'NumNuc'}.' valid nucleotides (A/T/G/C)',
								        '1' => 'Complete Deletion (Exclude Gap/Missing/Ambiguous sites. Allways used in LD analysis)');     

							  
		my $swsize = $$options{'WidthSW'};
		my $LDsingl = $usemuts_labels{$$options{'UseLDSinglets'}};
		my $WindowType = $windowtype_labels{$$options{'WindowType'}};
		my $UseMuts = $usemuts_labels{$$options{'UseMuts'}};
		my $CompDel = $completedeletion_labels{$$options{'CompleteDeletion'}};
		
		#don't add $alnStart when its 1 and correct coordinates starting with 0
		#$alnStart = 0 if $alnStart == 1;
		$$options{'StartPos'} = 1 if $$options{'StartPos'} == 0;
		
		my $pStart = ($$options{'StartPos'}+$alnStart)-1;
		my $pEnd = ($$options{'EndPos'}+$alnStart)-1;
		my $rsize = ($$options{'EndPos'} - $$options{'StartPos'});
		
		
		print "### Reestimation of Diversity Values\n";
		print "#Chromosome: ".$$options{'chr'}."\n";
		print "#Region: ".$pStart."..".$pEnd." ($rsize bp)\n";
		print "#Sliding-Window size: ".$$options{'WidthSW'}."bp\n";
		print "#Jump size: ".$$options{'JumpSW'}."bp\n";
		print "#Outgroup: ".$is_out."\n";
		print "#Singleton usage: ".$LDsingl."\n";
		print "#Window type: ".$WindowType."\n";
		print "#Use of Segregating sites/Mutations: ".$UseMuts."\n";
		print "#Gap/Missing/Ambigous site usage: ".$CompDel."\n";
		
		my $time_end = time - $time_start;
				
		foreach my $stat (keys %$results){
			print "\n#Statistic= $stat\n";
			
			foreach my $point (@{$$results{$stat}}){
				my $point_text = "";
				$point->{'startrefpos'} = ($point->{'startrefpos'} + $alnStart)-1 ;
				$point->{'endrefpos'} = ($point->{'endrefpos'} + $alnStart)-1 ;
				($rsize == $swsize) ? ($point_text = "Total value: ") : ($point_text = $point->{'startrefpos'}."..".($point->{'endrefpos'}+1) ) ;
				
				my $window_val = $point->{'val'};
				my $rounded = sprintf("%0.5f", $window_val);
				$rounded = "nan" if $stat eq 'k' && $rounded == 5;

				
				print $point_text."\t".$rounded."\n"; 
			}
			print "\n";
		}
		my $time_end2 = time - $time_start;
		print "\n\n### Time vs: $time_end s\n";
		print "### Time print: $time_end2 s\n";
		
		return 1;
	}
