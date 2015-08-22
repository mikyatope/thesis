use strict;
use warnings;
use Getopt::Long;


my (%options, %results, @cols)=();
my @order=('s','s_inter','eta','eta_e','pi','theta','tajd','fulidast','fulifast','k','fulid','fulif','faywuh','ldsites','d','dabs','dprime','dprimeabs','r2','h','hd','fufs');

die "Less than 22 stats in \@order" if (scalar(@order)!=22);

#default values (Fixed = can't not be changed by parameters)
$options{'RefSeq'}=1;
$options{'UseMuts'}=1;
$options{'WindowType'}=1;
$options{'Outgroup'}="last";
$options{'UseLDSinglets'}=0;
$options{'WidthSW'}=0;
$options{'JumpSW'}=0;
$options{'SeqChoice'}="all"; #fixed
$options{'StartPos'}=1; 
$options{'BlockDataFile'}="none"; #fixed
$options{'CompleteDeletion'}=0; #fixed
$options{'FixNum'}=1; #fixed, not used if CompleteDeletion=1
$options{'NumNuc'}=140; #fixed, not used if CD=1
$options{'IndivNames'}="a b c d"; #not used
$options{'RefPos'}=1; #fixed, the StartPos and EndPos in the reference Seq (0 in the aln)
$options{'EndPos'}=0; #if 0, the entire aln
$options{'SlidingWindow'}=1;
$options{'chr'}='3R';
$options{'PairLDFile'}="/tmp/LDpair.tsv";
$options{'seed'}=time;
$options{'trmod'} = "";

GetOptions(

           #inputfile
           'f|inputfile:s'        => \$options{'f'}, #phylip aln
           #inputfile
           'bdf|BlockDataFile:s'  => \$options{'BlockDataFile'}, 
           #seqchoice
  #         'sqch|SeqChoice:s'  => \$options{'SeqChoice'}, 

           #chr and track info
           'chr:s'                     => \$options{'chr'},
           'trmod|TrackModificator:s'  => \$options{'trmod'},

           #variscan options (read variscan doc)
           'm|UseMuts:i'          => \$options{'UseMuts'}, 
           'sw|SlidingWindow:i'   => \$options{'SlidingWindow'},
           'w|WidthSW:i'          => \$options{'WidthSW'},
           'j|JumpSW:i'           => \$options{'JumpSW'},
	       'wt|WindowType:i'      => \$options{'WindowType'},
	       'o|Outgroup:s'         => \$options{'Outgroup'},
	       'LDsin|UseLDSinglets:i'=> \$options{'UseLDSinglets'},
           'RefSeq:i'             => \$options{'RefSeq'}, 
           'EndPos:i'             => \$options{'EndPos'},
           'StartPos:i'           => \$options{'StartPos'},
           'CompleteDeletion:i'   => \$options{'CompleteDeletion'},
           'FixNum:i'             => \$options{'FixNum'},
           'NumNuc:i'             => \$options{'NumNuc'},

           #added by pablo
           'LDF|PairLDFile:s'     => \$options{'PairLDFile'},
	       'seed:i'		          => \$options{'seed'},
	   
	   #statistics
	   'S|s!'                 => \$options{'s'},
	   'S_inter|s_inter!'     => \$options{'s_inter'},
	   'Eta|eta!'             => \$options{'eta'},
       'Eta_e|eta_e!'         => \$options{'eta_e'},
	   'pi|Pi!'               => \$options{'pi'},
	   'theta|Theta!'         => \$options{'theta'},
	   'tajD!'                => \$options{'tajd'},
	   'fuliDast!'            => \$options{'fulidast'},
	   'fuliFast!'            => \$options{'fulifast'},
	   'K!'                   => \$options{'k'},
	   'fuliD!'               => \$options{'fulid'},
	   'fuliF!'               => \$options{'fulif'},
	   'FayWuH!'              => \$options{'faywuh'},
	   'LD_sites!'            => \$options{'ldsites'},
	   'D!'                   => \$options{'d'},
	   'Dabs!'                => \$options{'dabs'},
	   'Dprime!'              => \$options{'dprime'},
	   'Dprimeabs!'           => \$options{'dprimeabs'},
	   'r2!'                  => \$options{'r2'},
	   'h!'                   => \$options{'h'},
	   'hd!'                  => \$options{'hd'},
	   'FuFs!'                => \$options{'fufs'},
	   'missing!'		      => \$options{'missing'},
	   'netsites!'		      => \$options{'NetSites'}	
    );

my $variscanopt=&options2string(\%options, \@order, \@cols);

&Usage if (!$options{'f'} || scalar(@cols)==0);

unlink ($options{'PairLDFile'});

my $filename = "results/".$options{'chr'}.".".$options{'trmod'}.".".$options{'WidthSW'};
open CSV_FILE, ">".$filename.".csv" or die $!;
open WIG_FILE, ">".$filename.".wig" or die $!;

my $command='./variscan '.$options{"f"}.' '.$variscanopt;
print CSV_FILE "#Running $command\n";
print WIG_FILE "#Running $command\n";

(!&runAndSave(\@cols, $command, \%results)) && do{die "Error computing variscan\n";};
(!&printGBrowse(\%results, \%options)) && do{die;};

close CSV_FILE;
close WIG_FILE;

#####################################################################
############################FUNCTIONS################################
#####################################################################
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
	#discard perl-only options
	next if (grep(/$option/,@$order) || $option eq "inputfile" || $option eq "chr" || $option eq 'NetSites' || $option eq 'missing' || $option eq "trmod" );
	$string.="$option = $$options{$option}\|";
	
    }
    chop $string; $string.="\"";;
    return $string;
}
####################################################################
sub runAndSave($$$){
    my ($cols,$command, $results)=@_;

    (!open (COMMAND,"$command |")) && do{warn "$!\n";return 0;};
    while(<COMMAND>){
		chomp;
		next if ($_=~/^#/);
		$_=~s/^\s+//;	
		my @line=split(/\s+/,$_);
		if ($options{'missing'}){
			push(@{$$results{'missing'}},{'refpos'=>$line[2],'alnpos'=>$line[5],'val'=>$line[7]});   
		}
		if ($options{'NetSites'}){
			push(@{$$results{'NetSites'}},{'refpos'=>$line[2],'alnpos'=>$line[5],'val'=>$line[6]});   
		}
		
		for (my $i=8;$i<scalar(@line);$i++){
			push(@{$$results{$$cols[$i-8]}},{'startrefpos'=>$line[0],'endrefpos'=>$line[1], 'midrefpos'=>$line[2],'alnpos'=>$line[5],'val'=>$line[$i],'miss'=>$line[7],'net'=>$line[6] });
		}
#print "$_\n";
    }
    close(COMMAND);
    if (!defined($$results{$$cols[0]})){return 0;}
    return 1;
}
####################################################################
sub printGBrowse($$){
		my ($results, $options)=@_;
		my $chr = $$options{'chr'};
		my $ws = $$options{'WidthSW'};
		my $trmod = "";
		$trmod .= $$options{'trmod'};
		
		print CSV_FILE "#windowsize=$ws reference=$chr\n";
		print WIG_FILE "#windowsize=$ws reference=$chr\n";
		
		foreach my $stat (keys %$results){
		
			print CSV_FILE "###$stat\n";
			print WIG_FILE "###$stat\n";
			
			print CSV_FILE "chr\tstart\tstop\tvalue\tnet_sites\tmissing\n";
			
			print WIG_FILE 'track type=wiggle_0 name="'.$stat.'_'.$ws.'"  \n';
			print WIG_FILE "fixedStep chrom=$chr start=1 step=$ws span=$ws\n";
			
			foreach my $point (@{$$results{$stat}}){
				my $value = $point->{'val'};                                  #save value

				#transform values
				$value = "NA" if $point->{'val'} eq "nan";                    #convert "nan" to "0.0". Gbrowse bug, must change to "." when corrected.
				$value = sprintf('%.3f',0) if $point->{'val'} == 0    ;       #convert 0 to "0.000". Gbrowse bug with integers
				$value = sprintf('%.3f',$value) if ($value =~ /^[+-]?\d+$/);  #check for integers and add 3 decimals ".000". Gbrowse bug with integers

				#discard for less than 50% netsites
				$value = "NA" if $point->{'net'} < int($ws/2);
				#
				print CSV_FILE $chr."\t".$point->{'startrefpos'}."\t".$point->{'endrefpos'}."\t".$value."\t".$point->{'net'}."\t".$point->{'miss'}."\n"; 
			
				#discard for less than 50% netsites
				$value = "" if $value eq "NA";
				#
				print WIG_FILE $value."\n"; 

			}
			#print "\n";
		}

		return 1;
}

####################################################################
sub Usage(){
    print "perl VariScan.pl -f [alnfile.phy] [options] > output.gbrowse\n\nOptions:\n";

    die '
            ###script options
            f|inputfile: in phylip format [required]
            chr: chromosome to specify in custom track [default: 3R]

            ####variscan options (read variscan doc)
            m|UseMuts: Use total number of mutations (1) or the number of segregating sites (0) [default: 1],
            sw|SlidingWindow: Compute stat in SW [default: 1],  
            w|WidthSW: Window Length [default: 10000]
            j|JumpSW: Window Jump [default: 2000]  
            wt|WindowType: Window Type [default: 1] #see variascan doc 
            o|Outgroup: Outgroup sequences [default: last] 
            LDsin|UseLDSinglets: Use singletons in LD comp [default: 0]
            RefSeq: Refseq [default: 1]
            EndPos: Last position analysed [default: 0, the last]
            StartPos: First position analysed [default: 1]
            CompleteDeletion: CompleteDeletion(1) or not (0) [default: 1]
            FixNum: To specify when CompleteDeletion is 0, read doc [default: 0] 
            NumNuc: To specify when CompleteDeletion is 0, read doc [default: 0]
	    missing: prints missing column
	    NetSites: prints the nte number of sites analyzed
	    seed: seed passed on to variscan
	    
            ###statistics to compute (simply add -[stat] to compute it [required, at least one]):
            S|s
            S_inter|s_inter
            Eta|eta
            Eta_e|eta_e
            pi|Pi
            theta|Theta
            tajD
            fuliDast
            fuliFast
            K
            fuliD
            fuliF
            FayWuH
            LD_sites
            D
            Dabs
            Dprime
            Dprimeabs
            r2
            h
            hd
            FuFs
         ';
}
