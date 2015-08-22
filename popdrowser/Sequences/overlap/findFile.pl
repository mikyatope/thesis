my $start_cord = $ARGV[0];
my $end_cord = $ARGV[1];
my $region_size = $end_cord - $start_cord;
my $folder = "/home/share/dgrp.freeze1/overlap/2L/";

if ($region_size < 1000000) {
	$folder .= "no_out/1Mb/";
} elsif ($region_size < 3000000) {
	$folder .= "no_out/3Mb/";
} elsif ($region_size < 7000000) {
	$folder .= "no_out/7Mb/";
} else {
	#$folder ="whole"; #### Change!!!!	
}

print $folder."\n";
opendir(ODIR, $folder) or die("Cannot open directory $folder"); 
my @dirfiles = readdir(ODIR);
close ODIR;

my $previous_1 = 0;
foreach (@dirfiles) {
	next if $_ eq ".";
	next if $_ eq "..";
	$_ =~ m/.*\.(\d+)-(\d+)\.fa\.phy/;

	if ( $1 <= $start_cord && $2 >= $end_cord  ){ 
			print $_."\n"	if $1 > $previous_1;
			$previous_1 = $1;
	}
}