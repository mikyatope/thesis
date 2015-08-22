my @cards = (0..51);

shuffle( \@cards );    # permutes @array in place

print "@cards";

# randomly permutate @array in place
sub shuffle
{
	#fisher yates shuffle
    my $array = shift;
    my $i = @$array;
    while ( --$i )
    {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}

