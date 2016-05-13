# ignore, just a testing routine


open (IN2, 'pdbtm_unmasked_hval0_just_ss.fasta');
 
 while (my $line2=<IN2>){
	
	$count_insgesamt++;
 
	if ( substr( $line2, 0, 1 )  =~ 1|2|0|L|U|H|I ){
		$count++;
	}
	else {
		print $line2;
	}
 
 }

 print "$count\n";
 print $count_insgesamt;