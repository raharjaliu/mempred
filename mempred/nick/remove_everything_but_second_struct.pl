
#removes everything but the secondary structure information from a given fasta file 

open (IN, '../../data/sets/opm_unmasked_hval0.fasta');			#opm or pdbtm
open (my $write, '>', 'opm_unmasked_hval0_just_ss.fasta');

	$i = 0;
	
while (my $line=<IN>){
	
	$i++;
	if ($i=~3){
		print $write "$line";
		$i=0;
	}
	
}

close IN;
close $write;
