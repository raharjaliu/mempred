use Data::Dumper;

# searches given "*.just_ss.fasta" and returns profiles for helices with a length of 25

open (IN, 'opm_unmasked_hval0_just_ss_short_tester.fasta');
# open (my $write, '>', 'statistics_on_secondary_structure_opm.txt');
 
 $i=0;
 %helices=();
 %nonhelices=();
 while (my $line=<IN>){
	
	chomp($line);
	
	@split_string_into_chars = split("",$line);

	# my @c = $line =~ /H/g;
	$helix_number=0;
	$nonhelix_number=0;
	
	while($i<@split_string_into_chars){
	print "i: ".$i."\n";
		if($split_string_into_chars[$i] =~ /H/){

			while ($split_string_into_chars[$i] =~ /H/){
				$helix_string=$helix_string.'H';
				$helices{$helix_number}={$helix_string};
				$i++;
			}
			$helix_number++;
		}
		else{
			while ($split_string_into_chars[$i] !~ 'H'){
			# print "L";
			# print Dumper %nonhelices;
				$nonhelix_string=$nonhelix_string.$split_string_into_chars[$i];
				$nonhelices{$nonhelix_number}={$nonhelix_string};
				$i++;
			}
			$nonhelix_number++;
			
		}
		
	}

}

close IN;

print Dumper %helices;
print Dumper %nonhelices;



		# print "char: ". $char . "\n";
		# print Dumper $split_string_into_chars[$char];
		# print Dumper $split_string_into_chars[$char+1];
		# print "------------\n";
		# while($split_string_into_chars[$char] =~ H){
			# $counter++;
			# $helix_string=$helix_string+'H';
			# $helices{'helix'}={'$helix_string'}
			# next;
			
		# }
		# $index++;
		# print Dumper %helices;