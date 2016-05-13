
# writes statistics on the input file to two output files
# statistics_on_secondary_structure.txt contains statistics of the whole file and all sequences together
# statistics_on_secondary_structure_every_sequence.txt contains tab-delimited information about the different secondary structure elements (in percent) + length

open (IN2, 'opm_unmasked_hval0_just_ss.fasta');
open (my $write, '>', 'statistics_on_secondary_structure_opm.txt');
open (my $per_sequence, '>', 'statistics_on_secondary_structure_every_sequence_opm.txt');
 
 print $per_sequence "zeros\tones\ttwos\tU\tL\tH\tI\tLength of Sequence\n";
 
 $I_total=0;
 
 while (my $line=<IN2>){
	
	$zeros=0; $ones=0; $twos=0;
	$U=0; $L=0; $H=0; $I=0;
	
	chomp($line);
	
	@split_string_into_chars = split("",$line);
	
	$number_of_total_characters_in_string = @split_string_into_chars;
	$total_characters=$total_characters+$number_of_total_characters_in_string;
	
	foreach $char (@split_string_into_chars) {

		if($char =~ 0){
			$zeros_total++;
			$zeros++;
		}
		if($char =~ 1){
			$ones_total++;
			$ones++;
		}
		if($char =~ 2){
			$twos_total++;
			$twos++;
		}
		if($char =~ U){
			$U_total++;
			$U++;
		}
		if($char =~ L){
			$L_total++;
			$L++;
		}
		if($char =~ H){
			$H_total++;
			$H++;
		}
		if($char =~ h){		#in opm file, there are several lower-case H
			$H_total++;
			$H++;
		}
		if($char =~ I){
			$I_total++;
			$I++;
		}
	}

	
	$zeros_per_sequence = $zeros/$number_of_total_characters_in_string;
	$ones_per_sequence = $ones/$number_of_total_characters_in_string;
	$twos_per_sequence = $twos/$number_of_total_characters_in_string;
	$U_per_sequence = $U/$number_of_total_characters_in_string;
	$L_per_sequence = $L/$number_of_total_characters_in_string;
	$H_per_sequence = $H/$number_of_total_characters_in_string;
	$I_per_sequence = $I/$number_of_total_characters_in_string;
	
	$is_it_1=0;
	$is_it_1=$zeros_per_sequence+$ones_per_sequence+$twos_per_sequence+$U_per_sequence+$L_per_sequence+$H_per_sequence+$I_per_sequence;
	
	if ($is_it_1 eq 1){   #tests if the numbers add up to 1

	}
	else{
		print $line; print "\n";
		print $is_it_1;
		print "\n";
		print "$zeros/$number_of_total_characters_in_string". "+" . "$ones/$number_of_total_characters_in_string" . "+" . "$twos/$number_of_total_characters_in_string";
		print "\n\n";
	}
	
	
	 print $per_sequence "$zeros_per_sequence\t$ones_per_sequence\t$twos_per_sequence\t$U_per_sequence\t$L_per_sequence\t$H_per_sequence\t$I_per_sequence\t$number_of_total_characters_in_string\n";
 
 
	undef (@split_string_into_chars);
	$number_of_total_characters_in_string=0;
	undef ($char);
 
 
 }

$zeros_total_percentage=$zeros_total/$total_characters;
$ones_total_percentage=$ones_total/$total_characters;
$twos_total_percentage=$twos_total/$total_characters;
$U_total_percentage=$U_total/$total_characters;
$L_total_percentage=$L_total/$total_characters;
$H_total_percentage=$H_total/$total_characters;
$I_total_percentage=$I_total/$total_characters;
		
print $write "Secondary Structure Elements in total:\n" . "0:\t$zeros_total\t\t$zeros_total_percentage\n" . "1:\t$ones_total\t\t$ones_total_percentage\n" . "2:\t$twos_total\t\t$twos_total_percentage\n" . "U:\t$U_total\t\t$U_total_percentage\n" . "L:\t$L_total\t\t$L_total_percentage\n" . "H:\t$H_total\t\t$H_total_percentage\n" . "I:\t$I_total\t\t$I_total_percentage\n" . "Total:\t$total_characters\n";


 
 
  
 
	# if ( substr( $line2, 0, 1 )  =~ 1|2|0|L|U|H ){
		# $count++;
	# }
	# else {
		# print $line2;
	# }
	
	
 # print "$count\n";