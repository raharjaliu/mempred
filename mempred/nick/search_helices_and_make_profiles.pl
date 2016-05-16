
# searches "*.fasta" (commandline arg) and returns profiles for helices with a length of 25
# returns two text files: one for helices & one for non-helices


use Data::Dumper;
$commandline_arg=$ARGV[0];


open (IN, $commandline_arg);
open ($out_helix, '>', 'profiles_helix_'.$commandline_arg.'.txt');
open ($out_nonhelix, '>', 'profiles_nonhelix_'.$commandline_arg.'.txt');

 %id_seq_ss=();
 %helices=();
 %nonhelices=();
 %all=();
 
 $i=0;
 $j=0;
 $line ='';
 while (my $line=<IN>){		# parse data into hash
 
	chomp($line);
	push @{ $id_seq_ss{$j} }, $line;
	# $id_seq_ss{$j}={$line};
	$j++;
	
}

$counter=0;
$counter_2=0;
while ($counter<$j){		#main routine: searches for helices with a length > 25 and saves them and the corresponding sequence to a hash
	
	$stored_line_id='';
	$stored_line_id=$id_seq_ss{$counter}[0];
	$stored_line_seq=$id_seq_ss{$counter+1}[0];
	$stored_line_ss=$id_seq_ss{$counter+2}[0];
	$all{$counter_2}{$stored_line_id}{$stored_line_ss}={$stored_line_seq};
	$counter=$counter+3;
	

	@split_string_into_chars = split("",$stored_line_ss);
	@split_string_into_chars_seq = split("",$stored_line_seq);
	$helix_number=0;
	$nonhelix_number=0;
	$i=0;
		while($i<@split_string_into_chars){
			
			if($split_string_into_chars[$i] =~ m/(H|h)/){
				$helix_string='';
				$helix_string_seq='';

				while (($split_string_into_chars[$i] =~ m/(H|h)/)&& $i<@split_string_into_chars){		#helix
					$helix_string=$helix_string.'H';
					$helix_string_seq=$helix_string_seq.$split_string_into_chars_seq[$i];
					if (length($helix_string)>24){
						undef $helices{$stored_line_id}{'Helix '.$helix_number};
						$helices{$stored_line_id}{'Helix '.$helix_number}{$helix_string}=$helix_string_seq;
					}
					$i++;
				}
				$helix_number++;
			}
			else{
				$nonhelix_string='';
				$nonhelix_string_seq='';
				while (($split_string_into_chars[$i] !~ m/(h|H)/)&& $i<@split_string_into_chars){		#non-helix
					$nonhelix_string=$nonhelix_string.$split_string_into_chars[$i];
					$nonhelix_string_seq=$nonhelix_string_seq.$split_string_into_chars_seq[$i];
					if (length($nonhelix_string)>24){
						undef $nonhelices{$stored_line_id}{'Non-Helix '.$helix_number};
						$nonhelices{$stored_line_id}{'Non-Helix '.$nonhelix_number}{$nonhelix_string}=$nonhelix_string_seq;
					}
					$i++;
				}
				$nonhelix_number++;
			
			}
		
		}
	
	undef(@split_string_into_chars);
	$everythird++;
	$counter_2++;
}

%helices_25=();
%non_helices_25=();

print $out_helix "ID\tNumber_of_consecutive_helix\tHelix-sequence\tAA-Sequence\n";
print $out_nonhelix "ID\tNumber_of_consecutive_non-helix\tNon-helix-sequence\tAA-Sequence\n";

foreach my $id (keys %helices) {									#cuts helices at length 25
    foreach my $values (keys %{ $helices{$id} }) {
		foreach my $ss(keys %{ $helices{$id}{$values} }) {
			# $random_index = int(rand(length($ss)-25));
			# $ss_25 = substr($ss, $random_index, 25);
			$a = $helices{$id}{$values}{$ss};
			# $ss_25_seq_h = substr($a, $random_index, 25);
			# $helices_25{$id}{$values}{$ss_25}=$ss_25_seq_h;
			# print $out_helix "$id\t$ss_25\t$ss_25_seq_h\n";
			$rest=(int(length($ss)/25));
			$mod=(length($ss))%25;
			$plus=0;
			
			for ($j=0; $j<$rest; $j++){
				$random_index = int(rand($mod));
				$ss_25 = substr($ss, $random_index+($plus), 25);
				$ss_25_seq_h = substr($a, $random_index+($plus), 25);
				$helices_25{$id}{$values}{$j}{$ss_25}=$ss_25_seq_h;
				print $out_helix "$id\t$j\t$ss_25\t$ss_25_seq_h\n";
				$plus=25+$random_index+$plus;
				$mod=$mod-$random_index;
			}
			
		}
    }
}
foreach my $id (keys %nonhelices) {									#cuts non-helices at 25
    foreach my $values (keys %{ $nonhelices{$id} }) {
		foreach my $ss (keys %{ $nonhelices{$id}{$values} }) {
			# $random_index = int(rand(length($ss)-25));
			# $ss_25_nh = substr($ss, $random_index, 25);
			$b = $nonhelices{$id}{$values}{$ss};
			# $ss_25_seq_nh = substr($b, $random_index, 25);
			# $nonhelices_25{$id}{$values}{$ss_25_nh}=$ss_25_seq_nh;
			# print $out_nonhelix "$id\t$ss_25_nh\t$ss_25_seq_nh\n";
			
			$rest=(int(length($ss)/25));
			$mod=(length($ss))%25;
			$plus=0;
			for ($j=0; $j<$rest; $j++){
				$random_index = int(rand($mod));
				$ss_25_nh = substr($ss, $random_index+($plus), 25);
				$ss_25_seq_nh = substr($b, $random_index+($plus), 25);
				$helices_25{$id}{$values}{$j}{$ss_25_nh}=$ss_25_seq_nh;
				print $out_nonhelix "$id\t$j\t$ss_25_nh\t$ss_25_seq_nh\n";
				$plus=25+$random_index+$plus;
				$mod=$mod-$random_index;
			}
		}
    }
}



close IN;
close $out_helix;
close $out_nonhelix;