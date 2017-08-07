#!/usr/bin/perl 

use warnings;
use strict;

my ($dat, $flag, $file);

$dat = $ARGV[0];
$flag = $ARGV[1];
chomp($flag);

if (($ARGV[0] eq "") || ($ARGV[1] eq "")) {
	print "Usage: program filename flag\nPlease provide the required parameters and run again\n";
	exit;
}

my @arr = split("dat", $dat);
$file = $arr[0];


open(O1, ">".$file."dat.parse");
print  O1 ("Repeat Start\tRepeat End\tPeriod Size\tCopy No.\tAlignment Score\tConsensus\n");

system("cat ".$file."*txt.html > ".$file."tmp");

open(DAT, $dat) or die ("can not open trf data file");

while(my $line=<DAT>) {

my ($rep_start, $rep_end, $period_size, $copy_no, $pattern_size, $percent_match, $percent_indel, $align_score, $a_percent, $c_percent, $g_percent, $t_percent, $entropy, $consensus, $repeat);			
			chomp($line);
			if ($line =~ /^[0-9]/) {
			my @arr = split(" ", $line);

			$rep_start = $arr[0];
			$rep_end = $arr[1];
			$period_size = $arr[2];
			$copy_no = $arr[3];
			$pattern_size = $arr[4]; 
			$percent_match = $arr[5];
			$percent_indel = $arr[6];
			$align_score = $arr[7];
			$a_percent = $arr[8];
			$c_percent = $arr[9];
			$g_percent = $arr[10];
			$t_percent = $arr[11];
			$entropy = $arr[12];
			$consensus = $arr[13];
			$repeat = $arr[14];
			print O1 ("$rep_start\t$rep_end\t$period_size\t$copy_no\t$align_score\t$consensus\n");
			}

}


& flanking;
system("paste ".$file."dat.parse ".$file."txt.parse > ".$file."final.parse");
exit;

sub flanking {

my ($line_txt, $count, $start, $end, $left_start, $left_end, $right_start, $right_end, $left_seq, $right_seq);
$count = 0;	

open(O2, ">".$file."txt.parse" || die "can not open output file for alignment");
return if($flag==0);
print  O2 ("Left Flanking Sequence\tRight Flanking Sequence\n");

open(TXT, $file."tmp") or die ("can not open trf txt/html file");


while($line_txt=<TXT>) {

				chomp($line_txt);
				my @arr = split("[ :-]", $line_txt);

			if ($line_txt =~ "    Indices:") {
				$start = $arr[6]; 
				$end = $arr[8];
				$count = 1;
			}
			
			elsif  ($line_txt =~ "Left flanking sequence:") {
				$left_start = $arr[5]; 
				$left_end = $arr[9];
				chomp($line_txt = <TXT>);
					until ($line_txt eq "") {
						$left_seq .= $line_txt;
						chomp($line_txt = <TXT>);
						}
			}

			elsif  ($line_txt =~ "Right flanking sequence:") {
				$right_start = $arr[5]; 
				$right_end = $arr[9];
				chomp($line_txt = <TXT>);
					until ($line_txt eq "") {
						$right_seq .= $line_txt;
						chomp($line_txt = <TXT>);
						}
			}

			elsif((($line_txt =~ "Found at i:") || (eof)) && ($count == 1)){
				print  O2 ("$left_seq\t$right_seq\n");
    			$start = $end = $left_start = $left_end = $right_start = $right_end = $count = 0;
    			$left_seq = $right_seq = "";
			}


  }
}


__END__

USAGE:

 trfparser datfilename flag_value{0 or 1}
 
 Enter 0 as flag_value if you need to use only dat file, and 1 if you want to extract information from both dat and txt/html file.


AVAILABLE FIELDS:

By default, following information will be included in the .final.parse file:

.dat file: Repeat Start, Repeat End, Period Size, Copy No., Alignment Score, Consensus
.txt.html file: Left Flanking Sequence, Right Flanking Sequence

However, you can always modify the code to display other information, following is a list of all available fields and description:

.dat file

$rep_start			:		Indices of the repeat relative to the start of the sequence
$rep_end 	
$period_size		:		Period size of the repeat
$copy_no			:		Number of copies aligned with the consensus pattern
$pattern_size		:		Size of consensus pattern (may differ slightly from the period size)
$percent_match	:		Percent of matches between adjacent copies overall
$percent_indel	:		Percent of indels between adjacent copies overall
$align_score		:		Alignment score
$a_percent			:		Percent composition for each of the four nucleotides
$c_percent
$g_percent
$t_percent
$entropy			:		Entropy measure based on percent composition
$consensus		:		Consensus sequence
$repeat				:		Repeat sequence

.txt.html file

$start					:		Indices of the repeat relative to the start of the sequence
$end			
$left_start			:		Indices of the Left flanking sequence relative to the start of the sequence
$left_end		
$right_start		:		Indices of the Right flanking sequence relative to the start of the sequence
$right_end		
$left_seq				:		Left flanking sequence
$right_seq			:		Right flanking sequence

