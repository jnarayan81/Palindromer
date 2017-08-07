# $Id$ All messages
# Perl module for Palidromer PoliSubs;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2016 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CommonSubs  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jitendra <jnarayan81@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package PaliSubs;
use warnings::register;
#use Statistics::R;
#use diagnostics;

use Exporter qw(import);
 
our @EXPORT_OK = qw(yesORno deldir);


sub palindrome {
my($seqfile, $identity, $plot, $outfile) = @_;
use File::Temp qw(tempfile);
use File::Copy;
use File::Remove 'remove';

mkdir("$ENV{PWD}/$outfile") or die $!;
if ($plot eq "yes") { mkdir("$ENV{PWD}/$outfile/dotplotRES") or die $!;}
#mkdir('tmpRES') or die $!;
# It read one fasta sequence from multi fasta file and lastz it. Store the data in referece sequence fasta name/ID.

local $/ = "\n>";  # read by FASTA record
#my @myids=('scaffold1763size26795A','scaffold1763size26795','scaffold3384size10299','scaffold983size82572','scaffold2310size31877A','scaffold2310size31877','scaffold2830size13842','scaffold2691size24698','scaffold1182size56708','scaffold2502size16882','scaffold3486size20253','scaffold2330size18741','scaffold1662size56991','scaffold698size60806A','scaffold698size60806','scaffold1607size29713','scaffold1387size51584','scaffold1277size37928','scaffold1202size40069','scaffold1150size80769','scaffold1064size52433','scaffold158size147409A','scaffold158size147409','scaffold123size173649','scaffold50size285965');
my $seqCnt=0;
my $paliCnt=0;
my @allGC; my @allpalGC;
open FASTA, $seqfile;
while (<FASTA>) {
    chomp;
    my $seq = $_;
    my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
        $seq =~ s/^>*.+\n//;  # remove FASTA header
        $seq =~ s/\n//g;  # remove endlines

	$seqCnt++;
	my $GC=PaliSubs::process_contig($id, $seq);
	push @allGC, $GC;
        # Recursively check if string is palindrome my $aaa=palindrome_r($seq); print "$aaa\n";
	# Remove the file when the reference goes away with the UNLINK option
	my $tmp_fh = new File::Temp( UNLINK => 1 );
	#Write the sequence in tmp file
	print $tmp_fh ">$id\n$seq\n";
print "Looking for palindrome in $id :";
#next if length($seq) < 10000;
        #print "$id\n$seq\n";
	#if ( grep( /^$id$/, @myids ) ) {

	#LASTZ SETTINGS and RUN
	my $myLASTZ="lastz $tmp_fh $tmp_fh --chain --output=$outfile/$seqfile-$id.tmpal --format=general- --progress --ambiguous=iupac --identity=$identity --strand=minus";
        system ("$myLASTZ");
	#}
	#Check only of there is some reported palindrome
	if (-s "$outfile/$seqfile-$id.tmpal") {
   		#system ("dotter $tmp_fh $tmp_fh");
		copy("$tmp_fh", "$seqfile-sequence_$id.palfs");
		$paliCnt++;
		my $palGC=PaliSubs::process_contig($id, $seq);
		push @allpalGC, $palGC;
		#Plot DOTPLOT with R
		if ($plot eq "yes") { dotplotR ($tmp_fh, $identity, $id, $outfile);}
		#TRF SETTINGS and RUN
		system ("./trf/trf409.linux64 $seqfile-sequence_$id.palfs 2 5 7 80 10 50 2000 -l 6 -d");
		#File extension append the parameters at end ... Need to update it later
		#system ("perl trfparser_v1.pl $seqfile-sequence_$id.palfs.2.5.7.80.10.50.2000.dat , 1");
		#process_di_tri($id, $seq, 1);
		reformatTRF("$seqfile-sequence_$id.palfs.2.5.7.80.10.50.2000.dat" , 1, "$id");
		my @allfs_files = glob ('*.palfs');
		moveFiles(\@allfs_files,"$outfile");
		my @allfp_files = glob ('*.final.parse');
		moveFiles(\@allfp_files,"$outfile");
		remove( '*.html', '*.parse', '*.tmp','*.dat' );
	}

#last; # Terminate and check the first seq result
}

#contactinate all palindromic hits files
print "Deleting size zero files\n";
system("find $outfile -size 0 -delete"); # Delete size zero file


print "\nConcatinating ----------\n";
#system ("perl -please $outfile/*.tmpal > Palindrome.palfc");
system ("cd $outfile; for a in *.tmpal ; do cat $a >> Palindrome.palfc ; done; cd ..");

#my $outfile = shift;
#open OUT,">".$outfile or die "Could not open $outfile:$!\n";
#print OUT $_ while <>;

#system ("perl -e 'open OUT,">".shift;print OUT <>' output *.txt");

my @allfc_files = glob ('*.palfc');
moveFiles(\@allfc_files,"$outfile");

print  "\nDeleting .tmpal files\n";
#remove("$outfile/*.tmpal"); #delete all
close FASTA;

#Final STAT
printLines(50,'*');
my $stat=$paliCnt*100/$seqCnt;
my $sum = [map {$sum +=$_} @allpalGC]->[$#allpalGC];
my $avgPalGC=$sum/scalar(@allpalGC);
my $sumAll = [map {$sumAll +=$_} @allGC]->[$#allGC];
my $avgAllGC=$sumAll/scalar(@allGC);

print "FINAL STATS: 
Total sequence checked: $seqCnt 
Total palindrome sequence count: $paliCnt 
Percent of sequences are in palindrome: $stat%
Average GC of palindromic sequence: $avgPalGC% vs $avgAllGC% all\n";

printLines(50,'*');
}

sub printLines {
my ($num,$sign)=@_;
my $seplines = ($sign x $num)."\n";
print "\n$seplines";
}

sub dotplotR {
my ($fileD, $idt, $id, $outfile) = @_;
print "$fileD, $idt, $id, $outfile\n";

my $myDOTPLOT="lastz $fileD $fileD --chain --output=$outfile/dotplotRES/$id.dotplot --format=rdotplot --progress --ambiguous=iupac --identity=$idt --strand=minus";
system ("$myDOTPLOT");
system ("Rscript Rscripts/LASTZdotplot.R --title $id-dotplot -s 1000 -r 100 -o $outfile/dotplotRES/seeDotplot-$id $outfile/dotplotRES/$id.dotplot");

}

#subroutine to check tri di bases
sub process_di_tri {
  my($contig, $sequence, $extended) = @_;
  my $len = length($sequence);

  # count mono-nucleotides and total
  $x_counters{'total1'} += $len;
  my $newa = ($sequence =~ tr/a/a/);
  $x_counters{'a'} += $newa;
  my $newc = ($sequence =~ tr/c/c/);
  $x_counters{'c'} += $newc;
  my $newg = ($sequence =~ tr/g/g/);
  $x_counters{'g'} += $newg;
  my $newt = ($sequence =~ tr/t/t/);
  $x_counters{'t'} += $newt;
  $x_counters{'n'} += $len - ($newa + $newc + $newg + $newt);

  # count di- and tri-nucleotides
  for ($i = 0; $i < $len - 2; $i++) {
    print STDERR "$contig ... $i\n"
      if (($i % 500000) == 0) && ($i != 0);
    $x_counters{'total2'}++;
    $x_counters{'total3'}++;
    my $tri = substr($sequence, $i, 3);
    my $di = substr($tri, 0, 2);
    if ($di !~ /n/) {
      $x_counters{$di}++;
      }
    else {
      $x_counters{'nn'}++;
      }
    if ($tri !~ /n/) {
      $x_counters{$tri}++;
      }
    else {
      $x_counters{'nnn'}++;
      }
    } # end for ($i = 0; $i < $len - 2; $i++)
  if ($len > 1)		# take care of last di-nucleotide
    {
    $x_counters{'total2'}++;
    my $last2 = substr($sequence, -2, 2);
    if ($last2 !~ /n/) {
      $x_counters{$last2}++;
      }
    else {
      $x_counters{'nn'}++;
      }
    }

  # remove non-ACGT from counts unless requested
  unless ($extended) {
    $x_counters{'total1'} -= $x_counters{'n'};
#    $x_counters{'n'} = 0;
    $x_counters{'total2'} -= $x_counters{'nn'};
    $x_counters{'nn'} = 0;
    $x_counters{'total3'} -= $x_counters{'nnn'};
    $x_counters{'nnn'} = 0;
    }

  #if ($print) {
    #print STDOUT "\n$CONTIG=$contig\n\n" unless ($TOTALS > 1);
    #print_extended(%x_counters) unless ($TOTALS > 1);
    # update counters
print "$x_counters{'total1'} - $x_counters{'total2'} - $x_counters{'total3'} + $x_counters{'n'} - $x_counters{'nn'} -  $x_counters{'a'} - $x_counters{'c'} - $x_counters{'g'} - $x_counters{'t'}";
    #$x_totals{$_} += $x_counters{$_} foreach keys %x_counters;
    #%x_counters = map { ($_, 0) } keys %x_counters;
    #}
} # end do tri caclculation


###FUNCTION TO MOVE THE FILES WITH WILD
sub moveFiles {
    my ( $source_ref, $arc_dir ) = @_;
    my @old_files = @$source_ref;
   foreach my $old_file (@old_files)
         {
    #my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    #my $new_file = $arc_dir . $short_file_name;
    move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
   }
}


sub process_contig
  {
  my($contig, $sequence) = @_;
  my($len, $ACGTbases, $ATbases, $GCbases, $nonACGTbases);

  # Remove Contig name prefix?
  $contig =~ s/^.*([Cc]ontig)/$1/ if $SHORTEN_CONTIG_NAMES;
  $len = length($sequence);
  push @CONTIG_LENGTHS, $len;
  $Total_Bases += $len;
  $Max_Bases = $len if $Max_Bases < $len;
  $Min_Bases = $len if $Min_Bases > $len || $Min_Bases < 0;

  $ATbases = ($sequence =~ tr/aAtT/aAtT/);
  $GCbases = ($sequence =~ tr/cCgG/cCgG/);
  $ACGTbases = $ATbases + $GCbases;
  $nonACGTbases = $len - $ACGTbases;
  if ($ACGTbases)
    {
    $GC_per_cent = sprintf "%.1f", 100 * $GCbases / $ACGTbases;
    }
  else
    {
    $GC_per_cent = '-';
    }
  $Total_GC += $GCbases;
  $Total_ACGT += $ACGTbases;
  if ($nonACGTbases)
    {
    my $more_Max_Nons = ($Max_Nons < $MAX_PATTERN_MIN_RPT) ?
      $Max_Nons + 1 : $MAX_PATTERN_MIN_RPT;
    my @Nons = ($sequence =~ /[^acgtACGT]{$more_Max_Nons,}/g);
    foreach (@Nons)
      {
      my $l = length $_;
      $Max_Nons = $l if ($Max_Nons < $l);
      }
    $Total_Non_ACGT_Ends += length $1 if ($sequence =~ /^([^acgtACGT]+)/);
    if (substr($sequence, -1) =~ /[^acgtACGT]+$/)
      {
      my $rs = reverse $sequence;
      $Total_Non_ACGT_Ends += length $1 if ($rs =~ /^([^acgtACGT]+)[acgtACGT]/);
      }
    my $more_Max_Ns = ($Max_Ns < $MAX_PATTERN_MIN_RPT) ?
      $Max_Ns + 1 : $MAX_PATTERN_MIN_RPT;
    my @Ns = ($sequence =~ /[nN]{$more_Max_Ns,}/g);
    foreach (@Ns)
      {
      my $l = length $_;
      $Max_Ns = $l if ($Max_Ns < $l);
      }
    $Total_N_Ends += length $1 if ($sequence =~ /^([nN]+)/);
    if (uc substr($sequence, -1) eq 'N' && uc substr($sequence, 0, 1) ne 'N')
      {
      my $rs = uc reverse $sequence;
      $Total_N_Ends += length $1 if ($rs =~ /^(N+)/);
      }
    }
    #return ($contig, $len, $GC_per_cent, $nonACGTbases);
    return $GC_per_cent;
  } # end process_contig


sub flanking {
my ($file, $flag)=@_;
my ($line_txt, $count, $start, $end, $left_start, $left_end, $right_start, $right_end, $left_seq, $right_seq);
$count = 0;	

my $fileN="$file"."txt.parse";
open (my $fh2, '>', $fileN) or die "Could not open file $fileN $!";
return if ($flag==0);
print $fh2 ("Left Flanking Sequence\tRight Flanking Sequence\tsequence_id\n");
local $/ = "\n";
open(TXT, '<', $file."tmp") or die "Can not open trf txt/html file";
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
		print  $fh2 ("$left_seq\t$right_seq\t$id\n");
    		$start = $end = $left_start = $left_end = $right_start = $right_end = $count = 0;
    		$left_seq = $right_seq = "";
	}
  }
close TXT;
close $fh2;
}


sub reformatTRF {
my ($dat, $flag, $id)=@_;

if (($dat eq "") || ($flag eq "")) {
	print "Usage: program filename flag\nPlease provide the required parameters and run again\n";
	exit;
}
my @arr = split("dat", $dat); my $file = $arr[0];
my $filename="$file"."dat.parse";
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print  $fh ("Repeat Start\tRepeat End\tPeriod Size\tCopy No.\tAlignment Score\tConsensus\tsequence_id\n");

system("cat ".$file."*txt.html > ".$file."tmp");

local $/ = "\n";
open (DATA, '<', $dat) or die "Can not open TRF data file";
while(my $line=<DATA>) {
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
	print $fh ("$rep_start\t$rep_end\t$period_size\t$copy_no\t$align_score\t$consensus\t$id\n");
	}

}
close $fh;
close DATA;
flanking($file, $flag, $id);
system("paste ".$file."dat.parse ".$file."txt.parse > ".$file."final.parse");

=pod
USAGE:

 trfparser datfilename flag_value{0 or 1}
 
 Enter 0 as flag_value if you need to use only dat file, and 1 if you want to extract information from both dat and txt/html file.


AVAILABLE FIELDS:

By default, following information will be included in the .final.parse file:

.dat file: Repeat Start, Repeat End, Period Size, Copy No., Alignment Score, Consensus
.txt.html file: Left Flanking Sequence, Right Flanking Sequence

However, you can always modify the code to display other information, following is a list of all available fields and description:

.dat file

$rep_start	:		Indices of the repeat relative to the start of the sequence
$rep_end 	
$period_size	:		Period size of the repeat
$copy_no	:		Number of copies aligned with the consensus pattern
$pattern_size	:		Size of consensus pattern (may differ slightly from the period size)
$percent_match	:		Percent of matches between adjacent copies overall
$percent_indel	:		Percent of indels between adjacent copies overall
$align_score	:		Alignment score
$a_percent	:		Percent composition for each of the four nucleotides
$c_percent
$g_percent
$t_percent
$entropy	:		Entropy measure based on percent composition
$consensus	:		Consensus sequence
$repeat		:		Repeat sequence

.txt.html file

$start		:		Indices of the repeat relative to the start of the sequence
$end			
$left_start	:		Indices of the Left flanking sequence relative to the start of the sequence
$left_end		
$right_start	:		Indices of the Right flanking sequence relative to the start of the sequence
$right_end		
$left_seq	:		Left flanking sequence
$right_seq	:		Right flanking sequence

=cut

}


###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
  print LOG $message;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}


sub Who {
my $VERSION = shift;
print "\nPalindromer, version $VERSION by Jit and Nico\n"; 

}

# This function prints the script usage mode
sub printUsage {
my ($message) = @_;
if (defined $message) {
    print STDERR "\nERROR: $message\n";
  }

print << "End_Print_Usage";

Usage:

perl Palindromer.pl -f sampleSeq.fa -o TESTOUT2 -i 90

Mandatory parameters: 
  -f <infile> 	Provide the contigs/scaffods file.

  -o <outfolder> 	Provide the outfolder. 

  -i <num> 		Identify percentage for lastz.

  -p <yes or no>	Plot the R dotplot

End_Print_Usage

}


1;


