#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Temp qw(tempfile);
use File::Copy;
use Cwd qw();
use lib '.'; #find in root
use Getopt::Long;
use File::Remove 'remove';
use File::Path qw(make_path remove_tree);
use Fcntl qw( LOCK_EX LOCK_NB );
use PaliSubs;
use Cwd 'abs_path';

print <<'WELCOME';   
		 ___
               ((   ))
              ((AA-AA))   
	Palindromer v0.1 [Oct 19, 2016]  - [+]

Citation - automated detection of palindromes

License: Creative Commons Licence
Bug-reports and requests to: jitendra.narayanATunamur.be and nicolas.debortoliATunamur.be

-------------------------------------------------------------------
WELCOME

my (
	$infile,
	$outfile, 
	$identity,
	$plot,
	$clean,
	$repeats,
	$strand,
	$brutal,
	$logfile,
);

# Default option setting for palindromer tool
my $VERSION=0.1;
my $verbose=0; 		# Verbose set to 0;
my %options = ();
$brutal='no';

GetOptions(
	\%options,
    	'infile|f=s'    	=> \$infile,        	## Infile
    	'outfile|o=s' 		=> \$outfile,           ## Outfile
	'identity|i=i' 		=> \$identity, 		## identity percentage
	'plot|p=s' 		=> \$plot, 		## plot dotplot with R "yes or no"
	'clean|c=s' 		=> \$clean, 		## clean result yes or no
	'strand|s=s' 		=> \$strand, 		## plot dotplot with R "minus" or "plus" or "both"
	'repeats|r=s' 		=> \$repeats, 		## yes or no
	'brutal|b=s'		=> \$brutal,
    	'help|?|h!'     	=> sub { PaliSubs::printUsage($VERSION) },
   	'who|w!'     		=> sub { PaliSubs::Who($VERSION) },
	'verbose' 		=> \$verbose,
    	'logfile=s' 		=> \$logfile,		## logfile
	
) or die PaliSubs::printUsage();

if ((!$infile) or (!$outfile) or (!$identity) or (!$strand) or (!$clean) or (!$repeats)) { 
print "ERROR: You might forgot to provide right flags\n";
PaliSubs::printUsage(); exit; }

#Check if allready running same script anywhere.
flock(DATA,LOCK_EX|LOCK_NB)
  or  die "This script ($0) is already running.\n Wait the first instances to finish and try again later\n Sorry for inconvenience\n";

#Set up the path of all required tools - stored in the same folder
#$ENV{'PATH'} = "/bin:/usr/bin:/usr/bin/env:$ENV{PWD}/lastz-distrib-1.03.73/src:$ENV{PWD}/trf";

if ($plot) { 
	print "\nWARNINGS:You have selected doplot option, which might take more computation time\n\n";
	if ($plot =~ /^[0-9,.E]+$/ ) { print "Enter yes or no in -p <yes or no> option\n"; exit;}
	}

if ($identity !~ /^[0-9,.E]+$/ ) { print "Enter numeric value in $identity -i <number>\n"; exit;}

#Remove tmpRES folder
#remove_tree( "tmpRES");
#remove_tree( "$outfile"); # Remove older for palindrome
#also delete

my $path = abs_path("$outfile");
print "Outfile in $path\n";

#To check the palindrome
print "Cheking palindrome in $infile file\n\n";
PaliSubs::palindrome($infile, $identity, $plot, $path, $strand,$clean, $repeats, $brutal);

print "\nPlease check your $path folder for result :)\n\n";


__END__

In English, the term palindrome refers to a string of letters that have the same meaning written in both directions some classic English palindrome are kayak, civic, noon, and racecar.

As a set of paired sequences (one on each of the strands of a double strand of DNA), the palindromes recognized by restriction enzymes follow a slightly different set of rules. They are probably more properly referred to as palindromic sequences to distinguish them from language palindromes. Palindromic sequences are a short run of bases (typically 3 to 5 in length), follow by their complementary bases in reverse order. For example the recognition sequence for BamHI is GGATCC.

Note the first three bases GGA are followed by the complement of those three bases in reverse order: TCC. The complement to the whole six base strand is CCTAGG, read backwards (as it would be when reading from 5’ to 3’ on the complementary strand) is GGATCC, an exact match for the original strand.

This pattern makes it possible to reconstruct a palindromic sequence from one-half of one strand. For example, a six-base recognition sequence (e.g. TAGCTA) can be reconstructed from just knowing the first three bases on one strand:

Starting with the original sequence - TAG
Calculate the reverse complement of the sequence - ATC
Reverse the order of the reverse complement (CTA) and add it to the end of the forward strand - TAGCTA
Calculate the reverse complement of the whole forward strand to finish the reverse strand: ATCGAT
Final result: 
TAGCTA
ATCGAT
Having short stretches of DNA that read the same on both strands of double-stranded DNA allow restriction enzymes to cut both strands in the same place.

Note having palindromic sequences along short the short stretches over which restriction enzymes function is not common.

To prove this take any random stretch of DNA such as AGTCCGATCCGT
find its reverse complement: TCAGGCTAGGCA
flip it to the proper 5’-3’ orientation for the complementary strand: ACGGATCGCCACT
and you don’t get the same sequence back ACGGATCGCCACT ≠ AGTCCGATCCGT



#!/usr/bin/perl
#The previous line must be at the beginning of all Perl scripts

# Prompt the user for a palindrome
print "\nPlease enter a possible palindrome: ";

# Delete the "Enter" from the end of the inputted word
chomp ( $palindrome = <STDIN> );

# Take the palindrome and put each letter into an array
@palindrome = split( //, $palindrome );

# Reverse the array of letters so it is backwards
@backwards = reverse( @palindrome );

# Find the size of the palindrome
$sizeOfPalindrome = @palindrome;

# Loop through the two arrays and see if the word is
# the same backwards as it is forward
for( $i = 0; $i < $sizeOfPalindrome; $i++ ){
	# If any of the letters do not match then
	# the loop is broken
	if( @palindrome[$i] ne @backwards[$i] ){
		$pali = "";
		last;
	}
	else{
		$pali = "true";
	}
		
} 

# If the variable $pali from above is true then all of the
# words matched and the word entered is a palindrome
if( $pali ){
	print "\n\"$palindrome\" is a palindrome\n";
}
# Otherwise the word entered was not a palindrome
else{
	print "\n\"$palindrome\" is NOT a palindrome\n";
}



#! perl -slw
use strict;

my %invert; @invert{ qw[ A C G T ] } = qw[ T G C A ];

my $in = do{ local $/; <DATA> };
chomp $in;
print $in;

for my $p1 ( 1 .. length( $in ) -2 ) {
    next unless substr( $in, $p1, 1 ) eq $invert{ substr $in, $p1+1, 1 };
    my $pals = 0;
    for my $p2 ( 1 .. $p1 -1 ) {
        last unless substr( $in, $p1-$p2, 1 ) eq $invert{ substr $in, $p1+$p2+1, 1 };
        ++$pals;
    }
    if( $pals ) {
        printf "%s%s at %d\n", ' 'x($p1-$pals),
            substr( $in, $p1-$pals, ($pals+1)*2 ), $p1-$pals;
    }
}

__DATA__
AGAGGTCAGTCTGCATCGTATCGATCGTCGACGATCGATACGATGCAGACTGACGAGAG


[13:44:54.28] C:\test>1085446
AGAGGTCAGTCTGCATCGTATCGATCGTCGACGATCGATACGATGCAGACTGACGAGAG
           TGCA at 11
                   ATCGAT at 19
                     CGATCG at 21
    GTCAGTCTGCATCGTATCGATCGTCGACGATCGATACGATGCAGACTGAC at 4
                               CGATCG at 31
                                 ATCGAT at 33
                                           TGCA at 43



