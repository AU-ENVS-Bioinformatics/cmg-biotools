#!perl -w
# Authors: Tammi Camilla Vesth (tammi@cbs.dtu.dk)
# For license see /usr/biotools/CMG-biotools.license
#===============================================================
#===============================================================

#print "Running basicInfoPS.pl\n";

#===== USE =====
use strict;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Getopt::Long;


if ($#ARGV < 2){
	print "Usage: basicInfoPS.pl file.atStats.tab file.orf.fna file.CodonAaUsage\n";
	exit;
}	

my $atStats_file = $ARGV[0];
my $orfFna_file = $ARGV[1];
my $codonUsage_file = $ARGV[2];
my @atStatsArray ;
my $count = 50;

print "%!PS\n/Times-Roman findfont\n16 scalefont\nsetfont\n";

open (atStatsFile, $atStats_file);
while (<atStatsFile>) {
 	chomp;
	@atStatsArray = split('\t', $_);
}
close (atStatsFile);

open (codonUsageFile, $codonUsage_file);
my @codonUsageArray = <codonUsageFile>;
close (codonUsageFile);

open (orfFnaFile, $orfFna_file);
my @orfFnaArray = <orfFnaFile>;
chomp $orfFnaArray[0];
close (orfFnaFile);

my $name = $atStatsArray[0];
$name =~ s/\_/ /g;
$name =~ s/prodigal.orf.fna//g;
printf "newpath\n100 %d moveto\n", (700);
print "(Identifier:\t$name) show\n";

#printf "newpath\n100 %d moveto\n", (650);
#print "(Header name of orf.fna file:\t$orfFnaArray[0]) show\n";

for (my $i = 1; $i < ($#atStatsArray); $i+=2) {
	printf "newpath\n100 %d moveto\n", (650-$count);
	print "($atStatsArray[$i]\t$atStatsArray[$i+1]) show\n";
	$count += 50;
}

my $at = 0;
my $gc = 0;

foreach my $codonUsageLine (@codonUsageArray) {
	if ($codonUsageLine =~ m/bias/) {
		my @tmp = split ('\t', $codonUsageLine);
		printf "newpath\n100 %d moveto\n", (450);
		printf "(Bias in third position:\t%.4f) show\n", $tmp[1];
		printf "newpath\n100 %d moveto\n", (430);
		print "(Bias in third position:\t100% GC = +1 and 100% AT = -1) show\n";

	}	
}	

=cut
 This looks great!  
 A couple of small suggestions - can you switch the two star plots on top, so that it is Codon usage on the left, then amino acid usage?  
 I like the position-specific nucleotide usage - this is good.  
 In that blank space to the right, can you put some text?  
 I'm thinking perhaps a couple of lines of text, with the name of the organism and other relevant information given - perhaps something like this:
  
 Analysis Date: 8 March, 2011
 
 NCBI project ID: 61623 (1 replicon)
 Accession: AM114193
 Phylum: Euryarchaeota
 Name: uncultured methanogenic archaeon RC-I
 
 Genome Size: 3,179,916 bp 
 Genome AT content: 45.4%
 
 3117 predicted proteins (Prodigal v2.5)
 Ave. AT content: 43.6% +/- 0.0534
 Total Coding Bases: 2,682,015 bp
 Fraction coding: 84.3%
 
 Bias in third position = 0.36
 Amino acid bias = 0.30
  
 Those last two numbers are calculated like I mentioned in a previous email - 
 the 'bias in third position' is merely the (%G + %C) in third position (you plot this!), minus (%A + %T).  
 Thus, if something were 100% GC, then the number would be +1, whilst if something were 100% AT, the number would be -1.  
 The second number is (16 - #aa that make up 80% of the total)/20 - I estimated this to be ~10, so that's (16-10)/20, or 0.30
=put 
