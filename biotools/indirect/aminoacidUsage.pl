#!/usr/bin/perl
# Authors: Tammi Camilla Vesth (tammi@cbs.dtu.dk)
# For license see /usr/biotools/CMG-biotools.license
#=============================================================================================
=put
	Synopsis	perl aminoacidUsage.pl file.proteins.fsa > file.usage
	
	The code presented here is intended to analyze the codon and amino acid usage in a genome.
	Here the genome is presented as an open reading frame fasta formated file.
	Is is assumed that the translation startsite is the first position of the open reading frame.
	The sequence is not translated. Instead the amino acid counts are done by 
	
	=====	subrutine : aa_specific_cc	
	This subrutine counts the codon usage for each amino acid. 
	The counts are amino acid specific and by summing ud the codon counts for each amino acid the amino acid count is obtained.				
=						
=cut
#=============================================================================================
#	Use
#=============================================================================================
use strict;
use warnings;
use Carp;

use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqUtils;
use Bio::Tools::IUPAC;

#=============================================================================================
#	Input and error handling
#=============================================================================================
my ($file,$help,$out,$config, $seqin);

my $USAGE = "# Usage: aminoacidUsage.pl -f [file.orf.fsa] -o [file.aminoacidUsage]
# Script to calculate the amino acid usage for a FASTA formatted open reading frame file.
# -f must be a FASTA formatted file with amino acids.
# -o output is a tab delimeted file with one line per amino acid.\n";

GetOptions(
'f:s'  => \$file,
'o:s'  => \$out,
'c:s'  => \$config
);

print "# ERROR:	Output file name not defined\n\n",$USAGE and exit if (!defined $out );
print "# ERROR:	Input file not found or can not be opened($file)\n\n", $USAGE and exit if ( ! -e $file or -z $file);

#=============================================================================================
#	Define variables
#=============================================================================================
my $table = new Bio::SeqUtils;
my @BASES = $table->valid_aa(0);
my %all = $table->valid_aa(2);

#=============================================================================================
#	Read input file
#=============================================================================================
if (-e $file){
	$seqin = new Bio::SeqIO(-format => 'fasta', -file   => $file);
} else {
	print "# ERROR:	Could not open file ($file)\n\n", $USAGE and exit;
}

#=============================================================================================
#	Calculate amino acid usage
#=============================================================================================
my %composition;
my $total;
foreach my $base ( @BASES ) {
    $composition{$base} = 0;
}
while ( my $seq = $seqin->next_seq ) {
    if( $seq->alphabet ne 'protein' ) {
	confess("Must only provide amino acid sequences to aacomp...skipping this seq");
	next;
    }
    foreach my $base ( split(//,$seq->seq()) ) {
	$composition{uc $base}++;
	$total++;
    }
}

#=============================================================================================
#	Print amino acid usage
#=============================================================================================
my @print_aa_array = ("G", "A", "V", "L", "I", "F", "Y", "W", "H", "K", "R", "D", "E", "N", "Q", "S", "T", "M", "C", "P");

#=====	Print amino acid usage	
open (OUT, ">$out") or die $!;
foreach my $print_aa_key (@print_aa_array){
	printf OUT "aa\t$print_aa_key\t%.4f\n", ($composition{$print_aa_key}/$total)*100;
}
close (OUT);

=put
printf("%d aa\n",$total); 
printf("%5s %4s\n", 'aa', '#' );
my $ct = 0;
foreach my $base ( @BASES ) {
    printf(" %s %s %3d\n", $base, $all{$base}, $composition{$base} );
    $ct += $composition{$base};
}
printf( "%6s %s\n", '','-'x5);
printf( "%6s %3d\n", '',$ct);
=cut

__END__

