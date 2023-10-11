#!perl -w
# Authors: Tammi Camilla Vesth (tammi@cbs.dtu.dk)
# For license see /usr/biotools/CMG-biotools.license
#===============================================================
#  This code will calculate the AT content of a genome from a DNA fasta file with one or more entries.
#  It can be used on multiple files and then the results will be put in one file
#  The information can then be visualized in a histogrem as follows:
#  histo -c 5 output.file.tab
#  Plotting the AT contetn standard deviation:
#  histo -c 3 output.file.tab
#===============================================================

#ls -1 `pwd`/CodonUsage/Data/*orfs.fsa > test.list
use strict;

use Bio::SeqIO;
use Bio::Tools::SeqStats;
#use Statistics::Descriptive;

use Getopt::Long;
my $format = 'fasta';
my $file;
my $list;
my $help =0;
GetOptions('h|help|?' => \$help);

my @files = @ARGV;

my $USAGE = "usage: gccalc.pl [path to file(s)]\nFile format should be fasta\n";
if( $help ) { die $USAGE; }

foreach my $file (@files){
	chomp $file;
	my @values_array_gc = ();
	my @values_array_at = ();
	
	$file = shift unless $file;
	my $seqin;
	if( defined $file ) {
		print "Could not open file [$file]\n$USAGE" and exit unless -e $file;
		$seqin = new Bio::SeqIO(-format => $format, -file   => $file);
	} 
	else {
		print "$USAGE";
		exit;
	}


	my ($total_base, $total_gc, $total_at);

	while( my $seq = $seqin->next_seq ) {
		next if( $seq->length == 0 );
		if( $seq->alphabet eq 'protein' ) {
			warn("$file : gccalc does not work on amino acid sequences ...skipping this seq");
			exit;
		}
	
		my $seq_stats  =  Bio::Tools::SeqStats->new('-seq'=>$seq);
		my $hash_ref = $seq_stats->count_monomers();  # for DNA sequence
		
		$total_base += $seq->length;
		$total_gc += $hash_ref->{'G'} + $hash_ref->{'C'};
		$total_at += $hash_ref->{'A'} + $hash_ref->{'T'};
		push(@values_array_gc, ($hash_ref->{'G'} + $hash_ref->{'C'}) / $seq->length());
		push(@values_array_at, ($hash_ref->{'A'} + $hash_ref->{'T'}) / $seq->length());

	#	Make active to print numbers for each sequence
=put	
	print $seq->display_id , "\t";
	printf "GC %.5f\t", ($hash_ref->{'G'} + $hash_ref->{'C'}) /$seq->length();
	printf "AT %.5f\t", ($hash_ref->{'A'} + $hash_ref->{'T'}) /$seq->length();
	printf "%.5f\t", (($hash_ref->{'G'} + $hash_ref->{'C'}) /$seq->length()) + (($hash_ref->{'A'} + $hash_ref->{'T'}) /$seq->length());
	print "\n";
=cut 
	}	
	
	#printf "$file\tStDevAT:\t%.4f\tPerAT:\t%.4f\tTotalBases:\t%d\tValuesArray\t%d\n", standard_deviation(@values_array_at), ($total_at / $total_base)*100, $total_base, scalar @values_array_at;
	#printf "$file\tStDevAT:\t%.4f\tPerAT:\t%.4f\tTotalBases:\t%d\n", standard_deviation(@values_array_at), ($total_at / $total_base)*100, $total_base;
	printf "$file\tTotalBases:\t%d\tPerAT:\t%.2f\tStDevAT:\t%.2f\n", $total_base,($total_at / $total_base)*100,standard_deviation(@values_array_at);

}

#	Make active to print more infor for genome
=put
print "St.dev GC :\t", standard_deviation(@values_array_gc)."\n";
print "St.dev AT :\t", standard_deviation(@values_array_at)."\n";
printf "Total GC :\t%.4f out of %d bases\n", $total_gc / $total_base;
printf "Total AT :\t%.4f out of %d bases\n", $total_at / $total_base;
printf "Total bases :\t%d\n", $total_base;
=cut

sub standard_deviation {
	my(@numbers) = @_;
	#Prevent division by 0 error in case you get junk data
	return undef unless(scalar(@numbers));
	
	# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
		$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);
	
	# Step 2, find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	foreach my $num (@numbers) {
		$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);
	
	# Step 3, standard deviation is the square root of the
	# above mean
	my $std_dev = sqrt($mean2);
	return $std_dev;
}


__END__

=head1 NAME
 
 gccalc - GC content of nucleotide sequences
 
 =head1 SYNOPSIS
 
 gccalc [-f/--format FORMAT] [-h/--help] filename
 or
 gccalc [-f/--format FORMAT] < filename
 or
 gccalc [-f/--format FORMAT] -i filename
 
 =head1 DESCRIPTION
 
 This scripts prints out the GC content for every nucleotide sequence
 from the input file.
 
 =head1 OPTIONS
 
 The default sequence format is fasta.
 
 The sequence input can be provided using any of the three methods:
 
 =cut
