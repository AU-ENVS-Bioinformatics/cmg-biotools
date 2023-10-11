#!/usr/bin/perl
# Authors: Peter Fisher Hallin, Tammi Camilla Vesth (tammi@cbs.dtu.dk)
# For license see /usr/biotools/CMG-biotools.license
#=============================================================================================
=put
	Synopsis	perl CodonAaUsage.pl file.orf.fsa > file.usage
	
	The code presented here is intended to analyze the codon and amino acid usage in a genome.
	Here the genome is presented as an open reading frame fasta formated file.
	Is is assumed that the translation startsite is the first position of the open reading frame.
	The sequence is not translated. Instead the amino acid counts are done by 
	
	=====	subrutine :  codon_counts
	This subrutine counts the codon usage for the entire fasta file = whole genome.
	The method Bio::Tools::SeqStats->count_codons was used. 	
	This method counts codons in ONE sequence. 
	It was used in a way so that it sums up for all sequences = whole genome
						
	=====	subrutine : aa_specific_cc	
	This subrutine counts the codon usage for each amino acid. 
	The counts are amino acid specific and by summing ud the codon counts for each amino acid the amino acid count is obtained.				
=cut
#=============================================================================================
#	Use
#=============================================================================================
use strict;
use Carp;
use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqUtils;
use Bio::Tools::IUPAC;
use warnings;
use Bio::Perl;
use Bio::Tools::SeqStats;

#=============================================================================================
#	Obtain file information
#=============================================================================================

my $file;

if ($#ARGV == 0) {
	$file = $ARGV[0];
	chomp $file;
}

else { 
	print "perl CodonAaUsage.pl file.orf.fsa > file.usage\n";
	exit;
} 

#=============================================================================================
#	MAIN
#=============================================================================================
#=====	hash to determine the printing order of amino acids. ===== 
#=====	Nessesary since we order amino acids based on biochemical features =====
my @print_aa_array = ("G", "A", "V", "L", "I", "F", "Y", "W", "H", "K", "R", "D", "E", "N", "Q", "S", "T", "M", "C", "P");

#=====	hash to determine the printing order of codons sorted by last base, then seqound and then first position. =====
my @print_codon_array =("AAA", "CAA", "GAA", "TAA", "ACA", "CCA", "GCA", "TCA", "AGA", "CGA", "GGA", "TGA", "ATA", "CTA", "GTA", "TTA", 
			"AAC", "CAC", "GAC", "TAC", "ACC", "CCC", "GCC", "TCC", "AGC", "CGC", "GGC", "TGC", "ATC", "CTC", "GTC", "TTC", 
			"AAG", "CAG", "GAG", "TAG", "ACG", "CCG", "GCG", "TCG", "AGG", "CGG", "GGG", "TGG", "ATG", "CTG", "GTG", "TTG", 
			"AAT", "CAT", "GAT", "TAT", "ACT", "CCT", "GCT", "TCT", "AGT", "CGT", "GGT", "TGT", "ATT", "CTT", "GTT", "TTT");
my @print_pos_count = ("A","T","G","C");

my %temp_codon_counts = &codon_counts($file);
my %temp_results_hash = &aa_specific_cc(%temp_codon_counts);
my $temp_total_codon_count = &total_codon_count(%temp_codon_counts);
my %temp_pos_count_hash_of_hashes = &pos_codoncount(%temp_codon_counts);

#===================================================================
#=====	Print codon usage =====
foreach my $print_codon_key (@print_codon_array){
	#printf "$print_codon_key\t%.4f\n", $temp_results_hash{$print_codon_key};
	printf "codon\t$print_codon_key\t%.5f\t%d\n", ($temp_codon_counts{$print_codon_key}/$temp_total_codon_count)*100, $temp_codon_counts{$print_codon_key};

}
#=====	Print position specific nucleotide counts =====	
print "freq\tpos\tA\tA\tT\tT\tG\tG\tC\tC\n";
my $bias = 0;
foreach my $pos (sort keys %temp_pos_count_hash_of_hashes){
	print "freq\t$pos\t";
	foreach my $print_pos_key (@print_pos_count) {
		printf "$print_pos_key\t%.5f\t", ($temp_pos_count_hash_of_hashes{$pos}{$print_pos_key})/$temp_total_codon_count; 
		#=====	calculate bias in third position
		if ($pos =~ m/pos3/) {
			my $at = ($temp_pos_count_hash_of_hashes{"pos3"}{"A"}/$temp_total_codon_count)+($temp_pos_count_hash_of_hashes{"pos3"}{"T"}/$temp_total_codon_count);
			my $gc = ($temp_pos_count_hash_of_hashes{"pos3"}{"G"}/$temp_total_codon_count)+($temp_pos_count_hash_of_hashes{"pos3"}{"C"}/$temp_total_codon_count);
			$bias = $gc-$at;
		}	
	}	
	print "\n";		
}
#=====	Print bias on third position
printf "biasIn3pos\t%.4f\n", $bias;

#=====	Print amino acid usage	
foreach my $print_aa_key (@print_aa_array){
	printf "aa\t$print_aa_key\t%.4f\n", $temp_results_hash{$print_aa_key}*100;
}	



#=============================================================================================
#	SUBRUTINE : total codon count of fasta file/genomes
#=============================================================================================
sub total_codon_count
{
	my %tcc_codon_counts = @_;
	my $tcc_total_genome = 0;
	my %tcc_results_hash;
	#=====	Count codon usage in entire genome =====
	while (my ($tcc_k, $tcc_v) = each(%tcc_codon_counts) ) {
		$tcc_total_genome += $tcc_v;
	}
	return $tcc_total_genome;
}
#=============================================================================================
#	SUBRUTINE : position specific nucleotide usage
#=============================================================================================
sub pos_codoncount
{
	my %pcc_codon_counts = @_;
	my %pos_count_hash_of_hashes = (
	'pos1' => {'A' => 0, 'T' => 0,'C' => 0, 'G' => 0},
	'pos2' => {'A' => 0, 'T' => 0,'C' => 0, 'G' => 0},
	'pos3' => {'A' => 0, 'T' => 0,'C' => 0, 'G' => 0});
	
	
	#=====	Go trough each codon in codon counting hash =====
	foreach my $pcc_key (sort keys %pcc_codon_counts){
			
		#=====	Go trough each position (1,2,3) =====
		foreach my $pos (sort keys %pos_count_hash_of_hashes){
			my $substr = '';
		
			if ($pos =~ 'pos1' ){
				$substr = substr($pcc_key,0,1);
				$pos_count_hash_of_hashes{$pos}{$substr} = ($pcc_codon_counts{$pcc_key})+$pos_count_hash_of_hashes{$pos}{$substr};	
			}
			if ($pos =~ 'pos2' ){
				$substr = substr($pcc_key,1,1);
				$pos_count_hash_of_hashes{$pos}{$substr} = ($pcc_codon_counts{$pcc_key})+$pos_count_hash_of_hashes{$pos}{$substr};	
			}
			if ($pos =~ 'pos3' ){
				$substr = substr($pcc_key,2,1);
				$pos_count_hash_of_hashes{$pos}{$substr} = ($pcc_codon_counts{$pcc_key})+$pos_count_hash_of_hashes{$pos}{$substr};	
			}
		}
	}
	return %pos_count_hash_of_hashes;
}

#=============================================================================================
#	SUBRUTINE : aa specific codon counts
#=============================================================================================

sub aa_specific_cc
{
	my %aascc_codon_counts = @_;
	my $aascc_total_genome = 0;
	my %results_hash;
	#=====	Count codon usage in entire genome =====
	while (my ($aascc_k, $aascc_v) = each(%aascc_codon_counts) ) {
		$aascc_total_genome += $aascc_v;
	}
	
	#===== Define hash of hashes =====
	my %aascc_hash_of_hashes = (
	'A' => {'aa'=>0,'GCT'=>0, 'GCG'=>0, 'GCC'=>0, 'GCA'=>0},
	'C' => {'aa'=>0,'TGT'=>0,'TGC'=>0},
	'D' => {'aa'=>0,'GAT'=>0,'GAC'=>0},
	'E' => {'aa'=>0,'GAG'=>0,'GAA'=>0},
	'F' => {'aa'=>0,'TTT'=>0,'TTC'=>0},
	'G' => {'aa'=>0,'GGT'=>0,'GGG'=>0,'GGC'=>0,'GGA'=>0},
	'H' => {'aa'=>0,'CAT'=>0,'CAC'=>0},
	'K' => {'aa'=>0,'AAG'=>0,'AAA'=>0},
	'L' => {'aa'=>0,'TTG'=>0,'TTA'=>0,'CTT'=>0,'CTG'=>0,'CTC'=>0,'CTA'=>0},
	'I' => {'aa'=>0,'ATA'=>0,'ATC'=>0,'ATT'=>0,},
	'M' => {'aa'=>0,'ATG'=>0},
	'N' => {'aa'=>0,'AAT'=>0,'AAC'=>0},
	'P' => {'aa'=>0,'CCT'=>0,'CCG'=>0,'CCC'=>0,'CCA'=>0},
	'Q' => {'aa'=>0,'CAG'=>0,'CAA'=>0},
	'R' => {'aa'=>0,'CGT'=>0,'CGG'=>0,'CGC'=>0,'CGA'=>0,'AGG'=>0,'AGA'=>0},
	'S' => {'aa'=>0,'TCT'=>0,'TCG'=>0,'TCC'=>0,'TCA'=>0,'AGT'=>0,'AGC'=>0},
	'T' => {'aa'=>0,'ACT'=>0,'ACG'=>0,'ACC'=>0,'ACA'=>0},
	'V' => {'aa'=>0,'GTT'=>0,'GTA'=>0,'GTC'=>0,'GTG'=>0},
	'W' => {'aa'=>0,'TGG'=>0},
	'Y' => {'aa'=>0,'TAT'=>0,'TAC'=>0},	
	'stop' => {'aa'=>0,'TAG'=>0,'TGA'=>0,'TAA'=>0}
	);
	
	#=====	Put codon count values into the hash of hashes =====
	#=====	Goes trough each codon in codon_count =====
	foreach my $aascc_codon ( keys %aascc_codon_counts ){ 
		foreach my $k1 ( keys %aascc_hash_of_hashes ){ 
			if (exists $aascc_hash_of_hashes{$k1}{$aascc_codon}) {
				$aascc_hash_of_hashes{$k1}{$aascc_codon} = $aascc_codon_counts{$aascc_codon};
			}
		}	
	}
	
	my $aascc_total = 0;
	#=====	Go trough each amino acid in hash of hashes =====
	foreach my $k ( sort keys %aascc_hash_of_hashes ){		
		#=====	Calculate total codon count for amino acid =====
		foreach my $k2 ( sort keys %{$aascc_hash_of_hashes{$k}}){	
			$aascc_total+= $aascc_hash_of_hashes{$k}{$k2} unless $k2 =~ /aa/;	
		}
		$aascc_hash_of_hashes{$k}{aa} = ($aascc_total/$aascc_total_genome);
		$aascc_total = 0;
	}	
	my $aascc_total2 = 0;
	
	#=====	Print amino acid specific codon usage =====
	foreach my $f (sort keys %aascc_hash_of_hashes ){
		
		foreach my $f2 (sort keys %{$aascc_hash_of_hashes{$f}}){	
			$aascc_total2+= $aascc_hash_of_hashes{$f}{$f2} unless $f2 =~ /aa/;	
		}
		
		foreach my $f1 (sort keys %{$aascc_hash_of_hashes{$f}}){
			if ($aascc_total2 != 0) {
				($results_hash{$f1} = ($aascc_hash_of_hashes{$f}{$f1}/$aascc_total2)) unless $f1 =~ /aa/;
				#printf "$f1\t%.4f\n", $aascc_hash_of_hashes{$f}{$f1}/$aascc_total2 unless $f1 =~ /aa/;
			}
			else{
				($results_hash{$f1} =  0.0000) unless $f1 =~ /aa/;
				#print "$f1\t0.0000\n" unless $f1 =~ /aa/;
			}
		}	
		$aascc_total2 = 0;
	}	
	#=====	Save amino acid usage in result hash =====
	foreach my $m ( sort keys %aascc_hash_of_hashes ){
		($results_hash{$m} =  ($aascc_hash_of_hashes{$m}{aa})) ;
		#printf "$m\t%.4f\n" , $aascc_hash_of_hashes{$m}{aa};
	}
	return %results_hash;
}	

#=============================================================================================
#	SUBRUTINE : codon counts
#=============================================================================================
sub codon_counts
{	
	my $cc_file = $_[0];
	my $cc_dna_seq = new Bio::SeqIO(-format => "fasta", -file => $cc_file);
	my %cc_genome_hash = ();
	
	#=====	The reason for making this myself is that the BioPerl codon count only print the codons it finds =====
	#=====	For comparative purposes I need ALL codons for each file! =====
	my(%global_codon_counts) = (
	'TTT'=>0, 'TTG'=>0, 'TTC'=>0, 'TTA'=>0, 'TGT'=>0, 'TGG'=>0, 'TGC'=>0, 'TGA'=>0, 'TCT'=>0, 'TCG'=>0, 
	'TCC'=>0, 'TCA'=>0, 'TAT'=>0, 'TAG'=>0, 'TAC'=>0, 'TAA'=>0, 'GTT'=>0, 'GTG'=>0, 'GTC'=>0, 'GTA'=>0, 
	'GGT'=>0, 'GGG'=>0, 'GGC'=>0, 'GGA'=>0, 'GCT'=>0, 'GCG'=>0, 'GCC'=>0, 'GCA'=>0, 'GAT'=>0, 'GAG'=>0, 
	'GAC'=>0, 'GAA'=>0, 'CTT'=>0, 'CTG'=>0, 'CTC'=>0, 'CTA'=>0, 'CGT'=>0, 'CGG'=>0, 'CGC'=>0, 'CGA'=>0, 
	'CCT'=>0, 'CCG'=>0, 'CCC'=>0, 'CCA'=>0, 'CAT'=>0, 'CAG'=>0, 'CAC'=>0, 'CAA'=>0, 'ATT'=>0, 'ATG'=>0, 
	'ATC'=>0, 'ATA'=>0, 'AGT'=>0, 'AGG'=>0, 'AGC'=>0, 'AGA'=>0, 'ACT'=>0, 'ACG'=>0, 'ACC'=>0, 'ACA'=>0, 
	'AAT'=>0, 'AAG'=>0, 'AAC'=>0, 'AAA'=>0); 
		
	while (my $cc_seq = $cc_dna_seq->next_seq) {
		
		#======	Create a sequence object from each fasta entry =====	
		my $seqobj = Bio::PrimarySeq->new(-seq=>($cc_seq->seq), -alphabet=>'dna', -id=>'test');
		#======	Run count_codons and store results in one hash per sequence object =====	
		my $cc_hash_ref = Bio::Tools::SeqStats->count_codons($seqobj);
		
		#======	Go trough each codon in hash for each sequence =====
		foreach my $cc_dna_codon (sort keys %$cc_hash_ref) {
			#====== If codon already found add count to value of codon =====
			if (exists $cc_genome_hash{$cc_dna_codon}){
				$cc_genome_hash{$cc_dna_codon} += $cc_hash_ref->{$cc_dna_codon};
			}
			
			#========== If codon not already found, add codon and count as value =====
			else{ $cc_genome_hash{$cc_dna_codon} = $cc_hash_ref->{$cc_dna_codon}; }
		}		
	}
	while ( my ($cc_k, $cc_v) = each(%cc_genome_hash) ) {
		if (exists $global_codon_counts{$cc_k}) { $global_codon_counts{$cc_k} = $cc_v;}
	}
	return (%global_codon_counts);
}	
