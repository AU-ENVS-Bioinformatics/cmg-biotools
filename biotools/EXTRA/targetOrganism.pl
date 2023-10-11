#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Data::Dumper;
use Bio::Tools::RNAMotif; # or anything that gives start, end, strand info
use Bio::DB::GenBank;
use Getopt::Long;
use List::Util qw(first);

#========================================================================
# Print information to user
#========================================================================
my $ProgramName  = "targetOrganism.pl";

#print "#--------------------------------------------------\n";
#print "#------------ PROGRAM NAME: -----------------------------\n"; 
#print "#--------------------------------------------------\n";
#print "# $ProgramName \n";

#========================================================================
# Print usage
#========================================================================
sub print_usage { 
	print "#--------------------------------------------------\n";
	print "# USAGE:\n"; 
	print "#--------------------------------------------------\n";
	print "# perl $ProgramName -gi gi.list\n"; 
	print "# -gi, list of NCBI GI id number\n"; 
	print "#--------------------------------------------------\n";
	exit;
}

#========================================================================
# Get options
#========================================================================
my ($gilist, $out);

&GetOptions ("gi:s" =>  \$gilist);

unless (defined $gilist) {
	print "#--------------------------------------------------\n";
	print "# ERROR:\n";
	print "#--------------------------------------------------\n"; 
	print "# List of NCBI GI number not defined\n";
	print print_usage();  
}

#========================================================================
# Read file with GI numbers
#========================================================================
my (@gis, $gi);
open (FILE_gi, $gilist) or die "Couldn't open location file: $!";
while ($gi = <FILE_gi>) {
	chomp $gi;
    push(@gis, $gi);
}
close FILE_gi;

#========================================================================
# MAIN
#========================================================================
my $y = scalar(@gis);
my $x = 0;

open (rerunFILE, '>rerun_gi.txt');

for my $x_gi (@gis) {
	$x++;
	print "#------------------------------------------------------------------\n";
	print "# PROGRESS:\n"; 
	print "# Processing protein with GI number: $x_gi :: sequence number $x of $y\n";

	#my @results = ;
	#print Dumper(@results);
	if ((glob("*$x_gi*fsa"))[0] =~ m/$x_gi/) {
		print "#------------------------------------------------------------------\n";
		print "# WARNING:\n"; 
		print "# Files with GI number: $x_gi already exists. Delete file and re-run\n";	
		print "# NOTE: reoccuring GI numbers are saved in file rerun_gi.txt\n";			
		print rerunFILE $x_gi, "\n";

	}else {
		#print "Running targetgenes\n";
		&targetgenes($x_gi) ;
	}
}

close (rerunFILE); 

print "# FINISHED\n"; 
print "#------------------------------------------------------------------\n";
#========================================================================
# Get sequence from GenBank and retrieve corresponding genome identifyer
#========================================================================

sub targetgenes {
	my $gi = $_[0];

	my (@sorted, %HoH, %CDS, $proteinobj,$outfile);
	my $proteindb = Bio::DB::GenBank->new();
	if ($proteindb->get_Seq_by_gi($gi)) {
		$proteinobj = $proteindb->get_Seq_by_gi($gi);
		$outfile = "targetseqs_gi_".$gi.".fsa";
		open FILE, ">$outfile" or die $!;
		close FILE;
	}else {
		$outfile = "targetseqs_gi_".$gi.".fsa";
		open FILE, ">$outfile" or die $!;
		print FILE "GI not found!\n";
		close FILE;
	}
	
	($proteinobj->annotation)->get_Annotations('dblink') or next;
	my $genomeid = (($proteinobj->annotation)->get_Annotations('dblink'))[0]->primary_id;

	$HoH{$genomeid}{'protein_id'} = $proteinobj->id; 
	$HoH{$genomeid}{'gi'} = $gi ;

	for my $featobj ($proteinobj->get_SeqFeatures) { 
		if ($featobj->primary_tag eq 'CDS') {
			for my $tag ($featobj->get_all_tags){
				$HoH{$genomeid}{'protein_locus_tag'} = ($featobj->get_tag_values('locus_tag'))[0] if ($tag eq 'locus_tag');	
			}
		}
	}

	#========================================================================
	# Get genome from GenBank based on genome identifyer as found above
	#========================================================================
	my $genomedb = Bio::DB::GenBank->new();
	my $genomeobj = $genomedb->get_Seq_by_gi($genomeid);

	for my $featobj ($genomeobj->get_SeqFeatures) { 
		
		if ($featobj->primary_tag eq 'CDS') {
			$featobj->has_tag('locus_tag') or next;

			# Store info on target protein in genome GenBank file 
			if ($HoH{$genomeid}{'protein_locus_tag'} eq ($featobj->get_tag_values('locus_tag'))[0]) {
				$HoH{$genomeid}{'genome_start'} = $featobj->location->start;
				$HoH{$genomeid}{'genome_end'} = $featobj->location->end;
				$HoH{$genomeid}{'genome_locus_tag'} = ($featobj->get_tag_values('locus_tag'))[0];
				$HoH{$genomeid}{'genome_strand'} = 'complement' if ($HoH{$genomeid}{'genome_start'} > $HoH{$genomeid}{'genome_end'});
			}
		
			# Store information on other CDS's in genome GenBank file
			for my $tag ($featobj->get_all_tags){
				if ($tag eq 'translation'){
					$CDS{($featobj->get_tag_values('locus_tag'))[0]}{'translation'} = ($featobj->get_tag_values('translation'))[0];
				}
				$CDS{($featobj->get_tag_values('locus_tag'))[0]}{'strand'} = 'leading';
				if ($featobj->location->start > $featobj->location->end){				
					$CDS{($featobj->get_tag_values('locus_tag'))[0]}{'strand'} = 'complement' ;
					
				}
				$CDS{($featobj->get_tag_values('locus_tag'))[0]}{'start'} = $featobj->location->start;
				$CDS{($featobj->get_tag_values('locus_tag'))[0]}{'end'} = $featobj->location->end;
			}
		}
	}
	#my $gboutfile = "targetseqs_gi_".$gi.".gbk";
	#my $gb_out = Bio::SeqIO->new(-file   => ">$gboutfile", -format => "genbank");
	#$gb_out->write_seq($genomeobj);
	#print "# Wrote genome GenBank to file: $gboutfile\n";

	#========================================================================
	# Sort sequences based on positions and extract 5 up- and down-stream genes
	#========================================================================
	foreach my $cds ( sort { $CDS{$a}->{'start'} <=> $CDS{$b}->{'start'} } keys %CDS){ push (@sorted, $cds) }

	my $target_index = first { $sorted[$_] eq $HoH{$genomeid}{'protein_locus_tag'} } 0 .. $#sorted;
	my @returnedArray = @sorted[$target_index-5..$target_index+5];
	@returnedArray = grep {$_ ne ''} @returnedArray ;
	#print Dumper(@returnedArray);

	#========================================================================
	# Warn if GenBank file does not contain all the surounding genes
	#========================================================================
	print "# WARNING: GenBank file does not reach +/- 5 genes from target. Number of sequences = ", scalar(@returnedArray), "\n" if scalar(@returnedArray) < 11;

	my $seq_out = Bio::SeqIO->new(-file   => ">$outfile", -format => "fasta");

	print "# Wrote target sequences to file: $outfile\n";

	for my $locus_tag (@returnedArray) {
		if (exists $CDS{$locus_tag}){
			my $display_id = $locus_tag."_gi".$gi."_genomeid".$genomeid."_".($CDS{$locus_tag}->{'strand'}."_target_".$HoH{$genomeid}{'protein_locus_tag'});
			my $seqobj = Bio::Seq->new( -display_id => $display_id, -seq => ($CDS{$locus_tag}->{'translation'}));
			#print Dumper($seqobj);	
			$seq_out->write_seq($seqobj);
		}
	}

}

