#!/usr/bin/perl
# Read a fasta file and extract the DNA sequence data
# Translate it to protein and print it out in 25-character-long lines

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

my $file;

&GetOptions ("f:s" =>  \$file);

unless (defined $file) {
	print "# FASTA DNA file not defined \n";
	exit;
}

# this implicity uses the <> file stream
my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => $file);
my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => '>'.$file.".fsa");

my $prot_obj;

while (my $seq = $seqin->next_seq()) {
	my $len = $seq->length();
	#print ">", $seq->desc(), " ", $seq->id(), "_gene size kb $len\n";
	$prot_obj = $seq->translate(-frame => 0);

	$seqout->write_seq($prot_obj);
}
