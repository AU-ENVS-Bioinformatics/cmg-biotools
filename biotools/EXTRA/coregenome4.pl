#!/usr/bin/perl -w

# STUFF TO DO
#   Try to make it portable...
#      Need to do a nice test for required software
#   Need to rethink command-line options
#      Which means re-working the documentation
#   I have partial error checking, but only partial. Could make it better

use strict;
use Cwd;
use Digest::MD5 qw/ md5_hex /;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $ProgramName  = "coregenome.pl";
my $Version      = "4.0c";

# Checking Dependencies
my $FORMATDB	= `which formatdb 2>/dev/null`;
my $BLASTALL	= `which blastall 2>/dev/null`;
my $BLAT	= `which blat 2>/dev/null`;
my $MSUB	= `which msub 2>/dev/null`;
my $GZIP	= `which gzip 2>/dev/null`;
my $GUNZIP	= `which gunzip 2>/dev/null`;
my $PERL	= `which perl 2>/dev/null`;

#my $MCL	= `which mcl 2>/dev/null`;
my $MCL		= "which mcl 2>/dev/null";	# FIX

#my $R		= `which R 2>/dev/null`;
my $R		= `which R 2>/dev/null`;	# FIX

my $MSUB_A	= "nodes=1:ppn=4,mem=8gb,walltime=43200"; # Resource list for aligning
my $MSUB_G	= "nodes=1:ppn=2,mem=8gb,walltime=43200"; # Resource list for grouping
my $MSUB_P	= "nodes=1:ppn=2,mem=8gb,walltime=86400"; # Resource list for permuting

my ($Blast, $Blat, $DataDir, $Failure, $Mcl, $Msub, $Perms, $Permx, $Queue, $Slack, $Slc, $Wanthelp);
my $Size     = 1;
my $Numaxis  = 0;
my $RSpace   = 10;
my $A        = 50; # Default for SLC, alignment length
my $M        = 50; # Default for SLC, match identity
my $LOWER_A  = 20;
my $LOWER_M  = 50;

&GetOptions (
 "blast:s"	=>  \$Blast,	# Use blast alignments, optional string to give path of blastall
 "blat:s"	=>  \$Blat,	# Use blat alignments, optional string to give path of blat
 "mcl:s"	=>  \$Mcl,	# Use mcl clustering, optional string to give path of mcl
 "msub:s"	=>  \$Queue,	# Use of the queuing system. Can be negated
 "nomsub"	=>  \$Msub,	# Do not use the queue
 "numaxis!"	=>  \$Numaxis,  # Should a numeric x-axis be used instead of strain names?
 "o|out=s"      =>  \$DataDir,  # keep temporary files in this dir ( TRY AND CREATE UNLESS EXIST )
 "perms=i"	=>  \$Perms,	# Use permutations, and how many?
 "permx=i"	=>  \$Permx,	# Also permute the order of genomes on the x-axis (=> permuted core and pan curves)
 "r|rspace=i"	=>  \$RSpace,	# Add extra space on the right side of the plot used for long species names
 "size!"        =>  \$Size,     # Enables (or disables) the printing of genome size columns (default)
 "slc:s"	=>  \$Slc,	# Cutoff for SLC & force use on systems where MCL is available
 "slack=i"	=>  \$Slack,	# no. genomes a gene is allowed to be missing from and still be a core gene
 "h|help"       =>  \$Wanthelp  # print the help and exit
);

$Msub = ( defined $Msub ? undef : 1 );

print_help()  if ($Wanthelp);
print_usage() if(!scalar @ARGV && -t STDIN);

$BLASTALL = $Blast if ( defined $Blast && length $Blast );
$BLAT     = $Blat  if ( defined $Blat  && length $Blat  );
$MCL      = $Mcl   if ( defined $Mcl   && length $Mcl   );
for ( $BLASTALL, $BLAT, $FORMATDB, $MCL, $MSUB, $GZIP, $GUNZIP, $PERL, $R ) {
  chomp;
  $_ = undef unless (-x $_);
}
$BLASTALL = undef if ( defined $Blat  );
$BLAT     = undef if ( defined $Blast );
$MCL      = undef if ( defined $Slc   );
$Slc      = ( defined $Slc && length $Slc ? $Slc : "$A/$M" );
$Slack    = 0 unless defined $Slack;
$Perms    = $Permx if ( defined $Permx );
$Numaxis  = 1 if ( defined $Permx );
$Size     = 0 if ( defined $Permx );
$Queue    = ( defined $Queue && length $Queue ? "-q $Queue" : "" );
die "\nERROR: Invalid cutoff for '--slc'. See $ProgramName -h for more information.\n\n" unless ( $Slc =~ /^[0-9]+\/[0-9]+$/ );
die "\nERROR: Using both --slack and --perms is currently unsupported.\n\n" if ( $Slack && $Perms );
die "\nERROR: Permutations only possible if MCL is available.\nPlease make sure it is installed and supply a valid path to MCL if needed.\n\n" if ( $Perms && not defined $MCL );
($A , $M) = split /\//, $Slc;

# Store the args which change program functionality
my @CMD     = (( defined $BLASTALL ? $BLASTALL : $BLAT ), ( defined $MCL ? $MCL : "$A/$M" ));



die "\nERROR: $ProgramName cannot find one of 'blat' or 'blastall'.\nProvide path using one of '--blat' or '--blast'\n\n" unless ( defined $BLASTALL || defined $BLAT );
my $Host = `echo \$HOST`;
die "\nERROR: $ProgramName is not designed to run on (s/i)biology.\nExcecute from the CGE Cluster.\n\n" if ($Host =~ /^(interaction|sbiology|ibiology)$/);

(mkdir $DataDir || die "\n$ProgramName: You must specify a directory with '-o'. Use '-h' for details.\n\n") unless -d $DataDir;
mkdir "$DataDir/log" unless ( -d "$DataDir/log" );
$DataDir = cwd()."/".$DataDir unless ( $DataDir =~ /^\// ); # Is the dir a 





# --------------------------------------------------------------------
# %% MAIN PROGRAM %%
#

##############################
# Create config reading from stdin and prepare data files

my @config = check_jobmd5( );
@config = prepare_fasta(@config);



##############################
# Set up alignment and perform

chdir $DataDir;
@config = ( $BLAT ? align_blat(@config) : align_blastall(@config) );



##############################
# Perform clustering and collect jobs

@config = ( $MCL ? group_mcl(@config) : group_greedy(@config) );



##############################
# Create permuted data and Error Bars?

@config = permute(@config) if ( $Perms );



##############################
# Collect jobs

if ( $MSUB ) {
	my %jobids;
	foreach ( @config ) {
		$jobids{$_->{align_jobid}} = 1 if ( exists $_->{align_jobid} );
		$jobids{$_->{group_jobid}} = 1 if ( exists $_->{group_jobid} );
		$jobids{$_->{perms_jobid}} = 1 if ( exists $_->{perms_jobid} );
	}
	while ( %jobids ) {
		my @jobs = `showq`;
		foreach my $jobid ( keys %jobids ) {
			unless ( grep /^$jobid/, @jobs ) {
				warn "# jobid $jobid finished                                \n" ;
				delete $jobids{$jobid};
			}
		}
		print STDERR "# Waiting for ".(keys %jobids)." job(s) to finish...   \r";
		sleep 5 if ( keys %jobids );
	}
	warn "# All jobs finished!                                    \n";
}

sleep 30 unless (-s "$config[0]->{groupfile}"); # For some stupid reason, (the script might miss some jobids ?!)


##############################
# Perform core/pan calculus

warn "# Calculating Core Genome\n";
@config = core_genome(@config);
warn "# Calculating Pan-Genome\n";
@config = pan_genome(@config);
@config = total_clusters(@config);
@config = new_clusters(@config);
@config = error_bars(@config) if ( $Perms );



##############################
# Plot and terminate

plot(@config);
warn "#\n# ERROR: Some jobs were not finished correctly. Re-start to re-submit" if ( $Failure );
warn "# Job Finished!\n";

exit;




# --------------------------------------------------------------------
# %% SUBROUTINES %%
#

sub check_jobmd5 {
	my $cmd = "$DataDir/cmd";
	my $md5;
	if ( -s $cmd ) {
		warn "# Re-starting in directory: $DataDir. Existing files will not be re-created.\n";
		$md5 = md5("$cmd");
	} else {
		warn "# Using directory: $DataDir\n";
	}
	my @config = &parse_config;
	open CMD , ">$cmd.temp";
	print CMD join("\n", @CMD, map {$_->{source}} @config), "\n";
	close CMD;
	if ( defined $md5 && $md5 ne md5("$cmd.temp") ) {
		die "\nERROR: Job checksum does not match the one recorded for the existing directory!\nIf you change job-defining options you will need a new run.\nEither specify a new directory or remove the existing one.\n\n";
	} else {
		system "mv $cmd.temp $cmd";
	}
	return @config;
}

sub core_genome {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'core genome'} = 0;
		open GRP , "$config[$id]->{groupfile}" or die "\n$!\n$config[$id]->{groupfile}\nERROR: 'Group $id' cluster file missing. Restart to create.\n\n";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my %repr;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$repr{$1} = 1;
			}
			$config[$id]->{'core genome'}++ if (scalar ( keys %repr ) <= $id + 1 and
                                                            scalar ( keys %repr ) >= $id + 1 - $Slack);
		}
		close GRP;
	}
	return @config;
}

sub error_bars {
	warn "# Calculating Error-bars\n";
	my @config = @_;
	my @perm_config;
	my %sums = ('core genome' => [], 'pan genome' => [], 'new clusters' => []);
	open PAN,  ">$DataDir/pan.permutations" or die;
	open CORE, ">$DataDir/core.permutations" or die;
	for my $p ( 0..($Perms-1) ) {
		for my $id ( 0..$#config ) {
			$perm_config[$id]->{groupfile} = $config[$id]->{perm_files}->[$p];
			$perm_config[$id]->{checksum}  = $config[$id]->{perm_checksum}->[$p];
		}
		@perm_config = pan_genome(@perm_config);
		@perm_config = core_genome(@perm_config);
		@perm_config = new_clusters(@perm_config);
		print PAN join("\t",  map { $_->{'pan genome'}  } @perm_config), "\n";
		print CORE join("\t", map { $_->{'core genome'} } @perm_config), "\n";
		@{ $sums{ 'pan genome' } } = map { ( defined $sums{'pan genome'}->[$_] ? $sums{'pan genome'}->[$_] : 0 ) + $perm_config[$_]->{'pan genome'} } 0..$#config;
		@{ $sums{ 'core genome' } } = map { ( defined $sums{'core genome'}->[$_] ? $sums{'core genome'}->[$_] : 0 )  + $perm_config[$_]->{'core genome'} } 0..$#config;
		@{ $sums{ 'new clusters' } } = map { ( defined $sums{'new clusters'}->[$_] ? $sums{'new clusters'}->[$_] : 0 ) + $perm_config[$_]->{'new clusters'} } 0..$#config;
	}
	close PAN;
	close CORE;
	my $ptbl = "$DataDir/tbl.permutations";
	open TBL , ">$ptbl";
	my @keys = ("id", "description", "total genes", "total clusters", "new clusters",
                    "pan genome", "core genome");
	push @keys, "tag" if ( defined $config[0]->{tag} );
	print TBL join ("\t",@keys)."\n";
	foreach my $id ( 0 .. $#config ) {
		print TBL join("\t", map { exists $sums{$_} ? $sums{$_}->[$id] / $Perms : $config[$id]->{$_} } @keys), "\n";
	}
	close TBL;
	return @config;
}

sub fasta2fasta {
	my ($input,$output) = @_;
	my ($temp, $nprot, $skipped) = ("$output.temp.raw", 0, 0);
	my (%ids, $seq, $id);
	open RAW , "| sort > $temp" or die $!; 
	open FASTA , $input || die "\nERROR: Cannot open $input. Did you specify the correct source?\n\n";
	$id = <FASTA>;
	$id =~ s/>(.*)\n/$1/s;
	while (<FASTA>) {
		chomp;
		if (/^>(.*)/) {
			push @{ $ids{$seq} }, $id if ( $id );
			print RAW $seq, "\n" if ( $seq );
			($seq, $id) = ("", $1);
		}  else {
			$seq .= $_;
		}
		if ( eof ) {
			push @{ $ids{$seq} }, $id;
			print RAW $seq, "\n";
		}
	}
	close FASTA;
	close RAW;
	my $checksum = md5 ( $temp );
	open FASTA , "| $PERL -F\"\\t|\\n\" -ane '\$F[2] =~ s/([^\\n]{60})/\$1\\n/gi; chomp \$F[2]; printf \">%s %s\\n%s\\n\", \@F' > $output";
	open TAB , $temp or die $!;
	while (<TAB>) {
		chomp;
		if ( ! /^([A-Za-z]+)$/ ) {
			$skipped++;
			next;
		}
		$nprot++;
		my $id = shift @{ $ids{$1} };
		print FASTA "$checksum.$nprot\t$id /sourcefile=\"$input\"\t$1\n";
	}
	close TAB;
	close FASTA;
	unlink $temp;
	return ($nprot,$skipped,$checksum);
}

sub md5 {
	my @files = @_;
	my $md5 = Digest::MD5->new;
	for my $file (@files) {
		open(FILE, $file) or die "\nERROR: Can't open '$file': $!";
		binmode(FILE);
		while (<FILE>) {
			$md5->add($_);
		}
		close(FILE);
	}
	return $md5->hexdigest;
}

sub new_clusters {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'new clusters'} = 0;
		open GRP , "$config[$id]->{groupfile}" or die "\nERROR: 'Group $id' cluster file missing. Restart to create.\n\n";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my $is_new_family = 1;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$is_new_family = 0 if $1 ne $config[$id]->{checksum};
			}
			$config[$id]->{'new clusters'} += $is_new_family;
		}
		close GRP;
	}
	return @config;
}

sub pan_genome {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'pan genome'} = 0;
		open GRP , "$config[$id]->{groupfile}" or die "\nERROR: Cluster file missing, restart to create ($config[$id]->{groupfile})\n\n";
		while (<GRP>) {
			$config[$id]->{'pan genome'}++;
		}
		close GRP;
	}
	return @config;
}

sub parse_config {
	my @config;
	while ( <> ) {
		next if /^#/;
		chomp;
		my ($description,$source,$tag) = split /\t/;
		warn "# description=$description, source=$source\n";
		my $id = $#config + 1;
		( $config[$id]->{description} , $config[$id]->{source} , 
                  $config[$id]->{id} ,          $config[$id]->{tag} ,
		  $config[$id]->{groupfile} ,	$config[$id]->{mcl_logfile} ) = ( $description , $source , $id , $tag, "$DataDir/group_$id.dat", "$DataDir/log/$id.mcl.log" );
	}
	return @config;
}

sub permute {
	my @config = @_;
	my $seed = time() ^ ($$ + ($$ << 15));
	my %done;
	mkdir "$DataDir/permutations" unless (-d "$DataDir/permutations");
	if (-s "$DataDir/permutations/seed") {
		( $seed ) = map { chomp; $_ } `cat $DataDir/permutations/seed`;
		warn "# Re-starting permutation jobs (Seed retrieved as $seed)\n";
	} else {
		open SEED, ">$DataDir/permutations/seed";
		warn "# Starting permutation jobs (Seed set to $seed)\n";
		print SEED "$seed\n";
		close SEED;
	}
	srand $seed;
	my @perm_configs;
	for my $p ( 0..($Perms-1) ) { # Set number of perms
		my @i = 0..$#config;
		my @perm;
		while ( @i ) {
			push @perm, splice @i, int rand scalar(@i), 1; # Permute
		}
		for my $id ( 0..$#config ) {
			$perm_configs[$p]->[$id]->{align_jobid} = $config[$perm[$id]]->{align_jobid};
			$perm_configs[$p]->[$id]->{checksum}    = $config[$perm[$id]]->{checksum};
			$perm_configs[$p]->[$id]->{groupfile}   = "$DataDir/permutations/$id/group_". md5_hex(join("_", sort map { $perm[$_] } 0..$id)).".dat"; # Set some filename
			$perm_configs[$p]->[$id]->{mcl_logfile} = "/dev/null";
		}
	}
	for my $n ( 0..$#config ) {
		mkdir "$DataDir/permutations/$n" unless (-d "$DataDir/permutations/$n"); # make a dir for each pos...
		open RUN, ">$DataDir/permutations/$n/group.pl";
		for my $p ( 0..($Perms-1) ) {
			push @{ $config[$n]->{perm_files} }, $perm_configs[$p]->[$n]->{groupfile};
			push @{ $config[$n]->{perm_checksum} }, $perm_configs[$p]->[$n]->{checksum};
			next if ( exists $done{ $perm_configs[$p]->[$n]->{groupfile} });
			next if ( -s $perm_configs[$p]->[$n]->{groupfile} ); 		# Check if the clustering has already been done...
			write2script_mcl(-config => $perm_configs[$p], -genome=> $n, -fh => \*RUN);
			$done{$perm_configs[$p]->[$n]->{groupfile}} = 1;
		}
		close RUN;
		next if (-z "$DataDir/permutations/$n/group.pl" ); # If the script file is empty, then everything has been done already :-)
		my $ids;
		foreach my $g ( 0 .. $#config ) {
			$ids .= ":$config[$g]->{align_jobid}" if ( exists $config[$g]->{align_jobid} );
		}
		my ($fh, $file) = tempfile("XXXXXXXX", SUFFIX => ".sh", DIR => "$DataDir/permutations/$n/");
		print $fh "#!/usr/bin/sh\n\n$PERL $DataDir/permutations/$n/group.pl\n";
		close $fh;
		$config[$n]->{perms_jobid} = $MSUB ? `$MSUB -j oe -o $DataDir/log/perm.$n.log -r y $Queue -d $DataDir/permutations/$n -l $MSUB_P,depend=afterok$ids $file` : `$PERL $DataDir/permutations/$n/group.pl`;
		$config[$n]->{perms_jobid} =~ s/\n//g;
		if ( $config[$n]->{perms_jobid} ) {
			warn "# jobid $config[$n]->{perms_jobid} submitted using MCL (Permutating at position $n of $#config)\n";
		} else {
			warn "# job submission failed??! (Permutating at position $n)\n# --PLEASE RE-RUN SCRIPT!!\n";
			delete $config[$n]->{perms_jobid};
		}
	}
	return @config;
}

sub prepare_fasta {
	my @config = @_;
	for my $id ( 0 .. $#config ) {
		my ($nprot,$skipped,$checksum) = fasta2fasta( $config[$id]->{source} , "$DataDir/$id.fsa");
		mkdir "$DataDir/$checksum" unless (-d "$DataDir/$checksum");
		warn "# prepare_fasta($id): md5=$checksum, nprot=$nprot, skipped=$skipped\n";
		$config[$id]->{target} = "$DataDir/$checksum/$id.fsa";
		$config[$id]->{'total genes'} = $nprot;
		$config[$id]->{skipped} = $skipped;
		$config[$id]->{checksum} = $checksum;
		die "\nERROR: nprot == 0 for $config[$id]->{source}\n" if $nprot == 0;
	}
	return @config;
}

sub total_clusters {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		$config[$id]->{'total clusters'} = 0;
		open GRP , "$DataDir/group_$id.dat";
		while (<GRP>) {
			my ($group_id,$count,@members) = split /\t/;
			my $is_member = 0;
			foreach (@members) {
				next unless /^(.*)\.(.*)/;
				$is_member = 1 if $1 eq $config[$id]->{checksum};
			}
			$config[$id]->{'total clusters'} += $is_member;
		}
		close GRP;
	}
	return @config;
}





##############################
# Do BLASTALL Alignments

sub align_blastall {
	my @config = @_;
	my %jobids;
	foreach my $id ( 0 .. $#config ) {
		open RUN , ">$DataDir/$config[$id]->{checksum}/align.pl";
		print RUN "#!$PERL\nmy \%checksum = ( 0 => '$config[0]->{checksum}'";
		foreach my $a ( 1 .. $#config ) {
			print RUN ", $a => '$config[$a]->{checksum}'";
		}
		print RUN ");\n";
		print RUN <<EOS;

			chdir "$DataDir/$config[$id]->{checksum}"; # Ensure we are in the right dir
			open \$fsa_fh , ">$id.fsa";
			fasta("$DataDir/$id.fsa", \$fsa_fh);
			close \$fsa_fh;
			system "$FORMATDB -i $id.fsa -p T -t $config[$id]->{checksum}";
			my \$pipe = "$DataDir/$config[$id]->{checksum}/align.pipe";
			system "rm \$pipe" if (-e \$pipe);
			`mknod \$pipe p`;
			if ( !fork() ) {
				open \$pipe_fh , ">\$pipe" || die;
				foreach my \$a ( 0 .. $#config ) {
					next if ( -e "$DataDir/$config[$id]->{checksum}/\$checksum{\$a}-$config[$id]->{checksum}.align.gz" );
					fasta("$DataDir/\$a.fsa", \$pipe_fh);
				}
				close \$pipe_fh;
				exit;
			}
			my \%fh;
			foreach my \$checksum ( values %checksum ) {
				(open \$fh{\$checksum}, "| ".(-e "$GZIP" ? "$GZIP -c" : "cat")." > \$checksum-$config[$id]->{checksum}.align.temp.gz" || die) unless ( -e "\$checksum-$config[$id]->{checksum}.align.gz" );
			}
			open BLAST, "$BLASTALL -a 3 -e 1 -F 0 -m 8 -p blastp -i \$pipe -d $id.fsa |";
			while (<BLAST>) {
				my (\$qid, \$sid, \$pct, \$aln) = split /\\t|\\n/, \$_;
				my (\$qmd5, \$qn, \$qln) = split /\\./, \$qid;
				my (\$smd5, \$sn, \$sln) = split /\\./, \$sid;
				my \$alr = 100 * \$aln / (\$qln > \$sln ? \$qln : \$sln);
				print { \$fh{\$qmd5} } "\$qmd5.\$qn\\t\$smd5.\$sn\\t\$alr/\$pct\\n" if ( \$pct >= \$LOWER_M && \$alr >= \$LOWER_A );
			}
			while ( my (\$k, \$v) = each %fh ) {
				close \$v;
				system "mv \$k-$config[$id]->{checksum}.align.temp.gz \$k-$config[$id]->{checksum}.align.gz";
			}

			sub fasta {
				my (\$in_name, \$out_fh) = \@_;
				open FASTA , \$in_name;
				my (\$seq, \$qln) = ("", 0);
				while (<FASTA>) {
					if ( /^>/ ) {
						\$seq =~ s/^(>\\S+)/sprintf "%s.%d", \$1, \$qln/ie;
						print \$out_fh \$seq;
						(\$seq, \$qln) = ("", 0);
					} else {
						\$qln += length() - 1;
					}
					\$seq .= \$_;
					if ( eof ) {
						\$seq =~ s/^(>\\S+)/sprintf "%s.%d", \$1, \$qln/ie;
						print \$out_fh \$seq;
					}
				}
				close FASTA;
			}
			print (($id+1)."\n");

EOS
		close RUN;
		$config[$id]->{align_jobid} = $MSUB ? `$MSUB -j oe -o $DataDir/$config[$id]->{checksum}/blast.log -r y $Queue -d $DataDir/$config[$id]->{checksum} -l $MSUB_A $DataDir/$config[$id]->{checksum}/align.pl` : `$PERL $DataDir/$config[$id]->{checksum}/align.pl`;
		$config[$id]->{align_jobid} =~ s/\n//g;
		if ( $config[$id]->{align_jobid} ) {
			warn "# jobid $config[$id]->{align_jobid} submitted using blast (Genome=$id, md5=$config[$id]->{checksum})\n";
		} else {
			warn "# job submission failed??! (Genome=$id, md5=$config[$id]->{checksum})\n# --PLEASE RE-RUN SCRIPT!!\n";
			$Failure = 1;
			delete $config[$id]->{align_jobid};
		}
	}
	return @config;
}

sub align_blat {
	my @config = @_;
	my %jobids;
	foreach my $id ( 0..$#config ) {
		open RUN , ">$DataDir/$config[$id]->{checksum}/align.pl";
		print RUN "#!$PERL\nmy \%checksum = ( 0 => '$config[0]->{checksum}'";
		foreach my $a ( 1 .. $#config ) {
			print RUN ", $a => '$config[$a]->{checksum}'";
		}
		print RUN ");\n";
		print RUN <<EOS;

			print (($id+1)."\n");
			chdir "$DataDir/$config[$id]->{checksum}"; # Ensure we are in the right dir
			my \$pipe = "$DataDir/$config[$id]->{checksum}/align.pipe";
			system "rm \$pipe" if (-e \$pipe);
			`mknod \$pipe p`;
			my \$list = "blat.list";
			open LIST, ">\$list";
			foreach my \$a ( 0 .. $#config ) {
				print LIST "$DataDir/\$a.fsa\\n" unless ( -e "$DataDir/$config[$id]->{checksum}/\$checksum{\$a}-$config[$id]->{checksum}.align.gz" );
			}
			close LIST;
			exit unless ( -s \$list );
			if ( !fork() ) {
				my \%fh;
				foreach my \$checksum ( values %checksum ) {
					(open \$fh{\$checksum}, "| ".(-e "$GZIP" ? "$GZIP -c" : "cat")." > \$checksum-$config[$id]->{checksum}.align.temp.gz" || die) unless ( -e "\$checksum-$config[$id]->{checksum}.align.gz" );
				}
				open PSL, "tail -n 6 \$pipe | cut -f 1,2,10,11,14,15 |";
				while (<PSL>) {
					my (\$mat, \$mis, \$qid, \$qln, \$sid, \$sln) = split /\\t|\\n/, \$_;
					my (\$md5, \$n) = split /\\./, \$qid;
					my \$aln = \$mat + \$mis;
					my \$pct = \$mat / \$aln * 100;
					my \$alr = 100 * \$aln / (\$qln > \$sln ? \$qln : \$sln);
					print { \$fh{\$md5} } "\$qid\\t\$sid\\t\$alr/\$pct\\n" if ( \$pct >= \$LOWER_M && \$alr >= \$LOWER_A );
				}
				exit;
			}
			system "$BLAT -prot $DataDir/$id.fsa \$list \$pipe > blat.log";
			while (wait != -1) {};
			foreach my \$checksum ( values %checksum ) {
				system "mv \$checksum-$config[$id]->{checksum}.align.temp.gz \$checksum-$config[$id]->{checksum}.align.gz" if ( -e "\$checksum-$config[$id]->{checksum}.align.temp.gz" );
			}

EOS
		close RUN;
		$config[$id]->{align_jobid} = $MSUB ? `$MSUB -j oe -o $DataDir/$config[$id]->{checksum}/blat.log -r y $Queue -e -o  -d $DataDir/$config[$id]->{checksum} -l $MSUB_G $DataDir/$config[$id]->{checksum}/align.pl` : `$PERL $DataDir/$config[$id]->{checksum}/align.pl`;
		$config[$id]->{align_jobid} =~ s/\n//g;
		if ( $config[$id]->{align_jobid} ) {
			warn "# jobid $config[$id]->{align_jobid} submitted using blat (Genome=$id, md5=$config[$id]->{checksum})\n";
		} else {
			warn "# job submission failed??! (Genome=$id, md5=$config[$id]->{checksum})\n# --PLEASE RE-RUN SCRIPT!!\n";
			$Failure = 1;
			delete $config[$id]->{align_jobid};
		}
	}
	return @config;
}





##############################
# Do Clustering

sub group_greedy {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		next if ( -s "$config[$id]->{groupfile}" ); # Check if the clustering has already been done...
		open RUN, ">$DataDir/$config[$id]->{checksum}/group.pl";
		print RUN "#!$PERL\nmy \%checksum = ( 0 => '$config[0]->{checksum}'";
		foreach my $id ( 1 .. $#config ) {
			print RUN ", $id => '$config[$id]->{checksum}'";
		}
		print RUN ");\n";
		print RUN <<EOS;

			chdir "$DataDir/$config[$id]->{checksum}"; # Ensure we are in the right dir
			my \$pipe = "group.pipe";
			system "rm \$pipe" if (-e \$pipe);
			`mknod \$pipe p`;
			if ( !fork() ) {
				open PIPE , ">\$pipe";
				foreach my \$x (0 .. $id) {
					foreach my \$y (0 .. $id) {
						open ALIGN , (-e "$GUNZIP" ? "$GUNZIP -c" : "cat")." $DataDir/\$checksum{\$y}/\$checksum{\$x}-\$checksum{\$y}.align.gz |" || die;
						while (<ALIGN>) {
							my ( undef, undef, \$c ) = split /\\t|\\n/, \$_;
							my ( \$a, \$m ) = split /\\//, \$c;
							print PIPE \$_ if ( \$a >= $A && \$m >= $M );
						}
						close ALIGN;
					}
				}
				close PIPE;
				exit;
			}
			open GRP , \$pipe;
			open DAT , ">$config[$id]->{groupfile}.temp";
			my (%all, %COMBINE, %groupArray);
			while ( <GRP> ) {
				chomp;
				my (\$subject,\$query) = split("\t");
				\$query = \$subject unless \$query =~ /[A-Z0-9]+/gi ;
				( \$all{\$query}, \$all{\$subject}, \$COMBINE{\$query}{\$subject}, \$COMBINE{\$subject}{\$query} ) = ( 1,1,1,1 );
			}
			close GRP;
			&groupMainTree;
			foreach my \$group_no ( sort { \$a <=> \$b } ( keys %groupsArray )) {
				print DAT "\$group_no\\t".(\$#{\$groupsArray{\$group_no}}+1)."\\t".join("\\t", \@{\$groupsArray{\$group_no}})."\\n";
			}
			close DAT;
			system "mv $config[$id]->{groupfile}.temp $config[$id]->{groupfile}";

			sub groupMainTree {
				\$groupCount=1;
				foreach my \$entry (keys %all) {
					next if \$included{\$entry};
					\$included{\$entry} = 1;
					push \@{\$groupsArray{\$groupCount}}, \$entry;
					group(\$entry);
					\$groupCount++;
				}
			}

			sub group {
				my \$query = \$_[0];
				foreach my \$subject (keys \%{\$COMBINE{\$query}}) {
					next if \$included{\$subject};
					\$included{\$subject} = 1;
					push \@{\$groupsArray{\$groupCount}} , \$subject;
					group(\$subject);
				}
			}
			print (($id+1)."\n");

EOS
		close RUN;
		my $ids;
		foreach my $g ( 0 .. $id ) {
			$ids .= ":$config[$g]->{align_jobid}" if ( exists $config[$g]->{align_jobid} );
		}
		$config[$id]->{group_jobid} = $MSUB ? `$MSUB -j oe -o $DataDir/log/slc.$id.log -r y $Queue -d $DataDir/$config[$id]->{checksum} -l $MSUB_G,depend=afterok$ids $DataDir/$config[$id]->{checksum}/group.pl` : `$PERL $DataDir/$config[$id]->{checksum}/group.pl`;
		$config[$id]->{group_jobid} =~ s/\n//g;
		if ( $config[$id]->{group_jobid} ) {
			warn "# jobid $config[$id]->{group_jobid} submitted using single-linkage clustering (Building group $id of $#config)\n";
		} else {
			warn "# job submission failed??! (Genome=$id, md5=$config[$id]->{checksum})\n# --PLEASE RE-RUN SCRIPT!!\n";
			$Failure = 1;
			delete $config[$id]->{group_jobid};
		}
	}
	return @config;
}

sub group_mcl {
	my @config = @_;
	foreach my $id (0 .. $#config) {
		next if ( -s "$config[$id]->{groupfile}" ); # Check if the clustering has already been done...
		open RUN, ">$DataDir/$config[$id]->{checksum}/group.pl";
		write2script_mcl(-config => \@config, -genome=> $id, -fh => \*RUN);
				
		close RUN;
		my $ids;
		
		foreach my $g ( 0 .. $id ) {
			$ids .= ":$config[$g]->{align_jobid}" if ( exists $config[$g]->{align_jobid} );
			
		}

		$config[$id]->{group_jobid} = `$PERL $DataDir/$config[$id]->{checksum}/group.pl`;		
		#$config[$id]->{group_jobid} = $MSUB ? `$MSUB -j oe -o $config[$id]->{mcl_logfile} -r y $Queue -d $DataDir/$config[$id]->{checksum} -l $MSUB_G,depend=afterok$ids $DataDir/$config[$id]->{checksum}/group.pl` : `$PERL $DataDir/$config[$id]->{checksum}/group.pl`;
		#$config[$id]->{group_jobid} =~ s/\n//g;
		if ( $config[$id]->{group_jobid} ) {
			warn "# jobid $config[$id]->{group_jobid} submitted using MCL (Building group $id of $#config)\n";
		} else {
			warn "# job submission failed??! (Genome=$id, md5=$config[$id]->{checksum})\n# --PLEASE RE-RUN SCRIPT!!\n";
			$Failure = 1;
			delete $config[$id]->{group_jobid};
		}
	}
	return @config;
}

sub write2script_mcl {
	my %argv   = @_;
	my @config = @{ $argv{-config} };
	my $id     = $argv{-genome};
	my $fh     = $argv{-fh};
	print $fh "#!$PERL\nmy \%checksum = ( 0 => '$config[0]->{checksum}'";
	foreach my $id ( 1 .. $#config ) {
		print $fh ", $id => '$config[$id]->{checksum}'";
	}
	print $fh ");\n";
	print $fh <<EOS;

		chdir "$DataDir/$config[$id]->{checksum}"; # Ensure we are in the right dir
		my \$mcl_pipe = "$config[$id]->{groupfile}.mcl.pipe";
		system "rm \$mcl_pipe" if (-e \$mcl_pipe);
		`mknod \$mcl_pipe p`;
		my \$pipe = "$config[$id]->{groupfile}.group.pipe";
		system "rm \$pipe" if (-e \$pipe);
		`mknod \$pipe p`;
		if ( !fork() ) {
			system "$MCL \$pipe --abc -o \$mcl_pipe";
			exit;
		}
		if ( !fork() ) {
			open PIPE , ">> \$pipe";
			foreach my \$x (0 .. $id) {
				foreach my \$y (0 .. $id) {
					open ALIGN , (-e "$GUNZIP" ? "$GUNZIP -c" : "cat ")." $DataDir/\$checksum{\$y}/\$checksum{\$x}-\$checksum{\$y}.align.gz |" || die;
					while (<ALIGN>) {
						my ( \$qid, \$sid, \$c ) = split /\\t|\\n/, \$_;
						my ( \$a, \$m ) = split /\\//, \$c;
						print PIPE join("\t", \$qid, \$sid, \$a * \$m / 10000), "\n";
					}
					close ALIGN;
				}
			}
			close PIPE;
			exit;
		}
		open MCL , "cat \$mcl_pipe |";
		open OUT , ">$config[$id]->{groupfile}";
		while (<MCL>) {
			my \@a = split /\\t/, \$_;
			print OUT join("\\t", \$., scalar \@a, \$_);
		}
		close MCL;
		close OUT;
		unlink \$mcl_pipe;
		unlink \$pipe;
		print (($id+1)."\n");

EOS
}





##############################
# Draw Plot

sub plot {
	my @config = @_;
	my $tbl = "$DataDir/tbl";
	open TBL , ">$tbl";
	my @keys = ("id", "description", "total genes", "total clusters", "new clusters",
                    "pan genome", "core genome");
	push @keys, "tag" if ( defined $config[0]->{tag} );
	print TBL join ("\t",@keys)."\n";
	foreach my $id ( 0 .. $#config ) {
		print TBL join("\t", map { $config[$id]->{$_} } @keys), "\n";
	}
	close TBL;
	$tbl = ( $Permx ? "$DataDir/tbl.permutations" : $tbl ); # We create both tbl's, but only read from the one the user asked for
	open R , "| $R --vanilla > /dev/null" or return;
	warn "# Creating Plot\n";
	print R "
postscript('$DataDir/ps', title='$ProgramName v$Version');
data <- read.table('$tbl', sep='\t',dec='.', header=T, row.names=NULL);
if (0$Slack) data[1:($Slack), 'core.genome'] <- NA
layout( matrix(c(1,2), nrow= 2), heights=c(10,3))
op <- par(mar = c(0,2,2,1))
x<-rbind(data[,3], data[,4], data[,5])
data[which.max(data[,'core.genome']),'core.genome'] <- data[1,5]
ymax <- max(1.2*max(data[,6]), 2.0 * max(data[,3]))
col <- c(gray(0.8), gray(0.5), gray(0.2))
";
	if ( $Size ) {
		print R "
r <- barplot(x, beside=TRUE, plot=F)
rspace <- $RSpace * max(r) / 100
xlim <- c(min(r)-0.5, max(r)+0.5+rspace)
r <- barplot(x, beside=TRUE, ylim=c(0,ymax), xlim=xlim, col=col)
lines(r[3,1:(ncol(r)-sum(is.na(data[,'core.genome'])))], type='b', data[!is.na(data[,'core.genome']),'core.genome'], col='red', lwd=4)
lines(r[3,], type='b', data[,6], col='blue', lwd=4)
legend(1,ymax,c('Total Genes', 'Total Gene Clusters', 'New Gene Clusters','Core Genome','Pan Genome'),col=c(col,'red','blue'), lwd=c(4,4,4,4))
r.errorbars <- r[3,]
"	} else {
		print R "
r <- barplot(x[3,], width=1/1.2, plot=F)
rspace <- $RSpace * max(r) / 100
xlim <- c(min(r)-0.5, max(r)+0.5+rspace)
r <- barplot(x[3,], width=1/1.2, ylim=c(0,ymax), xlim=xlim, col=col[3])
lines(r[1:(length(r)-sum(is.na(data[,'core.genome'])))], type='b', data[!is.na(data[,'core.genome']),'core.genome'], col='red', lwd=4)
lines(r, type='b', data[,6], col='blue', lwd=4)
legend(1,ymax,c('New Gene Clusters','Core Genome','Pan Genome'),col=c(col[3],'red','blue'), lwd=c(4,4,4))
r <- t(r)
r.errorbars <- r
"	}
	if ( $Perms ) {
		print R "
pan.perms <- read.table('$DataDir/pan.permutations')
pan.perms.sd <- apply(pan.perms, 2, sd)
pan.perms.mean <- apply(pan.perms, 2, mean)
nz <- pan.perms.sd > 0
arrows(r.errorbars[nz], data[nz,6] - pan.perms.sd[nz], r.errorbars[nz], data[nz,6] + pan.perms.sd[nz], angle=90, code=3, length=0.1)
core.perms <- read.table('$DataDir/core.permutations')
core.perms.sd <- apply(core.perms, 2, sd)
core.perms.mean <- apply(core.perms, 2, mean)
nz <- core.perms.sd > 0
arrows(r.errorbars[nz], data[nz,'core.genome'] - core.perms.sd[nz], r.errorbars[nz], data[nz,'core.genome'] + core.perms.sd[nz], angle=90, code=3, length=0.1)
"	}
	print R "
if (0$Numaxis) {
	axis(1)
} else {
	op <- par(mar = c(1,2,0,1))
	plot.new()
	plot.window(ylim=c(0,10),xlim=xlim)
	for (i in 1:nrow(data)) {
		text(mean(r[,i]), 9, adj=0, data[i,2], cex=0.7, srt=-45)
	}
}
dev.off();
";
	close R;
	system "cat $DataDir/ps";
}





# --------------------------------------------------------------------
# %% Documentation %%
#

sub print_usage {
  print STDERR "Usage: perl $ProgramName [Optional Options...] -o [Dir] < [Config_file]\nSee $ProgramName -h for details.\n";
  exit;
}

sub print_help {
  my $blat  = ( defined $BLAT     ? $BLAT         : "'blat' not found!"     );
  my $blast = ( defined $BLASTALL ? $BLASTALL     : "'blastall' not found!" );
  my $mcl   = ( defined $MCL      ? $MCL          : "'mcl' not found!"      );
  my $msub  = ( defined $MSUB     ? "msub Found!" : "No queue found!"       );
  my $size  = ( defined $Size     ? "size"        : "nosize"                );
  my $queue = ( defined $Queue    ? "'$Queue'"    : "Do not set 'msub -q'"  );
  for ( $blat, $blast, $mcl, $msub, $size ) {
    chomp;
  }
  open LESS, "| less";
  print LESS <<EOH;

NAME
    $ProgramName - derive core-/pan genome for a list of genomes/samples.

SYNOPSIS
    perl $ProgramName [-o <dir>] [other options] < list > output.ps

DESCRIPTION
    Reads a list from STDIN or a file, defining the genomes/proteomes to be
    compared. The order of the list indicate in what order the genomes will
    appear in the plot. Each line must contain a description and data source,
    separated by tab. Example:

Campylobacter jejuni M1	/home/people/carsten/scripts/coregenome/testdata/38041.proteins.fsa
Campylobacter jejuni subsp. jejuni CF93-6	/home/people/carsten/scripts/coregenome/testdata/16265.proteins.fsa
Campylobacter jejuni subsp. jejuni HB93-13	/home/people/carsten/scripts/coregenome/testdata/16267.proteins.fsa
Campylobacter jejuni subsp. jejuni 84-25	/home/people/carsten/scripts/coregenome/testdata/16367.proteins.fsa
Campylobacter jejuni subsp. jejuni 81-176	/home/people/carsten/scripts/coregenome/testdata/16135.proteins.fsa
Campylobacter jejuni subsp. jejuni 260.94	/home/people/carsten/scripts/coregenome/testdata/16229.proteins.fsa
Campylobacter coli RM2228	/home/people/carsten/scripts/coregenome/testdata/12516.proteins.fsa
Campylobacter lari RM2100	/home/people/carsten/scripts/coregenome/testdata/12517.proteins.fsa
Campylobacter upsaliensis RM3195	/home/people/carsten/scripts/coregenome/testdata/12518.proteins.fsa
Campylobacter fetus subsp. fetus 82-40	/home/people/carsten/scripts/coregenome/testdata/16293.proteins.fsa

    Each proteome is aligned against all proteomes in the list and the naming
    is derived based in MD5 checksums of the raw proteins sequences. Based on
    the alignments, all proteins will be clustered into gene clusters.

    The program will for each proteome X calculate the follwing:

    - Number of proteins observed in X
    - Number of new gene clusters observed in genome in X compared to 0 .. (X-1)
    - Size of pan genome at genome X
    - Size of core genome at genome X

    Both the choice of alignment algorithm and clustering method are adjustable
    via options.

OPTIONS
    --blast [path of 'blastall' (Optional)]
	Set the alignments to be performed by BLAST.
	This is the most accurate approach, but it is much slower than BLAT. The
	script will search for the 'blastall' executable using 'which blast', but
	it is also possible to specify the path using the optional argument to
	this option. If found, BLAT is the default, BLAST the fall-back.
	BLAST used: $blast

    --blat [path of 'blat' (Optional)]
	Set the alignments to be performed by BLAT.
	BLAT is much faster and almost as accurate as BLAST. The script will
	search for the 'blat' executable using 'which blat', but it is also
	possible to specify the path using the optional argument to this option.
	If found, BLAT is the default, BLAST the fall-back.
	BLAT used: $blat

    --mcl [path of 'mcl' (Optional)]
	Set the gene clustering to be performed by MCL.
	This script features a build-in implementation of single-linkage
	clustering (SLC), but MCL is both faster and generally does a better
	job. Unlike SLC which uses an absolute cutoff, for MCL a score is
	calculated as the multiplum of the length of the alignment relative to
	the longest sequence and the reported percent identity; provided that
	these values are at least $LOWER_A % and $LOWER_M %, respectively.

	The script will search for the 'mcl' executable using 'which mcl', but
	it is also possible to specify the path using the optional argument to
	this option. If found, MCL will be the default, otherwise SLC will be
	used.
	MCL used: $mcl
	MCL is written by Stijn van Dongen; see:
	http://micans.org/mcl/man/mcl.html#references

    --msub <String> or --nomsub
	Set the script to use (or ignore) a queuing system.
	Using the queue jobs will run in parallel greatly speeding things up.
	Otherwise they will be executed one at a time. The 'msub' command
	itself needs to be findable with 'which msub'. The optional string
	argument allows the specification of a specific queue to use, if several
	exist. Default is to use the queue if available.
	Queuing system check: $msub
	Default queue: $queue
	
    --numaxis
	Draw a numeric x-axis instead of using strain names.
	NOTE: This option should be used with '--nosize' since the numeric axis
	will be extremely misleading in a plot with juxtaposed columns.

    -o [dir]
	Set the required output directory (the working directory).
	If it already exists, it will be checked for any files contained therein
	stemming from an earlier run, and any files found will not be remade. If
	it doesn't exist, it will be created.

    --perms <N>
	Add error bars from <N> permutations of the data.
	Performs <N> permutations of the data and uses the result to calculate
	error bars for the core and pan genome curves. NOTE: Permutations are
	extremely cpu-taxing as the clustering needs to be re-done for every
	single permutation of each genome. Use of MCL will be forces and even
	then only small values of <N> should be used. If MCL is unavailable,
	permutations are not possible.

    --permx <N>
	Add error bars and permute the order of genomes on the x-axis.
	Similar to '--perms' but the order of the genomes along the x-axis is
	permuted as well. Automatically forces the '--nosize' and '--numaxis'
	options. The size columns cannot meaningfully be drawn in this type of
	plot, and printing the strain names would be downright wrong.

    -r or --rspace [Integer]
	Add empty space to the right of the plot. Increase if the column names
	don't fit in the plotting window. Given as a relative value equal to a
	percentage of the uncompressed (original) plot; e.g. a value of 100
	would compress the chart to fill half of the page.
	Default: $RSpace

    --size, --nosize
	Enforce or suppress the printing of the Genome Size columns.
	Default: $size

    --slack <N>
	Number of genomes a gene is allowed to be missing from and still be
	considered part of the core genome. Notice that the core genome curve
	will not reach the rightmost genomes when using this option. This is
	natural, because the core genome of species 'n' now requires data from
	genome 'n+slack', which is unavailable for the rightmost species.
	Default: 0 (No slack)

    --slc [cutoff]
	Force use of single-linkage clustering (SLC) and adjust cutoff.
	Because MCL is both faster and generally does a better job than SLC, it
	is the default, if MCL is available. SLC would only be desirable if MCL
	is not installed, but use of SLC can be forced by specifying this option
	even if MCL is available.

	To form a linkage between two sequences, SLC requires a cutoff of the
	format "length/identity"; e.g. '50/50' would mean the program assumes
	two sequences being similar if the aligment identified spans 50 % or
	more of the longest of the sequence AND that 50 % of the residues are
	conserved within the alignment. Only the reciprocal match is considered.
	Default: $A/$M

VERSION
    Current Version: $ProgramName $Version

    Version Changes:

    Version 4.0c represents a major overhaul. It adapts the code to the CGE
    cluster while at the same time allowing it to run on other platforms (albeit
    without parallelization). It also incorporates MCL and allows for BLAT to
    substitute BLAST, both of which improves performance substantially.
    Exploiting the faster MCL algorithm, options for permutation of the data
    have been added.

    Version 3.0c has been updated to use the sbiology cueing system. This also
    changes the behaviour of the '-cpu' option which now de-facto controls the
    amount of simultaneous submissions to the cue. Thus this number should now
    generally be higher than for previous versions else "submission-lag" will be
    a factor. New '-cpu' default has been set.
    New options: '-a' and '-m' options to finetune the cutoff

    Version 2.8c includes smaller bugfixes and adjusts composition of temporary
    files to better integrate with the coregenes script.

    Version 2.7c added the 'rspace' option as well as the 'total gene families'
    column. Also replaced 'saco_convert' with 'tab2fasta' for preparing files
    for blast. 'saco_convert' has problems with long identifiers which can be
    disasterous.

    Version 2.6c mainly adjusts the composition of temporary files to better
    integrate with the coregenes script for extracting individual genes.

    Version 2.5c introduced the '-nosize' option.

    Version 2.4c introduced new options and fixed the bug which prevented the
    script from using newer versions of BLAST. As a consequence it should run
    faster.

AUTHOR
    Carsten Friis, carsten\@cbs.dtu.dk,
    Based on script made by Peter Fischer Hallin, pfh\@cbs.dtu.dk

EOH
  close LESS;
  exit;
}
