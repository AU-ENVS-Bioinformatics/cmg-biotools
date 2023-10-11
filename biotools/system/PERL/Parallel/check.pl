#!/usr/bin/perl
use ChildManager;

$|= 1;
my $nr_of_childprocesses = 15;
my $cm = new ChildManager($nr_of_childprocesses);
$cm->set_start_type('READ');

my $coco = 0;
my @fh;
my $index = 0;
my $active_children = 0;
foreach $entry (0..300){
        print "Mother: starting up child $entry\n";
	$fh[$index]=$cm->start ;# initiating child
	unless(defined $fh[$index]){ #Child process
	        for (my $i = int(rand(1000)); $i > 0; $i--) {
			print "Child: $entry\n"; }
		sleep(int(rand(10)));
		print "# Child output ends\n";
		exit;
	}
        # Parent
	$active_children++;
        $index = ($index+1) % $nr_of_childprocesses;
	next if $active_children < $nr_of_childprocesses;
print "Mother: reading childs output (index:$index, entry:$entry)\n";
	my $lcount = 0;
	my $line = '';
	my $thisfh = $fh[$index];
print "Mother: FILEHANDLE CLOSED (index:$index, entry:$entry)\n" if eof($thisfh);
print "Mother: FILEHANDLE UNDEF (index:$index, entry:$entry)\n" unless defined $thisfh;
        my $cline;
	while (defined ($line = <$thisfh>) and $line ne "# Child output ends\n" ) {
	   $cline = $line unless $lcount;
	   $lcount++; }
	close $thisfh;
	print $cline;
	$cline =~ m/(\d+)/o;
print "Mother: MISSING CHILD $coco\n" if $1 != $coco;
	$coco++;	
print "Mother: got $lcount lines\n";
print "Mother: Child did not give all output\n" unless defined $line;
print "Mother: children:",  join(' ', keys %{$cm->{'proc_tab'}}), "\n";
 #	$cm->wait_next_child('BLOCK');
#print "Mother: child died\n";
    	$active_children--;
}

while ($active_children > 0) {
        $index = ($index+1) % $nr_of_childprocesses;
	next unless defined $fh[$index];
print "Mother: reading childs output\n";
	my $lcount = 0;
	my $line = '';
	my $thisfh = $fh[$index];
	while (defined ($line = <$thisfh>) and $line ne "# Child output ends\n" ) {
	   print $line unless $lcount;
	   $lcount++; }
	close $fh[$index];
print "Mother: got $lcount lines\n";
print "Mother: Child did not give all output\n" unless defined $line;
	$active_children--;
}
$cm->wait_all_children;

print "End\n";
