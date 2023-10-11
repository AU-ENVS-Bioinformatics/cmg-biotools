#!/usr/bin/perl

# Authors: Peter Fisher Hallin
# For license see /usr/biotools/CMG-biotools.license

$program = "blastn";
$Identities = 0;

while (<>){
	if (/Query= +(.*)\n/) {
		#print "TEST: ",$_,"\n";
		
		$queryname = $1;
		chomp $queryname;
		$queryname =~ s/\s+/s/g;
		#my @split = split("_", $queryname);
		#$name = $split[0];
		#$qgs = $split[-2];
		#$qge = $split[-1];
		#print "TEST:" ,$name, "\t", $qgs, "\t", $qge, "\n";
		if ($queryname =~ /(.*\S)_+(\S.*)-+(\S.*)/) {
			$name = $1;
			$qgs = $2;
			$qge = $3;
		}
	}
	if (/\(+(.*) letters\)\n/) {
		$querylength = $1;
	}
	if (/^>+(.*?)\n/) {
		if ($querystart) {
			print $name, "\t", $program, "\t", $qgs, "\t", $qge, "\t", $sgs, "\t", $sge, "\t", $bits, "\t", $Escore, "\t", $Identities, "\t", $Positives, "\t", $querystart, "\t", $queryend, "\t", ($queryend > $querystart ? "+" : "-"), "\t", $sbjctstart, "\t", $sbjctend, "\t", ($sbjctend > $sbjctstart ? "+" : "-"), "\t", $queryname, "\t", $sbjctname, "\t", $strand1, "\t", $strand2, "\t", $sbjctlength, "\t", $querylength, "\t", $frame1, "\t", $frame2, "\n";
			$querystart = 0;
			$sbjctstart = 0;
			$sbjctname = $1;
			#print "TEST: ", $sbjctname, "\t", $program, "\t", $qgs, "\t", $qge, "\n";
			chomp $sbjctname;
		} else {	
			$sbjctname = $1;
			chomp $sbjctname;
		}
		$sbjctname =~ s/\s+/s/g;
		if ($sbjctname =~ /(.*\S)_+(\S.*)-+(\S.*)/) {
			$sgs = $2;
			$sge = $3;
		}
	} 
	if (/Length = +(.*)\n/) {
		$sbjctlength = $1;
	}
	if (/Score = +(.*) bits +(.*), Expect.* = +(.*)\n/) {
		if ($querystart) {
			print $name, "\t", $program, "\t", $qgs, "\t", $qge, "\t", $sgs, "\t", $sge, "\t", $bits, "\t", $Escore, "\t", $Identities, "\t", $Positives, "\t", $querystart, "\t", $queryend, "\t", ($queryend > $querystart ? "+" : "-"), "\t", $sbjctstart, "\t", $sbjctend, "\t", ($sbjctend > $sbjctstart ? "+" : "-"), "\t", $queryname, "\t", $sbjctname, "\t", $strand1, "\t", $strand2, "\t", $sbjctlength, "\t", $querylength, "\t", $frame1, "\t", $frame2, "\n";
			$querystart = 0;
			$sbjctstart = 0;
			$bits = $1; 
			$Escore = $3;
		} else {	
		$bits = $1; 
		$Escore = $3;
		}
	}
	if (/Identities = +(.*) \(+(.*)%\)\n/) {
		$Positives = $2; 
	}
	if (/Strand = +(.*) \/ +(.*)\n/) {
		$strand1 = $1;
		$strand2 = $2;
	}
	if (/Query: ([0-9]*) (.*) ([0-9]*)\n/){
		if ($querystart == 0) {
			$querystart = $1;
		}
		$queryend = $3;
	}
	if (/Sbjct: ([0-9]*) (.*) ([0-9]*)\n/) {
		if ($sbjctstart == 0) {
			$sbjctstart = $1;
		} 
		$sbjctend = $3;
	}
	if (/Reference:/ && $querystart) {
		print $name, "\t", $program, "\t", $qgs, "\t", $qge, "\t", $sgs, "\t", $sge, "\t", $bits, "\t", $Escore, "\t", $Identities, "\t", $Positives, "\t", $querystart, "\t", $queryend, "\t", ($queryend > $querystart ? "+" : "-"), "\t", $sbjctstart, "\t", $sbjctend, "\t", ($sbjctend > $sbjctstart ? "+" : "-"), "\t", $queryname, "\t", $sbjctname, "\t", $strand1, "\t", $strand2, "\t", $sbjctlength, "\t", $querylength, "\t", $frame1, "\t", $frame2, "\n";
		$querystart = 0;
		$sbjctstart = 0;
	}
}	

print $name, "\t", $program, "\t", $qgs, "\t", $qge, "\t", $sgs, "\t", $sge, "\t", $bits, "\t", $Escore, "\t", $Identities, "\t", $Positives, "\t", $querystart, "\t", $queryend, "\t", ($queryend > $querystart ? "+" : "-"), "\t", $sbjctstart, "\t", $sbjctend, "\t", ($sbjctend > $sbjctstart ? "+" : "-"), "\t", $queryname, "\t", $sbjctname, "\t", $strand1, "\t", $strand2, "\t", $sbjctlength, "\t", $querylength, "\t", $frame1, "\t", $frame2, "\n";


