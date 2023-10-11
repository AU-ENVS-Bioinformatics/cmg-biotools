#! /usr/bin/perl -w
# Authors: Peter Fisher Hallin
# For license see /usr/biotools/CMG-biotools.license

#use lib qw(/home/people/pfh/scripts/roseplot/PostScript-Simple-0.06p3/lib);
use PostScript::Simple;
use Getopt::Long;

####################################################################
# Different global variables
my $dia = 60;
my $pi = 3.14159265;
my $arbitraryFontHeight2mmScaleFactor = 10;
my @inputY;
my @labels;
my $donutEffect = 0.3;
####################################################################


if (defined($help)) {
	system "echo /home/people/pfh/scripts/roseplot/roseplot.help";
	die("\n");
}

&GetOptions 	(
				"axissize:f", \$scaleFontSize,
				"Tsize:f", \$headlineFontSize,
				"STsize:f", \$subHeadlineFontSize,
				"labelFontSize:f", \$labelFontSize,
				"axistitle:s", \$axistitle,
				"T:s", \$T,
				"ST:s", \$ST,
				"gridcircles:f", \$gridcircles,
				"fcolor:s", \$fcolor,
				"steps:f", \$steps,
				"output:s", \$output,
				"donut:f", \$donutEffect,
				"style:s", \$style, # 1: fill, # 2: outline # 3: diamond # 4: Box # 5: cross # 6: circles
				"cenX:f", \$cenX,
				"cenY:s", \$cenY,
				'help' => \$help,
				"Xcol:s", \$Xcol,
				"Ycol:s", \$Ycol,
				);

$scaleFontSize = 9 unless (defined($scaleFontSize));
$headlineFontSize = 20 unless (defined($headlineFontSize));
$subHeadlineFontSize = 15 unless (defined($subHeadlineFontSize));
$Xcol = 1 unless (defined($Xcol));
$Ycol = 2 unless (defined($Ycol));
$cenX = 100 unless (defined($cenX));
$cenY = 150 unless (defined($cenY));
$steps = 2 unless (defined($steps));
$color = "red" unless (defined($color));
$gridcircles = 5 unless (defined($gridcircles));
$T = "" unless (defined($T));
$ST = "" unless (defined($ST));
$axistitle = "" unless (defined($axistitle));
$style = "1" unless (defined($style));
$labelFontSize = "12" unless (defined($labelFontSize));
$fcolor = 'blue' unless defined ($fcolor);
#warn "Using axis title '$axistitle'.\n";
#warn "Using title '$T'.\n";
#warn "Using sub title '$ST'.\n";
#warn "Using $steps step(s) in circle perimeter labels.\n";
#warn "Using $gridcircles gridcircle(s).\n";
#warn "Using circle origo at ($cenX,$cenY).\n";
#warn "Using column $Xcol for obtaining X-data.\n";
#warn "Using column(s) $Ycol for obtaining Y-data.\n";

my @cen = ($cenX,$cenY);

if (defined($output)) {
	#warn "Using output filename: '$output'.\n";
} else {
	die ("Please specify output filename with option -output.\n");
}

if (defined($labelsString)) {
	$labelsString =~ s/\(|\)//g;
	@labels = split(/,/,$labelsString);
	warn "Detected ".($#labels+1)." elements in label vector.\n";
}

my $i = -1;
@cols = split(/,/,$Ycol);
while (defined($line=<>))  {
	chomp($line);
	next if ($line =~ /^#/);
	@split = split(/\t/,$line);
	$i++;
	if (defined($Xcol)) {
		push @inputX, $split[$Xcol-1]; 
	} else {
		push @inputX, $i;
	}
	for ($j = 0 ; $j <= $#cols ; $j++) {
		$inputY[$i][$j] = $split[$cols[$j]-1];
	}
}

my $yLbound = $inputY[0][0];
my $yUbound = $inputY[0][0];

for ($i = 0 ; $i <= $#inputX ; $i++) {
	for ($j = 0 ; $j <= $#cols ; $j++) {    
		my $value = $inputY[$i][$j];
		$yLbound = $value + 0 if ($yLbound >= $value);
		$yUbound = $value + 0 if ($yUbound <= $value);
	}
}


# create a new PostScript object
#$p = new PostScript::Simple(landscape=> 0,xsize => 297,ysize => 210,colour => 1,eps => 0,units => "mm");

$p = new PostScript::Simple(landscape=> 0 ,papersize => "A4",colour => 1,eps => 0,units => "mm");

# create a new page
$p->newpage;

# translates an (x,y) coordinate to polar
sub polar {
	$pi = 3.14159265;
	$x = $_[0];
	$localDia = $_[1];
	$angle = ($x/($#inputX+1))*(2*$pi);
	return (sin($angle)*$localDia+$cen[0],cos($angle)*$localDia+$cen[1]);
}

# create inner and outer thick circle
$p->setlinewidth( 0.4 );
$p->circle($cen[0],$cen[1],$dia);
$p->circle($cen[0],$cen[1],dia2donutDia(0));

# create radial grid lines
$p->setcolour(210,210,210);
$p->setlinewidth( 0.00001 );

@labels = @inputX;
if (@labels) {
	for ($i = 0 ; $i <= $#labels ; $i++) {
		$pi = 3.14159265;
		$angle = ($i/($#labels+1))*(2*$pi);
		($a,$b) = (sin($angle)*$dia+$cen[0],cos($angle)*$dia+$cen[1]);
		($A,$B) = (sin($angle)*dia2donutDia(0)+$cen[0],cos($angle)*dia2donutDia(0)+$cen[1]);
		$p->line($A,$B, $a,$b); 
	}
}
    
for ($i = 1  ; $i < $gridcircles ; $i++)  {
	$p->circle($cen[0],$cen[1],dia2donutDia($dia/($gridcircles)*$i));
}


for ($j = 0 ; $j <= $#cols ; $j++) {
	# draw mid polygon (rose)
	my @polygon;
	$p->setlinewidth( "0.37" );
	for ($i = 0 ; $i <= $#inputX ; $i++) {    
		$tempDia = (($inputY[$i][$j] - $yLbound)/($yUbound-$yLbound))*$dia;
		push @polygon, polar($i,dia2donutDia($tempDia)); 
	}
	$tempDia = (($inputY[0][$j] - $yLbound)/($yUbound-$yLbound))*$dia;
	push @polygon, polar(0,dia2donutDia($tempDia)); 
	@fcolorArray = split(/,/,$fcolor);
	@thisColor = split(/\s+/ , $fcolorArray[$j]);
	if ($#thisColor > 1) {
		$p->setcolour($thisColor[0],$thisColor[1],$thisColor[2]);
	} else {
		$p->setcolour($fcolorArray[$j]);
	}
	@styleArray = split(/,/,$style);
	# 1: fill, # 2: outline # 3: diamond # 4: Box # 5: cross
	if ($styleArray[$j] == 1) {
		$p->polygon({filled=>1}, @polygon); 
	} elsif ($styleArray[$j] == 2) {
		$p->polygon({filled=>0}, @polygon);
	} elsif ($styleArray[$j] == 3) {
		diamonds(@polygon);
	} elsif ($styleArray[$j] == 4) {
		boxes(@polygon);
	} elsif ($styleArray[$j] == 5) {
		cross(@polygon);
	} elsif ($styleArray[$j] == 6) {
		circles(@polygon);
	}
}

$p->setcolour(255,255,255);
$p->circle({filled => 1}, $cen[0],$cen[1],dia2donutDia(0));
$p->setcolour(0,0,0);
$p->setlinewidth( 0.4 );
$p->circle($cen[0],$cen[1],dia2donutDia(0));

sub circles {
	@inp = @_;
	for ( my $i = 2 ; $i <= $#inp ; $i += 2) {
		$p->circle({filled=>1},$inp[$i],$inp[$i+1],1);
	}
}

sub rotate {
	(my $angle,my $x,my $y) = @_;
	return ($x*cos($angle)-$y*sin($angle),$x*sin($angle)+$y*cos($angle));
}

sub boxes {
	@inp = @_;
	for ( my $i = 2 ; $i <= $#inp ; $i += 2) {
		my $angle = -(($i/2)/($#labels+1))*(2*$pi);
		my $size = 1;
		@P = (rotate($angle,-$size,$size),rotate($angle,$size,$size),rotate($angle,$size,-$size),rotate($angle,-$size,-$size));
		@P = ($P[0]+$inp[$i],$P[1]+$inp[$i+1],$P[2]+$inp[$i],$P[3]+$inp[$i+1],$P[4]+$inp[$i],$P[5]+$inp[$i+1],$P[6]+$inp[$i],$P[7]+$inp[$i+1]);
		$p->polygon({filled=>1},@P);
	}
}

sub cross {
	@inp = @_;
	for ( my $i = 2 ; $i <= $#inp ; $i += 2) {
		my $angle = -(($i/2)/($#labels+1))*(2*$pi);
		my $size = 1;
		@P = (rotate($angle,-$size,$size),rotate($angle,$size,$size),rotate($angle,$size,-$size),rotate($angle,-$size,-$size));
		@P = ($P[0]+$inp[$i],$P[1]+$inp[$i+1],$P[2]+$inp[$i],$P[3]+$inp[$i+1],$P[4]+$inp[$i],$P[5]+$inp[$i+1],$P[6]+$inp[$i],$P[7]+$inp[$i+1]);
		$p->line($P[0],$P[1],$P[4],$P[5]);
		$p->line($P[2],$P[3],$P[6],$P[7]);
	}
}
sub diamonds {
	@inp = @_;
	for ( my $i = 2 ; $i <= $#inp ; $i += 2) {
		my $angle = 45-(($i/2)/($#labels+1))*(2*$pi);
		my $size = 1;
		@P = (rotate($angle,-$size,$size),rotate($angle,$size,$size),rotate($angle,$size,-$size),rotate($angle,-$size,-$size));
		@P = ($P[0]+$inp[$i],$P[1]+$inp[$i+1],$P[2]+$inp[$i],$P[3]+$inp[$i+1],$P[4]+$inp[$i],$P[5]+$inp[$i+1],$P[6]+$inp[$i],$P[7]+$inp[$i+1]);
		$p->polygon({filled=>1},@P);
	}
}

sub dia2donutDia {
	return ($dia * $donutEffect ) + ( 1 - $donutEffect ) * $_[0];
}

# write circle labels
$p->setfont("Times-Roman", $labelFontSize);
$p->setlinewidth(0.1 );
for ($n = 0; $n <= $#labels; $n++) {
	my $f = $n/($#labels+1);
	$rest = $n % $steps;
	$dia2 = $dia+$dia*0.08*($rest+1);
	$p->setcolour(0,0,0);
	$p->circletext( {align => "outside"}, $cen[0], $cen[1], $dia2, -($f)*360+90 , $labels[$n]);
	($a,$b) = (sin($f*2*$pi)*$dia*1.02+$cen[0],cos($f*2*$pi)*$dia*1.02+$cen[1]);
	($c,$d) = (sin($f*2*$pi)*$dia2*0.98+$cen[0],cos($f*2*$pi)*$dia2*0.98+$cen[1]);
	$p->setcolour(100,100,100);
	$p->line($a,$b,$c,$d);
}

# write title and subtitle
$p->setfont("Times-Roman", $headlineFontSize);
$p->setcolour(0,0,0);
$p->text({align => "centre"},$cen[0],$dia*$steps*0.08+$cen[1]+1.3*$dia, $T);
$p->setfont("Times-Roman", $subHeadlineFontSize);
$p->setcolour(0,0,0);
$p->text({align => "centre"},$cen[0],$dia*$steps*0.08+$cen[1]+1.3*$dia-(3.1*$headlineFontSize/$arbitraryFontHeight2mmScaleFactor), $ST);

# write mid bound axis ticks
$p->setfont("Times-Roman", $scaleFontSize);
$gridcircles++;
for ($i = 0  ; $i < $gridcircles ; $i++)  {
	if ($i == 0 || $i == $gridcircles) {
		$p->setlinewidth( 0.4 );
	} else {
		$p->setlinewidth( 0.1 );
	}
	$p->line($dia*$steps*0.08+$dia+$cen[0]*1.18  , $cen[1]+dia2donutDia($dia/($gridcircles-1)*$i) , $dia*$steps*0.08+$dia+$cen[0]*1.2  , $cen[1]+dia2donutDia($dia/($gridcircles-1)*$i));
	$p->text({align => "right"},$dia*$steps*0.08+$dia+$cen[0]*1.17,$cen[1]+dia2donutDia($dia/($gridcircles-1)*$i)-$scaleFontSize/$arbitraryFontHeight2mmScaleFactor , sprintf("%0.2f",($i/($gridcircles-1))*($yUbound-$yLbound)+$yLbound));
}

$p->line($dia*$steps*0.08+$dia+$cen[0]*1.2,$cen[1]+dia2donutDia(0),$dia*$steps*0.08+$dia+$cen[0]*1.2,$cen[1]+$dia);

# put vertical axis text 
$p->text({align => "centre",rotate => "90"},$dia*$steps*0.08+$dia+$cen[0]*1.24,$cen[1]+dia2donutDia($dia/2) , $axistitle);

# write the output to a file
$p->output($output);
