#!/usr/bin/perl -w

# Authors: Peter Fisher Hallin
# For license see /usr/biotools/CMG-biotools.license
@_ = split(/[-_]/, $ARGV[0]);
$r = 0.1*substr($_[0], 0, 2);
$g = 0.1*substr($_[0], 2, 2);
$b = 0.1*substr($_[0], 4, 2);
printf "%1.2f %1.2f %1.2f\n", 0.3*$r, 0.3*$g, 0.3*$b unless $r == 1 && $g == 1 && $b == 1;
for ($i = 0; $i < $#_; $i++) {
  $r = 0.1*substr($_[$i], 0, 2);
  $rf = 0.1*substr($_[$i+1], 0, 2)-$r;
  $g = 0.1*substr($_[$i], 2, 2);
  $gf = 0.1*substr($_[$i+1], 2, 2)-$g;
  $b = 0.1*substr($_[$i], 4, 2);
  $bf = 0.1*substr($_[$i+1], 4, 2)-$b;
  $j = 0;
  while ($j < 0.895) {
    printf "%1.2f %1.2f %1.2f\n", 0.9*$r+$j*$rf, 0.9*$g+$j*$gf, 0.9*$b+$j*$bf;
    if ($rf > 0 || $gf > 0 || $bf > 0) {
      if ($j < 0.45) {
        $j += 0.1;
      }
      elsif ($j < 0.825) {
        $j += 0.05;
      }
      else {
        $j += 0.01;
      }
    }
    else {
      if ($j < 0.045) {
        $j += 0.01;
      }
      elsif ($j < 0.375) {
        $j += 0.05;
      }
      else {
        $j += 0.10;
      }
    }
  }
}
$r = 0.1*substr($_[$i], 0, 2);
$g = 0.1*substr($_[$i], 2, 2);
$b = 0.1*substr($_[$i], 4, 2);
printf "%1.2f %1.2f %1.2f\n", 0.9*$r, 0.9*$g, 0.9*$b;
printf "%1.2f %1.2f %1.2f\n", 0.3*$r, 0.3*$g, 0.3*$b unless $r == 1 && $g == 1 && $b == 1;
