# Authors: Peter Fisher Hallin
# For license see /usr/biotools/CMG-biotools.license
set terminal postscript portrait color 
set output 'test.bias.ps'
set border 3 front linetype -1 linewidth 1.000
set boxwidth 0.9 absolute
set style fill solid 1.00 noborder
set style data histograms

set grid noxtics nomxtics ytics nomytics noztics nomztics nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set key outside right top enhanced
set title "Position specific nucleotid usage" font "arial,20"
set ytics 0.1
set yrange [0 : 0.55]	
set xrange [0.5: 3.5]
plot  '< grep freq test.CodonAaUsage' using ($4):xtic(2) title 'A' linecolor rgb "#009900",\
'' using ($6) title 'T' linecolor rgb "#EE0000",\
'' using ($8) title 'G' linecolor rgb "#000011",\
'' using ($10) title 'C'  linecolor rgb "#000099"

