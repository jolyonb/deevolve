#!/bin/sh


outdir=logs
plotdir=plots
# Whats the plot going to be called?
plotname=plot_$1

p1=1
p2=2
WL=3
PL=4
SL=5
CL=13

################################################

filename=$outdir/run$1d.dat
echo $filename

gnuplot << EOF
 
set term postscript landscape color enhanced 20 solid
set output "$plotdir/$plotname.eps"

set contour base
set view map
unset key
set xlabel "w_0"
set ylabel "w_a"
unset surface

set style line 1 lc rgb "black" lw 4
set style line 2 lc rgb "black" lw 4
set style line 3 lc rgb "blue" lw 4
set style line 4 lc rgb "blue" lw 4
set style line 5 lc rgb "green" lw 4
set style line 6 lc rgb "green" lw 4
set style line 7 lc rgb "red" lw 4
set style line 8 lc rgb "red" lw 4


set style increment user
set cntrparam levels 2
splot "$filename" u $p1:$p2:$WL w l t "WMAP",\
"" u $p1:$p2:$PL w l t "Planck",\
"" u $p1:$p2:$SL w l t "SN1a",\
"" u $p1:$p2:$CL w l t "Combined"


EOF

echo "Done plotting"
epstopdf --autorotate=All $plotdir/$plotname.eps
echo "Conversion to PDF complete"
