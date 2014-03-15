#!/bin/sh

################################################
#
#   Plot script: multiple w_de(a)
#       User calls:
#     ./plot_compw.sh myrunID1 myrunID2 ...
#
#   This plots evolution of w_de for all these runs
#
################################################
################################################
# Put the name of the output directory and plot here

# Whats the directory containing the datafiles?
outdir=logs
# Whats the plot going to be called?
plotname=plot_compw
# What column number is w_de in the data?
wdeloc=13

################################################

for var in "$@"
do
    echo $var
    filename=$filename' '\"''$outdir'/run'$var'.dat'\"' u 2:'$wdeloc' w l lw 2 t '\"''$var''\"','
done

# Need to remove final character from filename
filename=${filename%?}

gnuplot << EOF

set term postscript landscape color enhanced 20
set output "$outdir/$plotname.eps"

wmin=-1.1
wmax=1.1

set log x
set xlabel "a"
set ylabel "w_{de}(i)"
set ytics -1.2,0.4
set format y "%2.1f"
set format x "10^{%L}"
set yr [wmin:wmax]
plot $filename

EOF

echo "Done plotting"
epstopdf --autorotate=All $outdir/$plotname.eps
echo "Conversion to PDF complete"
