# Put the name of the output directory and plot here

# Whats the directory containing the datafiles?
outdir=plots
# Whats the plot going to be called?
plotname=plot1

WMAP=2
SN=3
PLANCK=4
HUBBLE=5
DFGRS=6
SDSS=7
WZ=8
B=9
BB=10

gnuplot << EOF

set term postscript landscape color enhanced 20
set output "$outdir/$plotname.eps"
set yr[0:1.1]
set xl "w"
set yl "Likelihood (normalised)"
filename="out/chi2.data"
unset ytics
plot filename u 1:$WMAP w l lw 2 t "WMAP",\
"" u 1:$SN w l lw 2 t "SN",\
"" u 1:$PLANCK w l lw 2 t "Planck",\
"" u 1:$HUBBLE w l lw 2 t "Hubble",\
"" u 1:$DFGRS w l lw 2 t "6dfGRS",\
"" u 1:$SDSS w l lw 2 t "SDSS",\
"" u 1:$WZ w l lw 2 t "WiggleZ",\
"" u 1:$B w l lw 2 t "Boss9",\
"" u 1:$BB w l lw 2 t "Boss11"


EOF

echo "Done plotting"
epstopdf --autorotate=All $outdir/$plotname.eps
echo "Conversion to PDF complete"