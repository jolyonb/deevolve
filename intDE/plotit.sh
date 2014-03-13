#!/bin/sh

####################################################
#
# Plot script: w(a), Omega_i(a), H(a), a(\tau)
#   J. Pearson, Durham, March 2014
#
####################################################
# USEAGE
#
# To run for first time, need to change permissions
#  chmod +x plot-w-Om.sh
#
# To run normally, call with fileID
#
#  ./plot-w-Om.sh fileID
#
# NOTES
# 
# This assumes the output directory is "conf/"
#  - this can be modified in "outputfig" below
#
# The outputted figure will be named
# "conf/fileID_plot.eps" 
#  - i.e. in same dir, with same fileID prefix
#
####################################################

echo "Start plotting"

gnuplot << EOF


# -> Will plot Omega_i(a), w_i(a), H(a), a(\tau) 

# INPUT DATA FILE
inputdata="logs/run$1.dat"

# OUTPUT IMAGE FILE
outputfig="logs/run$1_plot.eps"

# Expected data format:
# column 1:  	tau
# column 2: 	a
# column 4: 	H
# column 9: 	Omega_M
# column 10: 	Omega_R	
# column 12: 	Omega_DE
# column 13: 	w_tot
# column 16: 	w_de

# If anything different, can change column numbers here...

tau_loc=1
a_loc=2
H_loc=4
OM_loc=9
OR_loc=10
ODE_loc=12
wtot_loc=13
wde_loc=16

# Do the plotting
set term postscript landscape color enhanced 20
set output outputfig

set multiplot
set size 0.5,0.5
set key spacing 1.5
set log x

# Plot a(\tau)

set origin 0,0.5
set xlabel "{/Symbol t}"
set ylabel "a({/Symbol t})"
set log
set format x "10^{%L}"
set format y "10^{%L}"
plot inputdata u tau_loc:a_loc  w l lw 2 notitle

# Plot H(a)

set origin 0.5,0.5
set xlabel "a"
set ylabel "H(a)"
set format x "10^{%L}"
set format y "10^{%L}"
plot inputdata u a_loc:3  w l lw 2 notitle

unset log
set format y "%2.1f"
set log x

# Plot Omega_i(a)

set origin 0,0
set xlabel "a"
set ylabel "Density fractions"
set format x "10^{%L}"
set ytics 0,0.2
plot inputdata u a_loc:OM_loc  w l lw 2 t "{/Symbol W}_{m}", "" u a_loc:OR_loc w l lw 2 t "{/Symbol W}_{r}", "" u a_loc:ODE_loc w l lw 2 t "{/Symbol W}_{DE}"

# Plot w_i(a)

set origin 0.5,0
set xlabel "a"
set ylabel "w_i"
set ytics -1.2,0.4
set format x "10^{%L}"
plot inputdata u a_loc:wtot_loc  w l lw 2 t "w_{tot}","" u a_loc:wde_loc  w l lw 2 t "w_{de}"


EOF

echo "Done plotting, converting to PDF"
epstopdf --autorotate=All logs/run$1_plot.eps
echo "Conversion to PDF complete"

