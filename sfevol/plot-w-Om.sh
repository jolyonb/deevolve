#!/bin/sh

####################################################
#
# Plot script: w(a), Omega_i(a), H(a), a(\tau)
#   J. Pearson, Durham, March 2014
#
#   J. Pearson, Durham, March 2014
#
####################################################
#
# USEAGE
#
# To run for first time, need to change permissions
#  chmod +x plot-w-Om.sh
#
# To run normally, call with dir & fileID as arguments
#   (note: dont need forward-slash on dir argument)
#
#  ./plot-w-Om.sh DIR FILEID
#
# NOTES
#
# The outputted figure will be named
# "DIR/FILEID_plot.eps"
#  - i.e. in same dir, with same fileID prefix
#
####################################################
echo "Start plotting"
gnuplot << EOF

# -> Will plot Omega_i(a), w_i(a), H(a), a(\tau)

# INPUT DATA FILE
inputdata="$1/$2_out.dat"

# OUTPUT IMAGE FILE
outputfig="$1/$2_plot.eps"

# change column numbers here...

tau_loc = 1
a_loc = 2
H_loc = 3
OM_loc = 7
OR_loc = 8
ODE_loc = 9
wtot_loc = 10
wde_loc = 11

# y-axis limits, for w_i & \Omega_i

wmin = -1.1
wmax = 1.1
Omin = -0.01
Omax = 1.01

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
set yr [Omin:Omax]
plot inputdata u a_loc:OM_loc  w l lw 2 t "{/Symbol W}_{m}", "" u a_loc:OR_loc w l lw 2 t "{/Symbol W}_{r}", "" u a_loc:ODE_loc w l lw 2 t "{/Symbol W}_{DE}"

# Plot w_i(a)

set origin 0.5,0
set xlabel "a"
set ylabel "w_i"
set ytics -1.2,0.4
set format x "10^{%L}"
set yr [wmin:wmax]
plot inputdata u a_loc:wtot_loc  w l lw 2 t "w_{tot}","" u a_loc:wde_loc  w l lw 2 t "w_{de}"


EOF

echo "Done plotting"
epstopdf --autorotate=All $1/$2_plot.eps
echo "Conversion to PDF complete"
