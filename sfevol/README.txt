CODE REQUIREMENTS
- GSL library
- C++ Boost (see http://www.boost.org). Just untar in /usr/local/


EVOLVER
----------------
To run evolver:

	> make clean
	> make
	> ./sfevol

User can specify simulation parameters in params.ini
The user can supply their own params file,

	> ./sfevol my_params

If a parameter isnt specified, the defaults will be used (see initalise.cpp)

PLOT OUTPUT
----------------
(1) plot-w-Om.sh
To plot a given runs output with gnuplot. 
Plots a(tau), H(a), w_i(a), Omega_i(a). 

First time:
	
	> chmod +x plot-w-Om.sh

Every subsequent time

	> ./plot-w-Om.sh DIR FILEID


(2) plot_compw.sh
This will plot w_de(a) for many runs. E.g.

	> ./plot_compw.sh run001 run002 run003

will plot w_de(a) for all three, on the same axes. 
The directory containing the data and resulting plot file name
must be specified in the script [SOMETHING TO MODIFY LATER]