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
To plot output
First time:
	
	> chmod +x plot-w-Om.sh

Every subsequent time

	> ./plot-w-Om.sh DIR FILEID

