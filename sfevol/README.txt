
To run evolver:

> make clean
> make
> ./sfevol

By default, outputs go into "conf/" (can be changed in code)


To plot output
First time:
> chmod +x plot-w-Om.sh

Every subsequent time
> ./plot-w-Om.sh FILEID

where FILEID is the field ID defined in main.cpp. The plot script outputs to "conf/" by default.

