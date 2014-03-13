
To run evolver:

> make clean
> make
> ./sfevol

By default, outputs go into "conf/" (can be changed in code)


To plot output
First time:
> chmod +x plot-w-Om.sh

Every subsequent time
> ./plot-w-Om.sh DIR FILEID


Installing boost on MacOS
- MacPorts
sudo port install boost
(will probably need 
cairomm
pix-buf2
poppler
.. install each via sudo port install X)

