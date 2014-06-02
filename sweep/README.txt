
1) generate *.inis
user must modify "usersinc.ini", with structure
PARAM_name PARAM_start PARAM_end PARAM_increment

- each parameter on a new line
- lines starting with "#" are deemed to be comments and are not read in
- that creates *.gni files, dumped into "ginis/"
> ./gini"

2) run evol code on each *.gni
> python evol.py

3) collect the chi2 information into a single file
> ./colchi2

All codes can be recompiled via
./install.sh
All codes can be run in correct order via
./runsweep.sh



