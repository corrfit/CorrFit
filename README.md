# CorrFit

Two install CorrFit two directories are needed:
* CorrFit,
* Spherical Harmonics,

as well as root installation.

First, Spherical Harmonics directory has to be compiled.
Afterwards, it is recommended to use these settings:

export MYSQLSTORAGE=1 

export C2MRANDOM=1

export FULLOPT=1

MYSQLSTORAGE - CorrFit will use MYSQL for storing calculations instead of local root file
C2MRANDOM - each thread randomly chooses which cell of the Chi2 map it is going to calculate (important when multithreading)
FULLOPT - if you do not plan to debug the program

And then it is enough to use:
* make
 
Additionally MySQL database is necessary for storing the calculated parameters.
 
