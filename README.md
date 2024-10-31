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
 
Additionally MySQL database is necessary for storing the calculated parameters:

1. Install MySQL
2. Open MySQL and run macro
CFStorage.structure.sql

Add to Interactions table the following records:
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (1,0,0,1);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (2,0,1,0);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (3,0,1,1);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (4,1,0,0);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (5,1,0,1);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (6,1,1,0);
INSERT INTO Interactions (interactionID,strong,coulomb,quantumStatistics) VALUES (7,1,1,1);

select * from Interactions; powinno dać tabelkę:
+---------------+--------+---------+-------------------+
| interactionID | strong | coulomb | quantumStatistics |
+---------------+--------+---------+-------------------+
|             1 |      0 |       0 |                 1 |
|             2 |      0 |       1 |                 0 |
|             3 |      0 |       1 |                 1 |
|             4 |      1 |       0 |                 0 |
|             5 |      1 |       0 |                 1 |
|             6 |      1 |       1 |                 0 |
|             7 |      1 |       1 |                 1 |
+---------------+--------+---------+-------------------+


3. Before running CorrFit, in the config file add MySQL info:
     mDBHost = sRPInstance->getPar("StorageDBHost");
      mDBUser = sRPInstance->getPar("StorageDBUser");
      mDBPass = sRPInstance->getPar("StorageDBPass");

Please add parameters: StorageDBHost StorageDBUser StorageDBPass to the parameter file

 
