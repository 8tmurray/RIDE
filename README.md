# loRIDE
logistic Retainment Interval Dose Escalation Design

The R programs in this directory facilitate implementing the loRIDE design and replicating the simulation studies reported in the manuscript.

ride-functions.R contains the R functions for implementing the loRIDE design

comparator-functions.R contains the R functions for implementing the comparator designs, i.e. CRM, BMA-CRM mTPI, Keyboard and BOIN

The remaining files 'source' the above files, so the above files just need to be placed in the working directory of the R session. 

dose-guides.R facilitates replication of Table 1 in the manuscript and Supplemental Table 1 in the Web Supplement.

classic-simulation.R facilitates replication of the classic simulation study in the manuscript.

contemporary-simulation.R facilitates replication of the contemporary simulation study in the manuscript.

The Results directory contains Classic and Contemporary directories that contain the .txt results files for the respective simulation. These results files may be loaded in R using code in the analysis section of the respective simulation.R program. Otherwise, these simulation.R programs will allow for full replication of these results files.
