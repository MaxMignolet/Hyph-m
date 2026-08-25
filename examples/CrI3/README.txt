The directories 01 through 05 contain the abinit input files. The directory 06_hyph-m contains all the necessary inputs for Hyph-m. There is no need to run Abinit to run directory 06.

If one wants to run directories 01 through 05, access to a cluster is required as well as ~300 GB of disk space. The inputs assume that the pseudopotentials are located under `$HOME/.abinit/pseudos`.

Directory 06 can be run on any laptop.

Directory list:

01_grs:
 - Ground state
02_dfpt: ddk, atomic pert, electric field pert, B-field pert
 - t0  -> t2:  d/dk perturbation
 - t3  -> t26: atomic perturbation
 - t27 -> t29: E-field perturbation
 - t30 -> t33: local B-field perturbation
 - The d/dk perturbation need to be performed before the other pert. Once the d/dk perturbation is computed, the rest can be launched in parallel.
03_freq_resp
 - computation of the Berry curvatures from the 1WFs
04_mrgddb
 - merges the ddb's from 02_dfpt and 03_freq_resp
05_anaddb
 - process the merged ddb in 04_mrgddb to obtain the berry curvatures and corrected stiffness matrices
06_hyph-m
 - computes the hybrid phonon-magnon modes by calling Hyph-m
