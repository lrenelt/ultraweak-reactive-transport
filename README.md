This repository contains the code to reproduce the results from 'An optimally stable approximation of reactive transport using discrete test and infinite trial spaces' submitted to FVCA10.

The code requires an existing installation of dune-PDELab (master as of 19-05-2023). An exemplary opts-file is contained in the repository (defaults.opts).

- To reproduce the data for the velocity and pressure field in Fig.1, execute 'solve_darcy_mixed.cc'. The result is given as a vtu-file.

- To reproduce the plots in Fig.2, run 'solve_transport_dg.cc' and solve_transport_normal_eq.cc'. Output is also given as vtu-files.

- To reproduce the convergence results from Fig.3, run 'conv_test_yasp.cc'. Resultsare given as csv-files.
