add_executable("conv_test_yasp" convergence_test_yasp.cc)
add_executable("solve_transport_normal_eq" solve_transport_normal_eq.cc)

# reference solutions
add_executable("solve_darcy_mixed" solve_darcy_mixed.cc)
add_executable("solve_transport_dg" solve_transport_dg.cc)

dune_symlink_to_source_files(FILES debug_parameters.ini)
dune_symlink_to_source_files(FILES conv_test_parameters.ini)
dune_symlink_to_source_files(FILES solver_config.ini)
