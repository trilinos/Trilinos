
tribits_extpkg_create_imported_all_libs_target_and_config_file( CUSOLVER
  INNER_FIND_PACKAGE_NAME  CUDAToolkit
  IMPORTED_TARGETS_FOR_ALL_LIBS  CUDA::cusolver )
# Above, the CUDA TPL should have already found CUDAToolkit so we just need to
# grab the target from it to form the CUSPARSE::all_libs target.
