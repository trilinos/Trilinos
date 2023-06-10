include_raw_cmake_build()

tribits_add_test(package1-prg  NOEXEPREFIX  NOEXESUFFIX
  NAME Prg  DIRECTORY ${PACKAGE_BINARY_DIR}/src  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Package1 Deps: tpl1" )
