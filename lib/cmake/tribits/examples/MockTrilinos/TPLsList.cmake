tribits_repository_define_tpls(
  MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/tpls/"    PT
  BLAS           cmake/TPLs/    PT
  LAPACK         cmake/TPLs/    PT
  Boost          cmake/TPLs/    ST
  Scotch         cmake/TPLs/    ST
  METIS          cmake/TPLs/    TS
  ParMETIS       cmake/TPLs/    ST
  CppUnit        cmake/TPLs/    ST
  ADOLC          cmake/TPLs/    ST
  ADIC           cmake/TPLs/    EX
  TVMET          cmake/TPLs/    ST
  #Zlib           cmake/TPLs/    PT  # Listed in zoltan/cmake/Dependencies.cmake!
  y12m           cmake/TPLs/    ST
  SuperLUDist    cmake/TPLs/    ST
  SuperLU        cmake/TPLs/    ST
  UMFPACK        cmake/TPLs/    ST
  AMD            cmake/TPLs/    TS
  PETSC          cmake/TPLs/    ST
  MUMPS          cmake/TPLs/    ST
  DUMMY          cmake/TPLs/    ST
  )
