#
# Define the list of TPLs, their find module names, and their classification
#
# TPL_NAME:
#
#   The name of the TPL used in the CMake cache variables TPL_ENABLE_${TPL_NAME}
#
# TPL_FINDMOD:
#
#   The name of the find module under Trilinos/cmake/TPLs that is used to get the
#   names of the TPLs.  If left empty '', it will just be set to FindTPL${TPL_NAME}.
#
#
# TPL_CLASSIFICATION:
#
#   PS: Primary Stable TPL
#
#     Primary Stable TPLs are those TPLs that a Trilinos developer must have
#     installed on their machine in order to be able to do Trilinos
#     development.  For example, we require that you have BLAS, LAPACK, and
#     MPI installed in order to do Trilinos development.  These are
#     fundamental dependencies that are needed in order to do precheckin
#     testing.
#
#   SS: Secondary Stable TPL
#
#     Secondary Stable TPLs are those TPLs that are not required in order to
#     be able to develop and test Trilinos before checkins but are none the
#     less offically supported.  Support for SS TPLs is tested as part of the
#     nightly testing process.
#
#   TS: Tertiary Stable TPL
#
#     Tertiary Stable TPLs are those TPLs that are supported TPLs but can not
#     be included in the set of SS TPLs because they may conflicit with other
#     SS Code.  For example, METIS is listed as a TS TPL because it conflicts
#     with ParMETIS which is declared as a SS TPL.
#
#   EX: Experimental TPL
#
#     Experimental TPLs are not offically supported.  They represent
#     experimental capabilities of Trilinos packages.  Support for EX TPLs is
#     never tested as part of the main nightly testing process.  However,
#     package developers are encouraged to set up their own nightly testing
#     for their EX TPLs for their packages.
#
# The default enable for all TPLs is empty "" reguardless of the category.
# The idea is that the enabling of the TPL will be done by the package and
# other enables that the user has to set.
#
# NOTE: The TPLs must be listed in the order of increasing dependencies (if
# such dependencies exist).
#

SET(Trilinos_TPLS_FINDMODS_CLASSIFICATIONS
  Peano           ""    EX
  CUDA            ""    EX
  Thrust          ""    EX
  Cusp            ""    EX
  TBB             ""    EX
  Pthread         ""    SS
  BinUtils        ""    SS
  ARPREC          ""    SS
  QD              ""    SS
  MPI             ""    PS
  BLAS            ""    PS
  LAPACK          ""    PS
  Boost           ""    SS
  Scotch          ""    SS
  METIS           ""    TS
  ParMETIS        ""    SS
  PaToH           ""    SS
  CppUnit         ""    SS
  ADOLC           ""    SS
  ADIC            ""    EX
  TVMET           ""    SS
  MF              ""    SS
  ExodusII        ""    SS
  Nemesis         ""    SS
  XDMF            ""    TS
  Netcdf          ""    SS
  y12m            ""    SS
  SuperLUDist     ""    SS
  SuperLUMT	  ""	SS
  SuperLU         ""    SS
  Zlib            ""    SS
  UMFPACK         ""    SS
  MA28            ""    TS
  AMD             ""    TS
  PETSC           ""    SS
  HYPRE           ""    EX
  BLACS           ""    SS
  SCALAPACK       ""    SS
  MUMPS           ""    SS
  PARDISO_MKL     ""    EX
  Oski            ""    SS
  TAUCS           ""    SS
  ForUQTK         ""    EX
  Dakota          ""    EX
  HIPS            ""    EX
  HDF5            ""    EX
  MATLAB          ""    EX
  CASK            ""    EX
  SPARSKIT        ""    SS
  QT              ""    SS
  gtest           ""    EX
  BoostLib        ""    SS
  OpenNURBS       ""    EX
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
