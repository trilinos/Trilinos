#
# Define the list of TPLs and their classification
#
#
# TPL classifications are:
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
#     SS Code.  For example, METIS is listed as a TS package because it
#     conflicts with ParMETIS which is declared as a SS TPL.
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
#
# NOTE: The TPLs must be listed in the order of increasing dependencies (if
# such dependencies exist).
#

SET(Trilinos_TPLS_AND_CLASSIFICATIONS
  TBB            EX
  Pthread        SS
  MPI            PS
  BLAS           PS
  LAPACK         PS
  Boost          SS
  Scotch         SS
  METIS          TS
  ParMETIS       SS
  PaToH          SS
  CppUnit        SS
  ADOLC          SS
  TVMET          SS
  MF             SS
  ExodusII       SS
  Nemesis        SS
  ZoltanTpl      TS
  y12m           SS
  SuperLUDist    SS
  SuperLU        SS
  Zlib		 SS
  UMFPACK        SS
  MA28           TS
  AMD            TS
  PETSC          SS
  HYPRE          EX
  BLACS          SS
  SCALAPACK      SS
  MUMPS          SS
  Oski           SS
  TAUCS          SS
  ForUQTK	 EX
  Dakota	 EX
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
