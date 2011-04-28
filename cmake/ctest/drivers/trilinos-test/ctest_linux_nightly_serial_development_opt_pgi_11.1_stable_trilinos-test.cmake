
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.trilinos-test.pgi11.cmake")

#
# Set the options specific to this build case
#

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME SERIAL_OPT_DEV1)
SET(CTEST_TEST_TYPE EXPERIMENTAL)
#SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)

SET(Trilinos_PACKAGES FEApp Mesquite Didasko MOOCHO TrilinosCouplings Moertel NOX RBGen Stratimikos ML Komplex Ifpack Trios Pamgen Amesos Galeri AztecOO Claps Pliris OptiPack EpetraExt Triutils GlobiPack Shards Zoltan Kokkos RTOp ThreadPool TrilinosFramework)

SET( EXTRA_EXCLUDE_PACKAGES Rythmos Piro)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
  "-DTrilinos_DATA_DIR:STRING=$ENV{TRILINOSDATADIRECTORY}"
  "-DTrilinos_ENABLE_Rythmos:BOOL=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
