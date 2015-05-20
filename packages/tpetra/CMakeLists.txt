############################################################
# Define the package
############################################################

# NOTE: If the TBB (Intel's Threading Building Blocks) TPL is enabled,
# it's probably not possible to build Tpetra with shadowing warnings.
#TRIBITS_PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS CLEANED )
TRIBITS_PACKAGE( Tpetra CLEANED )

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)


# mfh 14 Apr 2015: We don't have a test build for Kokkos or Tpetra
# with C++11 disabled, so I am REQUIRING C++11.  If you don't like
# this, set me up a test machine with C++11 OFF.  I don't have time to
# do it myself.
#
# I left a special case for CUDA because I'm not sure whether TriBITS
# counts C++11 as enabled in that case; please feel free to revise if
# you know better.
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
ASSERT_DEFINED(TPL_ENABLE_CUDA)
IF (NOT ${PROJECT_NAME}_ENABLE_CXX11 AND NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "As of Trilinos 12.0, C++11 is REQUIRED when building Kokkos and its downstream packages, which include Tpetra.  You MUST enable C++11 when building these packages.  Trilinos should enable this by default, if your C++11 compiler supports C++11 (that is, if it's not too old).  If you set Trilinos_ENABLE_CXX11=OFF, then you are disabling C++11.  DON'T DO THAT if you want to build Kokkos or its downstream packages, including Tpetra.  If you aren't happy about that, then you need to disable Kokkos explicitly, by setting Trilinos_ENABLE_Kokkos=OFF.  This WILL disable Tpetra and its downstream packages.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Kokkos_Refactor
  TPETRA_ENABLE_KOKKOS_REFACTOR
  "Enable the \"Kokkos refactor\" version of Tpetra.  In order to enable this, you must also enable the following subpackages that do not currently build by default: KokkosCore, KokkosContainers, KokkosAlgorithms (all in Kokkos package), TeuchosKokkosCompat, TeuchosKokkosComm (both in Teuchos package), and TpetraKernels (in Tpetra package).  This option will be enabled by default in the next Trilinos release.  Once that happens, you will no longer need to enable all of those subpackages explicitly."
  ON
  )

############################################################
# Process subpackages
############################################################

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# Identify what files to exclude when making a source tarball.
#
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  core/doc/CodingGuidelines
  core/doc/TpetraDesign
  )

#
# Do standard package postprocessing
#
TRIBITS_PACKAGE_POSTPROCESS()
