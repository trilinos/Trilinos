
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Moertel)

#
# B) Set up package-specific options
#

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Experimental
#  HAVE_KOKKOSCLASSIC_EXPERIMENTAL
#  "Enable experimental."
#  OFF )

IF(Moertel_ENABLE_Mesh)

  # Exodus is optional, it can be missing or come in through SEACAS
  IF(${PACKAGE_NAME}_ENABLE_SEACASExodus)
    SET(HAVE_MOERTEL_EXODUS ${${PACKAGE_NAME}_ENABLE_SEACASExodus})
    SET(HAVE_MOERTEL_SEACAS ON)
    INCLUDE_DIRECTORIES(${SEACASExodus_INCLUDE_DIRS})
  ENDIF()

  # Nemesis is optional, it can be missing or come in through SEACAS
  # It is only needed when MPI is enabled, and if Exodus is present
  IF(TPL_ENABLE_MPI AND HAVE_MOERTEL_EXODUS)
    IF(${PACKAGE_NAME}_ENABLE_SEACASNemesis)
      SET(HAVE_MOERTEL_NEMESIS ${${PACKAGE_NAME}_ENABLE_SEACASNemesis})
    ENDIF()
  ENDIF()

ENDIF()




#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
