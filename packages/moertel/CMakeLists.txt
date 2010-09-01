
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Moertel)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Experimental
#  HAVE_KOKKOS_EXPERIMENTAL
#  "Enable experimental."
#  OFF )

IF(TPL_ENABLE_ExodusII)
ADD_OPTION_AND_DEFINE(TPL_ENABLE_ExodusII
  HAVE_MOERTEL_EXODUS
  "Enable ExodusII reader in Moertel requires ExodusII TPL)."
  OFF )
  INCLUDE_DIRECTORIES(${TPL_ExodusII_INCLUDE_DIRS})
ENDIF()

IF(TPL_ENABLE_MPI)
  IF(TPL_ENABLE_Nemesis)
  ADD_OPTION_AND_DEFINE(TPL_ENABLE_Nemesis
    HAVE_MOERTEL_NEMESIS
    "Enable ExodusII/Nemesis reader in Moertel requires Nemesis TPL)."
    OFF )
    INCLUDE_DIRECTORIES(${TPL_Nemesis_INCLUDE_DIRS})
  ENDIF()
ENDIF()


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src-mesh)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
