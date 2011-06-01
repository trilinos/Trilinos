
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

# Exodus is optional, it can be missing, come in through a TPL, or Trios
IF(TPL_ENABLE_ExodusII)
ADD_OPTION_AND_DEFINE(TPL_ENABLE_ExodusII
  HAVE_MOERTEL_EXODUS
  "Enable ExodusII reader in Moertel requires ExodusII TPL."
  OFF )
  INCLUDE_DIRECTORIES(${TPL_ExodusII_INCLUDE_DIRS})
  #ELSEIF(Trilinos_ENABLE_Trios)
  #ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_Trios
ELSEIF(Trios_ENABLE_EXODUS)
  ADD_OPTION_AND_DEFINE(Trios_ENABLE_EXODUS
  HAVE_MOERTEL_EXODUS
  "Enable ExodusII reader in Moertel requires Trios)."
  OFF )
  INCLUDE_DIRECTORIES(${Trios_INCLUDE_DIRS})
ENDIF()

# Nemesis is optional, it can be missing, come in through a TPL, or Trios
# It is only needed when MPI is enabled, and if Exodus is present
IF(TPL_ENABLE_MPI AND HAVE_MOERTEL_EXODUS)
  IF(TPL_ENABLE_Nemesis)
  ADD_OPTION_AND_DEFINE(TPL_ENABLE_Nemesis
    HAVE_MOERTEL_NEMESIS
    "Enable ExodusII/Nemesis reader in Moertel requires Nemesis TPL."
    OFF )
    INCLUDE_DIRECTORIES(${TPL_Nemesis_INCLUDE_DIRS})
	#ELSEIF(Trilinos_ENABLE_Trios)
	#ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_Trios
  ELSEIF(Trios_ENABLE_NEMESIS)
  ADD_OPTION_AND_DEFINE(Trios_ENABLE_NEMESIS
    HAVE_MOERTEL_NEMESIS
    "Enable ExodusII reader in Moertel requires Trios)."
    OFF )
  ENDIF()
ENDIF()




#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
