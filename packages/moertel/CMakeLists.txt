
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

IF(Moertel_ENABLE_Mesh)

  # Exodus is optional, it can be missing, come in through a TPL, or SEACAS
  IF(TPL_ENABLE_ExodusII)
    ADD_OPTION_AND_DEFINE(TPL_ENABLE_ExodusII
      HAVE_MOERTEL_EXODUS
      "Enable ExodusII reader in Moertel requires ExodusII TPL."
      OFF )
    INCLUDE_DIRECTORIES(${TPL_ExodusII_INCLUDE_DIRS})
    #ELSEIF(Trilinos_ENABLE_SEACAS)
    #ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_SEACAS
  ELSEIF(${PACKAGE_NAME}_ENABLE_SEACASExodus)
    SET(HAVE_MOERTEL_EXODUS ${${PACKAGE_NAME}_ENABLE_SEACASExodus})
    SET(HAVE_MOERTEL_SEACAS ON)
    INCLUDE_DIRECTORIES(${SEACASExodus_INCLUDE_DIRS})
  ENDIF()

  # Nemesis is optional, it can be missing, come in through a TPL, or SEACAS
  # It is only needed when MPI is enabled, and if Exodus is present
  IF(TPL_ENABLE_MPI AND HAVE_MOERTEL_EXODUS)
    IF(TPL_ENABLE_Nemesis)
      ADD_OPTION_AND_DEFINE(TPL_ENABLE_Nemesis
        HAVE_MOERTEL_NEMESIS
        "Enable ExodusII/Nemesis reader in Moertel requires Nemesis TPL."
        OFF )
      INCLUDE_DIRECTORIES(${TPL_Nemesis_INCLUDE_DIRS})
	  #ELSEIF(Trilinos_ENABLE_SEACAS)
	  #ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_SEACAS
    ELSEIF(${PACKAGE_NAME}_ENABLE_SEACASNemesis)
      SET(HAVE_MOERTEL_NEMESIS ${${PACKAGE_NAME}_ENABLE_SEACASNemesis})
    ENDIF()
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
