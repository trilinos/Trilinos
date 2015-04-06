

#
# A) Define the package
#

TRIBITS_PACKAGE(Moertel)

#
# B) Set up package-specific options
#

# For some reason, cannot set the default to ${${PROJECT_NAME}_ENABLE_SEACAS}
TRIBITS_ADD_OPTION_AND_DEFINE( 
  Moertel_ENABLE_Exodus
  MOERTEL_HAVE_EXODUS
  "Enable Exodus mesh support in Moertel"
#  ${${PROJECT_NAME}_ENABLE_SEACAS}
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE( 
  Moertel_ENABLE_Morkon
  MOERTEL_HAVE_MORKON
  "Enable new Kokkos MorKon subpackage support"
  OFF 
)

IF(Moertel_ENABLE_Exodus)

    INCLUDE_DIRECTORIES(${SEACAS_INCLUDE_DIRS})

ENDIF()


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
IF(Moertel_ENABLE_Morkon)
  ADD_SUBDIRECTORY(morkon)
ENDIF()

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
