INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( SEACAS ENABLE_SHADOWING_WARNINGS )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_Netcdf)
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXODUS
  HAVE_SEACAS_EXODUS
  "Enable support for Exodus in ${PACKAGE_NAME}."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_NEMESIS
  HAVE_SEACAS_NEMESIS
  "Enable support for Nemesis in ${PACKAGE_NAME}."
  ON )
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus and Nemesis support not available since TPL_ENABLE_Netcdf was not set.")
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IOSS
  HAVE_SEACAS_IOSS
  "Enable support for Ioss in ${PACKAGE_NAME}."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_FEM_IO
  HAVE_SEACAS_FEM_IO
  "Enable support for the fem_io (helper/interface library for connecting Ioss to application) in ${PACKAGE_NAME}."
  OFF )

#
# C) Add the libraries, tests, and examples
#

IF (${PACKAGE_NAME}_ENABLE_EXODUS)
ADD_SUBDIRECTORY(exodus)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_NEMESIS)
ADD_SUBDIRECTORY(nemesis)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_IOSS)
ADD_SUBDIRECTORY(ioss)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_FEM_IO)
ADD_SUBDIRECTORY(fem_io)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
