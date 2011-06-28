INCLUDE(AddOptionAndDefine)

#
# B) Set up package-specific options
#

IF (${PACKAGE_NAME}_ENABLE_Netcdf)
<<<<<<< HEAD
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

  ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXODUS
    HAVE_SEACAS_EXODUS
    "Enable support for Exodus in ${PACKAGE_NAME}."
    OFF )

  ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_NEMESIS
    HAVE_SEACAS_NEMESIS
    "Enable support for Nemesis in ${PACKAGE_NAME}."
    OFF )
=======
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
>>>>>>> SierraToTrilinosBranch
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_IOSS
  HAVE_SEACAS_IOSS
  "Enable support for Ioss in ${PACKAGE_NAME}."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_CHACO
  HAVE_SEACAS_CHACO
  "Enable support for chaco in ${PACKAGE_NAME}."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_APREPRO_LIB
  HAVE_SEACAS_APREPRO_LIB
  "Enable support for aprepro library in ${PACKAGE_NAME}."
  ON )

#
# C) Add the libraries, tests, and examples
#

<<<<<<< HEAD
=======
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EXODUS)
>>>>>>> SierraToTrilinosBranch
IF (${PACKAGE_NAME}_ENABLE_EXODUS)
ADD_SUBDIRECTORY(exodus)
ENDIF()

<<<<<<< HEAD
=======
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NEMESIS)
>>>>>>> SierraToTrilinosBranch
IF (${PACKAGE_NAME}_ENABLE_NEMESIS)
ADD_SUBDIRECTORY(nemesis)
ENDIF()

<<<<<<< HEAD
=======
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_IOSS)
>>>>>>> SierraToTrilinosBranch
IF (${PACKAGE_NAME}_ENABLE_IOSS)
ADD_SUBDIRECTORY(ioss)
ENDIF()

<<<<<<< HEAD
=======
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_CHACO)
>>>>>>> SierraToTrilinosBranch
IF (${PACKAGE_NAME}_ENABLE_CHACO)
ADD_SUBDIRECTORY(chaco)
ENDIF()

<<<<<<< HEAD
=======
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_APREPRO_LIB)
>>>>>>> SierraToTrilinosBranch
IF (${PACKAGE_NAME}_ENABLE_APREPRO_LIB)
ADD_SUBDIRECTORY(aprepro_lib)
ENDIF()
