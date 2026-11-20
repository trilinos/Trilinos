INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Teko)

#
# B) Set up package-specific options
#

IF(${PACKAGE_NAME}_ENABLE_MueLu)
  GLOBAL_SET(HAVE_MUELU ON)
ELSE()
  GLOBAL_SET(HAVE_MUELU OFF)
ENDIF()

# enable development mode
IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  SET(Teko_ENABLE_DEV_MODE ON)
ELSE()
  SET(Teko_ENABLE_DEV_MODE OFF)
ENDIF()

# enable development mode
IF(${PROJECT_NAME}_ENABLE_DEBUG_MSGS)
  SET(Teko_DEBUG_OFF OFF)
ELSE()
  SET(Teko_DEBUG_OFF ON)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(tests)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

TRIBITS_EXCLUDE_FILES(tests/data/nsjac.mm tests/data/nsjac_1.mm)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
