

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(Pamgen)


#
# B) Set up package-specific options
#
IF(${PACKAGE_NAME}_ENABLE_Boost)
  SET(HAVE_PAMGEN_BOOST ON )
ENDIF()



#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)


TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
