INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(MueLu)

#
# B) Set up package-specific options
#

PACKAGE_ADD_DEBUG_OPTION()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PROFILING
  HAVE_${PACKAGE_NAME_UC}_PROFILING
  "Enable detailed timing and memory profiling of MueLu."
  OFF )
IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
  SET(HAVE_MUELU_EPETRA_AND_EPETRAEXT ON)
  #TODO? Also Set PACKAGE_NAME_ENABLE_Epetra_and_EpetraExt ?
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(gallery) # iff EXAMPLE or TEST ENABLED
ADD_SUBDIRECTORY(xpetra)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(utils)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)
PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
