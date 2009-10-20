
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Tifpack )

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_${PACKAGE_NAME}_Trilinos
  "Building as a Trilinos package"
  ON
  )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


