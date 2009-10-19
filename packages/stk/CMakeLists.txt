
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

# ENABLE_SHADOWING_WARNINGS  not possible, because of TBB
PACKAGE( STK )

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_STK_Trilinos
  "Building as a Trilinos package"
  ON
  )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(stk_util)

ADD_SUBDIRECTORY(stk_mesh)


#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()

