
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Teko)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(Teko_ENABLE_BLOCK_SMOOTHERS
#  Teko_ENABLE_ML_SMOOTHERS
#  "Enable ml block smoothers."
#  OFF )

# Enable ML block Smoothers
IF(ML_ENABLE_TekoSmoothers)
   SET(Teko_ENABLE_ML_SMOOTHERS ON)
ELSE()
   SET(Teko_ENABLE_ML_SMOOTHERS OFF)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(tests)
PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

PACKAGE_ARCH_EXCLUDE_FILES(
  tests/data/nsjac.mm
  tests/data/nsjac_1.mm
  src/mlutils
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
