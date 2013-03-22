INCLUDE(TribitsPackageMacros)
#INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(MeshingGenie)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

#TRIBITS_PACKAGE_POSTPROCESS()
