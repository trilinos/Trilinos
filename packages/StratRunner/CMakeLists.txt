INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(StratRunner)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)


PACKAGE_POSTPROCESS()
