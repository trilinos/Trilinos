#
# A) Define the package
#

TRIBITS_PACKAGE(Zoltan2Sphynx ENABLE_SHADOWING_WARNINGS)


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
