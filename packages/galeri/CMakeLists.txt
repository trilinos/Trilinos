
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Galeri)

#
# B) Set up package-specific options
#

# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#this dir only had "EXTRA_DIST" files in it
#ADD_SUBDIRECTORY(src-fem)  

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example example-fem)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
