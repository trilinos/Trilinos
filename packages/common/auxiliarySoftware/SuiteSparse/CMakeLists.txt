

#
# A) Define the package
#

TRIBITS_PACKAGE(TrilinosSS)

#
# B) Set up package-specific options
#

# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
