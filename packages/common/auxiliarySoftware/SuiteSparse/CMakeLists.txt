

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
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
