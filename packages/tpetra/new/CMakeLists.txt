INCLUDE(CombinedOption)
INCLUDE(TribitsETISupport)

#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(New)

# This subpackage inherits explicit template instantiation (ETI)
# support from TpetraCore, on which it depends.  Thus, we don't need
# to add an ETI option or enable ETI support here.

#
# Add libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
