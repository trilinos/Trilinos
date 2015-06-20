
#
# Define the subpackage
#

TRIBITS_SUBPACKAGE(KokkosCompat)

#
# Set up subpackage-specific configuration macros
#

#-----------------------------------------------------------------------------

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
