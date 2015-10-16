INCLUDE(TribitsETISupport)

#
# Define the subpackage
#
TRIBITS_SUBPACKAGE(Kernels)

# Define ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION,
# which defaults to ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Add this subpackage to the list of subpackages that do explicit
# template instantiation (ETI).
TRIBITS_ADD_ETI_SUPPORT()

# If building in debug mode, define the HAVE_TPETRAKERNELS_DEBUG macro.
TRIBITS_ADD_DEBUG_OPTION()

#
# Set up subpackage-specific configuration options
#

#
# "Optimization level" for TpetraKernels computational kernels.  The
# higher the level, the more code variants get generated, and thus the
# longer the compile times.  However, more code variants mean both
# better performance overall, and more uniform performance for corner
# cases.  Values of current interest (24 Apr 2014) are 0, 1, and 2.
#
TRIBITS_ADD_OPTION_AND_DEFINE( KokkosLinAlg_Opt_Level
  KOKKOSLINALG_OPT_LEVEL
  "Optimization level for TpetraKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems."
  "1"
  )

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(perf_test)
TRIBITS_ADD_TEST_DIRECTORIES(unit_test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()

