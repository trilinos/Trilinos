#
# A) Define the package
#
tribits_package( SimpleCxx  ENABLE_SHADOWING_WARNINGS  CLEANED )

#
# B) Platform-specific checks
#
include(CheckFor__int64)
check_for___int64(HAVE_SIMPLECXX___INT64)

#
# C) Set up package-specific options
#
tribits_add_debug_option()
tribits_add_show_deprecated_warnings_option()

#
# D) Add the libraries, tests, and examples
#
add_subdirectory(src)
tribits_add_test_directories(test)
#tribits_add_example_directories(example)

#
# E) Do standard post processing
#
tribits_package_postprocess()
