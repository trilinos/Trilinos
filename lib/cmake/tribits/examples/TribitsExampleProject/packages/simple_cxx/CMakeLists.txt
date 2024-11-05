#
# A) Define the package
#
tribits_package( SimpleCxx  ENABLE_SHADOWING_WARNINGS  CLEANED )

#
# B) Platform-specific checks
#
include(CheckFor__int64)
check_for___int64(HAVE_SIMPLECXX___INT64)
tribits_pkg_export_cache_var(HAVE_SIMPLECXX___INT64)

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

# Set a variable that will be used in downstream packages

if (SimpleCxx_ENABLE_SimpleTpl)
  set(simpletplText "simpletpl ")
else()
  set(simpletplText)
endif()

global_set(EXPECTED_SIMPLECXX_AND_DEPS
  "SimpleCxx ${simpletplText}headeronlytpl")
tribits_pkg_export_cache_var(EXPECTED_SIMPLECXX_AND_DEPS)

#
# E) Do standard post processing
#
tribits_package_postprocess()
