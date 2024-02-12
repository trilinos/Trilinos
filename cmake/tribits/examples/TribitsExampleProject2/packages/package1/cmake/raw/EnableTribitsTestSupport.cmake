set(Package1_USE_TRIBITS_TEST_FUNCTIONS  OFF  CACHE  BOOL
  "Use TriBITS testing functions")
set(Package1_TRIBITS_DIR  ""  CACHE  PATH
  "Path to TriBITS implementation base dir (e.g. TriBITS/tribits)")
if (Package1_USE_TRIBITS_TEST_FUNCTIONS  AND  Package1_TRIBITS_DIR)
  # Pull in and turn on TriBITS testing support
  include("${Package1_TRIBITS_DIR}/core/test_support/TribitsAddTest.cmake")
  include("${Package1_TRIBITS_DIR}/core/test_support/TribitsAddAdvancedTest.cmake")
  set(Package1_ENABLE_TESTS ON)
endif()
