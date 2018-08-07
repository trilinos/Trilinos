INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_BUILDS.cmake")

# Disable SEACAS test that fails on mutrino (#2815)
ATDM_SET_ENABLE(SEACASExodus_exodus_unit_tests_DISABLE ON)
