TRIBITS_SUBPACKAGE(BDDC)

# BDDC package is deprecated in February 2022; to be removed by next release
MESSAGE(WARNING "***** WARNING:  Package BDDC is deprecated in February 2022; it will be removed by next release *****")

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()
      
ADD_SUBDIRECTORY(src)

IF (Trilinos_ENABLE_Gtest)
  TRIBITS_ADD_TEST_DIRECTORIES(test)
  TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
ELSE()
  MESSAGE(STATUS "BDDC disables tests/examples as Trilinos disables Gtest")   
ENDIF()



TRIBITS_SUBPACKAGE_POSTPROCESS()
