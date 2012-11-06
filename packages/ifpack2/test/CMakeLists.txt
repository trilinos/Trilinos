TRIBITS_ADD_LIBRARY(
    ifpack2_tpetra_testinst
    SOURCES TestInstantiations.cpp
    TESTONLY
)

TRIBITS_ADD_TEST_DIRECTORIES(belos)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)
#IF (Ifpack2_ENABLE_Anasazi) 
#  TRIBITS_ADD_TEST_DIRECTORIES(anasazi)
#ENDIF()
