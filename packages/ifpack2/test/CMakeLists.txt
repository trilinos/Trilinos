TRIBITS_ADD_LIBRARY(
    ifpack2_tpetra_testinst
    SOURCES TestInstantiations.cpp
    TESTONLY
)

TRIBITS_ADD_TEST_DIRECTORIES(belos)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)
TRIBITS_ADD_TEST_DIRECTORIES(vanka)

# mfh 08 Aug 2014: I found this already disabled.  However, Ifpack2
# still has an optional test dependency on Anasazi.  This is a bad
# idea because eigensolvers normally drive linear solvers and
# preconditioners, not the other way around.  Thus, if we want to keep
# this test, it should go in Anasazi, not in Ifpack2.  Please do not
# reenable the optional test dependency of Ifpack2 on Anasazi.

#IF (Ifpack2_ENABLE_Anasazi) 
#  TRIBITS_ADD_TEST_DIRECTORIES(anasazi)
#ENDIF()
