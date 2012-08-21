INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddTest)

#
#  When adding test directories if you want the test to be distributed
#  simply add the subdirectory like normal. If you do not want the test
#  to be distributed then you need to test to add it like:
#  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/<subdir>)
#    ADD_SUBDIRECTORY(<subdir>)
#  ENDIF()
#
#  This will prevent configure errors from not being able to find the
#  test's directory when it has been removed from distributions
#


#IF(PERL_EXECUTABLE)
#  SET(ZOLTAN_TEST_DRIVER CMND ${PERL_EXECUTABLE} ARGS ../ctest_zoltan.pl)
#ELSE()
#  MESSAGE(STATUS "Skipping Zoltan tests because PERL_EXECUTABLE is not found.")
#  RETURN()
#ENDIF()

#IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ch_brack2_3)
#  ADD_SUBDIRECTORY(ch_brack2_3)
#ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LoggerTest
  SOURCES LoggerTest.cpp
  COMM serial mpi
  NOEXEPREFIX
)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TimerTest
  SOURCES TimerTest.c
  COMM serial mpi
  NOEXEPREFIX
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SignalTest
  SOURCES SignalTest.c
  COMM serial mpi
  NOEXEPREFIX
)


IF (TRIOS_ENABLE_TRACING)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TraceTest
    SOURCES TraceTest.c
    COMM serial mpi
    NOEXEPREFIX
  )
ENDIF (TRIOS_ENABLE_TRACING)
