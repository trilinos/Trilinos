
# @HEADER
#  ************************************************************************
#
#                    Trios: Trilinos I/O Support
#                  Copyright 2011 Sandia Corporation
#
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#
#  1. Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
#  3. Neither the name of the Corporation nor the names of the
#  contributors may be used to endorse or promote products derived from
#  this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
#  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
#
#  *************************************************************************
# @HEADER


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
