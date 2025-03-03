# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

if(HAVE_STK_Trilinos)
ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)
endif()

FILE(GLOB SOURCES CONFIGURE_DEPENDS *.cpp
                                    command_line/*.cpp
                                    diag/*.cpp
                                    environment/*.cpp
                                    ngp/*.cpp
                                    parallel/*.cpp
                                    registry/*.cpp
                                    util/*.cpp)

IF (NOT STK_HAS_MPI)
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestDistributedIndex.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestGenerateParallelConsistentIDs.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestGenerateUniqueConsistentIDs.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestIndexGapFinder.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestLoc.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestIreduceReplacement.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestMPIDatatypeGenerator.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestMPITagManager.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestParallelVectorConcat.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestImprobeAll.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestImprobeSelf.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/parallel/UnitTestMPIFinalizationCallback.cpp")
ENDIF()

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_util_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NUM_MPI_PROCS 1-6
    NOEXEPREFIX
  )
else()
  add_executable(stk_util_utest ${SOURCES})
  target_link_libraries(stk_util_utest stk_unit_test_utils)
  target_link_libraries(stk_util_utest stk_util_registry)
  target_link_libraries(stk_util_utest stk_util_diag)
  target_link_libraries(stk_util_utest stk_util_env)
  target_link_libraries(stk_util_utest stk_util_parallel)
  target_link_libraries(stk_util_utest stk_util_command_line)
  target_link_libraries(stk_util_utest stk_unit_main)

  INSTALL(TARGETS stk_util_utest DESTINATION ${STK_INSTALL_BINDIR})

  FOREACH(NUMPROCS 1 2 3 4 6 8)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS})
    add_test(NAME "stk_util_utest_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_util_utest 
             --gtest_filter=*-UnitTestDemangle.UnitTest:FPExceptions.*)
    set_property(TEST "stk_util_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fp_np${NUMPROCS})
    add_test(NAME "stk_util_utest_floating_point_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fp_np${NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_util_utest 
             --gtest_filter=FPExceptions.*)
    set_property(TEST "stk_util_utest_floating_point_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})

  ENDFOREACH()

  FOREACH(NUMPROCS 1)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/demangle_np${NUMPROCS})
    add_test(NAME "stk_util_demangle_utest_np${NUMPROCS}"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/demangle_np${NUMPROCS} 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/stk_util_utest 
             --gtest_filter=UnitTestDemangle.UnitTest)
    set_property(TEST "stk_util_demangle_utest_np${NUMPROCS}" PROPERTY PROCESSORS ${NUMPROCS})

  ENDFOREACH()

endif()

