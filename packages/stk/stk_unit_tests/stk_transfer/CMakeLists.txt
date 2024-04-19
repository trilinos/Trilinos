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

FILE(GLOB SOURCES *.cpp)

if(NOT STK_ENABLE_STKMiddle_mesh)
  LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/UnitTestConservativeTransfer.cpp")
endif()

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE(stk_transfer_utest
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    COMM serial mpi
    NOEXEPREFIX
  )

  TRIBITS_ADD_TEST(stk_transfer_utest
    ARGS "--gtest_filter=GeomXferImpl.coarseSearchImpl_emptyDomain"
    COMM serial mpi
    NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
    NOEXEPREFIX
  )
else()
  add_executable(stk_transfer_utest ${SOURCES})
  target_link_libraries(stk_transfer_utest stk_transfer)
  target_link_libraries(stk_transfer_utest stk_unit_test_utils)
  if (STK_ENABLE_STKMiddle_mesh)
    target_link_libraries(stk_transfer_utest stk_transfer_fixtures)
  endif()
  target_link_libraries(stk_transfer_utest stk_util_parallel)
  target_link_libraries(stk_transfer_utest stk_unit_main)

  target_include_directories(stk_transfer_utest PUBLIC
    ${STK_TOPLEVEL_SOURCE_DIR}/stk_unit_test_utils
  )

  add_test("stk_transfer_utest" sh -c "mpiexec --np 2 ${CMAKE_CURRENT_BINARY_DIR}/stk_transfer_utest --gtest_filter=GeomXferImpl.coarseSearchImpl_emptyDomain")
endif()

