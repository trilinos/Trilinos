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

SET(SOURCES "")
FILE(GLOB SOURCES *.cpp)
#disable because of dependency on data file not in Trilinos
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/readAttributes.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/howToWriteMeshWithInternalSidesets.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/QueryExoVars.cpp")
LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/RenamedInputFields.cpp")

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(stk_io_doc_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    NOEXEPREFIX
  )
else()
  add_executable(stk_io_doc_tests ${SOURCES})
  target_link_libraries(stk_io_doc_tests PUBLIC stk_io_util)
  target_link_libraries(stk_io_doc_tests PUBLIC stk_io)
  target_link_libraries(stk_io_doc_tests PUBLIC stk_unit_test_utils)
  target_link_libraries(stk_io_doc_tests PUBLIC stk_mesh_fixtures)
  target_link_libraries(stk_io_doc_tests PUBLIC stk_unit_main)

  add_test(NAME "stk_io_doc_tests" COMMAND stk_io_doc_tests)
endif()

if(TARGET stk_io_doc_tests)
  target_include_directories(stk_io_doc_tests PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
  )
endif()
