// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "gtest/gtest.h"
#include "stk_util/parallel/BroadcastArg.hpp"  // for BroadcastArg
#include "stk_util/parallel/Parallel.hpp"      // for MPI_Comm_rank, MPI_COMM_WORLD
#include "stk_util/stk_config.h"               // for STK_HAS_MPI
#include <algorithm>                           // for copy
#include <iostream>                            // for operator<<, basic_ostream::operator<<, bas...
#include <memory>                              // for allocator_traits<>::value_type
#include <string>                              // for string, basic_string, char_traits
#include <vector>                              // for vector


#if defined(STK_HAS_MPI)  // means that MPI is available

namespace {

    void makeArg(const std::vector<std::string>& strings, int& argc, char**& argv) {
        argc = strings.size();
        if(argc == 0) {
            argv=nullptr;
        } else {
            argv = new char*[argc];
        }
        for(int i=0; i<argc; ++i) {
            argv[i] = new char[strings[i].size()+1];
            std::copy(strings[i].begin(), strings[i].end(), argv[i]);
            argv[i][strings[i].size()] = '\0';
        }
    
    }
    void delArg(int& argc, char**& argv) {
        for(int i=0; i<argc; ++i) {
            delete [] argv[i];
        }        
        delete [] argv;
        argv=nullptr;
        argc=0;
    }


TEST(BroadcastArg, emptyList)
{
    int parallel_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &parallel_rank);
    int argc=0;
    char** argv = nullptr;
    if(parallel_rank==0) {
       makeArg({}, argc, argv);
    }
    stk::BroadcastArg bcast(MPI_COMM_WORLD, argc, argv);
    delArg(argc, argv);
    EXPECT_EQ(bcast.m_argc, 0);
    EXPECT_EQ(bcast.m_argv, nullptr);    
}

TEST(BroadcastArg, oneEntry)
{
    int parallel_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &parallel_rank);
    int argc=0;
    char** argv = nullptr;
    if(parallel_rank==0) {
       makeArg({"hello"}, argc, argv);
    }
    stk::BroadcastArg bcast(MPI_COMM_WORLD, argc, argv);

    std::cout<<parallel_rank<<": AAA: "<<bcast.m_argv[0]<<std::endl;


    delArg(argc, argv);
    EXPECT_EQ(bcast.m_argc, 1);
    EXPECT_EQ(std::string(bcast.m_argv[0]),"hello");
}

TEST(BroadcastArg, twoEntry)
{
    int parallel_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &parallel_rank);
    int argc=0;
    char** argv = nullptr;
    if(parallel_rank==0) {
        makeArg({"hello","world"}, argc, argv);
    }
    stk::BroadcastArg bcast(MPI_COMM_WORLD, argc, argv);

    delArg(argc, argv);
    EXPECT_EQ(bcast.m_argc, 2);
    EXPECT_EQ(std::string(bcast.m_argv[0]),"hello");
    EXPECT_EQ(std::string(bcast.m_argv[1]),"world");
}

TEST(BroadcastArg, emptyStrings)
{
    int parallel_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &parallel_rank);
    int argc=0;
    char** argv = nullptr;
    if(parallel_rank==0) {
        makeArg({"", "hello", "", "world", ""}, argc, argv);
    }
    stk::BroadcastArg bcast(MPI_COMM_WORLD, argc, argv);

    delArg(argc, argv);
    EXPECT_EQ(bcast.m_argc, 5);
    EXPECT_EQ(std::string(bcast.m_argv[0]),"");
    EXPECT_EQ(std::string(bcast.m_argv[1]),"hello");
    EXPECT_EQ(std::string(bcast.m_argv[2]),"");
    EXPECT_EQ(std::string(bcast.m_argv[3]),"world");
    EXPECT_EQ(std::string(bcast.m_argv[4]),"");
}


} 

#endif // STK_HAS_MPI

