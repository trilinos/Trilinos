// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <gtest/gtest.h>
#include <sstream>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc

// macro for short term solution to not break builds for Trilinos.
// Dependency on unitTestUtils here is the main worry.

#if defined(STK_BUILT_IN_SIERRA)

#include <exampleMeshes/StkMeshFromGeneratedMesh.h>
#include <stk_mesh/base/Comm.hpp>

namespace
{

std::string getGeneratedMeshString(const int xdim, const int ydim, const int zdim)
{
    std::ostringstream oss;
    oss << "generated: " << xdim << "x" << ydim << "x" << zdim;
    return oss.str();
}

//DocTest1
TEST( CommMeshCounts, Serial )
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numprocs = stk::parallel_machine_size(communicator);
    if ( numprocs == 1 )
    {
        const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2);
        unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

        std::vector<size_t> comm_mesh_counts;
        stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts);

        size_t goldNumElements = 10*20*2;
        EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);
    }
}
//DocTest2
TEST( CommMeshCounts, Parallel )
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numprocs = stk::parallel_machine_size(communicator);

    const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2*numprocs);
    unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    std::vector<size_t> comm_mesh_counts;
    stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts);

    size_t goldNumElements = 10*20*2*numprocs;
    EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);
}
//DocTest3
TEST( CommMeshCountsWithStats, Parallel )
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numprocs = stk::parallel_machine_size(communicator);

    const std::string generatedMeshSpec = getGeneratedMeshString(10,20,2*numprocs);
    unitTestUtils::exampleMeshes::StkMeshCreator stkMesh(generatedMeshSpec, communicator);

    std::vector<size_t> comm_mesh_counts;
    std::vector<size_t> min_counts;
    std::vector<size_t> max_counts;

    stk::mesh::comm_mesh_counts(*stkMesh.getBulkData(), comm_mesh_counts, min_counts, max_counts);

    size_t goldNumElements = 10*20*2*numprocs;
    EXPECT_EQ(goldNumElements, comm_mesh_counts[stk::topology::ELEMENT_RANK]);

    size_t goldMinNumElements = 10*20*2;
    EXPECT_EQ(goldMinNumElements, min_counts[stk::topology::ELEMENT_RANK]);

    size_t goldMaxNumElements = goldMinNumElements;
    EXPECT_EQ(goldMaxNumElements, max_counts[stk::topology::ELEMENT_RANK]);
}
//EndDocTest
}
#endif
