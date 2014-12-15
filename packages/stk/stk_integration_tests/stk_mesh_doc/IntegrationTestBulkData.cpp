// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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



#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <algorithm>                    // for sort
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>  // for test_shift_ring
#include <unit_tests/Setup8Quad4ProcMesh.hpp>
#include "unit_tests/UnitTestMeshUtils.hpp"
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include <stk_mesh/base/FieldBLAS.hpp>  // for stk::mesh::field_fill
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>
#include <unit_tests/BulkDataTester.hpp>

namespace
{

//BEGIN_DOC1
TEST(BulkData_test, use_entity_ids_for_resolving_sharing)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::unit_test::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if(stkMeshBulkData.parallel_size() == 2)
    {
        std::string exodusFileName = getOption("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_bulk_data(stkMeshBulkData);
            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
            exodusFileReader.create_input_mesh();
            exodusFileReader.populate_bulk_data();
        }
    }

    stkMeshBulkData.set_use_entity_ids_for_resolving_sharing(false);
    EXPECT_FALSE(stkMeshBulkData.use_entity_ids_for_resolving_sharing());

    stkMeshBulkData.set_use_entity_ids_for_resolving_sharing(true);
    EXPECT_TRUE(stkMeshBulkData.use_entity_ids_for_resolving_sharing());
}

TEST(BulkData_test, testTwoDimProblemForSharingOfDifferentEdgesWithSameNodesFourProc)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    const int spatialDim = 2;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::unit_test::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if ( stkMeshBulkData.parallel_size() == 4 )
    {
        std::string exodusFileName = getOption("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_bulk_data(stkMeshBulkData);
            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
            exodusFileReader.create_input_mesh();
            // With in populate_bulk_data, the option set_use_entity_ids_for_resolving_sharing is set to true
            exodusFileReader.populate_bulk_data();
        }

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(30u, globalCounts[stk::topology::EDGE_RANK]);
    }
}
//END_DOC1

TEST(BulkData_test, test3DProblemSharingOfDifferentFacesWithSameNodesTwoProc)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::unit_test::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if ( stkMeshBulkData.parallel_size() == 2 )
    {
        std::string exodusFileName = getOption("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_bulk_data(stkMeshBulkData);
            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
            exodusFileReader.create_input_mesh();
            exodusFileReader.populate_bulk_data();
        }

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(2u, globalCounts[stk::topology::FACE_RANK]);
    }
}

TEST(BulkData_test, test3DProblemSharingOfDifferentFacesWithSameNodesOneProc)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::mesh::unit_test::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);
    if ( stkMeshBulkData.parallel_size() == 1 )
    {
        std::string exodusFileName = getOption("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_bulk_data(stkMeshBulkData);
            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
            exodusFileReader.create_input_mesh();
            exodusFileReader.populate_bulk_data();
        }

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(2u, globalCounts[stk::topology::FACE_RANK]);
    }
}

}
// empty namespace
