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



#include "../stk_unit_tests/stk_mesh/Setup8Quad4ProcMesh.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include <algorithm>                    // for sort
#include <exception>                    // for exception
#include <gtest/gtest.h>
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stddef.h>                     // for size_t
#include <stdexcept>                    // for logic_error, runtime_error
#include <stdlib.h>                     // for exit
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/FieldBLAS.hpp>  // for stk::mesh::field_fill
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc

namespace
{
using stk::unit_test_util::build_mesh;

//BEGIN_DOC1
TEST(BulkData_test, use_entity_ids_for_resolving_sharing)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if(stkMeshBulkData.parallel_size() == 2)
    {
        std::string exodusFileName = stk::unit_test_util::get_option("-i", "mesh.exo");

        stk::io::fill_mesh(exodusFileName, stkMeshBulkData);
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
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if ( stkMeshBulkData.parallel_size() == 4 )
    {
        std::string exodusFileName = stk::unit_test_util::get_option("-i", "mesh.exo");

        stk::io::fill_mesh(exodusFileName, stkMeshBulkData);

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(15u, globalCounts[stk::topology::EDGE_RANK]);
    }
}
//END_DOC1

TEST(BulkData_test, test3DProblemSharingOfDifferentFacesWithSameNodesTwoProc)
{
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

    if ( stkMeshBulkData.parallel_size() == 2 )
    {
        std::string exodusFileName = stk::unit_test_util::get_option("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
            stk::io::fill_mesh_preexisting(exodusFileReader, exodusFileName, stkMeshBulkData);
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
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);
    if ( stkMeshBulkData.parallel_size() == 1 )
    {
        std::string exodusFileName = stk::unit_test_util::get_option("-i", "mesh.exo");

        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);
            exodusFileReader.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
            stk::io::fill_mesh_preexisting(exodusFileReader, exodusFileName, stkMeshBulkData);
        }

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(1u, globalCounts[stk::topology::FACE_RANK]);
    }
}

TEST(IntegrationTest, PartChangeGenerated)
{
    //demonstrates expected behavior, no shells so part changes work fine. see ticket 13216
    MPI_Comm communicator = MPI_COMM_WORLD;

    const int spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
    stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;
    stk::mesh::MetaData& stkMeshMetaData = stkMeshBulkData.mesh_meta_data();
    if (stkMeshBulkData.parallel_size() != 4) return;
    stk::mesh::Part& partToAdd = stkMeshMetaData.declare_part("urp_part", stk::topology::ELEM_RANK);
    stk::mesh::PartVector add_parts;
    add_parts.push_back(&partToAdd);

    {
        stk::io::StkMeshIoBroker exodusFileReader(communicator);
        exodusFileReader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
        stk::io::fill_mesh_preexisting(exodusFileReader, "generated:1x1x4", stkMeshBulkData);
        stkMeshBulkData.modification_begin();
        stk::mesh::EntityVector entities;
        stk::mesh::Selector select_owned(stkMeshMetaData.locally_owned_part());
        {
            const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, select_owned);
            for(size_t i = 0; i < buckets.size(); ++i)
            {
                for (unsigned j = 0; j < buckets[i]->size(); j++) {
                    stk::mesh::Entity entity = (*buckets[i])[j];
                    entities.push_back(entity);
                }
            }
        }
        for (size_t i = 0; i < entities.size(); ++i) {
            stkMeshBulkData.change_entity_parts(entities[i], add_parts);
        }
        EXPECT_NO_THROW(stkMeshBulkData.modification_end());
    }
}

TEST(IntegrationTest, ShellPartChangeCylinder)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 4) { GTEST_SKIP(); }

  const std::string exodusFileName = "cyl_3block.g";

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;
  stk::mesh::MetaData& stkMeshMetaData = stkMeshBulkData.mesh_meta_data();
  stk::mesh::PartVector add_parts {&stkMeshMetaData.declare_part("urp_part", stk::topology::ELEM_RANK)};

  stk::io::StkMeshIoBroker exodusFileReader(communicator);
  exodusFileReader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  stk::io::fill_mesh_preexisting(exodusFileReader, exodusFileName, stkMeshBulkData);
  stk::mesh::EntityVector entities;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::ELEM_RANK, stkMeshMetaData.locally_owned_part(), entities);
  stkMeshBulkData.modification_begin();
  stkMeshBulkData.change_entity_parts(entities, add_parts);
  EXPECT_NO_THROW(stkMeshBulkData.modification_end());
}

// now using graph to create faces during mesh read, split coincident test cases fail
TEST(IntegrationTest, ShellPartChange2Hexes2Shells)
{
    //demonstrates failing when reading a mesh with shells and changing parts, see ticket 13216
    MPI_Comm communicator = MPI_COMM_WORLD;
    const std::string exodusFileName = "ALefLRA.e";

    const int spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
    stk::mesh::BulkData& stkMeshBulkData = *bulkPtr;
    stk::mesh::MetaData& stkMeshMetaData = stkMeshBulkData.mesh_meta_data();
    if (stkMeshBulkData.parallel_size() != 4) return;
    stk::mesh::Part& partToAdd = stkMeshMetaData.declare_part("urp_part", stk::topology::ELEM_RANK);
    stk::mesh::PartVector add_parts;
    add_parts.push_back(&partToAdd);


    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
    stk::io::fill_mesh_preexisting(exodusFileReader, exodusFileName, stkMeshBulkData);
    stk::mesh::EntityVector entities;
    stk::mesh::Selector select_owned(stkMeshMetaData.locally_owned_part());

    {
        const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEM_RANK, select_owned);
        for(size_t i = 0; i < buckets.size(); ++i)
        {
            for (unsigned j = 0; j < buckets[i]->size(); j++) {
                stk::mesh::Entity entity = (*buckets[i])[j];
                entities.push_back(entity);
            }
        }
    }
    stkMeshBulkData.modification_begin();
    for (size_t i = 0; i < entities.size(); ++i) {
        stkMeshBulkData.change_entity_parts(entities[i], add_parts);
    }
    EXPECT_NO_THROW(stkMeshBulkData.modification_end());
}

}
// empty namespace
