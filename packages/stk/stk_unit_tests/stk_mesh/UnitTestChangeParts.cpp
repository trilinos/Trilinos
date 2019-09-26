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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_THROW, etc
#include <stdexcept>                    // for runtime_error
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_unit_test_utils/TextMesh.hpp>
#include <string>                       // for string

#include "mpi.h"                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc

namespace
{

TEST(UnitTestChangeParts, test_throw_on_internal_part_change)
{
    const int spatialDim = 3;
    stk::mesh::MetaData metaData(spatialDim);
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::mesh::BulkData bulkData(metaData, communicator);

    std::string generatedMeshSpec = "generated:1x1x4";
    stk::io::fill_mesh(generatedMeshSpec, bulkData);

    stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, 1);

    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts;

    addParts.push_back(&metaData.locally_owned_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    addParts.clear();
    addParts.push_back(&metaData.globally_shared_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    addParts.clear();
    removeParts.push_back(&metaData.locally_owned_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

    removeParts.clear();
    removeParts.push_back(&metaData.globally_shared_part());
    EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);
}

TEST(UnitTestChangeParts, test_batch_part_change)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int p_size = stk::parallel_machine_size( pm );

    if (p_size != 1) {
      return;
    }

    const int spatialDim = 3;
    stk::mesh::MetaData metaData(spatialDim);
    stk::mesh::BulkData bulkData(metaData, pm, stk::mesh::BulkData::NO_AUTO_AURA);

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    stk::mesh::Part& part = metaData.declare_part_with_topology("new_part", stk::topology::NODE);
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, bulkData);

    stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1u);
    EXPECT_TRUE(bulkData.is_valid(elem1));

    stk::mesh::EntityVector nodes(bulkData.begin_nodes(elem1), bulkData.begin_nodes(elem1)+bulkData.num_nodes(elem1));
    EXPECT_EQ(8u, nodes.size());

    for(stk::mesh::Entity node : nodes) {
        EXPECT_FALSE(bulkData.bucket(node).member(part));
    }

    stk::mesh::PartVector add_parts(1, &part);
    bulkData.batch_change_entity_parts(nodes, add_parts, {});

    for(stk::mesh::Entity node : nodes) {
        EXPECT_TRUE(bulkData.bucket(node).member(part));
    }

    stk::mesh::PartVector remove_parts(1, &part);
    bulkData.batch_change_entity_parts(nodes, {}, remove_parts);

    for(stk::mesh::Entity node : nodes) {
        EXPECT_FALSE(bulkData.bucket(node).member(part));
    }
}

TEST(UnitTestChangeParts, test_superset_and_subset_part_change)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int p_size = stk::parallel_machine_size( pm );

    if (p_size != 1) {
      return;
    }

    const int spatialDim = 3;
    stk::mesh::MetaData metaData(spatialDim);
    stk::mesh::BulkData bulkData(metaData, pm, stk::mesh::BulkData::NO_AUTO_AURA);

    stk::mesh::Part& supersetPart = metaData.declare_part_with_topology("parent", stk::topology::NODE);
    stk::mesh::Part& subsetPart1  = metaData.declare_part_with_topology("child 1", stk::topology::NODE);
    stk::mesh::Part& subsetPart2  = metaData.declare_part_with_topology("child 2", stk::topology::NODE);

    metaData.declare_part_subset(supersetPart, subsetPart1);
    metaData.declare_part_subset(supersetPart, subsetPart2);

    stk::mesh::FieldBase &field = metaData.declare_field< stk::mesh::Field<double> >(stk::topology::NODE_RANK, "sam");
    stk::mesh::put_field_on_mesh(field, supersetPart, (stk::mesh::FieldTraits<stk::mesh::Field<double>>::data_type*) nullptr);

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, bulkData);

    stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, 1u);
    EXPECT_TRUE(bulkData.is_valid(node1));

    stk::mesh::Entity node2 = bulkData.get_entity(stk::topology::NODE_RANK, 2u);
    EXPECT_TRUE(bulkData.is_valid(node2));

    EXPECT_FALSE(bulkData.bucket(node1).member(supersetPart));
    EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart1));
    EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart2));

    EXPECT_FALSE(bulkData.bucket(node2).member(supersetPart));
    EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
    EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart2));

    bulkData.modification_begin();
    stk::mesh::PartVector add_parts(1, &subsetPart1);
    bulkData.change_entity_parts(node1, add_parts, {});
    bulkData.modification_end();

    EXPECT_TRUE(bulkData.bucket(node1).member(supersetPart));
    EXPECT_TRUE(bulkData.bucket(node1).member(subsetPart1));
    EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart2));

    bulkData.modification_begin();
    add_parts[0] = &supersetPart;
    bulkData.change_entity_parts(node2, add_parts, {});
    bulkData.modification_end();

    EXPECT_TRUE(bulkData.bucket(node2).member(supersetPart));
    EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
    EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart2));

    stk::mesh::PartVector partVector{&subsetPart1, &subsetPart2};
    stk::mesh::Selector selector = supersetPart & !stk::mesh::selectUnion(partVector);
    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(selector, bulkData.buckets(stk::topology::NODE_RANK), nodes);

    EXPECT_EQ(1u, nodes.size());
    EXPECT_EQ(2u, bulkData.identifier(nodes[0]));

    bulkData.modification_begin();
    add_parts[0] = &subsetPart2;
    bulkData.change_entity_parts(node2, add_parts, {});
    bulkData.modification_end();

    EXPECT_TRUE(bulkData.bucket(node2).member(supersetPart));
    EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
    EXPECT_TRUE(bulkData.bucket(node2).member(subsetPart2));

    double* node1Data = (double*) stk::mesh::field_data(field, node1);
    double* node2Data = (double*) stk::mesh::field_data(field, node1);

    EXPECT_TRUE(node1Data != nullptr);
    EXPECT_TRUE(node2Data != nullptr);


    for(unsigned i=3; i<8u; ++i) {
        stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, i);
        EXPECT_TRUE(bulkData.is_valid(node));
        double* nodeData = (double*) stk::mesh::field_data(field, node);
        EXPECT_TRUE(nodeData == nullptr);
    }
}

}
