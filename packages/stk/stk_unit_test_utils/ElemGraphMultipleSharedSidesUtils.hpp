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
#ifndef _ElemGraphMultipleSharedSidesUtils_hpp_
#define _ElemGraphMultipleSharedSidesUtils_hpp_

#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"

namespace {

class TwoElemMultipleSharedSideTester : public ::testing::Test
{
public:
    TwoElemMultipleSharedSideTester(size_t spatialDim, const std::vector<stk::mesh::EntityId> nodeIDs[2], const std::vector<stk::mesh::EntityId> &sharedNodeIds, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    : meta(spatialDim),
      skinPart(meta.declare_part_with_topology("skin", get_side_topology())),
      activePart(meta.declare_part("active")),
      bulkData(meta, MPI_COMM_WORLD, auraOption)
    {
        if (bulkData.parallel_size() <= 2)
        {
            make_mesh_2_elems_connected_through_multiple_sides(nodeIDs, sharedNodeIds);
        }
    }

    stk::topology get_side_topology()
    {
        return get_element_topology().side_topology();
    }

    stk::topology get_element_topology()
    {
        if(meta.spatial_dimension() == 3)
            return stk::topology::HEX_8;
        else
            return stk::topology::QUAD_4_2D;
    }

    void make_mesh_2_elems_connected_through_multiple_sides(const std::vector<stk::mesh::EntityId> nodeIDs[2], const std::vector<stk::mesh::EntityId> &sharedNodeIds)
    {
        bulkData.modification_begin();
        create_elements(bulkData, nodeIDs, activePart);
        add_shared_nodes(bulkData, sharedNodeIds);
        bulkData.modification_end();
    }

    void create_element(stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityId> nodeIds[2], stk::mesh::Part& activePart, stk::mesh::EntityId id)
    {
        stk::mesh::Part& elemPart = bulkData.mesh_meta_data().get_topology_root_part(get_element_topology());
        stk::mesh::Entity elem = stk::mesh::declare_element(bulkData, elemPart, id, nodeIds[id-1]);
        bulkData.change_entity_parts(elem, {&activePart});
    }

    bool i_should_create_elem_1(stk::mesh::BulkData& bulkData)
    {
        return bulkData.parallel_size() == 1 || bulkData.parallel_rank() == 0;
    }

    bool i_should_create_elem_2(stk::mesh::BulkData& bulkData)
    {
        return bulkData.parallel_size() == 1 || bulkData.parallel_rank() == 1;
    }

    void create_elements(stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityId> nodeIds[2], stk::mesh::Part& activePart)
    {
        if(i_should_create_elem_1(bulkData))
            create_element(bulkData, nodeIds, activePart, 1);
        if(i_should_create_elem_2(bulkData))
            create_element(bulkData, nodeIds, activePart, 2);
    }

    void add_shared_nodes(stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityId>& sharedNodeIds)
    {
        int otherProc = 1 - bulkData.parallel_rank();
        for(stk::mesh::EntityId nodeID : sharedNodeIds) {
            stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, nodeID);
            bulkData.add_node_sharing(node, otherProc);
        }
    }

protected:
    stk::mesh::MetaData meta;
    stk::mesh::Part& skinPart;
    stk::mesh::Part& activePart;
    stk::mesh::BulkData bulkData;
};

static const std::vector<stk::mesh::EntityId> twoElemTwoSharedSideNodeIDs[2] = {
         {1, 2, 3, 4, 5, 6, 7, 8},
         {3, 2, 9, 4, 7, 6, 10, 8} };
static const std::vector<stk::mesh::EntityId> twoElemTwoSharedSideSharedNodeIds = {2, 3, 4, 6, 7, 8};

class TwoElemTwoSharedSideTester : public TwoElemMultipleSharedSideTester
{
public:
    TwoElemTwoSharedSideTester() :
        TwoElemMultipleSharedSideTester(3, twoElemTwoSharedSideNodeIDs, twoElemTwoSharedSideSharedNodeIds, stk::mesh::BulkData::AUTO_AURA)
    {}
};



void remove_part_if_owned(stk::mesh::BulkData& bulkData, stk::mesh::Entity entity, stk::mesh::Part& part)
{
    if (bulkData.is_valid(entity) && bulkData.bucket(entity).owned())
    {
        bulkData.change_entity_parts(entity, {}, {&part});
    }
}

void remove_element_from_part(stk::mesh::BulkData& bulkData, stk::mesh::EntityId elemId, stk::mesh::Part& activePart)
{
    stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
    bulkData.modification_begin();
    remove_part_if_owned(bulkData, elem, activePart);
    bulkData.modification_end();
}

void test_skinned_mesh(stk::mesh::BulkData& bulkData, size_t expectedNumFacesPerElement)
{
    ElemGraphTestUtils::test_num_faces_per_element(bulkData, {expectedNumFacesPerElement, expectedNumFacesPerElement});
    std::vector<size_t> global_mesh_counts;
    stk::mesh::comm_mesh_counts(bulkData, global_mesh_counts);
    EXPECT_EQ(expectedNumFacesPerElement*global_mesh_counts[stk::topology::ELEMENT_RANK], global_mesh_counts[bulkData.mesh_meta_data().side_rank()]);
}

void test_total_sides_and_sides_per_element(stk::mesh::BulkData& bulkData, size_t expectedTotalFaces, const std::vector<size_t>& expectedFacesPerElement)
{
    ElemGraphTestUtils::test_num_faces_per_element(bulkData, expectedFacesPerElement);
    std::vector<size_t> global_mesh_counts;
    stk::mesh::comm_mesh_counts(bulkData, global_mesh_counts);
    EXPECT_EQ(expectedTotalFaces, global_mesh_counts[bulkData.mesh_meta_data().side_rank()]);
}


stk::mesh::EntityVector get_killed_elements(stk::mesh::BulkData& bulkData)
{
    stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
    if (bulkData.is_valid(elem2) && bulkData.bucket(elem2).owned())
    {
        return {elem2};
    }
    return stk::mesh::EntityVector();
}

void test_element_death_with_multiple_shared_sides(stk::mesh::BulkData& bulkData, stk::mesh::Part& activePart, stk::mesh::Part& skinPart)
{
    stk::mesh::ElemElemGraph elem_elem_graph(bulkData, activePart);
    remove_element_from_part(bulkData, 2, activePart);
    process_killed_elements(bulkData, elem_elem_graph, get_killed_elements(bulkData), activePart, {&activePart, &skinPart});
    test_total_sides_and_sides_per_element(bulkData, 2u, {2u, 2u});
}


stk::mesh::EntityId get_connected_elem_id(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph& elem_elem_graph, stk::mesh::Entity localElem, size_t i)
{
    if(bulkData.parallel_size() > 1)
    {
        return elem_elem_graph.get_connected_remote_id_and_via_side(localElem, i).id;
    }
    else
    {
        return bulkData.identifier(elem_elem_graph.get_connected_element_and_via_side(localElem, i).element);
    }
}

bool is_connected_elem_local(int numProcs, stk::mesh::ElemElemGraph& elem_elem_graph, stk::mesh::Entity localElem, size_t i)
{
    if(numProcs > 1)
    {
        return !elem_elem_graph.is_connected_elem_locally_owned(localElem, i);
    }
    else
    {
        return elem_elem_graph.is_connected_elem_locally_owned(localElem, i);
    }
}

void test_local_or_remote_connections(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph& elem_elem_graph, stk::mesh::EntityId expectedRemoteId, stk::mesh::Entity localElem, size_t i)
{
    EXPECT_TRUE(is_connected_elem_local(bulkData.parallel_size(), elem_elem_graph, localElem, i));
    EXPECT_EQ(expectedRemoteId, get_connected_elem_id(bulkData, elem_elem_graph, localElem, i));
}

void test_elements_connected_n_times(stk::mesh::BulkData &bulkData, stk::mesh::ElemElemGraph& elem_elem_graph, stk::mesh::Entity localElem, stk::mesh::EntityId remoteId, size_t numTimesConnected)
{
    ASSERT_EQ(numTimesConnected, elem_elem_graph.get_num_connected_elems(localElem));
    for(size_t i=0; i<numTimesConnected; ++i)
    {
        test_local_or_remote_connections(bulkData, elem_elem_graph, remoteId, localElem, i);
    }
}

void test_elems_kissing_n_times(stk::mesh::BulkData& bulkData, stk::mesh::Part& activePart, size_t numKisses)
{
    stk::mesh::ElemElemGraph elem_elem_graph(bulkData, activePart);
    stk::mesh::EntityId id = 1 + bulkData.parallel_rank();
    stk::mesh::EntityId remoteId = 2u-bulkData.parallel_rank();
    test_elements_connected_n_times(bulkData, elem_elem_graph, bulkData.get_entity(stk::topology::ELEM_RANK, id), remoteId, numKisses);
}



static const std::vector<stk::mesh::EntityId> twoElemThreeSharedSideNodeIDs[2] = {
         {1, 2, 3, 4, 5, 6, 7, 8},
         {3, 2, 1, 4, 7, 6, 10, 8} };
static const std::vector<stk::mesh::EntityId> twoElemThreeSharedSideSharedNodeIds = {1, 2, 3, 4, 6, 7, 8};

class TwoElemThreeSharedSideTester : public TwoElemMultipleSharedSideTester
{
public:
    TwoElemThreeSharedSideTester() :
        TwoElemMultipleSharedSideTester(3, twoElemThreeSharedSideNodeIDs, twoElemThreeSharedSideSharedNodeIds, stk::mesh::BulkData::AUTO_AURA)
    {}
};

class TwoElemThreeSharedSideNoAuraTester : public TwoElemMultipleSharedSideTester
{
public:
    TwoElemThreeSharedSideNoAuraTester() :
        TwoElemMultipleSharedSideTester(3, twoElemThreeSharedSideNodeIDs, twoElemThreeSharedSideSharedNodeIds, stk::mesh::BulkData::NO_AUTO_AURA)
    {}
};

static std::vector<stk::mesh::EntityId> twoQuadTwoSharedSideNodeIDs[2] = {
         {1, 2, 4, 3},
         {2, 5, 3, 4} };
static std::vector<stk::mesh::EntityId> twoQuadTwoSharedSideSharedNodeIds = {2, 3, 4};


class TwoElem2dTwoSharedSideTester : public TwoElemMultipleSharedSideTester
{
public:
    TwoElem2dTwoSharedSideTester() :
        TwoElemMultipleSharedSideTester(2, twoQuadTwoSharedSideNodeIDs, twoQuadTwoSharedSideSharedNodeIds, stk::mesh::BulkData::NO_AUTO_AURA)
    {}
};


}

#endif

