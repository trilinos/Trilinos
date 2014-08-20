
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/ConnectivityMap.hpp>  // for ConnectivityMap
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc

namespace {
TEST(stkMeshHowTo, useCustomConnectivityMap)
{
    const unsigned spatialDimension = 2;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());
    stk::mesh::Part &tri_part = metaData.declare_part_with_topology("tri_part", stk::topology::TRIANGLE_3_2D);
    stk::mesh::Part &edge_part = metaData.declare_part_with_topology("edge_part", stk::topology::LINE_2);

    metaData.commit();

    stk::mesh::ConnectivityMap custom_connectivity = stk::mesh::ConnectivityMap::none();
    //Now set which connectivities we want enabled:
    custom_connectivity(stk::topology::ELEM_RANK, stk::topology::NODE_RANK) = stk::mesh::ConnectivityMap::fixed();
    custom_connectivity(stk::topology::ELEM_RANK, stk::topology::EDGE_RANK) = stk::mesh::ConnectivityMap::fixed();
    custom_connectivity(stk::topology::EDGE_RANK, stk::topology::NODE_RANK) = stk::mesh::ConnectivityMap::fixed();
    //for most non-trivial cases, you must at least enable upward node -> elem connections.
    custom_connectivity(stk::topology::NODE_RANK, stk::topology::ELEM_RANK) = stk::mesh::ConnectivityMap::dynamic();

    //Now verify that node->edge, and edge->elem connections are disabled, but
    //elem->node and edge->node connections are allowed:
    EXPECT_FALSE(custom_connectivity.valid(stk::topology::NODE_RANK, stk::topology::EDGE_RANK));
    EXPECT_FALSE(custom_connectivity.valid(stk::topology::EDGE_RANK, stk::topology::ELEM_RANK));
    EXPECT_TRUE(custom_connectivity.valid(stk::topology::ELEM_RANK, stk::topology::NODE_RANK));
    EXPECT_TRUE(custom_connectivity.valid(stk::topology::EDGE_RANK, stk::topology::NODE_RANK));

    bool add_fmwk_data = false;
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD, add_fmwk_data, &custom_connectivity);
    mesh.modification_begin();

    //set up 1 element (3-node triangle) with elem->node and edge->node connections
    stk::mesh::EntityId elemId = 1;
    stk::mesh::EntityId elemNodeIds[] = {1, 2, 3};
    stk::mesh::EntityId elemEdgeIds[] = {6, 7, 8};
    stk::mesh::Entity elemNodes[3];
    stk::mesh::Entity elemEdges[3];
    stk::mesh::Entity elem = mesh.declare_entity(stk::topology::ELEM_RANK, elemId, tri_part);
    elemNodes[0] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[0]);
    elemNodes[1] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[1]);
    elemNodes[2] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[2]);

    elemEdges[0] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[0], edge_part);
    elemEdges[1] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[1], edge_part);
    elemEdges[2] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[2], edge_part);

    //downward element -> node connectivity
    mesh.declare_relation(elem, elemNodes[0], 0);
    mesh.declare_relation(elem, elemNodes[1], 1);
    mesh.declare_relation(elem, elemNodes[2], 2);

    //downward edge -> node connectivity
    mesh.declare_relation(elemEdges[0], elemNodes[0], 0); mesh.declare_relation(elemEdges[0], elemNodes[1], 1);
    mesh.declare_relation(elemEdges[1], elemNodes[1], 0); mesh.declare_relation(elemEdges[1], elemNodes[2], 1);
    mesh.declare_relation(elemEdges[2], elemNodes[2], 0); mesh.declare_relation(elemEdges[2], elemNodes[0], 1);
    mesh.modification_end();

    //now test upward connectivity which stk-mesh would have created automatically if we didn't disable it.
    unsigned expectedNumElemsPerEdge = 0;//this would be 1 if we didn't disable
    unsigned expectedNumEdgesPerNode = 0;//this would be 2 if we didn't disable
    EXPECT_EQ(expectedNumElemsPerEdge, mesh.num_elements(elemEdges[0]));
    EXPECT_EQ(expectedNumEdgesPerNode, mesh.num_edges(elemNodes[0]));
} }
