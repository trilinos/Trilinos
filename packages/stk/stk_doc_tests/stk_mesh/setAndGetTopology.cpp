#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityId
namespace stk { namespace mesh { class Part; } }

namespace {
TEST(stkMeshHowTo, setAndGetTopology)
{
    const unsigned spatialDimension = 3;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());

    //There are two methods of attaching topology to a part:

    //method 1: declare a part with a specified topology
    stk::mesh::Part &tetPart = metaData.declare_part_with_topology("tetElementPart", stk::topology::TET_4);

    stk::mesh::Part &hexPart = metaData.declare_part("hexElementPart");
    //method 2: set a topology on an existing part:
    stk::mesh::set_topology(hexPart, stk::topology::HEX_8);

    metaData.commit();
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD);
    mesh.modification_begin();
    stk::mesh::EntityId elem1Id = 1, elem2Id = 2;
    stk::mesh::Entity elem1=mesh.declare_entity(stk::topology::ELEMENT_RANK, elem1Id, tetPart);
    stk::mesh::Entity elem2=mesh.declare_entity(stk::topology::ELEMENT_RANK, elem2Id, hexPart);

    //No common nodes between elements - so topologically disconnected elem1 and elem2
    for(unsigned node_ord = 0 ; node_ord < 4; ++node_ord)
    {
      stk::mesh::Entity new_node = mesh.declare_entity(stk::topology::NODE_RANK, node_ord+100*elem1Id);
      mesh.declare_relation( elem1 , new_node , node_ord);
    }

    for(unsigned node_ord = 0 ; node_ord < 8; ++node_ord)
    {
      stk::mesh::Entity new_node2 = mesh.declare_entity(stk::topology::NODE_RANK, node_ord+100*elem2Id);
      mesh.declare_relation( elem2 , new_node2 , node_ord);
    }

    mesh.modification_end();
    //Note: this test is only about setting/getting topology, so we didn't bother to
    //create nodes to connect to the above elements, etc.

    stk::topology elem1_topology = mesh.bucket(elem1).topology();
    stk::topology elem2_topology = mesh.bucket(elem2).topology();

    EXPECT_EQ(stk::topology::TET_4, elem1_topology);
    EXPECT_EQ(stk::topology::HEX_8, elem2_topology);
}
}
