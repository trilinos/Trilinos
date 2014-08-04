#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <vector>

namespace {

//Mapping
TEST(stk_topology_how_to, map_topologies_to_ranks )
{
    stk::topology topology = stk::topology::INVALID_TOPOLOGY;
    EXPECT_EQ(stk::topology::INVALID_RANK, topology.rank());

    std::vector<stk::topology> node_rank_topologies;
    node_rank_topologies.push_back(stk::topology::NODE);

    std::vector<stk::topology> edge_rank_topologies;
    edge_rank_topologies.push_back(stk::topology::LINE_2);
    edge_rank_topologies.push_back(stk::topology::LINE_3);

    std::vector<stk::topology> face_rank_topologies;
    face_rank_topologies.push_back(stk::topology::TRI_3);
    face_rank_topologies.push_back(stk::topology::TRIANGLE_3);
    face_rank_topologies.push_back(stk::topology::TRI_4);
    face_rank_topologies.push_back(stk::topology::TRIANGLE_4);
    face_rank_topologies.push_back(stk::topology::TRI_6);
    face_rank_topologies.push_back(stk::topology::TRIANGLE_6);
    face_rank_topologies.push_back(stk::topology::QUAD_4);
    face_rank_topologies.push_back(stk::topology::QUADRILATERAL_4);
    face_rank_topologies.push_back(stk::topology::QUAD_8);
    face_rank_topologies.push_back(stk::topology::QUADRILATERAL_8);
    face_rank_topologies.push_back(stk::topology::QUAD_9);
    face_rank_topologies.push_back(stk::topology::QUADRILATERAL_9);

    std::vector<stk::topology> element_rank_topologies;
    element_rank_topologies.push_back(stk::topology::PARTICLE);
    element_rank_topologies.push_back(stk::topology::LINE_2_1D);
    element_rank_topologies.push_back(stk::topology::LINE_3_1D);
    element_rank_topologies.push_back(stk::topology::BEAM_2);
    element_rank_topologies.push_back(stk::topology::BEAM_3);
    element_rank_topologies.push_back(stk::topology::SHELL_LINE_2);
    element_rank_topologies.push_back(stk::topology::SHELL_LINE_3);

    element_rank_topologies.push_back(stk::topology::TRI_3_2D);
    element_rank_topologies.push_back(stk::topology::TRIANGLE_3_2D);
    element_rank_topologies.push_back(stk::topology::TRI_4_2D);
    element_rank_topologies.push_back(stk::topology::TRIANGLE_4_2D);
    element_rank_topologies.push_back(stk::topology::TRI_6_2D);
    element_rank_topologies.push_back(stk::topology::TRIANGLE_6_2D);
    element_rank_topologies.push_back(stk::topology::QUAD_4_2D);
    element_rank_topologies.push_back(stk::topology::QUADRILATERAL_4_2D);
    element_rank_topologies.push_back(stk::topology::QUAD_8_2D);
    element_rank_topologies.push_back(stk::topology::QUADRILATERAL_8_2D);
    element_rank_topologies.push_back(stk::topology::QUAD_9_2D);
    element_rank_topologies.push_back(stk::topology::QUADRILATERAL_9_2D);

    element_rank_topologies.push_back(stk::topology::SHELL_TRI_3);
    element_rank_topologies.push_back(stk::topology::SHELL_TRIANGLE_3);
    element_rank_topologies.push_back(stk::topology::SHELL_TRI_4);
    element_rank_topologies.push_back(stk::topology::SHELL_TRIANGLE_4);
    element_rank_topologies.push_back(stk::topology::SHELL_TRI_6);
    element_rank_topologies.push_back(stk::topology::SHELL_TRIANGLE_6);

    element_rank_topologies.push_back(stk::topology::SHELL_QUAD_4);
    element_rank_topologies.push_back(stk::topology::SHELL_QUADRILATERAL_4);
    element_rank_topologies.push_back(stk::topology::SHELL_QUAD_8);
    element_rank_topologies.push_back(stk::topology::SHELL_QUADRILATERAL_8);
    element_rank_topologies.push_back(stk::topology::SHELL_QUAD_9);
    element_rank_topologies.push_back(stk::topology::SHELL_QUADRILATERAL_9);

    element_rank_topologies.push_back(stk::topology::TET_4);
    element_rank_topologies.push_back(stk::topology::TETRAHEDRON_4);
    element_rank_topologies.push_back(stk::topology::TET_8);
    element_rank_topologies.push_back(stk::topology::TETRAHEDRON_8);
    element_rank_topologies.push_back(stk::topology::TET_10);
    element_rank_topologies.push_back(stk::topology::TETRAHEDRON_10);
    element_rank_topologies.push_back(stk::topology::TET_11);
    element_rank_topologies.push_back(stk::topology::TETRAHEDRON_11);

    element_rank_topologies.push_back(stk::topology::PYRAMID_5);
    element_rank_topologies.push_back(stk::topology::PYRAMID_13);
    element_rank_topologies.push_back(stk::topology::PYRAMID_14);
    element_rank_topologies.push_back(stk::topology::WEDGE_6);
    element_rank_topologies.push_back(stk::topology::WEDGE_15);
    element_rank_topologies.push_back(stk::topology::WEDGE_18);
    element_rank_topologies.push_back(stk::topology::QUADRILATERAL_9_2D);
    element_rank_topologies.push_back(stk::topology::QUADRILATERAL_9_2D);

    element_rank_topologies.push_back(stk::topology::HEX_8);
    element_rank_topologies.push_back(stk::topology::HEXAHEDRON_8);
    element_rank_topologies.push_back(stk::topology::HEX_20);
    element_rank_topologies.push_back(stk::topology::HEXAHEDRON_20);
    element_rank_topologies.push_back(stk::topology::HEX_27);
    element_rank_topologies.push_back(stk::topology::HEXAHEDRON_27);

    unsigned num_nodes_in_super_element = 10;
    element_rank_topologies.
      push_back(stk::create_superelement_topology(num_nodes_in_super_element));

//END-MAPPING

    // add a topology of invalid_rank
    unsigned zeroNodes = 0;
    element_rank_topologies.push_back(stk::create_superelement_topology(zeroNodes));

    ASSERT_EQ(1u, node_rank_topologies.size());
    ASSERT_EQ(2u, edge_rank_topologies.size());
    ASSERT_EQ(12u, face_rank_topologies.size());
    ASSERT_EQ(55u, element_rank_topologies.size());

    for (size_t i=0;i<node_rank_topologies.size();i++)
    {
        EXPECT_TRUE(node_rank_topologies[i].is_valid());
        EXPECT_EQ(stk::topology::NODE_RANK, node_rank_topologies[i].rank());
    }

    for (size_t i=0;i<edge_rank_topologies.size();i++)
    {
        EXPECT_TRUE(edge_rank_topologies[i].is_valid());
        EXPECT_EQ(stk::topology::EDGE_RANK, edge_rank_topologies[i].rank());
    }

    for (size_t i=0;i<face_rank_topologies.size();i++)
    {
        EXPECT_TRUE(face_rank_topologies[i].is_valid());
        EXPECT_EQ(stk::topology::FACE_RANK, face_rank_topologies[i].rank());
    }

    for (size_t i=0;i<element_rank_topologies.size()-1;i++)
    {
        EXPECT_TRUE(element_rank_topologies[i].is_valid());
        EXPECT_EQ(stk::topology::ELEMENT_RANK, element_rank_topologies[i].rank());
    }

    EXPECT_FALSE(element_rank_topologies.back().is_valid());
    EXPECT_EQ(stk::topology::INVALID_RANK, element_rank_topologies.back().rank());
}

}

