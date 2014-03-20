#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


namespace {

TEST(stk_topology_how_to, map_topologies_to_ranks )
{
    stk::topology topology = stk::topology::INVALID_TOPOLOGY;
    EXPECT_EQ(stk::topology::INVALID_RANK, topology.rank());

    std::vector<stk::topology> node_rank_topologies;
    node_rank_topologies.push_back(stk::topology::NODE);

    ASSERT_EQ(2u, node_rank_topologies.size());

    std::vector<stk::topology> edge_rank_topologies;
    edge_rank_topologies.push_back(stk::topology::LINE_2);
    edge_rank_topologies.push_back(stk::topology::LINE_3);

    ASSERT_EQ(2u, edge_rank_topologies.size());

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

    ASSERT_EQ(12u, face_rank_topologies.size());

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

    // zero noded super element
    unsigned num_nodes_in_super_element = 10;
    element_rank_topologies.push_back(stk::create_superelement_topology(num_nodes_in_super_element));

    // add a topology of invalid_rank
    unsigned zeroNodes = 0;
    element_rank_topologies.push_back(stk::create_superelement_topology(zeroNodes));

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

TEST(stk_topology_understanding, zero_dim_element)
{
    stk::topology sphere = stk::topology::PARTICLE;

    EXPECT_TRUE(sphere.is_valid());
    EXPECT_FALSE(sphere.has_homogeneous_faces());
    EXPECT_FALSE(sphere.is_shell());

    EXPECT_TRUE(sphere.rank() != stk::topology::NODE_RANK);
    EXPECT_TRUE(sphere.rank() != stk::topology::EDGE_RANK);
    EXPECT_TRUE(sphere.rank() != stk::topology::FACE_RANK);
    EXPECT_TRUE(sphere.rank() != stk::topology::CONSTRAINT_RANK);
    EXPECT_TRUE(sphere.rank() == stk::topology::ELEMENT_RANK);

    EXPECT_EQ(sphere.side_rank(), stk::topology::NODE_RANK);

    EXPECT_EQ(sphere.dimension(),1u);
    EXPECT_EQ(sphere.num_nodes(),1u);
    EXPECT_EQ(sphere.num_vertices(),1u);
    EXPECT_EQ(sphere.num_edges(),0u);
    EXPECT_EQ(sphere.num_faces(),0u);
    EXPECT_EQ(sphere.num_permutations(),1u);
    EXPECT_EQ(sphere.num_positive_permutations(),1u);

    EXPECT_FALSE(sphere.defined_on_spatial_dimension(0));

    EXPECT_TRUE(sphere.defined_on_spatial_dimension(1));
    EXPECT_TRUE(sphere.defined_on_spatial_dimension(2));
    EXPECT_TRUE(sphere.defined_on_spatial_dimension(3));

    EXPECT_EQ(sphere.base(),stk::topology::PARTICLE);
}

TEST(stk_topology_understanding, constraint_ranks)
{
    // that's it folks, there is no topology that provides a constraint rank
}

void checkForValidOffsets(const unsigned numNodes, const unsigned *offsets, const unsigned *goldValuesOffsets)
{
    for (unsigned i=0;i<numNodes;i++)
    {
        EXPECT_EQ(goldValuesOffsets[i], offsets[i]);
    }
}

void checkPermutedNodeIds(const unsigned numNodes, const unsigned *offsets, const unsigned *nodeIds, const unsigned *goldIds)
{
    for (unsigned i=0;i<numNodes;i++)
    {
        EXPECT_EQ(goldIds[offsets[i]], nodeIds[i]);
    }
}

void checkNodeOrderingAndOffsetsForFaces(const stk::topology &element, const unsigned *elementNodes, const unsigned *goldValuesOffsets)
{
    ASSERT_TRUE(element.has_homogeneous_faces());
    unsigned numNodesPerFace = element.face_topology(0).num_nodes();

    unsigned *offsets = new unsigned[numNodesPerFace];
    unsigned *nodeIds = new unsigned[numNodesPerFace];

    for (unsigned index=0; index<element.num_faces();index++)
    {
        element.face_node_ordinals(index, offsets);
        element.face_nodes(elementNodes, index, nodeIds);

        checkForValidOffsets(numNodesPerFace, offsets, &goldValuesOffsets[numNodesPerFace*index]);
        checkPermutedNodeIds(numNodesPerFace, offsets, nodeIds, elementNodes);
    }

    delete [] nodeIds; nodeIds = 0;
    delete [] offsets; offsets = 0;
}

void checkNodeOrderingAndOffsetsForEdges(const stk::topology &element, const unsigned *elementNodes, const unsigned *goldValuesOffsets)
{
    unsigned numNodesPerEdge = element.edge_topology().num_nodes();
    unsigned *offsets = new unsigned[numNodesPerEdge];
    unsigned *nodeIds = new unsigned[numNodesPerEdge];

    for (unsigned index=0; index<element.num_edges();index++)
    {
        element.edge_node_ordinals(index, offsets);
        element.edge_nodes(elementNodes, index, nodeIds);

        checkForValidOffsets(numNodesPerEdge, offsets, &goldValuesOffsets[numNodesPerEdge*index]);
        checkPermutedNodeIds(numNodesPerEdge, offsets, nodeIds, elementNodes);
    }

    delete [] nodeIds; nodeIds = 0;
    delete [] offsets; offsets = 0;
}

void checkNodeOrderingAndOffsetsForPermutations(const stk::topology &element, const unsigned *elementNodes, const unsigned *goldValuesOffsets)
{
    unsigned numNodes = element.num_nodes();
    unsigned *offsets = new unsigned[numNodes];
    unsigned *nodeIds = new unsigned[numNodes];

    for (unsigned index=0; index<element.num_permutations();index++)
    {
        element.permutation_node_ordinals(index, offsets);
        element.permutation_nodes(elementNodes, index, nodeIds);

        checkForValidOffsets(numNodes, offsets, &goldValuesOffsets[numNodes*index]);
        checkPermutedNodeIds(numNodes, offsets, nodeIds, elementNodes);
    }

    delete [] nodeIds; nodeIds = 0;
    delete [] offsets; offsets = 0;
}

TEST(stk_topology_understanding, one_dim_higher_order_element)
{
    stk::topology secondOrderBeam = stk::topology::BEAM_3;

    EXPECT_TRUE(secondOrderBeam.is_valid());
    EXPECT_FALSE(secondOrderBeam.has_homogeneous_faces());
    EXPECT_FALSE(secondOrderBeam.is_shell());

    EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::NODE_RANK);
    EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::EDGE_RANK);
    EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::FACE_RANK);
    EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::CONSTRAINT_RANK);
    EXPECT_TRUE(secondOrderBeam.rank() == stk::topology::ELEMENT_RANK);

    EXPECT_TRUE(secondOrderBeam.side_rank() == stk::topology::EDGE_RANK);

    EXPECT_EQ(2u, secondOrderBeam.dimension());
    EXPECT_EQ(3u, secondOrderBeam.num_nodes());
    EXPECT_EQ(2u, secondOrderBeam.num_vertices());
    EXPECT_EQ(1u, secondOrderBeam.num_edges());

    EXPECT_EQ(0u, secondOrderBeam.num_faces());
    EXPECT_EQ(1u, secondOrderBeam.num_positive_permutations());
    EXPECT_EQ(2u, secondOrderBeam.num_permutations());

    EXPECT_FALSE(secondOrderBeam.defined_on_spatial_dimension(0));
    EXPECT_FALSE(secondOrderBeam.defined_on_spatial_dimension(1));

    EXPECT_TRUE(secondOrderBeam.defined_on_spatial_dimension(2));
    EXPECT_TRUE(secondOrderBeam.defined_on_spatial_dimension(3));

    EXPECT_TRUE(secondOrderBeam.base() == stk::topology::BEAM_2);

    unsigned beamNodes[3] = { 10, 20, 14 }; // 10 *-------*-------* 20
                                            //            14
    {
        unsigned goldValuesOffsets[3] = { 0, 1, 2 };
        checkNodeOrderingAndOffsetsForEdges(secondOrderBeam, beamNodes, goldValuesOffsets);
    }

    {
        unsigned goldValuesOffsets[6] = {
                0, 1, 2,
                1, 0, 2
        };

        checkNodeOrderingAndOffsetsForPermutations(secondOrderBeam, beamNodes, goldValuesOffsets);
    }
}

TEST(stk_topology_understanding, two_dim_higher_order_element)
{
    stk::topology secondOrderTriShell = stk::topology::SHELL_TRIANGLE_6;
    EXPECT_TRUE(secondOrderTriShell == stk::topology::SHELL_TRI_6);

    EXPECT_TRUE(secondOrderTriShell.is_valid());
    EXPECT_TRUE(secondOrderTriShell.has_homogeneous_faces());
    EXPECT_TRUE(secondOrderTriShell.is_shell());

    EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::NODE_RANK);
    EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::EDGE_RANK);
    EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::FACE_RANK);
    EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::CONSTRAINT_RANK);
    EXPECT_TRUE(secondOrderTriShell.rank() == stk::topology::ELEMENT_RANK);

    EXPECT_TRUE(secondOrderTriShell.side_rank() == stk::topology::FACE_RANK);

    EXPECT_EQ(3u, secondOrderTriShell.dimension());
    EXPECT_EQ(6u, secondOrderTriShell.num_nodes());
    EXPECT_EQ(3u, secondOrderTriShell.num_vertices());
    EXPECT_EQ(3u, secondOrderTriShell.num_edges());
    EXPECT_EQ(2u, secondOrderTriShell.num_faces());

    // permutations are the number of ways the number of vertices can be permuted
    EXPECT_EQ(6u, secondOrderTriShell.num_permutations());
    // positive permutations are ones that the normal is maintained
    EXPECT_EQ(3u, secondOrderTriShell.num_positive_permutations());

    EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(0));
    EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(1));
    EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(2));

    EXPECT_TRUE(secondOrderTriShell.defined_on_spatial_dimension(3));

    EXPECT_TRUE(secondOrderTriShell.base() == stk::topology::SHELL_TRI_3);
    EXPECT_TRUE(secondOrderTriShell.base() == stk::topology::SHELL_TRIANGLE_3);

    unsigned shellNodes[6] = { 10, 11, 12, 100, 101, 102 }; // first 3 are vertex nodes (picture?)

    {
        unsigned goldValuesEdgeOffsets[9] = {
                0, 1, 3,
                1, 2, 4,
                2, 0, 5
        };

        checkNodeOrderingAndOffsetsForEdges(secondOrderTriShell, shellNodes, goldValuesEdgeOffsets);
    }

    {
        unsigned goldValuesFaceNodeOffsets[12] = {
                0, 1, 2, 3, 4, 5,
                0, 2, 1, 5, 4, 3
        };

        checkNodeOrderingAndOffsetsForFaces(secondOrderTriShell, shellNodes, goldValuesFaceNodeOffsets);
    }

    {
        unsigned goldValueOffsetsPerm[36] = {
                0, 1, 2, 3, 4, 5,
                2, 0, 1, 5, 3, 4,
                1, 2, 0, 4, 5, 3,
                0, 2, 1, 5, 4, 3,
                2, 1, 0, 4, 3, 5,
                1, 0, 2, 3, 5, 4
        };

        checkNodeOrderingAndOffsetsForPermutations(secondOrderTriShell, shellNodes, goldValueOffsetsPerm);
    }
}



TEST(stk_topology_understanding, three_dim_linear_element)
{
    stk::topology hex8 = stk::topology::HEX_8;
    EXPECT_TRUE(hex8 == stk::topology::HEXAHEDRON_8);

    EXPECT_TRUE(hex8.is_valid());
    EXPECT_TRUE(hex8.has_homogeneous_faces());
    EXPECT_FALSE(hex8.is_shell());

    EXPECT_TRUE(hex8.rank() != stk::topology::NODE_RANK);
    EXPECT_TRUE(hex8.rank() != stk::topology::EDGE_RANK);
    EXPECT_TRUE(hex8.rank() != stk::topology::FACE_RANK);
    EXPECT_TRUE(hex8.rank() != stk::topology::CONSTRAINT_RANK);
    EXPECT_TRUE(hex8.rank() == stk::topology::ELEMENT_RANK);

    EXPECT_TRUE(hex8.side_rank() == stk::topology::FACE_RANK);

    EXPECT_EQ(3u, hex8.dimension());
    EXPECT_EQ(8u, hex8.num_nodes());
    EXPECT_EQ(8u, hex8.num_vertices());
    EXPECT_EQ(12u, hex8.num_edges());
    EXPECT_EQ(6u, hex8.num_faces());

    // permutations are the number of ways the number of vertices can be permuted
    EXPECT_EQ(1u, hex8.num_permutations());
    // positive permutations are ones that the normal is maintained
    EXPECT_EQ(1u, hex8.num_positive_permutations());

    EXPECT_FALSE(hex8.defined_on_spatial_dimension(0));
    EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
    EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));

    EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

    EXPECT_TRUE(hex8.base() == stk::topology::HEX_8);

    unsigned hex8Nodes[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

    {
        stk::topology goldEdgeTopology = stk::topology::LINE_2;
        EXPECT_EQ(goldEdgeTopology, hex8.edge_topology());

        unsigned goldNumNodesPerEdge = 2;
        ASSERT_EQ(goldNumNodesPerEdge, hex8.edge_topology().num_nodes());
        unsigned goldValuesEdgeOffsets[24] = {
                0, 1,
                1, 2,
                2, 3,
                3, 0,
                4, 5,
                5, 6,
                6, 7,
                7, 4,
                0, 4,
                1, 5,
                2, 6,
                3, 7
        };

        checkNodeOrderingAndOffsetsForEdges(hex8, hex8Nodes, goldValuesEdgeOffsets);
    }

    {
        stk::topology goldFaceTopology = stk::topology::QUAD_4;
        unsigned goldNumNodesPerFace = 4;
        for (unsigned faceIndex=0;faceIndex<hex8.num_faces();faceIndex++)
        {
            EXPECT_EQ(goldFaceTopology, hex8.face_topology(faceIndex));
            ASSERT_EQ(goldNumNodesPerFace, hex8.face_topology(faceIndex).num_nodes());
        }

        unsigned goldValuesFaceOffsets[24] = {
                0, 1, 5, 4,
                1, 2, 6, 5,
                2, 3, 7, 6,
                0, 4, 7, 3,
                0, 3, 2, 1,
                4, 5, 6, 7
        };

        checkNodeOrderingAndOffsetsForFaces(hex8, hex8Nodes, goldValuesFaceOffsets);
    }

    {
        unsigned goldValueOffsetsPerm[8] = {
               0, 1, 2, 3, 4, 5, 6, 7
        };

        checkNodeOrderingAndOffsetsForPermutations(hex8, hex8Nodes, goldValueOffsetsPerm);
    }
}


TEST(stk_topology_understanding, equivalent_elements)
{
    std::pair<bool, unsigned> areElementsEquivalent;

    {
        unsigned hex1[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
        unsigned hex2[8] = { 0, 1, 2, 3, 4, 7, 6, 5 };
        unsigned hex3[8] = { 4, 5, 6, 7, 0, 1, 2, 3 };

        stk::topology hex8 = stk::topology::HEX_8;

        areElementsEquivalent = hex8.equivalent(hex1, hex2);
        EXPECT_FALSE(areElementsEquivalent.first);

        areElementsEquivalent = hex8.equivalent(hex1, hex3);
        // should be true, but is not implemented
        EXPECT_FALSE(areElementsEquivalent.first);
    }

    {
        unsigned triangle_1[3] = {0, 1, 2};
        unsigned triangle_2[3] = {0, 2, 1};

        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;

        areElementsEquivalent = triangular_shell.equivalent(triangle_1, triangle_2);

        EXPECT_TRUE(areElementsEquivalent.first);

        unsigned permutation_index = areElementsEquivalent.second;
        unsigned goldValue = 3;
        EXPECT_EQ(goldValue, permutation_index); // From previous unit test, this is the 4th permutation
    }

}

void verifyPermutationsForTriangle(stk::topology triangular_shell, unsigned* triangle_1_node_ids, unsigned* gold_triangle_1_permutations)
{
    ASSERT_TRUE(stk::topology::SHELL_TRIANGLE_3 == triangular_shell);
    unsigned triangle_1_permutation[3];
    for (unsigned i=0;i<triangular_shell.num_permutations();i++)
    {
        triangular_shell.permutation_nodes(triangle_1_node_ids, i, triangle_1_permutation);
        EXPECT_TRUE(gold_triangle_1_permutations[3*i+0] == triangle_1_permutation[0] &&
                    gold_triangle_1_permutations[3*i+1] == triangle_1_permutation[1] &&
                    gold_triangle_1_permutations[3*i+2] == triangle_1_permutation[2]);
    }
}

TEST(stk_topology_understanding, lexicographical_smallest_permutation)
{
    {
        unsigned triangle_node_ids[3] = {10, 8, 12};

        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;

        unsigned gold_triangle_permutations[18]= {
                10, 8, 12,
                12, 10, 8,
                8, 12, 10, // lexicographical smallest permutation by node ids if considering only positive permutations
                10, 12, 8,
                12, 8, 10,
                8, 10, 12  // lexicographical smallest permutation by node ids if considering all permutations
        };

        verifyPermutationsForTriangle(triangular_shell, triangle_node_ids, gold_triangle_permutations);

        bool usePositivePermutationsOnly = false;
        unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation(triangle_node_ids, usePositivePermutationsOnly);
        unsigned gold_lexicographical_smallest_permutation_index = 5;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);

        usePositivePermutationsOnly = true;
        permutation_index = triangular_shell.lexicographical_smallest_permutation(triangle_node_ids, usePositivePermutationsOnly);
        gold_lexicographical_smallest_permutation_index = 2;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
}

TEST(stk_topology_understanding, sub_topology)
{
    stk::topology hex20 = stk::topology::HEX_20;
    unsigned hex20Nodes[20] = {
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9, 10, 11,
            12, 13, 14, 15,
            16, 17, 18, 19
    };

    unsigned numFaces = hex20.num_sub_topology(stk::topology::FACE_RANK);
    EXPECT_EQ(6u, numFaces);

    unsigned faceIndex=2;
    stk::topology top = hex20.sub_topology(stk::topology::FACE_RANK, faceIndex);
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, top);

    unsigned nodeIdsFace[8];
    hex20.sub_topology_nodes(hex20Nodes, stk::topology::FACE_RANK, faceIndex, nodeIdsFace);

    unsigned goldIdsFace[8] = { 2, 3, 7, 6, 10, 15, 18, 14 };
    for (unsigned i=0;i<hex20.face_topology(faceIndex).num_nodes();i++)
    {
        EXPECT_EQ(goldIdsFace[i], nodeIdsFace[i]);
    }
}

TEST(stk_topology_understanding, sides)
{
    stk::topology hex20 = stk::topology::HEX_20;
    EXPECT_EQ(6u, hex20.num_sides());

    stk::topology quad8 = stk::topology::SHELL_QUADRILATERAL_8;
    EXPECT_EQ(2u, quad8.num_sides());

    stk::topology wedge = stk::topology::WEDGE_15;
    EXPECT_EQ(5u, wedge.num_sides());
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(0));
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(1));
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(2));
    EXPECT_EQ(stk::topology::TRIANGLE_6, wedge.side_topology(3));
    EXPECT_EQ(stk::topology::TRIANGLE_6, wedge.side_topology(4));

}

TEST(stk_topology_understanding, superelements)
{
    unsigned eightNodes=8;
    stk::topology validSuperElement = stk::create_superelement_topology(eightNodes);
    EXPECT_TRUE(validSuperElement.is_superelement());
    EXPECT_TRUE(stk::topology::ELEMENT_RANK == validSuperElement.rank());
    EXPECT_EQ(eightNodes, validSuperElement.num_nodes());
    EXPECT_EQ(0u, validSuperElement.num_edges());
    EXPECT_EQ(0u, validSuperElement.num_faces());
    EXPECT_EQ(0u, validSuperElement.num_permutations());
    EXPECT_EQ(0u, validSuperElement.num_sides());
    EXPECT_EQ(0u, validSuperElement.dimension());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.face_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.edge_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.base());
    EXPECT_FALSE(validSuperElement.has_homogeneous_faces());
    EXPECT_FALSE(validSuperElement.is_shell());

    unsigned zeroNodes=0;
    stk::topology invalidSuperElement = stk::create_superelement_topology(zeroNodes);
    EXPECT_FALSE(invalidSuperElement.is_superelement());
    EXPECT_TRUE(stk::topology::INVALID_RANK == invalidSuperElement.rank());
    EXPECT_EQ(zeroNodes, invalidSuperElement.num_nodes());
    EXPECT_EQ(0u, invalidSuperElement.num_edges());
    EXPECT_EQ(0u, invalidSuperElement.num_faces());
    EXPECT_EQ(0u, invalidSuperElement.num_permutations());
    EXPECT_EQ(0u, invalidSuperElement.num_sides());
    EXPECT_EQ(0u, invalidSuperElement.dimension());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.face_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.edge_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.base());
    EXPECT_FALSE(invalidSuperElement.has_homogeneous_faces());
    EXPECT_FALSE(invalidSuperElement.is_shell());
}

}

// stk_mesh should get away from generic entity, the entity keyword altogether, and use words
// like node, face, edge, and element.
// then figure out how this impacts stk::topology... (Martin: no need for node_rank and node_topology)
// this is a bleed through of internal workings to the user which doesn't make sense.

// remove constraint rank from topology

// Permutations for any continuum element are NOT DONE. Does anyone need them?

// consider adding a picture for each test for documentation?
