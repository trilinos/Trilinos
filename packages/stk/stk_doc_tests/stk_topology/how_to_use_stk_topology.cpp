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

#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <vector>

#ifndef __IBMCPP__

namespace {

void verifyPermutationsForTriangle(unsigned* triangle_1_node_ids, unsigned* gold_triangle_1_permutations)
{
    stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;
    unsigned triangle_1_permutation[3];
    for (unsigned i=0; i<triangular_shell.num_permutations(); i++) {
        triangular_shell.permutation_nodes(triangle_1_node_ids, i, triangle_1_permutation);

        EXPECT_TRUE(gold_triangle_1_permutations[3*i+0] == triangle_1_permutation[0] &&
                    gold_triangle_1_permutations[3*i+1] == triangle_1_permutation[1] &&
                    gold_triangle_1_permutations[3*i+2] == triangle_1_permutation[2]);
    }
}

//BEGIN_lexicographical
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

        verifyPermutationsForTriangle(triangle_node_ids, gold_triangle_permutations);

        bool usePositivePermutationsOnly = false;
        unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation((unsigned*)triangle_node_ids, usePositivePermutationsOnly);
        unsigned gold_lexicographical_smallest_permutation_index = 5;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);

        usePositivePermutationsOnly = true;
        permutation_index = triangular_shell.lexicographical_smallest_permutation((unsigned*)triangle_node_ids, usePositivePermutationsOnly);
        gold_lexicographical_smallest_permutation_index = 2;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
}
//END_lexicographical

//BEGIN_Preserves_polarity_lexicographical
TEST(stk_topology_understanding, lexicographical_smallest_permutation_preserve_polarity)
{
    {
        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;
        unsigned shell_node_ids[3] = {10, 8, 12};
        {
            unsigned triangle_node_ids[3] = {12, 10, 8};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
            unsigned expected_positive_permutation = 2;

            EXPECT_EQ(expected_positive_permutation, permutation_index);
            EXPECT_LT(expected_positive_permutation, triangular_shell.num_positive_permutations());
        }
        {
            unsigned triangle_node_ids[3] = {12, 8, 10};

            unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)triangle_node_ids, (unsigned*)shell_node_ids);
            unsigned expected_negative_permutation = 5;

            EXPECT_EQ(expected_negative_permutation, permutation_index);
            EXPECT_GE(expected_negative_permutation, triangular_shell.num_positive_permutations());
        }
    }
}

TEST(stk_topology_understanding, quad_lexicographical_smallest_permutation_preserve_polarity)
{
    {
        stk::topology quad_shell = stk::topology::SHELL_QUAD_4;
        unsigned shell_node_ids[4] = {1, 2, 3, 4};
        {
            unsigned quad_node_ids[4] = {1, 2, 3, 4};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
            unsigned expected_positive_permutation = 0;

            EXPECT_EQ(expected_positive_permutation, permutation_index);
            EXPECT_LT(expected_positive_permutation, quad_shell.num_positive_permutations());
        }

        {
            unsigned quad_node_ids[4] = {1, 4, 3, 2};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
            unsigned expected_negative_permutation = 4;

            EXPECT_EQ(expected_negative_permutation, permutation_index);
            EXPECT_GE(expected_negative_permutation, quad_shell.num_positive_permutations());
        }

        {
            unsigned quad_node_ids[4] = {4, 2, 3, 1};

            unsigned permutation_index = quad_shell.lexicographical_smallest_permutation_preserve_polarity((unsigned*)quad_node_ids, (unsigned*)shell_node_ids);
            unsigned expected_invalid_permutation = 8;

            EXPECT_EQ(expected_invalid_permutation, permutation_index);
            EXPECT_EQ(expected_invalid_permutation, quad_shell.num_permutations());
        }
    }
}
//END_Preserves_polarity_lexicographical

//SubTopology
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
    hex20.sub_topology_nodes((unsigned*)hex20Nodes, stk::topology::FACE_RANK, faceIndex, (unsigned*)nodeIdsFace);

    unsigned goldIdsFace[8] = { 2, 3, 7, 6, 10, 15, 18, 14 };
    for (unsigned i=0;i<hex20.face_topology(faceIndex).num_nodes();i++)
    {
        EXPECT_EQ(goldIdsFace[i], nodeIdsFace[i]);
    }
}

//Sides
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

//Superelements
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
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.face_topology(0));
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.edge_topology(0));
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
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.face_topology(0));
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.edge_topology(0));
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.base());
    EXPECT_FALSE(invalidSuperElement.has_homogeneous_faces());
    EXPECT_FALSE(invalidSuperElement.is_shell());
}
//Done

//beginCheckForPositivePolarity
TEST(stk_topology_how_to, check_for_positive_polarity)
{
    stk::topology quad4Topology = stk::topology::QUAD_4;

    ASSERT_EQ(8u, quad4Topology.num_permutations());
    ASSERT_EQ(4u, quad4Topology.num_positive_permutations());

    EXPECT_TRUE( quad4Topology.is_positive_polarity(0));
    EXPECT_TRUE( quad4Topology.is_positive_polarity(1));
    EXPECT_TRUE( quad4Topology.is_positive_polarity(2));
    EXPECT_TRUE( quad4Topology.is_positive_polarity(3));
    EXPECT_TRUE(!quad4Topology.is_positive_polarity(4));
    EXPECT_TRUE(!quad4Topology.is_positive_polarity(5));
    EXPECT_TRUE(!quad4Topology.is_positive_polarity(6));
    EXPECT_TRUE(!quad4Topology.is_positive_polarity(7));

    //or, print it and examine the output:
    stk::verbose_print_topology(std::cout, quad4Topology);
}
//endCheckForPositivePolarity

}

#endif // __IBMCPP__
