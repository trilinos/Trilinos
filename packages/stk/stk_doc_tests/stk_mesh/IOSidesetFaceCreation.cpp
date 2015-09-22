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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include <stddef.h>                     // for size_t, NULL

namespace stk { namespace mesh { class BulkData; } }

namespace
{


//BEGIN2hex1sideset
bool is_positive_permutation(stk::mesh::BulkData & mesh,
                                          stk::mesh::Entity face,
                                          stk::mesh::Entity hex,
                                          unsigned face_ordinal)
{
    stk::topology faceTopology = mesh.bucket(face).topology();
    stk::mesh::EntityVector face_node_ids(mesh.num_nodes(face));
    for (unsigned face_node_count=0; face_node_count < mesh.num_nodes(face); ++face_node_count) {
        face_node_ids[face_node_count] = mesh.begin_nodes(face)[face_node_count];
    }
    stk::topology connectedElemTopology = mesh.bucket(hex).topology();
    stk::mesh::EntityVector element_side_nodes(mesh.num_nodes(face));
    connectedElemTopology.face_nodes(mesh.begin_nodes(hex), face_ordinal, &element_side_nodes[0]);
    std::pair<bool, unsigned> permutation = faceTopology.equivalent(face_node_ids, element_side_nodes);

    bool is_a_valid_permutation = permutation.first;
    EXPECT_TRUE(is_a_valid_permutation);
    bool is_positive_permutation = permutation.second < faceTopology.num_positive_permutations();
    return is_positive_permutation;
}

TEST(StkMeshHowTo, StkIO2Hex1SidesetFaceCreation)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        //  -------  |S  -------             -------  |F  -------
        //  |     |  |I  |     |             |     |  |A  |     |
        //  |HEX1 5<-|D  4 HEX2| --STK-IO--> |HEX1 5<-|C->4 HEX2|
        //  |     |  |E  |     |             |     |  |E  |     |
        //  -------  |S  -------             -------  |   -------
        //           |E                               |----> face is put into
        //           |T                                       part surface_1
        //                                            |---> orientation points outward
        //                                                   from Hex1 face5

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        stkMeshIoBroker.add_mesh_database("ALA.e", stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
        stk::mesh::EntityVector all_faces;
        stk::mesh::get_entities(mesh, stk::topology::FACE_RANK, all_faces);
        std::sort(all_faces.begin(),all_faces.end());
        unsigned expected_num_faces = 1;
        ASSERT_EQ(expected_num_faces, all_faces.size());
        size_t face_index = 0;
        stk::mesh::Entity face = all_faces[face_index];
        unsigned expected_connected_elements = 2;
        ASSERT_EQ(expected_connected_elements, mesh.num_elements(face));

        EXPECT_TRUE(mesh.bucket(face).member(*mesh.mesh_meta_data().get_part("surface_1")));

        const stk::mesh::Entity * connected_elements = mesh.begin_elements(face);
        const stk::mesh::ConnectivityOrdinal * which_side_of_element = mesh.begin_element_ordinals(face);

        {
            int element_count = 0;
            stk::mesh::Entity hex_2 = connected_elements[element_count];
            EXPECT_EQ(2u, mesh.identifier(hex_2));
            unsigned expected_face_ordinal = 4;
            EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
            EXPECT_FALSE(is_positive_permutation(
                    mesh, face, hex_2, expected_face_ordinal));
        }

        {
            int element_count = 1;
            stk::mesh::Entity hex_1 = connected_elements[element_count];
            EXPECT_EQ(1u, mesh.identifier(hex_1));
            unsigned expected_face_ordinal = 5;
            EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
            EXPECT_TRUE(is_positive_permutation(
                    mesh, face, hex_1, expected_face_ordinal));
        }

    }
}
//END2hex1sideset

//BEGIN2hex2shell3sideset
TEST(StkMeshHowTo, StkIO2Hex2Shell3SidesetFaceCreation)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        //  -------  |S |S| |S|  |S |S  -------
        //  |     |  |I |H| |H|  |I |I  |     |
        //  |HEX1 5<-|D |E| |E0<-|D |D->4 HEX2|
        //  |     |  |E |L| |L|  |E |E  |     | |
        //  -------  |S |L| |L|  |S |S  ------- |
        //           |E  3   4   |E |E          |
        //           |T          |T |T          STK
        //                                      IO
        //                                       |
        //                                       V
        //
        //  -------  |F  |S|  |S|             |F  -------
        //  |     |  |A--|H|->1H|             |A  |     |
        //  |HEX1 5<-|C->1E|  |E0<------------|C->4 HEX2|
        //  |     |  |E  |L0<-|L|-------------|E  |     |
        //  -------  |   |L|  |L|             |   -------
        //           |    3    4              |
        //           |---> orientation        |-->orientation
        //           |---> in surface_1 part  |-->in surface_2 and
        //                                         surface_3 parts


        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        stkMeshIoBroker.add_mesh_database("ALefLRA.e", stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
        stk::mesh::EntityVector all_faces;
        stk::mesh::get_entities(mesh, stk::topology::FACE_RANK, all_faces);
        std::sort(all_faces.begin(),all_faces.end());
        unsigned expected_num_faces = 2;
        ASSERT_EQ(expected_num_faces, all_faces.size());

        size_t face_index = 0;
        {
            stk::mesh::Entity face = all_faces[face_index];
            unsigned expected_connected_elements = 3;
            ASSERT_EQ(expected_connected_elements, mesh.num_elements(face));

            EXPECT_TRUE(mesh.bucket(face).member(*mesh.mesh_meta_data().get_part("surface_1")));

            const stk::mesh::Entity * connected_elements = mesh.begin_elements(face);
            const stk::mesh::ConnectivityOrdinal * which_side_of_element = mesh.begin_element_ordinals(face);

            {
                int element_count = 0;
                stk::mesh::Entity shell_3 = connected_elements[element_count];
                EXPECT_EQ(3u, mesh.identifier(shell_3));
                unsigned expected_face_ordinal = 1;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_FALSE(is_positive_permutation(
                        mesh, face, shell_3, expected_face_ordinal));
            }
            {
                int element_count = 1;
                stk::mesh::Entity shell_4 = connected_elements[element_count];
                EXPECT_EQ(4u, mesh.identifier(shell_4));
                unsigned expected_face_ordinal = 1;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_FALSE(is_positive_permutation(
                        mesh, face, shell_4, expected_face_ordinal));
            }
            {
                int element_count = 2;
                stk::mesh::Entity hex_1 = connected_elements[element_count];
                EXPECT_EQ(1u, mesh.identifier(hex_1));
                unsigned expected_face_ordinal = 5;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_TRUE(is_positive_permutation(
                        mesh, face, hex_1, expected_face_ordinal));
            }
        }

        face_index = 1;
        {
            stk::mesh::Entity face = all_faces[face_index];
            unsigned expected_connected_elements = 3;
            ASSERT_EQ(expected_connected_elements, mesh.num_elements(face));

            EXPECT_TRUE(mesh.bucket(face).member(*mesh.mesh_meta_data().get_part("surface_2")));
            EXPECT_TRUE(mesh.bucket(face).member(*mesh.mesh_meta_data().get_part("surface_3")));

            const stk::mesh::Entity * connected_elements = mesh.begin_elements(face);
            const stk::mesh::ConnectivityOrdinal * which_side_of_element = mesh.begin_element_ordinals(face);

            {
                int element_count = 0;
                stk::mesh::Entity shell_3 = connected_elements[element_count];
                EXPECT_EQ(3u, mesh.identifier(shell_3));
                unsigned expected_face_ordinal = 0;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_TRUE(is_positive_permutation(
                        mesh, face, shell_3, expected_face_ordinal));
            }
            {
                int element_count = 1;
                stk::mesh::Entity shell_4 = connected_elements[element_count];
                EXPECT_EQ(4u, mesh.identifier(shell_4));
                unsigned expected_face_ordinal = 0;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_TRUE(is_positive_permutation(
                        mesh, face, shell_4, expected_face_ordinal));
            }
            {
                int element_count = 2;
                stk::mesh::Entity hex_2 = connected_elements[element_count];
                EXPECT_EQ(2u, mesh.identifier(hex_2));
                unsigned expected_face_ordinal = 4;
                EXPECT_EQ(expected_face_ordinal, which_side_of_element[element_count]);
                EXPECT_FALSE(is_positive_permutation(
                        mesh, face, hex_2, expected_face_ordinal));
            }
        }
    }
}
//END2hex2shell3sideset


} //end empty namespace
