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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include <vector>                       // for vector
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH

namespace stk { namespace mesh { class BulkData; } }

namespace
{


//BEGIN2hex
TEST(StkMeshHowTo, CreateFacesTwoHexes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    //  -----------
    //  |    |    |
    //  |HEX1|HEX2|
    //  |    |    |
    //  -----------
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database("AA.e", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::create_all_sides(mesh, mesh.mesh_meta_data().universal_part());

    //  ------  F  ------
    //  |    |  A  |    |
    //  |HEX1|<-C->|HEX2|   Also external faces!
    //  |    |  E  |    |
    //  ------  !  ------

    unsigned first_bucket = 0;
    unsigned first_element_in_bucket = 0;
    stk::mesh::Entity first_element = (*mesh.buckets(stk::topology::ELEMENT_RANK)[first_bucket])[first_element_in_bucket];
    stk::mesh::Entity internal_face = mesh.begin_faces(first_element)[5];

    unsigned num_elements_connected_to_single_face = 2;
    EXPECT_EQ(num_elements_connected_to_single_face, mesh.num_elements(internal_face));

    unsigned num_expected_external_faces = 10u;
    unsigned num_expected_internal_faces = 1u;
    unsigned num_expected_faces = num_expected_external_faces + num_expected_internal_faces;
    stk::mesh::Selector all_entities = mesh.mesh_meta_data().universal_part();
    std::vector<size_t> entity_counts;
    stk::mesh::count_entities(all_entities, mesh, entity_counts);
    EXPECT_EQ(num_expected_faces, entity_counts[stk::topology::FACE_RANK]);
  }
}
//END2hex

//BEGINshell
TEST(StkMeshHowTo, CreateFacesSingleShell)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    //  S
    //  H
    //  E
    //  L
    //  L
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database("e.e", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::create_all_sides(mesh, mesh.mesh_meta_data().universal_part());

    //  F  S  F
    //  A  H  A
    //  C->E<-C
    //  E  L  E
    //  1  L  2

    unsigned first_bucket = 0;
    unsigned first_element_in_bucket = 0;
    stk::mesh::Entity first_element = (*mesh.buckets(stk::topology::ELEMENT_RANK)[first_bucket])[first_element_in_bucket];
    stk::mesh::Entity face_one = mesh.begin_faces(first_element)[0];
    unsigned num_elements_connected_to_face_one = 1;
    EXPECT_EQ(num_elements_connected_to_face_one, mesh.num_elements(face_one));

    stk::mesh::Entity face_two = mesh.begin_faces(first_element)[1];
    unsigned num_elements_connected_to_face_two = 1;
    EXPECT_EQ(num_elements_connected_to_face_two, mesh.num_elements(face_two));

    EXPECT_NE(face_one, face_two);

    unsigned num_expected_faces = 2u;
    stk::mesh::Selector all_entities = mesh.mesh_meta_data().universal_part();
    std::vector<size_t> entity_counts;
    stk::mesh::count_entities(all_entities, mesh, entity_counts);
    EXPECT_EQ(num_expected_faces, entity_counts[stk::topology::FACE_RANK]);
  }
}
//ENDshell

//BEGINhexshellhex
TEST(StkMeshHowTo, CreateFacesTwoHexesInternalShell)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    //  ------S------
    //  |    |H|    |
    //  |HEX1|E|HEX2|
    //  |    |L|    |
    //  ------L------
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database("AeA.e", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::create_all_sides(mesh, mesh.mesh_meta_data().universal_part());

    //  ------  F  S  F  ------
    //  |    |  A  H  A  |    |
    //  |HEX1|<-C->E<-C->|HEX2|   Also external faces!
    //  |    |  E  L  E  |    |
    //  ------  1  L  2  ------

    unsigned first_bucket = 0;
    unsigned first_element_in_bucket = 0;
    stk::mesh::Entity first_element = (*mesh.buckets(stk::topology::ELEMENT_RANK)[first_bucket])[first_element_in_bucket];
    stk::mesh::Entity internal_face_one = mesh.begin_faces(first_element)[5];
    unsigned num_elements_connected_to_face_one = 2;
    EXPECT_EQ(num_elements_connected_to_face_one, mesh.num_elements(internal_face_one));

    unsigned second_element_in_bucket = 1;
    stk::mesh::Entity second_element = (*mesh.buckets(stk::topology::ELEMENT_RANK)[first_bucket])[second_element_in_bucket];
    stk::mesh::Entity internal_face_two = mesh.begin_faces(second_element)[4];
    unsigned num_elements_connected_to_face_two = 2;
    EXPECT_EQ(num_elements_connected_to_face_two, mesh.num_elements(internal_face_two));

    EXPECT_NE(internal_face_one, internal_face_two);

    unsigned num_expected_external_faces = 10u;
    unsigned num_expected_internal_faces = 2u;
    unsigned num_expected_faces = num_expected_external_faces + num_expected_internal_faces;
    stk::mesh::Selector all_entities = mesh.mesh_meta_data().universal_part();
    std::vector<size_t> entity_counts;
    stk::mesh::count_entities(all_entities, mesh, entity_counts);
    EXPECT_EQ(num_expected_faces, entity_counts[stk::topology::FACE_RANK]);
  }
}
//ENDhexshellhex


} //end empty namespace
