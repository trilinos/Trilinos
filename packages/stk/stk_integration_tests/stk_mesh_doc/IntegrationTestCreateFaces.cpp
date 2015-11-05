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

#include <stddef.h>                     // for size_t
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp>  // for create_faces
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/GearsFixture.hpp>  // for GearsFixture, etc
#include <stk_mesh/fixtures/TetFixture.hpp>  // for TetFixture
#include <stk_mesh/fixtures/degenerate_mesh.hpp>  // for VectorFieldType, etc
#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/FaceTestingUtils.hpp>

namespace
{


bool check_face_elem_connectivity(stk::mesh::BulkData& mesh, int count1=-1, int count2=-1, int count3=-1, int count4=-1, int count5=-1, int count6=-1, bool debug=false) {

    stk::mesh::FieldBase const * coord = mesh.mesh_meta_data().coordinate_field();
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    bool extra_face_not_accounted_for = false;
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool all_x_equal_half = true;
            stk::mesh::Entity const * node = mesh.begin_nodes(face);
            for (unsigned node_count = 0; node_count < mesh.num_nodes(face); ++node_count) {
                double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, node[node_count]));
                if (xyz[0] != 0.5) {
                    all_x_equal_half = false;
                    break;
                }
            }
            if (all_x_equal_half) {
                if (count1 == static_cast<int>(mesh.num_elements(face))) {
                    count1 = -1;
                }
                else if (count2 == static_cast<int>(mesh.num_elements(face))) {
                    count2 = -1;
                }
                else if (count3 == static_cast<int>(mesh.num_elements(face))) {
                    count3 = -1;
                }
                else if (count4 == static_cast<int>(mesh.num_elements(face))) {
                    count4 = -1;
                }
                else if (count5 == static_cast<int>(mesh.num_elements(face))) {
                    count5 = -1;
                }
                else if (count6 == static_cast<int>(mesh.num_elements(face))) {
                    count6 = -1;
                }
                else {
                    extra_face_not_accounted_for = true;
                }
                if (debug) {
                    std::cout << "num_elements:" << mesh.num_elements(face) << std::endl;
                    stk::mesh::Entity const * elements = mesh.begin_elements(face);
                    for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                        std::cout << "elem_count " << elem_count << " id " << mesh.entity_key(elements[elem_count])  << " for face_count " << face_count << " for bucket count " << bucket_count << std::endl;
                    }
                    std::cout.flush();
                }
            }
        }
    }
    if (!debug && count6 == -1 && count5 == -1 &&  count4 == -1 && count3 == -1 && count2 == -1 && count1 == -1 && !extra_face_not_accounted_for) {
        return true;
    }
    else {
        if (debug) {
            return false;
        }
        else {
            return check_face_elem_connectivity(mesh, count1, count2, count3, count4, count5, count6, true);
        }
    }
}

bool read_file_check_face_elem_connectivity(std::string filename, bool create_faces, int count1=-1, int count2=-1, int count3=-1, int count4=-1, int count5=-1, int count6=-1) {

    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return check_face_elem_connectivity(mesh, count1, count2, count3, count4, count5, count6);

}

void read_file_dump_mesh(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    std::cout << "filename: " << filename << std::endl;
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    mesh.dump_all_mesh_info(std::cout);

}


//The Magical Alphabet of Hexes, Shells & Sidesets
//
// A = hex in block A
// B = hex in block B
// e = shell in block E
// f = shell in block F
// L = sideset associated with the side on the left
// R = "        "           "   "  "     "   " right
// D = sideset containing 2 sides, one associated to left and one to right
// X = sideset associated with all sides on this surface
// S = sideset associated with the surface itself (doesn't exist, unpronounceable)
// H = a sideset composed of sides to the right, can be placed multiply in the word
// J = two hexes in block A connected to the same 8 nodes
// K = two hexes in block B connected to the same 8 nodes
//
// .e = the language of our Patron Saint Exodus
//
// RR = pronounced like a pirate
// RRR = roll the R

TEST(StkIo, NoCreateFaces)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u,  read_file_count_sides("AA.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADDA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADDB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeDA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeDB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADe.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeLB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADeRB.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("A.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AeA.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AeB.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("Ae.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AefA.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AefB.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("Aef.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("AL.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeDA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeDB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeDfRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeDfRB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALe.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeLB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeRB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALRA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALRB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReDA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReDB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARe.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReLB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReRB.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("e.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("eL.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeXfRA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADReA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ADReB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALReA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALReB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARReA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARReB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("Re.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ReL.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALefRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ARefLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AeDfA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALJ.e", false));
    }
}

TEST(StkIo, NoCreateFacesFullyConnected)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_FALSE(read_file_fully_connected_stk("AA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADDA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADDB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeDA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeDB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADe.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeLA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeLB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADeRB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("A.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AeA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AeB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("Ae.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AefA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AefB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("Aef.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AL.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeDA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeDB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeDfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeDfRB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALe.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeLA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeLB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeRB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALRB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReDA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReDB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARe.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReLA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReLB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AReRB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("e.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("eL.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALeXfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADReA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ADReB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALReA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALReB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARReA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARReB.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("Re.e", false));
        EXPECT_TRUE(read_file_fully_connected_stk("ReL.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALefRA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ARefLA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("AeDfA.e", false));
        EXPECT_FALSE(read_file_fully_connected_stk("ALJ.e", false));
    }
}

TEST(StkIo, NoCreateFacesSharedFacesDifferentElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADDA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADDB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeDA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeDB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADe.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeLA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeLB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeRA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ADeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("A.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AeA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AeB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("Ae.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AefA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AefB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("Aef.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("AL.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeDA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeDB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeDfRA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeDfRB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALe.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeLA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeLB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeRA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeRB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALRA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALRB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ARA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ARB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReDA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReDB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ARe.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReLA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReLB.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReRA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AReRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("e.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("eL.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALeXfRA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADReA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ADReB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALReA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALReB.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ARReA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ARReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("Re.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements_stk("ReL.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ALefRA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("ARefLA.e", false));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements_stk("AeDfA.e", false));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements_stk("ALJ.e", false));
    }
}

TEST(StkIo, NoCreateFacesSharedFacesSameElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("A.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AeA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AeB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("Ae.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AefA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AefB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("Aef.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeDfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeDfRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AReRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("e.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("eL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALeXfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ADReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("Re.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ReL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALefRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ARefLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("AeDfA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements_stk("ALJ.e", false));
    }
}

TEST(StkIo, CreateFacesCheckFaceElemConnectivity)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_TRUE(read_file_check_face_elem_connectivity("A.e",       true, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("e.e",       true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Ae.e",      true, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AL.e",      true, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("eL.e",      true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Re.e",      true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ReL.e",      true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AA.e",      true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AB.e",      true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADA.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADB.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     true, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeA.e",     true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeB.e",     true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Aef.e",     true, 3, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALB.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     true, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARA.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARB.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARe.e",     true, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDA.e",    true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDB.e",    true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefA.e",    true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefB.e",    true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRA.e",    true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRB.e",    true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRB.e", true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeXfRA.e", true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReA.e",   true, 3, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReB.e",   true, 3, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e",   true, 3, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReB.e",   true, 3, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e",   true, 1, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReB.e",   true, 1, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALefRA.e",   true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARefLA.e",   true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeDfA.e",   true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALJ.e",   true, 3));
    }
}

TEST(StkIo, FixSidesetIssues) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   false, 2, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", false, 3, 3));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     false, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     false, 2, 1));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e",   false, 2, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeLRA.e",   false, 2, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e",   false, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   false, 2));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALefRA.e",  false, 1, 5));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARefLA.e",  false, 3, 3));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeDfA.e",   false, 3, 3));
//        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALJ.e",     false, 3));
    }
}

TEST(StkIo, FixSidesetIssuesNoShells) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRA.e",    false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALJ.e",     false, 3));
    }
}

TEST(StkIo, NoCreateFacesCheckFaceElemConnectivity)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_TRUE(read_file_check_face_elem_connectivity("A.e",       false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("e.e",       false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AL.e",      false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("eL.e",      false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Re.e", false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ReL.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AA.e",      false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AB.e",      false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Ae.e",      false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADB.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeA.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeB.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Aef.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALB.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARB.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARe.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDA.e",    false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDB.e",    false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefA.e",    false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefB.e",    false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRA.e",    false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRB.e",    false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRB.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRB.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeXfRA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReA.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReB.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReB.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReB.e", false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALefRA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARefLA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeDfA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALJ.e", false, 3));
    }
}

TEST(StkIo, confirm_face_connectivity_AeA)
{
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database("AeA.e", stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

    stk::mesh::create_faces(mesh);

    //want to specifically check this connectivity
    //
    //  -----   F   S   F   -----
    //  |   |   A   H   A   |   |
    //  |HEX|-5>C<0-E-1>C<4-|HEX|
    //  | A |   E   L   E   | B |
    //  -----   D   L   C   -----
    //
    // -0> means the 0th ordinal face connection, -1> means the 1st, etc.
    //
    // since this is the current behavior and want to test it in STK as well

    bool face1okay = false;
    bool face2okay = false;
    typedef std::vector<stk::mesh::EntityId>  EntityIdVector;
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            stk::mesh::EntityKey face_key = mesh.entity_key(face);
            stk::mesh::Entity const * elements = mesh.begin_elements(face);
            for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                for (unsigned other_elem_count = elem_count;
                        other_elem_count < mesh.num_elements(face); ++other_elem_count) {
                    if ((elem_count != other_elem_count) &&
                            (elements[elem_count] != elements[other_elem_count])) {
                        stk::mesh::Entity elem1 = elements[elem_count];
                        stk::mesh::Entity const * face_array1 = mesh.begin_faces(elem1);
                        stk::mesh::Entity elem2 = elements[other_elem_count];
                        stk::mesh::Entity const * face_array2 = mesh.begin_faces(elem2);
                        if (mesh.bucket(elem1).topology() == stk::topology::HEX_8 && mesh.entity_key(face_array1[5]) == face_key &&
                                mesh.bucket(elem2).topology() == stk::topology::SHELL_QUAD_4 && mesh.entity_key(face_array2[0]) == face_key) {
                            face1okay = true;

                            stk::topology faceTopology = mesh.bucket(face).topology();
                            EntityIdVector face_on_shell_node_ids(faceTopology.num_nodes());
                            EntityIdVector face_on_hex_node_ids(faceTopology.num_nodes());

                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_hex_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array1[5])[count]).id();
                            }
                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_shell_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array2[0])[count]).id();
                            }

                            unsigned permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(face_on_shell_node_ids, face_on_hex_node_ids);
                            EXPECT_LT(permutation, faceTopology.num_positive_permutations());
                        }
                        if (mesh.bucket(elem1).topology() == stk::topology::SHELL_QUAD_4 &&  mesh.entity_key(face_array1[0]) == face_key &&
                                mesh.bucket(elem2).topology() == stk::topology::HEX_8 &&  mesh.entity_key(face_array2[5]) == face_key) {
                            face1okay = true;
                            stk::topology faceTopology = mesh.bucket(face).topology();
                            EntityIdVector face_on_shell_node_ids(faceTopology.num_nodes());
                            EntityIdVector face_on_hex_node_ids(faceTopology.num_nodes());

                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_hex_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array2[5])[count]).id();
                            }
                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_shell_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array1[0])[count]).id();
                            }

                            unsigned permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(face_on_shell_node_ids, face_on_hex_node_ids);
                            EXPECT_LT(permutation, faceTopology.num_positive_permutations());
                        }
                        if (mesh.bucket(elem1).topology() == stk::topology::HEX_8 &&  mesh.entity_key(face_array1[4]) == face_key &&
                                mesh.bucket(elem2).topology() == stk::topology::SHELL_QUAD_4 &&  mesh.entity_key(face_array2[1]) == face_key) {
                            face2okay = true;
                            stk::topology faceTopology = mesh.bucket(face).topology();
                            EntityIdVector face_on_shell_node_ids(faceTopology.num_nodes());
                            EntityIdVector face_on_hex_node_ids(faceTopology.num_nodes());

                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_hex_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array1[4])[count]).id();
                            }
                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_shell_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array2[1])[count]).id();
                            }

                            unsigned permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(face_on_shell_node_ids, face_on_hex_node_ids);
                            EXPECT_LT(permutation, faceTopology.num_positive_permutations());
                        }
                        if (mesh.bucket(elem1).topology() == stk::topology::SHELL_QUAD_4 &&  mesh.entity_key(face_array1[1]) == face_key &&
                                mesh.bucket(elem2).topology() == stk::topology::HEX_8 &&  mesh.entity_key(face_array2[4]) == face_key) {
                            face2okay = true;

                            stk::topology faceTopology = mesh.bucket(face).topology();
                            EntityIdVector face_on_shell_node_ids(faceTopology.num_nodes());
                            EntityIdVector face_on_hex_node_ids(faceTopology.num_nodes());

                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_hex_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array2[4])[count]).id();
                            }
                            for (unsigned count=0; count < faceTopology.num_nodes(); ++count) {
                                face_on_shell_node_ids[count] = mesh.entity_key(mesh.begin_nodes(face_array1[1])[count]).id();
                            }

                            unsigned permutation = faceTopology.lexicographical_smallest_permutation_preserve_polarity(face_on_shell_node_ids, face_on_hex_node_ids);
                            EXPECT_LT(permutation, faceTopology.num_positive_permutations());

                        }


                    }
                }
            }
        }
    }
    EXPECT_TRUE(face1okay);
    EXPECT_TRUE(face2okay);
}

TEST(StkIo, DISABLED_dump_mesh)
{
//    read_file_dump_mesh("ALeDfRA.e", true);
//    read_file_dump_mesh("ALeXfRA.e", true);
    read_file_dump_mesh("ADe.e", false);
}

} // empty namespace
