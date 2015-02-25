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
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
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

namespace
{

bool fully_connected_elements_to_faces(stk::mesh::BulkData& mesh)
{
    bool fully_connected = true;
    stk::mesh::BucketVector const & elem_buckets = mesh.buckets(stk::topology::ELEMENT_RANK);
    for (size_t bucket_count=0, bucket_end=elem_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *elem_buckets[bucket_count];
        const unsigned num_expected_faces = bucket.topology().num_faces();
        for (size_t elem_count=0, elem_end=bucket.size(); elem_count < elem_end; ++elem_count) {
            stk::mesh::Entity elem = bucket[elem_count];
            if (num_expected_faces != mesh.num_faces(elem)) {
                fully_connected = false;
                break;
            }
        }
    }
    return fully_connected;
}

unsigned count_shared_faces_between_different_elements(stk::mesh::BulkData& mesh) {
    unsigned shared_face_count = 0;
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool is_face_shared = false;
            stk::mesh::Entity const * elements = mesh.begin_elements(face);
            for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                for (unsigned other_elem_count = elem_count;
                        other_elem_count < mesh.num_elements(face); ++other_elem_count) {
                    if ((elem_count != other_elem_count) &&
                            (elements[elem_count] != elements[other_elem_count])) {
                        is_face_shared = true;
                        break;
                    }
                }
            }
            if (is_face_shared) {
                ++shared_face_count;
            }
        }
    }
    return shared_face_count;
}

unsigned count_shared_faces_between_same_element(stk::mesh::BulkData& mesh) {
    unsigned shared_face_count = 0;
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool is_face_shared = false;
            stk::mesh::Entity const * elements = mesh.begin_elements(face);
            for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                for (unsigned other_elem_count = elem_count;
                        other_elem_count < mesh.num_elements(face); ++other_elem_count) {
                    if ((elem_count != other_elem_count) &&
                            (elements[elem_count] == elements[other_elem_count])) {
                        is_face_shared = true;
                        break;
                    }
                }
            }
            if (is_face_shared) {
                ++shared_face_count;
            }
        }
    }
    return shared_face_count;
}

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

unsigned read_file_count_sides(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    std::vector<unsigned> countVec;
    stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, countVec);
    return countVec[2];

}

unsigned read_file_fully_connected(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return fully_connected_elements_to_faces(mesh);

}

unsigned read_file_shared_faces_different_elements(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return count_shared_faces_between_different_elements(mesh);
}

unsigned read_file_shared_faces_same_elements(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return count_shared_faces_between_same_element(mesh);

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
//
// .e = the language of our Patron Saint Exodus
//
// RR = pronounced like a pirate
// RRR = roll the R

TEST(StkIo, CreateFacesCountFaces)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(11u,  read_file_count_sides("AA.e", true));
        EXPECT_EQ(11u,  read_file_count_sides("AB.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ADA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ADB.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ADDA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ADDB.e", true));
        EXPECT_EQ(14u,  read_file_count_sides("ADeDA.e", true));
        EXPECT_EQ(14u,  read_file_count_sides("ADeDB.e", true));
        EXPECT_EQ(8u,   read_file_count_sides("ADe.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ADeLA.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ADeLB.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ADeRA.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ADeRB.e", true));
        EXPECT_EQ(6u,   read_file_count_sides("A.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AeA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AeB.e", true));
        EXPECT_EQ(7u,   read_file_count_sides("Ae.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AefA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AefB.e", true));
        EXPECT_EQ(7u,   read_file_count_sides("Aef.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALB.e", true));
        EXPECT_EQ(6u,   read_file_count_sides("AL.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ALeDA.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("ALeDB.e", true));
        EXPECT_EQ(14u,  read_file_count_sides("ALeDfRA.e", true));
        EXPECT_EQ(14u,  read_file_count_sides("ALeDfRB.e", true));
        EXPECT_EQ(7u,   read_file_count_sides("ALe.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALeLA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALeLB.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALeRA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALeRB.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALRA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("ALRB.e", true));
        EXPECT_EQ(11u,  read_file_count_sides("ARA.e", true));
        EXPECT_EQ(11u,  read_file_count_sides("ARB.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("AReDA.e", true));
        EXPECT_EQ(13u,  read_file_count_sides("AReDB.e", true));
        EXPECT_EQ(7u,   read_file_count_sides("ARe.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AReLA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AReLB.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AReRA.e", true));
        EXPECT_EQ(12u,  read_file_count_sides("AReRB.e", true));
        EXPECT_EQ(2u,   read_file_count_sides("e.e", true));
        EXPECT_EQ(2u,   read_file_count_sides("eL.e", true));
        EXPECT_EQ(12u,   read_file_count_sides("AXA.e", true));
        EXPECT_EQ(16u,   read_file_count_sides("ALeXfRA.e", true));
        EXPECT_EQ(11u,   read_file_count_sides("ASA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ALeSfRA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ADHeHA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ADHeHB.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ADReA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ADReB.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ALHeHA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ALHeHB.e", true));
        EXPECT_EQ(14u,   read_file_count_sides("ALHeHfRA.e", true));
        EXPECT_EQ(14u,   read_file_count_sides("ALHeHfRB.e", true));
        EXPECT_EQ(12u,   read_file_count_sides("ARHeHA.e", true));
        EXPECT_EQ(12u,   read_file_count_sides("ARHeHB.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ALReA.e", true));
        EXPECT_EQ(13u,   read_file_count_sides("ALReB.e", true));
        EXPECT_EQ(12u,   read_file_count_sides("ARReA.e", true));
        EXPECT_EQ(12u,   read_file_count_sides("ARReB.e", true));
        EXPECT_EQ(2u,   read_file_count_sides("Re.e", true));
        EXPECT_EQ(2u,   read_file_count_sides("ReL.e", true));
        EXPECT_EQ(14u,   read_file_count_sides("AXeA.e", true));
    }
}

TEST(StkIo, CreateFacesFullyConnected)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_TRUE(read_file_fully_connected("AA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADDA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADDB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeDA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeDB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADe.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeLA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeLB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADeRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("A.e", true));
        EXPECT_TRUE(read_file_fully_connected("AeA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AeB.e", true));
        EXPECT_TRUE(read_file_fully_connected("Ae.e", true));
        EXPECT_TRUE(read_file_fully_connected("AefA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AefB.e", true));
        EXPECT_TRUE(read_file_fully_connected("Aef.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALB.e", true));
        EXPECT_TRUE(read_file_fully_connected("AL.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeDA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeDB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeDfRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeDfRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALe.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeLA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeLB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARB.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReDA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReDB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARe.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReLA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReLB.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("AReRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("e.e", true));
        EXPECT_TRUE(read_file_fully_connected("eL.e", true));
        EXPECT_TRUE(read_file_fully_connected("AXA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeXfRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ASA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALeSfRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADHeHA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADHeHB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADReA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ADReB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALHeHA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALHeHB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALHeHfRA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALHeHfRB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARHeHA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARHeHB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALReA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ALReB.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARReA.e", true));
        EXPECT_TRUE(read_file_fully_connected("ARReB.e", true));
        EXPECT_TRUE(read_file_fully_connected("Re.e", true));
        EXPECT_TRUE(read_file_fully_connected("ReL.e", true));
        EXPECT_TRUE(read_file_fully_connected("AXeA.e", true));
    }
}

TEST(StkIo, CreateFacesSharedFacesDifferentElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("AA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("AB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADe.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeLA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeLB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeRA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("A.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AeA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AeB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("Ae.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AefA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AefB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("Aef.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AL.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeDfRA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeDfRB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALe.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeLA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeLB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeRA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ARA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ARB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReDA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReDB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ARe.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AReLA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AReLB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AReRA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("AReRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("e.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("eL.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AXA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeXfRA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ASA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ALeSfRA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADHeHA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADHeHB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADReA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADReB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALHeHA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALHeHB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALHeHfRA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALHeHfRB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ARHeHA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ARHeHB.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALReA.e", true));
        EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALReB.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ARReA.e", true));
        EXPECT_EQ(2u, read_file_shared_faces_different_elements("ARReB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("Re.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ReL.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AXeA.e", true));
    }
}

TEST(StkIo, CreateFacesSharedFacesSameElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADe.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("A.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AeA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AeB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Ae.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AefA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AefB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Aef.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AL.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDfRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDfRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALe.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReDA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReDB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARe.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReLA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReLB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("e.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("eL.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AXA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeXfRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ASA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeSfRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADHeHA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADHeHB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADReA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADReB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHfRA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHfRB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARHeHA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARHeHB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALReA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALReB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARReA.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARReB.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Re.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ReL.e", true));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AXeA.e", true));
    }
}

TEST(StkIo, NoCreateFaces)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u,  read_file_count_sides("AA.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("AB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADDA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADDB.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ADeDA.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ADeDB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADe.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADeLA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADeLB.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADeRA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADeRB.e", false));
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
        EXPECT_EQ(3u,  read_file_count_sides("ALeDA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ALeDB.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ALeDfRA.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ALeDfRB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ALe.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeLB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALeRB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALRB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARB.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("AReDA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("AReDB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARe.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReLA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReLB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReRA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AReRB.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("e.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("eL.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("AXA.e", false));
        EXPECT_EQ(6u,  read_file_count_sides("ALeXfRA.e", false));
        EXPECT_EQ(0u,  read_file_count_sides("ASA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ALeSfRA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADHeHA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ADHeHB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADReA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ADReB.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ALHeHA.e", false));
        EXPECT_EQ(3u,  read_file_count_sides("ALHeHB.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ALHeHfRA.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("ALHeHfRB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ARHeHA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ARHeHB.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALReA.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ALReB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARReA.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("ARReB.e", false));
        EXPECT_EQ(1u,  read_file_count_sides("Re.e", false));
        EXPECT_EQ(2u,  read_file_count_sides("ReL.e", false));
        EXPECT_EQ(4u,  read_file_count_sides("AXeA.e", false));
    }
}

TEST(StkIo, NoCreateFacesFullyConnected)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_FALSE(read_file_fully_connected("AA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADDA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADDB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeDA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeDB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADe.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeLA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeLB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADeRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("A.e", false));
        EXPECT_FALSE(read_file_fully_connected("AeA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AeB.e", false));
        EXPECT_FALSE(read_file_fully_connected("Ae.e", false));
        EXPECT_FALSE(read_file_fully_connected("AefA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AefB.e", false));
        EXPECT_FALSE(read_file_fully_connected("Aef.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALB.e", false));
        EXPECT_FALSE(read_file_fully_connected("AL.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeDA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeDB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeDfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeDfRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALe.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeLA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeLB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARB.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReDA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReDB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARe.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReLA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReLB.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("AReRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("e.e", false));
        EXPECT_FALSE(read_file_fully_connected("eL.e", false));
        EXPECT_FALSE(read_file_fully_connected("AXA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeXfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ASA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALeSfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADHeHA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADHeHB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADReA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ADReB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALHeHA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALHeHB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALHeHfRA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALHeHfRB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARHeHA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARHeHB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALReA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ALReB.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARReA.e", false));
        EXPECT_FALSE(read_file_fully_connected("ARReB.e", false));
        EXPECT_FALSE(read_file_fully_connected("Re.e", false));
        EXPECT_TRUE(read_file_fully_connected("ReL.e", false));
        EXPECT_FALSE(read_file_fully_connected("AXeA.e", false));
    }
}

TEST(StkIo, NoCreateFacesSharedFacesDifferentElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("A.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AeA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AeB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("Ae.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AefA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AefB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("Aef.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDfRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AReRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("e.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("eL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AXA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeXfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ASA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeSfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALHeHfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALHeHfRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("Re.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("ReL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_different_elements("AXeA.e", false));
    }
}

TEST(StkIo, NoCreateFacesSharedFacesSameElements)
{
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("A.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AeA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AeB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Ae.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AefA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AefB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Aef.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDfRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReDA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReDB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARe.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReLA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReLB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AReRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("e.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("eL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AXA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeXfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ASA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeSfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHfRA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALHeHfRB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARHeHA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARHeHB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARReA.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ARReB.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("Re.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("ReL.e", false));
        EXPECT_EQ(0u, read_file_shared_faces_same_elements("AXeA.e", false));
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
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADA.e",     true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADB.e",     true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     true, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeA.e",     true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeB.e",     true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Aef.e",     true, 3, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALB.e",     true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     true, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARA.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARB.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARe.e",     true, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDA.e",    true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDB.e",    true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefA.e",    true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefB.e",    true, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRA.e",    true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRB.e",    true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDA.e",   true, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDB.e",   true, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", true, 2, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRB.e", true, 2, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AXA.e",     true, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeXfRA.e", true, 1, 1, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ASA.e",     true, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeSfRA.e", true, 3, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADHeHA.e",  true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADHeHB.e",  true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHA.e",  true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHB.e",  true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHfRA.e",true, 3, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHfRB.e",true, 3, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARHeHA.e",  true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARHeHB.e",  true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReB.e",   true, 2, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReB.e",   true, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AXeA.e",      true, 1, 1, 1, 1));
    }
}

TEST(StkIo, DISABLED_FixSidesetIssues) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", false, 3, 3));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     false, 2, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeLRA.e",   false, 2, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e",   false, 2));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   false, 2));
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
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADA.e",     false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADB.e",     false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADe.e",     false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeA.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AeB.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("Aef.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALA.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALB.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALe.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARA.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARB.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARe.e",     false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDA.e",    false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADDB.e",    false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefA.e",    false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AefB.e",    false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRA.e",    false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALRB.e",    false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDA.e",   false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeDB.e",   false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLA.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeLB.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRA.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADeRB.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDA.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDB.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLA.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeLB.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRA.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeRB.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDA.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReDB.e",   false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLA.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReLB.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRA.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AReRB.e",   false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRA.e", false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeDfRB.e", false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AXA.e",     false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeXfRA.e", false, 1, 1, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ASA.e",     false));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALeSfRA.e", false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADHeHA.e", false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADHeHB.e", false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReA.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ADReB.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHA.e", false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHB.e", false, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHfRA.e", false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALHeHfRB.e", false, 1, 1, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARHeHA.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARHeHB.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReA.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ALReB.e", false, 1, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReA.e", false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("ARReB.e", false, 1));
        EXPECT_TRUE(read_file_check_face_elem_connectivity("AXeA.e", false, 1, 1, 1, 1));
    }
}

TEST(StkIo, confirm_face_connectivity_AeA)
{
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
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
    read_file_dump_mesh("AeA.e", true);
}

} // empty namespace
