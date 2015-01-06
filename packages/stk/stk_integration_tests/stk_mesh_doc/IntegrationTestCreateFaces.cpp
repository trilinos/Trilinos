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
            mesh.begin_elements(face);
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
            mesh.begin_elements(face);
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

unsigned read_file_count_sides(std::string filename, bool create_faces) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
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
    return 0u;
}

unsigned read_file_fully_connected(std::string filename, bool create_faces) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
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
    return false;
}

unsigned read_file_shared_faces_different_elements(std::string filename, bool create_faces) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
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
    return 0u;
}

unsigned read_file_shared_faces_same_elements(std::string filename, bool create_faces) {
    const int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numprocs == 1) {
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
    return 0u;
}

TEST(StkIo, CreateFacesCountFaces)
{
    EXPECT_EQ(11u,  read_file_count_sides("AA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("AB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADDA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ADDB.e", true));
    EXPECT_EQ(6u,   read_file_count_sides("AD.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADeDA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ADeDB.e", true));
    EXPECT_EQ(7u,   read_file_count_sides("ADe.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADeLA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ADeLB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ADeRA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ADeRB.e", true));
    EXPECT_EQ(6u,   read_file_count_sides("A.e", true));
    EXPECT_EQ(11u,  read_file_count_sides("AeA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("AeB.e", true));
    EXPECT_EQ(6u,   read_file_count_sides("Ae.e", true));
    EXPECT_EQ(11u,  read_file_count_sides("AefA.e", true));
    EXPECT_EQ(11u,  read_file_count_sides("AefB.e", true));
    EXPECT_EQ(6u,   read_file_count_sides("Aef.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ALA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ALB.e", true));
    EXPECT_EQ(6u,   read_file_count_sides("AL.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ALeDA.e", true));
    EXPECT_EQ(14u,  read_file_count_sides("ALeDB.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ALeDeRA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ALeDfRA.e", true));
    EXPECT_EQ(14u,  read_file_count_sides("ALeDfRB.e", true));
    EXPECT_EQ(7u,   read_file_count_sides("ALe.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("ALeLA.e", true));
    EXPECT_EQ(14u,  read_file_count_sides("ALeLB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ALeRA.e", true));
    EXPECT_EQ(14u,  read_file_count_sides("ALeRB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ALRA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ALRB.e", true));
    EXPECT_EQ(11u,  read_file_count_sides("ARA.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("ARB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("AReDA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("AReDB.e", true));
    EXPECT_EQ(7u,   read_file_count_sides("ARe.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("AReLA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("AReLB.e", true));
    EXPECT_EQ(12u,  read_file_count_sides("AReRA.e", true));
    EXPECT_EQ(13u,  read_file_count_sides("AReRB.e", true));
    EXPECT_EQ(2u,   read_file_count_sides("eD.e", true));
    EXPECT_EQ(1u,   read_file_count_sides("e.e", true));
    EXPECT_EQ(2u,   read_file_count_sides("eL.e", true));
}

TEST(StkIo, CreateFacesFullyConnected)
{
    EXPECT_TRUE(read_file_fully_connected("AA.e", true));
    EXPECT_TRUE(read_file_fully_connected("AB.e", true));
    EXPECT_TRUE(read_file_fully_connected("ADA.e", true));
    EXPECT_TRUE(read_file_fully_connected("ADB.e", true));
    EXPECT_TRUE(read_file_fully_connected("ADDA.e", true));
    EXPECT_TRUE(read_file_fully_connected("ADDB.e", true));
    EXPECT_TRUE(read_file_fully_connected("AD.e", true));
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
    EXPECT_TRUE(read_file_fully_connected("ALeDeRA.e", true));
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
    EXPECT_TRUE(read_file_fully_connected("eD.e", true));
    EXPECT_TRUE(read_file_fully_connected("e.e", true));
    EXPECT_TRUE(read_file_fully_connected("eL.e", true));
}

TEST(StkIo, CreateFacesSharedFacesDifferentElements)
{
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADDA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADDB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AD.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeDA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeDB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADe.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeLA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeLB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ADeRB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("A.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AeA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AeB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("Ae.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AefA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AefB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("Aef.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AL.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDeRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDfRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeDfRB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALe.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeLA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeLB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ALeRA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeRB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALRB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ARA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ARB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReDA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReDB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("ARe.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReLA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReLB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_different_elements("AReRB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("eD.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("e.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("eL.e", true));
}

TEST(StkIo, CreateFacesSharedFacesSameElements)
{
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AD.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeDB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADe.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeLB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADeRB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("A.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("AeA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("AeB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("Ae.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("AefA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("AefB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("Aef.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AL.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDB.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDeRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("ALeDfRA.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("ALeDfRB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("ALe.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLA.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeLB.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("ALeRA.e", true));
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
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("eD.e", true));
    EXPECT_EQ(1u, read_file_shared_faces_same_elements("e.e", true));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("eL.e", true));
}

TEST(StkIo, NoCreateFaces)
{
    EXPECT_EQ(0u,  read_file_count_sides("AA.e", false));
    EXPECT_EQ(0u,  read_file_count_sides("AB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADDA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADDB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AD.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADeDA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADeDB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADe.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADeLA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADeLB.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ADeRA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ADeRB.e", false));
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
    EXPECT_EQ(3u,  read_file_count_sides("ALeDeRA.e", false));
    EXPECT_EQ(3u,  read_file_count_sides("ALeDfRA.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALeDfRB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ALe.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALeLA.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALeLB.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALeRA.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALeRB.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("ALRA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ALRB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ARA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ARB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AReDA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AReDB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("ARe.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AReLA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AReLB.e", false));
    EXPECT_EQ(2u,  read_file_count_sides("AReRA.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("AReRB.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("eD.e", false));
    EXPECT_EQ(0u,  read_file_count_sides("e.e", false));
    EXPECT_EQ(1u,  read_file_count_sides("eL.e", false));
}
TEST(StkIo, NoCreateFacesFullyConnected)
{
    EXPECT_FALSE(read_file_fully_connected("AA.e", false));
    EXPECT_FALSE(read_file_fully_connected("AB.e", false));
    EXPECT_FALSE(read_file_fully_connected("ADA.e", false));
    EXPECT_FALSE(read_file_fully_connected("ADB.e", false));
    EXPECT_FALSE(read_file_fully_connected("ADDA.e", false));
    EXPECT_FALSE(read_file_fully_connected("ADDB.e", false));
    EXPECT_FALSE(read_file_fully_connected("AD.e", false));
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
    EXPECT_FALSE(read_file_fully_connected("ALeDeRA.e", false));
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
    EXPECT_FALSE(read_file_fully_connected("eD.e", false));
    EXPECT_FALSE(read_file_fully_connected("e.e", false));
    EXPECT_FALSE(read_file_fully_connected("eL.e", false));
}

TEST(StkIo, NoCreateFacesSharedFacesDifferentElements)
{
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ADDB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("AD.e", false));
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
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("ALeDeRA.e", false));
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
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("eD.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("e.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_different_elements("eL.e", false));
}

TEST(StkIo, NoCreateFacesSharedFacesSameElements)
{
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDA.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ADDB.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("AD.e", false));
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
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("ALeDeRA.e", false));
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
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("eD.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("e.e", false));
    EXPECT_EQ(0u, read_file_shared_faces_same_elements("eL.e", false));
}

} // empty namespace
