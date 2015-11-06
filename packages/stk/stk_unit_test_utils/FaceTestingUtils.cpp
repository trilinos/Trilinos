// Copyright (c) 2015, Sandia Corporation.
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

#include "FaceTestingUtils.hpp"
#include "ioUtils.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh)
{
    std::vector<unsigned> countVec;
    stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, countVec);
    return countVec[mesh.mesh_meta_data().side_rank()];
}

unsigned read_file_create_faces_count_sides(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    stk::mesh::create_faces(mesh);
    return count_sides_in_mesh(mesh);
}

unsigned read_file_count_sides(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    return count_sides_in_mesh(mesh);
}

bool is_face_fully_connected(const stk::mesh::BulkData& mesh, stk::mesh::MeshIndex faceMeshIndex)
{
    const unsigned num_expected_faces = faceMeshIndex.bucket->topology().num_faces();
    if (num_expected_faces != mesh.num_faces(stk::mesh::impl::get_entity(faceMeshIndex)))
        return false;
    return true;
}

bool fully_connected_elements_to_faces(const stk::mesh::BulkData& mesh)
{
    bool fully_connected = true;
    stk::mesh::impl::for_each_entity_run(mesh, stk::topology::ELEMENT_RANK,
        [&fully_connected](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& meshIndex)
        {
          fully_connected &= is_face_fully_connected(mesh,meshIndex);
        }
    );
    return fully_connected;
}

unsigned read_file_create_faces_fully_connected_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    stk::mesh::create_faces(mesh);
    return fully_connected_elements_to_faces(mesh);
}

unsigned read_file_fully_connected_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    return fully_connected_elements_to_faces(mesh);
}

bool is_face_shared_between_different_elements(const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
{
    stk::mesh::Entity const * elements = mesh.begin_elements(face);
    for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count)
        for (unsigned other_elem_count = elem_count; other_elem_count < mesh.num_elements(face); ++other_elem_count)
            if ((elem_count != other_elem_count) && (elements[elem_count] != elements[other_elem_count]))
                return true;
    return false;
}

unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& mesh)
{
    unsigned shared_face_count = 0;
    stk::mesh::impl::for_each_entity_run(mesh, stk::topology::FACE_RANK,
        [&shared_face_count](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& meshIndex)
        {
          if (is_face_shared_between_different_elements(mesh,stk::mesh::impl::get_entity(meshIndex)))
              ++shared_face_count;
        }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    stk::mesh::create_faces(mesh);
    return count_shared_faces_between_different_elements(mesh);
}

unsigned read_file_shared_faces_different_elements_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    return count_shared_faces_between_different_elements(mesh);
}

bool is_face_shared_between_same_element(const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
{
    stk::mesh::Entity const * elements = mesh.begin_elements(face);
    for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count)
        for (unsigned other_elem_count = elem_count; other_elem_count < mesh.num_elements(face); ++other_elem_count)
            if ((elem_count != other_elem_count) && (elements[elem_count] == elements[other_elem_count]))
                return true;
    return false;
}

unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& mesh)
{
    unsigned shared_face_count = 0;
    stk::mesh::impl::for_each_entity_run(mesh, stk::topology::FACE_RANK,
      [&shared_face_count](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& meshIndex)
      {
        stk::mesh::Entity face = stk::mesh::impl::get_entity(meshIndex);
        if (is_face_shared_between_same_element(mesh,face)) ++shared_face_count;
      }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    stk::mesh::create_faces(mesh);
    return count_shared_faces_between_same_element(mesh);
}

unsigned read_file_shared_faces_same_elements_stk(std::string filename)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    return count_shared_faces_between_same_element(mesh);
}

bool is_node_not_at_x_equal_half(const stk::mesh::BulkData& mesh, stk::mesh::Entity node)
{
    double *xyz = static_cast<double *>(stk::mesh::field_data(*mesh.mesh_meta_data().coordinate_field(), node));
    return (xyz[0] != 0.5);
}

bool is_face_at_x_equal_half(const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
{
    stk::mesh::Entity const * node = mesh.begin_nodes(face);
    for (unsigned node_count = 0; node_count < mesh.num_nodes(face); ++node_count)
        if (is_node_not_at_x_equal_half(mesh,node[node_count]))
            return false;
    return true;
}

stk::mesh::EntityVector get_faces_at_x_equal_half(const stk::mesh::BulkData& mesh)
{
    stk::mesh::EntityVector faces_at_x_equal_half;
    stk::mesh::impl::for_each_entity_run(mesh, stk::topology::FACE_RANK,
      [&faces_at_x_equal_half](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& meshIndex)
      {
        if (is_face_at_x_equal_half(mesh,stk::mesh::impl::get_entity(meshIndex)))
            faces_at_x_equal_half.push_back(stk::mesh::impl::get_entity(meshIndex));
      }
    );
    return faces_at_x_equal_half;
}

bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::vector<int>& countsIn, bool debug)
{
    std::vector<int> counts(6,-1);
    std::copy(countsIn.begin(),countsIn.end(),counts.begin());
    stk::mesh::EntityVector faces_at_x_equal_half = get_faces_at_x_equal_half(mesh);
    bool extra_face_not_accounted_for = false;
    for(stk::mesh::Entity face : faces_at_x_equal_half) {
        if (is_face_at_x_equal_half(mesh,face)) {
            if (counts[0] == static_cast<int>(mesh.num_elements(face))) {
                counts[0] = -1;
            }
            else if (counts[1] == static_cast<int>(mesh.num_elements(face))) {
                counts[1] = -1;
            }
            else if (counts[2] == static_cast<int>(mesh.num_elements(face))) {
                counts[2] = -1;
            }
            else if (counts[3] == static_cast<int>(mesh.num_elements(face))) {
                counts[3] = -1;
            }
            else if (counts[4] == static_cast<int>(mesh.num_elements(face))) {
                counts[4] = -1;
            }
            else if (counts[5] == static_cast<int>(mesh.num_elements(face))) {
                counts[5] = -1;
            }
            else {
                extra_face_not_accounted_for = true;
            }
            if (debug) {
                std::cout << "num_elements:" << mesh.num_elements(face) << std::endl;
                stk::mesh::Entity const * elements = mesh.begin_elements(face);
                for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                    std::cout << "elem_count " << elem_count << " id " << mesh.entity_key(elements[elem_count])  << " for face " << mesh.entity_key(face) << std::endl;
                }
                std::cout.flush();
            }
        }
    }
    if (!debug && counts[5] == -1 && counts[4] == -1 &&  counts[3] == -1 && counts[2] == -1 && counts[1] == -1 && counts[0] == -1 && !extra_face_not_accounted_for) {
        return true;
    }
    else {
        if (debug) {
            return false;
        }
        else {
            return check_face_elem_connectivity(mesh, counts, true);
        }
    }
}

bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::vector<int>& counts)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    stk::mesh::create_faces(mesh);
    return check_face_elem_connectivity(mesh, counts);

}

bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::vector<int>& counts)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);
    return check_face_elem_connectivity(mesh, counts);

}
