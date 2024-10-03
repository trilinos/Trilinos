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

#include "FaceTestingUtils.hpp"
#include "ioUtils.hpp"
#include "BuildMesh.hpp"
#include <set>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <string>
#include <vector>

using stk::unit_test_util::build_mesh;

unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh)
{
    std::vector<size_t> countVec;
    stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, countVec);
    return countVec[mesh.mesh_meta_data().side_rank()];
}

unsigned read_file_create_faces_count_sides(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return count_sides_in_mesh(*mesh);
}

unsigned read_file_count_sides(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return count_sides_in_mesh(*mesh);
}

bool is_face_fully_connected(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
{
    const unsigned num_expected_faces = mesh.bucket(elem).topology().num_faces();
    if (num_expected_faces != mesh.num_faces(elem))
        return false;
    return true;
}

bool fully_connected_elements_to_faces(const stk::mesh::BulkData& bulk)
{
    bool fully_connected = true;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::ELEMENT_RANK,
        [&fully_connected](const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
        {
          fully_connected &= is_face_fully_connected(mesh, elem);
        }
    );
    return fully_connected;
}

unsigned read_file_create_faces_fully_connected_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return fully_connected_elements_to_faces(*mesh);
}

unsigned read_file_fully_connected_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return fully_connected_elements_to_faces(*mesh);
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

unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& bulk)
{
    unsigned shared_face_count = 0;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::FACE_RANK,
        [&shared_face_count](const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
        {
          if (is_face_shared_between_different_elements(mesh, face))
              ++shared_face_count;
        }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return count_shared_faces_between_different_elements(*mesh);
}

unsigned read_file_shared_faces_different_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return count_shared_faces_between_different_elements(*mesh);
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

unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& bulk)
{
    unsigned shared_face_count = 0;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::FACE_RANK,
      [&shared_face_count](const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
      {
        if (is_face_shared_between_same_element(mesh,face))
            ++shared_face_count;
      }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return count_shared_faces_between_same_element(*mesh);
}

unsigned read_file_shared_faces_same_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return count_shared_faces_between_same_element(*mesh);
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

stk::mesh::EntityVector get_faces_at_x_equal_half(const stk::mesh::BulkData& bulk)
{
    stk::mesh::EntityVector faces_at_x_equal_half;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::FACE_RANK,
      [&faces_at_x_equal_half](const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
      {
        if (is_face_at_x_equal_half(mesh, face))
            faces_at_x_equal_half.push_back(face);
      }
    );
    return faces_at_x_equal_half;
}

std::set<unsigned> get_face_connectivity_at_x_equal_half(const stk::mesh::BulkData& mesh)
{
    std::set<unsigned> faces;
    for(stk::mesh::Entity face : get_faces_at_x_equal_half(mesh))
        faces.insert(mesh.num_elements(face));
    return faces;
}

std::ostream& operator<<(std::ostream& os, const std::set<unsigned>& data)
{
    os << "{ " << stk::util::join(data.begin(),data.end(),", ") << " }";
    return os;
}

bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::set<unsigned>& gold_faces)
{
    std::set<unsigned> current_faces = get_face_connectivity_at_x_equal_half(mesh);
    if (current_faces == gold_faces) {
        return true;
    }
    std::cout << "gold_faces = " << gold_faces << ", current_faces = " << current_faces << std::endl;
    return false;
}

bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return check_face_elem_connectivity(*mesh, counts);

}

bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return check_face_elem_connectivity(*mesh, counts);

}


namespace simple_fields {

unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh)
{
    std::vector<size_t> countVec;
    stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, countVec);
    return countVec[mesh.mesh_meta_data().side_rank()];
}

unsigned read_file_create_faces_count_sides(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return ::count_sides_in_mesh(*mesh);
}

unsigned read_file_count_sides(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return ::count_sides_in_mesh(*mesh);
}

bool fully_connected_elements_to_faces(const stk::mesh::BulkData& bulk)
{
    bool fully_connected = true;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::ELEMENT_RANK,
        [&fully_connected](const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
        {
          fully_connected &= is_face_fully_connected(mesh, elem);
        }
    );
    return fully_connected;
}

unsigned read_file_create_faces_fully_connected_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return ::fully_connected_elements_to_faces(*mesh);
}

unsigned read_file_fully_connected_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return ::fully_connected_elements_to_faces(*mesh);
}

unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& bulk)
{
    unsigned shared_face_count = 0;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::FACE_RANK,
        [&shared_face_count](const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
        {
          if (::is_face_shared_between_different_elements(mesh, face))
              ++shared_face_count;
        }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return ::count_shared_faces_between_different_elements(*mesh);
}

unsigned read_file_shared_faces_different_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return ::count_shared_faces_between_different_elements(*mesh);
}

unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& bulk)
{
    unsigned shared_face_count = 0;
    stk::mesh::for_each_entity_run_no_threads(bulk, stk::topology::FACE_RANK,
      [&shared_face_count](const stk::mesh::BulkData& mesh, stk::mesh::Entity face)
      {
        if (::is_face_shared_between_same_element(mesh,face))
            ++shared_face_count;
      }
    );
    return shared_face_count;
}

unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return ::count_shared_faces_between_same_element(*mesh);
}

unsigned read_file_shared_faces_same_elements_stk(std::string filename)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return ::count_shared_faces_between_same_element(*mesh);
}

bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::set<unsigned>& gold_faces)
{
    std::set<unsigned> current_faces = ::get_face_connectivity_at_x_equal_half(mesh);
    if (current_faces == gold_faces) {
        return true;
    }
    std::cout << "gold_faces = " << gold_faces << ", current_faces = " << current_faces << std::endl;
    return false;
}

bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    stk::mesh::create_all_sides(*mesh, mesh->mesh_meta_data().universal_part(), {}, false);
    return ::check_face_elem_connectivity(*mesh, counts);

}

bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts)
{
    std::shared_ptr<stk::mesh::BulkData> mesh = build_mesh(MPI_COMM_WORLD);
    stk::io::fill_mesh(filename, *mesh);
    return ::check_face_elem_connectivity(*mesh, counts);

}

} // namespace simple_fields


namespace stk
{
namespace unit_test_util
{

stk::mesh::Entity declare_element_side_with_nodes(stk::mesh::BulkData &mesh,
                                                  stk::mesh::Entity elem,
                                                  const stk::mesh::EntityVector &nodes,
                                                  stk::mesh::EntityId globalId,
                                                  stk::mesh::Part &part)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(mesh, elem, mesh.mesh_meta_data().side_rank(), nodes);
    return mesh.declare_element_side(elem, ordinalAndPermutation.first, stk::mesh::ConstPartVector{&part});
}

stk::mesh::Entity declare_element_to_edge_with_nodes(stk::mesh::BulkData &mesh, stk::mesh::Entity elem, const stk::mesh::EntityVector &sub_topology_nodes,
        stk::mesh::EntityId global_sub_topology_id, stk::mesh::Part &part)
{
    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
            get_ordinal_and_permutation(mesh, elem, stk::topology::EDGE_RANK, sub_topology_nodes);

    if((ordinalAndPermutation.first == stk::mesh::INVALID_CONNECTIVITY_ORDINAL) || (ordinalAndPermutation.second
            == stk::mesh::Permutation::INVALID_PERMUTATION))
    {
        stk::mesh::Entity invalid;
        invalid = stk::mesh::Entity::InvalidEntity;
        return invalid;
    }

    stk::mesh::Entity side = mesh.get_entity(stk::topology::EDGE_RANK, global_sub_topology_id);
    if(!mesh.is_valid(side))
    {
        side = mesh.declare_edge(global_sub_topology_id, stk::mesh::ConstPartVector{&part});
        for(unsigned i = 0; i < sub_topology_nodes.size(); ++i)
            mesh.declare_relation(side, sub_topology_nodes[i], i);
    }
    else
    {
        const stk::mesh::Entity* sideNodes = mesh.begin_nodes(side);
        unsigned numNodes = mesh.num_nodes(side);
        STK_ThrowRequireMsg(sub_topology_nodes.size() == numNodes,
                        "declare_element_to_sub_topology_with_nodes ERROR, side already exists with different number of nodes");
        for(unsigned i = 0; i < numNodes; ++i)
        {
            STK_ThrowRequireMsg(sub_topology_nodes[i] == sideNodes[i],
                            "declare_element_to_sub_topology_with_nodes ERROR, side already exists with different node connectivity");
        }
    }

    mesh.declare_relation(elem, side, ordinalAndPermutation.first, ordinalAndPermutation.second);
    return side;
}


stk::mesh::Part *get_surface_part_with_id(const stk::mesh::MetaData &meta, int id)
{
    const stk::mesh::PartVector &all_parts = meta.get_parts();

    for(auto part : all_parts)
    {
        if((part->primary_entity_rank() == meta.side_rank()) && (part->id() == id))
        {
            return part;
        }
    }

    return nullptr;
}

namespace simple_fields {

stk::mesh::Entity declare_element_side_with_nodes(stk::mesh::BulkData &mesh,
                                                  stk::mesh::Entity elem,
                                                  const stk::mesh::EntityVector &nodes,
                                                  stk::mesh::EntityId globalId,
                                                  stk::mesh::Part &part)
{
  return stk::unit_test_util::declare_element_side_with_nodes(mesh, elem, nodes, globalId, part);
}

stk::mesh::Entity declare_element_to_edge_with_nodes(stk::mesh::BulkData &mesh,
                                                     stk::mesh::Entity elem,
                                                     const stk::mesh::EntityVector &sub_topology_nodes,
                                                     stk::mesh::EntityId global_sub_topology_id,
                                                     stk::mesh::Part &part)
{
  return stk::unit_test_util::declare_element_to_edge_with_nodes(mesh, elem, sub_topology_nodes, global_sub_topology_id, part);
}

stk::mesh::Part *get_surface_part_with_id(const stk::mesh::MetaData &meta, int id)
{
  return stk::unit_test_util::get_surface_part_with_id(meta, id);
}

} // namespace simple_fields

}}
