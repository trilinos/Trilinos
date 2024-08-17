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

#ifndef FACETESTINGUTILS_HPP_
#define FACETESTINGUTILS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <set>                       // for set
#include <string>                    // for string
#include "stk_mesh/base/Entity.hpp"  // for Entity
#include "stk_mesh/base/Types.hpp"   // for EntityId, EntityVector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_count_sides(std::string filename);

unsigned read_file_count_sides(std::string filename);

bool fully_connected_elements_to_faces(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_fully_connected_stk(std::string filename);

unsigned read_file_fully_connected_stk(std::string filename);

unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename);

unsigned read_file_shared_faces_different_elements_stk(std::string filename);

unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& mesh);

unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename);

unsigned read_file_shared_faces_same_elements_stk(std::string filename);

bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::set<unsigned>& counts);

bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned count_sides_in_mesh(const stk::mesh::BulkData& mesh);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_create_faces_count_sides(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_count_sides(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
bool fully_connected_elements_to_faces(const stk::mesh::BulkData& mesh);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_create_faces_fully_connected_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_fully_connected_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned count_shared_faces_between_different_elements(const stk::mesh::BulkData& mesh);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_create_faces_shared_faces_different_elements_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_shared_faces_different_elements_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned count_shared_faces_between_same_element(const stk::mesh::BulkData& mesh);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_create_faces_shared_faces_same_elements_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
unsigned read_file_shared_faces_same_elements_stk(std::string filename);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
bool check_face_elem_connectivity(const stk::mesh::BulkData& mesh, const std::set<unsigned>& counts);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
bool read_file_create_faces_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
bool read_file_check_face_elem_connectivity_stk(std::string filename, const std::set<unsigned>& counts);

} // namespace simple_fields

namespace stk
{
namespace unit_test_util
{
stk::mesh::Entity declare_element_side_with_nodes(stk::mesh::BulkData &mesh,
                                                  stk::mesh::Entity elem,
                                                  const stk::mesh::EntityVector &nodes,
                                                  stk::mesh::EntityId globalId,
                                                  stk::mesh::Part &part);

stk::mesh::Entity declare_element_to_edge_with_nodes(stk::mesh::BulkData &mesh,
                                                     stk::mesh::Entity elem,
                                                     const stk::mesh::EntityVector &sub_topology_nodes,
                                                     stk::mesh::EntityId global_sub_topology_id,
                                                     stk::mesh::Part &part);

stk::mesh::Part *get_surface_part_with_id(const stk::mesh::MetaData &meta, int id);

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Entity declare_element_side_with_nodes(stk::mesh::BulkData &mesh,
                                                  stk::mesh::Entity elem,
                                                  const stk::mesh::EntityVector &nodes,
                                                  stk::mesh::EntityId globalId,
                                                  stk::mesh::Part &part);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Entity declare_element_to_edge_with_nodes(stk::mesh::BulkData &mesh,
                                                     stk::mesh::Entity elem,
                                                     const stk::mesh::EntityVector &sub_topology_nodes,
                                                     stk::mesh::EntityId global_sub_topology_id,
                                                     stk::mesh::Part &part);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
stk::mesh::Part *get_surface_part_with_id(const stk::mesh::MetaData &meta, int id);

} // namespace simple_fields

}
}

#endif // FACETESTINGUTILS_HPP_
