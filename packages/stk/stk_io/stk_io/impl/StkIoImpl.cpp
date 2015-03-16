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

#include <stk_io/impl/StkIoImpl.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk { namespace io { namespace impl {

stk::mesh::Entity get_or_create_face_at_element_side(stk::mesh::BulkData & bulk,
                                                     stk::mesh::Entity elem,
                                                     int side_ordinal,
                                                     int new_face_global_id,
                                                     stk::mesh::Part & part)
{
    stk::mesh::Entity new_face = stk::mesh::Entity();
    unsigned elem_num_faces = bulk.num_faces(elem);
    const stk::mesh::Entity * elem_faces = bulk.begin_faces(elem);
    const stk::mesh::ConnectivityOrdinal * elem_ord_it = bulk.begin_face_ordinals(elem);
    for (unsigned i=0 ; i<elem_num_faces ; ++i) {
        if (elem_ord_it[i] == static_cast<unsigned>(side_ordinal)) {
            new_face = elem_faces[i];
            break;
        }
    }
    if (!bulk.is_valid(new_face)) {
        new_face = stk::mesh::declare_element_side(bulk, new_face_global_id, elem, side_ordinal, &part);
    } else {
        stk::mesh::PartVector add_parts(1, &part);
        bulk.change_entity_parts(new_face, add_parts );
    }
    return new_face;
}

void create_shell_status(const std::vector<stk::topology> & elements_touching_surface, stk::topology original_element_topology, std::vector<ShellStatus> & element_shell_status) {
    element_shell_status.resize(elements_touching_surface.size(),NO_SHELLS);
    const bool original_element_is_shell = original_element_topology.is_shell();
    bool got_shells = original_element_is_shell;
    if (!original_element_is_shell) {
        for (size_t i=0 ; i<elements_touching_surface.size() ; ++i) {
            if (elements_touching_surface[i].is_shell()) {
                got_shells = true;
                break;
            }
        }
    }
    if (got_shells) {
        for (size_t i=0 ; i<elements_touching_surface.size() ; ++i) {
            if (original_element_is_shell == elements_touching_surface[i].is_shell() ) {
                element_shell_status[i] = YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS;
            }
            else {
                element_shell_status[i] = YES_SHELLS_ONE_SHELL_ONE_SOLID;
            }
        }
    }
}



void connect_face_to_other_elements(stk::mesh::BulkData & bulk,
                                      stk::mesh::Entity face,
                                      stk::mesh::Entity elem_with_face,
                                      int elem_with_face_side_ordinal
)
{
    stk::topology elem_topology = bulk.bucket(elem_with_face).topology();
    stk::topology side_topology = elem_topology.face_topology(elem_with_face_side_ordinal);
    int num_side_nodes = side_topology.num_nodes();
    std::vector<stk::mesh::Entity> side_nodes(num_side_nodes);
    elem_topology.face_nodes(bulk.begin_nodes(elem_with_face),elem_with_face_side_ordinal,&side_nodes[0]);
    std::vector<stk::mesh::Entity> common_elements;
    stk::mesh::impl::find_elements_these_nodes_have_in_common(bulk,num_side_nodes,&side_nodes[0],common_elements);

    std::vector<stk::topology> element_topology_touching_surface_vector(common_elements.size());
    for (size_t i=0 ; i<common_elements.size() ; ++i) {
        element_topology_touching_surface_vector[i] = bulk.bucket(common_elements[i]).topology();
    }
    std::vector<ShellStatus> element_shell_status;
    create_shell_status(element_topology_touching_surface_vector, elem_topology, element_shell_status);
    for (size_t count=0 ; count<common_elements.size() ; ++count) {
        if (common_elements[count] != elem_with_face) {
            stk::mesh::Entity other_elem = common_elements[count];
            stk::topology other_elem_topology = element_topology_touching_surface_vector[count];
            for (unsigned other_elem_side = 0; other_elem_side < other_elem_topology.num_faces() ; ++other_elem_side) {
                stk::topology other_elem_side_topology = other_elem_topology.face_topology(other_elem_side);
                std::vector<stk::mesh::Entity> other_elem_side_nodes(other_elem_side_topology.num_nodes());
                other_elem_topology.face_nodes(bulk.begin_nodes(other_elem),other_elem_side,&other_elem_side_nodes[0]);
                if (should_face_be_connected_to_element_side(side_nodes,other_elem_side_nodes,other_elem_side_topology,element_shell_status[count])) {
                    stk::mesh::declare_element_side(bulk, other_elem, face, other_elem_side);
                    break;
                }
            }
        }
    }
}

} } } // namespace impl io stk
























