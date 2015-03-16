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

#ifndef STK_IO_IMPL_HPP
#define STK_IO_IMPL_HPP

#include <vector>
#include <stk_mesh/base/Types.hpp>
#include "stk_topology/topology.hpp"
#include "stk_mesh/base/Entity.hpp"

namespace stk { namespace mesh { class BulkData; } }


namespace stk { namespace io { namespace impl {

stk::mesh::Entity get_or_create_face_at_element_side(stk::mesh::BulkData & bulk,
                                                     stk::mesh::Entity elem,
                                                     int side_ordinal,
                                                     int new_face_global_id,
                                                     stk::mesh::Part & part);

void connect_face_to_other_elements(stk::mesh::BulkData & bulk,
                                    stk::mesh::Entity face,
                                    stk::mesh::Entity elem_with_face,
                                    int elem_with_face_side_ordinal);


enum ShellStatus {
    NO_SHELLS = 25,
    YES_SHELLS_ONE_SHELL_ONE_SOLID = 32,
    YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS = 46
};

void create_shell_status(const std::vector<stk::topology> & elements_touching_surface, stk::topology original_element_topology, std::vector<ShellStatus> & element_shell_status);

template<typename ENTITY_ID>
bool should_face_be_connected_to_element_side(std::vector<ENTITY_ID> & face_nodes,
                                              std::vector<ENTITY_ID> & element_side_nodes,
                                              stk::topology element_side_topology,
                                              ShellStatus  shell_status)
{
    bool should_connect = false;
    const std::pair<bool, unsigned> equiv_result = element_side_topology.equivalent(face_nodes, element_side_nodes);
    const bool nodes_match = equiv_result.first;
    if (nodes_match) {
       if (NO_SHELLS == shell_status) {
           should_connect = true;
       }
       else {
           const unsigned permutation_of_element_side = equiv_result.second;
           const bool element_side_polarity_matches_face_nodes = permutation_of_element_side < element_side_topology.num_positive_permutations();
           if (YES_SHELLS_ONE_SHELL_ONE_SOLID == shell_status) {
               should_connect = !element_side_polarity_matches_face_nodes;
           }
           else { // YES_SHELLS_BOTH_SHELS_OR_BOTH_SOLIDS
               should_connect = element_side_polarity_matches_face_nodes;
           }
       }
    }
    return should_connect;
}

} } } // namespace impl io stk































#endif // STK_IO_IMPL_HPP
