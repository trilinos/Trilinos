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

#ifndef MESH_DIAGNOSTICS_HPP
#define MESH_DIAGNOSTICS_HPP

#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <string>
#include <vector>
#include <map>

namespace stk { namespace mesh {

typedef std::map<stk::mesh::EntityId, std::vector<std::pair<stk::mesh::EntityId, int>>> SplitCoincidentInfo;

class BulkData;

stk::mesh::SplitCoincidentInfo get_split_coincident_elements(stk::mesh::BulkData& bulkData);
std::vector<std::string> get_messages_for_split_coincident_elements(const stk::mesh::BulkData& bulkData, const stk::mesh::SplitCoincidentInfo & splitCoincidentElements);

std::vector<stk::mesh::EntityKeyProc> get_non_unique_key_procs(const stk::mesh::BulkData& bulkData);
std::vector<std::string> get_non_unique_key_messages(const stk::mesh::BulkData& bulkData, const std::vector<stk::mesh::EntityKeyProc> &badKeyProcs);

std::vector<stk::mesh::Entity> get_orphaned_owned_sides(const stk::mesh::BulkData& bulkData);
std::vector<stk::mesh::Entity> get_orphaned_sides_with_attached_element_on_different_proc(const stk::mesh::BulkData& bulkData);
std::vector<std::string> get_messages_for_orphaned_owned_sides(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& keys);

std::vector<stk::mesh::Entity> get_solo_sides_without_element_on_different_proc(const stk::mesh::BulkData& bulkData);
std::vector<std::string> get_messages_for_solo_sides(const stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& entities);


void throw_if_any_proc_has_false(MPI_Comm comm, bool is_all_ok_locally);

} }

#endif
