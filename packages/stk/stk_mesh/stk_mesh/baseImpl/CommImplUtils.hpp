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

#ifndef STK_MESH_BASEIMPL_COMMIMPLUTILS_HPP
#define STK_MESH_BASEIMPL_COMMIMPLUTILS_HPP

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------
//These functions are not part of the public API of stk-mesh.
//They are intended for use internally in the implementation of
//stk-mesh capabilities.
//----------------------------------------------------------------------
//

bool ghost_id_is_found_in_comm_data(const PairIterEntityComm& comm_data,
                                    int entity_owner,
                                    int ghost_id);

bool all_ghost_ids_are_found_in_comm_data(const PairIterEntityComm& comm_data,
                                          int entity_owner,
                                          const std::vector<int>& recvd_ghost_ids);

void comm_shared_procs(PairIterEntityComm commInfo,
                       std::vector<int>& sharingProcs);

void fill_sorted_procs(PairIterEntityComm ec, std::vector<int>& procs);

unsigned fill_procs_shared_then_ghosted(PairIterEntityComm ec, std::vector<int>& procs);

void fill_ghosting_procs(const PairIterEntityComm& ec, unsigned ghost_id, std::vector<int>& procs);

inline
bool is_comm_ordered(const PairIterEntityComm& ec)
{
  int n = ec.size();
  for (int i=1; i<n; ++i) {
    if (!(ec[i-1] < ec[i])) {
      return false;
    }
  }
  return true;
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // STK_MESH_BASEIMPL_COMMIMPLUTILS_HPP

