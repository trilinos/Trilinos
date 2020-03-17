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

#ifndef STK_MESH_BASEIMPL_MESHCOMMIMPLUTILS_HPP
#define STK_MESH_BASEIMPL_MESHCOMMIMPLUTILS_HPP

//----------------------------------------------------------------------

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>

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

void pack_induced_memberships_for_entities_less_than_element_rank(
          const BulkData& bulk_data,
          stk::CommSparse& comm,
          const EntityCommListInfoVector & entityCommListInfoVec );

void append_parts_from_sharer_to_owner(const BulkData& bulk,
                                       Entity entity,
                                       OrdinalVector &partOrdinals);

void pack_induced_memberships( const BulkData& bulk_data,
                               stk::CommSparse & comm ,
                               const EntityCommListInfoVector & entityCommListInfoVec );

void pack_part_memberships( const BulkData& meshbulk, stk::CommSparse & comm ,
                            const std::vector<EntityProc> & send_list );

void unpack_induced_parts_from_sharers(OrdinalVector& induced_parts,
                                       const EntityCommInfoVector& entity_comm_info,
                                       stk::CommSparse& comm,
                                       EntityKey expected_key);

void pack_and_send_induced_parts_from_sharers_to_owners(const BulkData& bulkData, stk::CommSparse& comm, EntityCommListInfoVector& entity_comm_list);

void filter_out_unneeded_induced_parts(const BulkData& bulkData,
                                       Entity entity,
                                       const OrdinalVector& induced_parts,
                                       OrdinalVector& remove_parts);
} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // STK_MESH_BASEIMPL_MESHCOMMIMPLUTILS_HPP

