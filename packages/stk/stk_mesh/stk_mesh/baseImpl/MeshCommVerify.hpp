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

#ifndef STK_MESH_BASEIMPL_MESHCOMMVERIFY_HPP
#define STK_MESH_BASEIMPL_MESHCOMMVERIFY_HPP

//----------------------------------------------------------------------

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/EntityCommListInfo.hpp>

#include <vector>
#include <functional>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
class BulkData;

namespace impl {

//----------------------------------------------------------------------
//These functions are not part of the public API of stk-mesh.
//They are intended for use internally in the implementation of
//stk-mesh capabilities.
//----------------------------------------------------------------------
//

void unpack_not_owned_verify_compare_closure_relations( const BulkData & mesh,
                                               Entity           entity,
                                               std::vector<Relation> const& recv_relations,
                                               bool& bad_rel);

void unpack_not_owned_verify_compare_parts(const BulkData &  mesh,
                                           Entity            entity,
                                           PartVector const& recv_parts,
                                           bool&             bad_part);

bool verify_parallel_attributes_for_bucket(const Bucket& bucket,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           std::ostream& error_log);

void pack_owned_verify(const BulkData& mesh,
                       const EntityCommDatabase& commDB,
                       const EntityCommListInfoVector& commList,
                       CommSparse& commSparse);

void unpack_not_owned_verify_compare_comm_info( const BulkData& mesh,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                                CommBuffer&            buf,
                                                Entity                 entity,
                                                EntityKey &            recv_entity_key,
                                                int       &            recv_owner_rank,
                                                unsigned  &            recv_comm_count,
                                                PartVector&    recv_parts,
                                                std::vector<Relation>& recv_relations,
                                                std::vector<int>    &  recv_comm,
                                                bool&                  bad_comm);

void unpack_not_owned_verify_report_errors(const BulkData& mesh,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           Entity entity,
                                           bool bad_key,
                                           bool bad_own,
                                           bool bad_part,
                                           bool bad_rel,
                                           bool bad_comm,
                                           EntityKey            recv_entity_key,
                                           int                  recv_owner_rank,
                                           PartVector const&    recv_parts,
                                           std::vector<Relation> const& recv_relations,
                                           std::vector<int>    const&  recv_comm,
                                           std::ostream & error_log);

bool unpack_not_owned_verify(const BulkData& mesh,
                             const EntityCommListInfoVector& commList,
                      const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                             CommSparse& commSparse,
                             std::ostream& error_log);

bool verify_parallel_attributes(const BulkData& mesh,
                                const EntityCommDatabase& commDB,
                                const EntityCommListInfoVector& commList,
                                const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                std::ostream & error_log );

bool comm_mesh_verify_parallel_consistency(const BulkData& mesh,
                                           const EntityCommDatabase& commDB,
                                           const EntityCommListInfoVector& commList,
                                           const std::function<PairIterEntityComm(Entity)>& getEntityComm,
                                           std::ostream & error_log );

void check_matching_parts_count(unsigned partsCount, int rank, int commSize, MPI_Comm comm);

void check_matching_parts(const PartVector& parts, unsigned partsCount, int rank, int commSize, MPI_Comm comm);

void check_matching_parts_across_procs(const PartVector& parts, MPI_Comm comm);

void check_matching_selectors_and_parts_across_procs(const Selector& selector,
                                                     const PartVector& add_parts,
                                                     const PartVector& remove_parts,
                                                     MPI_Comm comm);

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // STK_MESH_BASEIMPL_MESHCOMMVERIFY_HPP

