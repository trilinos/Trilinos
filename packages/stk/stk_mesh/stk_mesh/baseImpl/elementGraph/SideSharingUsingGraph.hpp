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

#ifndef SIDESHARINGUSINGGRAPH_HPP_
#define SIDESHARINGUSINGGRAPH_HPP_
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
#include <vector>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk { namespace mesh { class BulkData; } }

namespace stk { namespace mesh {

struct SideSharingData
{
    stk::mesh::impl::IdViaSidePair elementAndSide;
    stk::mesh::Entity side;
    int sharingProc;
    std::vector<int> allSharingProcs;
    int owningProc;
    stk::mesh::EntityId chosenSideId;
    stk::mesh::OrdinalVector partOrdinals;
    SideSharingData() : elementAndSide({0,-1}), side(stk::mesh::Entity()), sharingProc(-1), owningProc(-1), chosenSideId(stk::mesh::InvalidEntityId) {}
    SideSharingData(const stk::mesh::impl::IdViaSidePair& sidePair, stk::mesh::Entity sideIn, int sharing_proc, const std::vector<int> &all_sharing_procs, int owning_proc, stk::mesh::EntityId chosen_id)
    : elementAndSide(sidePair), side(sideIn), sharingProc(sharing_proc), allSharingProcs(all_sharing_procs), owningProc(owning_proc), chosenSideId(chosen_id) {}
};

stk::mesh::EntityVector fill_shared_entities_that_need_fixing(const stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank);

void fill_sharing_data(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph &graph, const stk::mesh::EntityVector& sidesThatNeedFixing, std::vector<SideSharingData>& sideSharingDataThisProc, std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides);

void allocate_and_send(stk::CommSparse& comm, const std::vector<SideSharingData>& sideSharingDataThisProc, const std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides);

void unpack_data(stk::CommSparse& comm, int my_proc_id, int num_procs, std::vector<SideSharingData>& sideSharingDataThisProc);

}}



#endif /* SIDESHARINGUSINGGRAPH_HPP_ */
