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

#include <stk_mesh/baseImpl/CommImplUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

bool ghost_id_is_found_in_comm_data(const PairIterEntityComm& comm_data,
                                    int entity_owner,
                                    int ghost_id)
{
  bool found_ghost_id = false;
  for (size_t i = 0; i < comm_data.size(); ++i) {
    if ((comm_data[i].ghost_id == static_cast<unsigned>(ghost_id)) &&
        (comm_data[i].proc == entity_owner)) {
          found_ghost_id = true;
          break;
    }
  }   
  return found_ghost_id;
}   

bool all_ghost_ids_are_found_in_comm_data(const PairIterEntityComm& comm_data,
                                          int entity_owner,
                                          const std::vector<int>& recvd_ghost_ids)
{
  bool found_all_ghost_ids = true;
  for (int ghost_id : recvd_ghost_ids) {
    if (!ghost_id_is_found_in_comm_data(comm_data, entity_owner, ghost_id)) {
      found_all_ghost_ids = false;
      break;
    }   
  }   
  return found_all_ghost_ids;
}

void comm_shared_procs(PairIterEntityComm commInfo,
                       std::vector<int>& sharingProcs)
{
  sharingProcs.clear();
  for(; !commInfo.empty(); ++commInfo) {
    if (commInfo->ghost_id == BulkData::SHARED) {
      sharingProcs.push_back(commInfo->proc);
    }
    else {
      break;
    }
  }
}

void fill_sorted_procs(PairIterEntityComm ec, std::vector<int>& procs)
{
  procs.clear();
  procs.reserve(ec.size());
  int prevProc = -1;
  bool isSorted = true;
  for(; !ec.empty(); ++ec) {
    isSorted = isSorted && (prevProc < ec->proc);
    procs.push_back( ec->proc );
    prevProc = ec->proc;
  }
  if (!isSorted) {
    stk::util::sort_and_unique(procs);
  }
}

unsigned fill_procs_shared_then_ghosted(PairIterEntityComm ec, std::vector<int>& procs)
{
//ok this is a bit tricky so we actually need a comment...
//ec contains pairs of 'ghost_id,proc' indicating which procs an entity is
//communicated with, and why. (And remember, 'shared' is one of the ghost_ids...)
//
//There are lots of possible combinations, such as:
// - an entity can be shared and ghosted to the same proc.
// - an entity can be ghosted twice to the same proc (two different ghost-ids)
//
//One thing we can count on is that the pairs are sorted by ghost_id, meaning all
//the shared come first, and then the rest of the ghosts.
//
//We want the procs list to have shared-procs first, sorted and unique.
//We want the ghost procs to come after shared, and we only want ghost procs
//that are NOT also shared.
//I.e., we want procs overall to be unique, and in two sorted sections
//for shared-procs and ghosted-procs.
//
//See unit tests in stk_mesh/UnitTestCommImplUtils.cpp.
//
  procs.clear();
  procs.reserve(ec.size());
  
  while(!ec.empty() && ec->ghost_id == BulkData::SHARED) {
    procs.push_back( ec->proc );
    ++ec;
  }

  const unsigned numShared = procs.size();

  //At this point we can safely assume the 'numShared' procs
  //are sorted and unique.
  //So now we will add any ghost procs that are not also shared.
 
  while(!ec.empty()) {
    if (std::find(procs.begin(), procs.end(), ec->proc) == procs.end()) {
      procs.push_back(ec->proc);
    }
    ++ec;
  }

  const unsigned numGhosts = procs.size() - numShared;
  if (numGhosts >= 2) {
    std::vector<int>::iterator beginGhosts = procs.begin()+numShared;
    std::sort(beginGhosts, procs.end());
  }
 
  return numShared;
}

void fill_ghosting_procs(const PairIterEntityComm& ec, unsigned ghostID, std::vector<int>& procs)
{
  procs.clear();
  const int n = ec.size(); 
  for (int i=0; i<n; ++i) {
    if ((ghostID == ec[i].ghost_id) || (ec[i].ghost_id == (BulkData::SYMM_INFO+ghostID))) {
      procs.push_back( ec[i].proc );
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

