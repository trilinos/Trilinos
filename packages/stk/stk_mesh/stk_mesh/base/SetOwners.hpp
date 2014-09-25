// Copyright (c) 2013, Sandia Corporation.
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

#ifndef stk_mesh_SetOwners_hpp
#define stk_mesh_SetOwners_hpp

#include <set>
#include <map>

#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

typedef std::less<int> LowestRankSharingProcOwns;
typedef std::greater<int> HighestRankSharingProcOwns;

/** Sets the owner for shared entities according to the template parameter OwnershipRule.
 * OwnershipRule is used as a the comparison operator in a std::set.
 * The default behavior of stk::mesh is to give ownership to the highest-rank sharing proc.
*/
template<class OwnershipRule>
void set_owners(BulkData& mesh)
{
  typedef std::set<int,OwnershipRule> ProcSet ;

  const int local_proc = mesh.parallel_rank();

  std::vector<EntityProc> entity_new_owners;

  const EntityCommListInfoVector& entity_comm = mesh.comm_list();

  for ( size_t i=0; i<entity_comm.size(); ++i) {
    Entity const entity = entity_comm[i].entity;;

    const PairIterEntityComm sharing = mesh.entity_comm_map_sharing(entity_comm[i].key);

    if ( ! sharing.empty() && entity_comm[i].owner == local_proc ) {
      ProcSet proc_set ;

      proc_set.insert( local_proc );

      for ( size_t j = 0 ; j < sharing.size() ; ++j ) {
        proc_set.insert( sharing[j].proc );
      }

      const int new_owner_proc = *proc_set.begin();

      entity_new_owners.push_back(std::make_pair( entity, new_owner_proc ) );
    }
  }

  mesh.modification_begin();

  mesh.change_entity_owner( entity_new_owners );

  mesh.modification_end();
}

}//namespace mesh
}//namespace stk

#endif // stk_mesh_SetOwner_hpp
