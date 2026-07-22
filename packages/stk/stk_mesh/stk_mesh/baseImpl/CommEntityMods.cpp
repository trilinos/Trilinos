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

#include <stk_mesh/baseImpl/CommEntityMods.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk {
namespace mesh {

namespace impl {

CommEntityMods::CommEntityMods(const BulkData& bulkData,
                               const EntityCommDatabase& commDB,
                               const EntityCommListInfoVector& commList)
 : m_bulkData(bulkData),
   m_commSparse(bulkData.parallel()),
   m_commList(commList),
   m_commDB(commDB)
{
}

CommEntityMods::~CommEntityMods()
{
}

void CommEntityMods::communicate(PackOption packOption)
{
  m_packOption = packOption;

  pack();

  m_commSparse.communicate();

  unpack();
}

void CommEntityMods::pack()
{
  pack_entity_mods();

  m_commSparse.allocate_buffers();

  if (need_to_send()) {
    pack_entity_mods();
  }
}

bool CommEntityMods::need_to_send()
{
  bool needToSend = false;

  for (int procNumber=0; procNumber < m_commSparse.parallel_size(); ++procNumber) {
    if (m_commSparse.send_buffer(procNumber).capacity() > 0) {
      needToSend = true;
      break;
    }
  }

  return needToSend;
}

void CommEntityMods::pack_entity_mods()
{
  const bool packShared = m_packOption == PACK_SHARED || m_packOption == PACK_ALL;
  const bool packGhosted = m_packOption == PACK_GHOSTED || m_packOption == PACK_ALL;

  for ( EntityCommListInfoVector::const_iterator
        i = m_commList.begin() ; i != m_commList.end() ; ++i ) { 
    if (i->entity_comm != -1) {
      Entity entity = i->entity;
      EntityState status = m_bulkData.is_valid(entity) ? m_bulkData.state(entity) : Deleted;

      if ( status == Modified || status == Deleted ) { 
        int owned_closure_int = m_bulkData.owned_closure(entity) ? 1 : 0;

        for ( PairIterEntityComm ec = m_commDB.comm(i->entity_comm); ! ec.empty() ; ++ec )
        {   
          if ( ( packGhosted && ec->ghost_id > BulkData::SHARED ) || ( packShared && ec->ghost_id == BulkData::SHARED ) )
          {
            m_commSparse.send_buffer( ec->proc )
                .pack<EntityKey>( i->key )
                .pack<unsigned>( ec->ghost_id )
                .pack<EntityState>( status )
                .pack<int>(owned_closure_int);
  
            const bool promotingGhostToShared =
              packGhosted && owned_closure_int==1 && !m_bulkData.bucket(entity).owned();
            if (promotingGhostToShared) {
              m_commSparse.send_buffer(m_commSparse.parallel_rank())
                  .pack<EntityKey>( i->key )
                  .pack<unsigned>( ec->ghost_id )
                  .pack<EntityState>( status )
                  .pack<int>(owned_closure_int);
            }
          }   
        }   
      }    
    }        
  }
}

void make_unique(std::vector<EntityParallelState>& pllStates)
{
  unsigned idx = 0;
  for(unsigned i=1; i<pllStates.size(); ++i) {
    if (pllStates[idx].comm_info.key == pllStates[i].comm_info.key &&
        pllStates[idx].from_proc == pllStates[i].from_proc)
    {
      pllStates[idx].state = std::max(pllStates[idx].state, pllStates[i].state);
      pllStates[idx].remote_owned_closure = (pllStates[idx].remote_owned_closure || pllStates[i].remote_owned_closure);
      pllStates[i].comm_info.key = EntityKey();
    }
    else {
      idx = i;
    }
  }

  pllStates.erase(std::remove_if(pllStates.begin(), pllStates.end(), [&](const EntityParallelState& pllState) { return pllState.comm_info.key == EntityKey(); }),
                  pllStates.end());
}

void CommEntityMods::unpack()
{
  for ( int procNumber = 0 ; procNumber < m_commSparse.parallel_size() ; ++procNumber ) { 
    CommBuffer & buf = m_commSparse.recv_buffer( procNumber );
    EntityKey key; 
    EntityState state;
    unsigned ghostId;
    int remote_owned_closure_int;
    bool remote_owned_closure;

    while ( buf.remaining() ) { 

      buf.unpack<EntityKey>( key )
          .unpack<unsigned>( ghostId )
          .unpack<EntityState>( state )
          .unpack<int>( remote_owned_closure_int);
      remote_owned_closure = ((remote_owned_closure_int==1)?true:false);

      EntityCommListInfoVector::const_iterator iter = std::lower_bound(m_commList.begin(), m_commList.end(), key);
      if (iter != m_commList.end() && iter->key == key) {
        EntityCommListInfo info = *iter;
        if (ghostId == BulkData::SHARED) {
          m_sharedMods.emplace_back(EntityParallelState{procNumber, state, info, remote_owned_closure});
        }
        else {
          int remoteProc = (procNumber == m_bulkData.parallel_rank()) ? m_bulkData.parallel_owner_rank(info.entity) : procNumber;
          m_ghostedMods.emplace_back(EntityParallelState{remoteProc, state, info, remote_owned_closure});
        }
      }
    }    
  }

  std::sort(m_sharedMods.begin(), m_sharedMods.end());
  std::sort(m_ghostedMods.begin(), m_ghostedMods.end());
  make_unique(m_ghostedMods);
}

}}} // end namepsace stk mesh impl

