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

#include <stk_mesh/baseImpl/MeshModLogObserver.hpp>
#include <stk_mesh/baseImpl/PrintEntityState.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace stk {
namespace mesh {
namespace impl {

MeshModLogObserver::MeshModLogObserver(const BulkData& mesh, const EntityKey& keyToWatch, std::ostream& ostrm,
                                       PatchType patchType,
                                       unsigned startAtModCycle)
 : ModificationObserver(stk::mesh::STK_INTERNAL),
   m_mesh(mesh),
   m_keys(),
   m_os(ostrm),
   m_patchType(patchType),
   m_startAtModCycle(startAtModCycle)
{
  m_keys.insert(keyToWatch);
}

void MeshModLogObserver::modification_begin_notification()
{
  if ((m_mesh.synchronized_count()+1) < m_startAtModCycle) { return; }

  m_os << "modification_begin mod-cycle=" << m_mesh.synchronized_count() << std::endl;
  update_patch();
  print_patch();
}

void MeshModLogObserver::entity_added(Entity entity)
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  if (match(entity)) {
    m_os << "P" << m_mesh.parallel_rank() << " mod-cycle=" << m_mesh.synchronized_count() << ", declare_entity " << m_mesh.entity_key(entity) << std::endl;
  }
}

void MeshModLogObserver::entity_deleted(Entity entity)
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  if (match(entity)) {
    m_os << "P" << m_mesh.parallel_rank() << " mod-cycle=" << m_mesh.synchronized_count() << ", destroy_entity " << m_mesh.entity_key(entity) << std::endl;
  }
}

void MeshModLogObserver::relation_destroyed(Entity from, Entity to, ConnectivityOrdinal ordinal)
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  if (match(from) || match(to)) {
    m_os << "P" << m_mesh.parallel_rank() << " mod-cycle=" << m_mesh.synchronized_count() << ", destroy_relation "
         << m_mesh.entity_key(from) << " -> " << m_mesh.entity_key(to) << " ordinal=" << ordinal << std::endl;
  }
}

void MeshModLogObserver::relation_declared(Entity from, Entity to, ConnectivityOrdinal ordinal)
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  if (match(from) || match(to)) {
    m_os << "P" << m_mesh.parallel_rank() << " mod-cycle=" << m_mesh.synchronized_count() << ", declare_relation "
         << m_mesh.entity_key(from) << " -> " << m_mesh.entity_key(to) << " ordinal=" << ordinal << std::endl;
  }
}

void MeshModLogObserver::started_modification_end_notification()
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  m_os << "start modification_end mod-cycle=" << m_mesh.synchronized_count() << std::endl;
  update_patch();
  print_patch();
}

void MeshModLogObserver::finished_modification_end_notification()
{
  if (m_mesh.synchronized_count() < m_startAtModCycle) { return; }

  m_os << "finish modification_end mod-cycle=" << m_mesh.synchronized_count() << std::endl;
  update_patch();
  print_patch();
}

void MeshModLogObserver::add_key_to_watch(const EntityKey& keyToWatch)
{
  m_keys.insert(keyToWatch);
}

bool MeshModLogObserver::match(const Entity entity) const
{
  EntityKey key = m_mesh.entity_key(entity);
  return match(key);
}

bool MeshModLogObserver::match(const EntityKey& key) const
{
  return m_mesh.synchronized_count() >= m_startAtModCycle && m_patchKeys.find(key) != m_patchKeys.end();
}

void MeshModLogObserver::update_patch()
{
  m_patchKeys.clear();
  m_patchKeys.insert(m_keys.begin(), m_keys.end());
  m_entityPatch.clear();

  for(const EntityKey& key : m_keys) {
    Entity entity = m_mesh.get_entity(key);
    if (!m_mesh.is_valid(entity)) { continue; }

    m_entityPatch.push_back(entity);
  }

  if (m_patchType == NO_PATCH) { return; }


  for(const EntityKey& key : m_keys) {
    Entity entity = m_mesh.get_entity(key);
    if (!m_mesh.is_valid(entity)) { continue; }

    StoreEntity storeEntity(m_mesh);
  
    if (m_patchType == UPWARD_PATCH) {
      VisitUpwardClosure(m_mesh, entity, storeEntity);
    }
    else if (m_patchType == UP_DOWN_PATCH) {
      VisitUpwardClosure(m_mesh, entity, storeEntity);
      EntityVector upward;
      storeEntity.store_visited_entities_in_vec(upward);
      VisitClosure(m_mesh, upward.begin(), upward.end(), storeEntity);
    }
    else {
      VisitAuraClosure(m_mesh, entity, storeEntity);
    }
  
    storeEntity.store_visited_entities_in_vec(m_entityPatch);
  }

  stk::util::sort_and_unique(m_entityPatch, EntityLess(m_mesh));

  for(Entity ent : m_entityPatch) {
    m_patchKeys.insert(m_mesh.entity_key(ent));
  }
}

void MeshModLogObserver::print_patch()
{
  static const std::vector<EntityRank> reverseRanks = {stk::topology::ELEM_RANK, stk::topology::FACE_RANK, stk::topology::EDGE_RANK, stk::topology::NODE_RANK};
  static const std::vector<std::string> revRankStrings = {"elems", "faces", "edges", "nodes"};

  for(unsigned i=0; i<m_entityPatch.size(); ++i) {
    unsigned reverseIdx = m_entityPatch.size() - i - 1;
    Entity ent = m_entityPatch[reverseIdx];
    const bool shared = m_mesh.bucket(ent).shared();
    const bool inAura = m_mesh.bucket(ent).in_aura();
    const bool customRecvGhost = m_mesh.in_receive_custom_ghost(m_mesh.entity_key(ent));
    m_os << "P" << m_mesh.parallel_rank() << " " << m_mesh.entity_key(ent)
         << " {P" << m_mesh.parallel_owner_rank(ent) << (shared?",shrd":"") << (inAura?",aura":"") << (customRecvGhost?",CRG":"") << m_mesh.state(ent)<<"}, ";
    for(unsigned r=0; r<reverseRanks.size(); ++r) {
      const EntityRank rank = reverseRanks[r];
      if (rank == m_mesh.entity_rank(ent)) { continue; }

      const unsigned numConn = m_mesh.num_connectivity(ent, rank);
      const Entity* conn = m_mesh.begin(ent, rank);
      m_os << revRankStrings[r] << "{";
      for(unsigned n=0; n<numConn; ++n) {
        m_os << m_mesh.identifier(conn[n]);
        if (n < (numConn-1)) { m_os << " "; }
      }
      m_os << "} ";
    }
    m_os << std::endl;
  }
}

} //namespace impl
} //namespace mesh
} //namespace stk

