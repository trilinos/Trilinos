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

#ifndef STK_MeshModLogObserver_HPP
#define STK_MeshModLogObserver_HPP

#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <iostream>

namespace stk {
namespace mesh {

class BulkData;

namespace impl {

enum PatchType : char {
  NO_PATCH,
  UPWARD_PATCH,
  UP_DOWN_PATCH,
  AURA_CLOSURE_PATCH
};

class MeshModLogObserver : public ModificationObserver
{
public:
  MeshModLogObserver(const BulkData& mesh, const EntityKey& keyToWatch, std::ostream& ostrm,
                     PatchType patchType = NO_PATCH,
                     unsigned startAtModCycle = 1);

  void modification_begin_notification() override;

  void entity_added(Entity entity) override;

  void entity_deleted(Entity entity) override;

  void relation_destroyed(Entity from, Entity to, ConnectivityOrdinal ordinal) override;
  
  void relation_declared(Entity from, Entity to, ConnectivityOrdinal ordinal) override;

  void started_modification_end_notification() override;
 
  void finished_modification_end_notification() override;

  virtual ~MeshModLogObserver() {}

  void add_key_to_watch(const EntityKey& keyToWatch);

  stk::mesh::EntityVector get_entity_patch() const { return m_entityPatch; }

  bool match(const Entity entity) const;
  bool match(const EntityKey& key) const;

  unsigned start_mod_cycle() const { return m_startAtModCycle; }

private:
  void update_patch();

  void print_patch();

  const BulkData& m_mesh;
  std::set<EntityKey> m_keys;
  std::ostream& m_os;
  PatchType m_patchType;
  unsigned m_startAtModCycle;
  std::set<EntityKey> m_patchKeys;
  EntityVector m_entityPatch;
};

} //namespace impl
} //namespace mesh
} //namespace stk

#endif

