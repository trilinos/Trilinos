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

#ifndef STK_MESH_ENTITYPROCMAPPING_HPP
#define STK_MESH_ENTITYPROCMAPPING_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/util/SortAndUnique.hpp>

#include <vector>
#include <set>

namespace stk {
namespace mesh {

struct EntityAndProcs
{
  EntityAndProcs(const EntityProc& entityProc)
  : entity(entityProc.first),
    proc(entityProc.second)
  {}

  EntityAndProcs(Entity ent, int p)
  : entity(ent),
    proc(p)
  {}

  Entity entity;
  int proc;
  std::vector<int> procs;
};

inline
bool is_valid(Entity entity)
{
  return entity.local_offset() != 0;
}

class EntityProcMapping {
public:
  EntityProcMapping(unsigned sizeOfEntityIndexSpace)
  : entityOffsets(sizeOfEntityIndexSpace, -1),
    entitiesAndProcs()
  {}

  void addEntityProc(Entity entity, int proc)
  {
    int offset = entityOffsets[entity.local_offset()];
    if (offset < 0) {
      entityOffsets[entity.local_offset()] = entitiesAndProcs.size();
      entitiesAndProcs.emplace_back(entity, proc);
    }
    else {
      EntityAndProcs& entityAndProcs = entitiesAndProcs[offset];
      entityAndProcs.entity = entity;
      if (entityAndProcs.proc < 0) {
        stk::util::insert_keep_sorted_and_unique(proc, entityAndProcs.procs);
      }
      else if (entityAndProcs.proc != proc) {
        entityAndProcs.procs.reserve(2);
        entityAndProcs.procs.push_back(proc);
        stk::util::insert_keep_sorted_and_unique(entityAndProcs.proc, entityAndProcs.procs);
        entityAndProcs.proc = -1;
      }
    }
  }

  void addEntityProc(const EntityProc& entityProc)
  {
    addEntityProc(entityProc.first, entityProc.second);
  }

  void eraseEntityProc(Entity entity, int proc)
  {
    int offset = entityOffsets[entity.local_offset()];
    if (offset < 0) {
      return;
    }
    else {
      EntityAndProcs& eap = entitiesAndProcs[offset];
      if (eap.proc == proc) {
        eap.proc = -1;
        eap.entity = Entity();
        entityOffsets[entity.local_offset()] = -1;
      }
      else {
        eap.proc = -1;
        std::vector<int>::iterator iter = std::find(eap.procs.begin(),
                                                    eap.procs.end(),
                                                    proc);
        if (iter != eap.procs.end()) {
          eap.procs.erase(iter);
          if (eap.procs.empty()) {
            eap.entity = Entity();
          }
        }
      }
    }
  }

  void eraseEntityProc(const EntityProc& entityProc)
  {
    eraseEntityProc(entityProc.first, entityProc.second);
  }

  size_t get_num_entities() const { return entitiesAndProcs.size(); }

  bool find(Entity entity, int proc) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset >= 0) {
      const EntityAndProcs& eap = entitiesAndProcs[offset];
      return (eap.entity == entity) &&
         ((eap.proc == proc) ||
          (std::find(eap.procs.begin(), eap.procs.end(), proc) != eap.procs.end())
         );
    }
    return false;
  }

  bool find(const EntityProc& entityProc) const
  {
    return find(entityProc.first, entityProc.second);
  }

  size_t get_num_procs(Entity entity) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset >= 0) {
      return (entitiesAndProcs[offset].entity == entity) &&
             entitiesAndProcs[offset].proc < 0 ?
                 entitiesAndProcs[offset].procs.size() : 1;
    }
    return 0;
  }

  template<typename SetType>
  void fill_set(SetType& entityProcSet)
  {
    entityProcSet.clear();
    for(const EntityAndProcs& entProcs : entitiesAndProcs) {
      if (is_valid(entProcs.entity) && entProcs.proc >= 0) {
        entityProcSet.insert(EntityProc(entProcs.entity, entProcs.proc));
      }
      else if (is_valid(entProcs.entity)) {
        for(int p : entProcs.procs) {
          entityProcSet.insert(EntityProc(entProcs.entity, p));
        }
      }
    }
  }

  template<typename VecType>
  void fill_vec(VecType& entityProcVec)
  {
    entityProcVec.clear();
    size_t lengthEstimate = static_cast<size_t>(std::floor(1.2*entitiesAndProcs.size()));
    entityProcVec.reserve(lengthEstimate);
    for(const EntityAndProcs& entProcs : entitiesAndProcs) {
      if (is_valid(entProcs.entity) && entProcs.proc >= 0) {
        entityProcVec.emplace_back(EntityProc(entProcs.entity, entProcs.proc));
      }
      else if (is_valid(entProcs.entity)) {
        for(int p : entProcs.procs) {
          entityProcVec.emplace_back(EntityProc(entProcs.entity, p));
        }
      }
    }
  }

private:
  std::vector<int> entityOffsets;
  std::vector<EntityAndProcs> entitiesAndProcs;
};

}
}

#endif
