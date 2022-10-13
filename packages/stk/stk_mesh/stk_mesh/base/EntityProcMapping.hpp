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

  size_t get_num_procs() const
  {
    return proc<0 ? procs.size() : 1;
  }

  bool find_proc(int p) const
  {
    return ((proc == p) ||
            ((proc < 0) && (std::find(procs.begin(), procs.end(), p) != procs.end())));
  }

  void add_proc(int p)
  {
    if (proc < 0) {
      stk::util::insert_keep_sorted_and_unique(p, procs);
    }
    else if (proc != p) {
      procs.reserve(2);
      procs.push_back(proc);
      stk::util::insert_keep_sorted_and_unique(p, procs);
      proc = -1;
    }
  }

  void erase_proc(int p)
  {
    if (proc == p) {
      proc = -1;
    }
    else if (proc < 0) {
      std::vector<int>::iterator iter = std::find(procs.begin(), procs.end(), p);
      if (iter != procs.end()) {
        procs.erase(iter);
        if (procs.empty()) {
          proc = -1;
        }
        else if (procs.size() == 1) {
          proc = procs[0];
          procs.clear();
        }
      }
    }
  }

  Entity entity;
  int proc;
  std::vector<int> procs;
};

class EntityProcMapping {
public:
  EntityProcMapping(unsigned sizeOfEntityIndexSpace = 1024)
  : entityOffsets(sizeOfEntityIndexSpace, -1),
    entitiesAndProcs()
  {}

  void reset(unsigned sizeOfEntityIndexSpace)
  {
    for(int& n : entityOffsets) {
      if (n != -1) {
        n = -1;
      }
    }
    entityOffsets.resize(sizeOfEntityIndexSpace, -1);
    entitiesAndProcs.clear();
  }
 
  void addEntityProc(Entity entity, int proc)
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset < 0) {
      entityOffsets[entity.local_offset()] = entitiesAndProcs.size();
      entitiesAndProcs.emplace_back(entity, proc);
    }
    else {
      entitiesAndProcs[offset].add_proc(proc);
    }
  }

  void addEntityProc(const EntityProc& entityProc)
  {
    addEntityProc(entityProc.first, entityProc.second);
  }

  void eraseEntityProc(Entity entity, int proc)
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset < 0) {
      return;
    }
    else {
      entitiesAndProcs[offset].erase_proc(proc);
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
    return (offset >= 0) ? entitiesAndProcs[offset].find_proc(proc) : false;
  }

  bool find(const EntityProc& entityProc) const
  {
    return find(entityProc.first, entityProc.second);
  }

  bool find(Entity entity) const
  {
    return (entityOffsets[entity.local_offset()] >= 0);
  }

  EntityAndProcs* find_entity_procs(Entity entity)
  {
    const int offset = entityOffsets[entity.local_offset()];
    return  offset >= 0 ? &entitiesAndProcs[offset] : nullptr;
  }

  const EntityAndProcs* find_entity_procs(Entity entity) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    return  offset >= 0 ? &entitiesAndProcs[offset] : nullptr;
  }

  size_t get_num_procs(Entity entity) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset >= 0 && entitiesAndProcs[offset].entity == entity) {
      return entitiesAndProcs[offset].get_num_procs();
    }
    return 0;
  }

  template<class Alg>
  void visit_entity_procs(const Alg& alg)
  {
    for(const EntityAndProcs& entProcs : entitiesAndProcs) {
      if (entProcs.proc >= 0) {
        alg(entProcs.entity, entProcs.proc);
      }
      else {
        for(int p : entProcs.procs) {
          alg(entProcs.entity, p);
        }
      }
    }
  }

  template<typename SetType>
  void fill_set(SetType& entityProcSet)
  {
    entityProcSet.clear();
    visit_entity_procs([&entityProcSet](Entity ent, int proc){entityProcSet.insert(EntityProc(ent,proc));});
  }

  template<typename VecType>
  void fill_vec(VecType& entityProcVec)
  {
    size_t lengthEstimate = static_cast<size_t>(std::floor(1.2*entitiesAndProcs.size()));
    entityProcVec.reserve(lengthEstimate);
    entityProcVec.clear();
    visit_entity_procs([&entityProcVec](Entity ent, int proc){entityProcVec.push_back(EntityProc(ent,proc));});
  }

private:
  std::vector<int> entityOffsets;
  std::vector<EntityAndProcs> entitiesAndProcs;
};

}
}

#endif
