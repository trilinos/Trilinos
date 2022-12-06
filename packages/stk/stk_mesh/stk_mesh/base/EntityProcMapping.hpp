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
#include <stk_util/util/MCSR.hpp>

#include <vector>
#include <set>

namespace stk {
namespace mesh {

constexpr unsigned s_initialCapacity = 128;
constexpr int s_invalidProc = -1;

class EntityProcMapping
{
public:
  EntityProcMapping(unsigned sizeOfEntityIndexSpace = 1024)
  : entityOffsets(sizeOfEntityIndexSpace, -1),
    entities(),
    entityProcs(0, s_invalidProc)
  {}

  void reset(unsigned sizeOfEntityIndexSpace)
  {
    std::fill(entityOffsets.begin(), entityOffsets.end(), -1);
    entityOffsets.resize(sizeOfEntityIndexSpace, -1);
    entities.clear();
    entityProcs.clear(s_initialCapacity);
  }
 
  void deallocate()
  {
    std::vector<int> tmpInts;
    entityOffsets.swap(tmpInts);
    EntityVector tmpEntities;
    entities.swap(tmpEntities);
    entityProcs.clear(0);
  }

  void addEntityProc(Entity entity, int proc)
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset < 0) {
      unsigned newOffset = entityProcs.num_rows();
      entityOffsets[entity.local_offset()] = newOffset;
      entities.push_back(entity);
      ThrowAssert(entities.size() == newOffset+1);
      entityProcs.add_row();
      entityProcs.add_item(newOffset, proc);
    }
    else {
      entityProcs.add_item(offset, proc);
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
      entityProcs.remove_item(offset, proc);
    }
  }

  void eraseEntityProc(const EntityProc& entityProc)
  {
    eraseEntityProc(entityProc.first, entityProc.second);
  }

  size_t get_num_entities() const { return entityProcs.num_rows(); }

  bool find(Entity entity, int proc) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset >= 0) {
      const int* beg = entityProcs.begin(offset);
      const int* end = entityProcs.end(offset);
      return std::find(beg, end, proc) != end;
    }
    return false;
  }

  bool find(const EntityProc& entityProc) const
  {
    return find(entityProc.first, entityProc.second);
  }

  bool find(Entity entity) const
  {
    return (entityOffsets[entity.local_offset()] >= 0);
  }

  size_t get_num_procs(Entity entity) const
  {
    const int offset = entityOffsets[entity.local_offset()];
    if (offset >= 0) {
      return entityProcs.size(offset);
    }
    return 0;
  }

  template<class Alg>
  void visit_entity_procs(const Alg& alg)
  {
    ThrowAssert(entities.size() == entityProcs.num_rows());
    for(unsigned i=0; i<entities.size(); ++i) {
      const unsigned numProcs = entityProcs.size(i);
      const int* procs = entityProcs.begin(i);
      for(unsigned p=0; p<numProcs; ++p) {
        alg(entities[i], procs[p]);
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
    entityProcVec.reserve(entityProcs.total_num_items());
    entityProcVec.clear();
    visit_entity_procs([&entityProcVec](Entity ent, int proc){entityProcVec.push_back(EntityProc(ent,proc));});
  }

private:
  std::vector<int> entityOffsets;
  stk::mesh::EntityVector entities;
  stk::util::MCSR<int> entityProcs;
};

}
}

#endif
