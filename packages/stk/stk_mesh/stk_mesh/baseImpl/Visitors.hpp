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

#ifndef stk_mesh_base_impl_Visitors_hpp
#define stk_mesh_base_impl_Visitors_hpp

//----------------------------------------------------------------------

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, class DESIRED_ENTITY>
void VisitClosureGeneral(
        const BulkData & mesh,
        Entity inputEntity,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    if (desired_entity(inputEntity)) {
      do_this(inputEntity);
      const EntityRank inputEntityRank = mesh.entity_rank(inputEntity);
      for (EntityRank rank = stk::topology::NODE_RANK ; rank < inputEntityRank ; ++rank) {
          unsigned num_entities_of_rank = mesh.num_connectivity(inputEntity,rank);
          if (num_entities_of_rank > 0) {
            const bool dontRecurse = rank == stk::topology::NODE_RANK ||
                                     inputEntityRank <= stk::topology::ELEM_RANK;
            const Entity * entities = mesh.begin(inputEntity,rank);

            for (unsigned i=0 ; i<num_entities_of_rank ; ++i) {
                if (dontRecurse) {
                  if (desired_entity(entities[i])) {
                    do_this(entities[i]);
                  }
                }
                else {
                    VisitClosureGeneral(mesh,entities[i],do_this,desired_entity);
                }
            }
          }
      }
    }
}

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitClosureGeneral(
        const stk::mesh::BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    for (FORWARD_ITERATOR entity_iterator = start ; entity_iterator != finish ; ++entity_iterator)
    {
        VisitClosureGeneral<DO_THIS_FOR_ENTITY_IN_CLOSURE,DESIRED_ENTITY>(mesh,*entity_iterator,do_this,desired_entity);
    }
}

template <typename VECTOR>
struct StoreInVector {
    StoreInVector(VECTOR & vec_in) : ev(vec_in) {}
    void operator()(stk::mesh::Entity entity) {
      ev.push_back(entity);
    }
    VECTOR & ev;
};

template <typename SET>
struct StoreInSet {
    StoreInSet(SET & set_in) : es(set_in) {}
    void operator()(stk::mesh::Entity entity) {
      es.insert(entity);
    }
    SET & es;
};

struct AlwaysVisit {
    bool operator()(Entity entity) { return true; }
};

struct OnlyVisitOnce {
    OnlyVisitOnce(const BulkData& mesh_in)
    : mesh(mesh_in), already_visited(mesh.get_size_of_entity_index_space(), false) {}
    bool operator()(Entity entity) {
        if (mesh.is_valid(entity) && !already_visited[entity.local_offset()]) {
            already_visited[entity.local_offset()] = true;
            return true;
        }
        return false;
    }
    const BulkData& mesh;
    std::vector<bool> already_visited;
};

struct OnlyVisitUnchanged
{
    OnlyVisitUnchanged(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity){
        if (mesh.state(entity) == Unchanged) {
            return true;
        }
        return false;
    }
    BulkData & mesh;
};

struct OnlyVisitLocallyOwnedOnce {
    OnlyVisitLocallyOwnedOnce(const BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity)
    {
        return ovo(entity) && mesh.bucket(entity).owned();
    }
    const BulkData& mesh;
    OnlyVisitOnce ovo;
};

struct OnlyVisitSharedOnce {
    OnlyVisitSharedOnce(const BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity)
    {
        if (ovo(entity) && !mesh.in_shared(mesh.entity_key(entity))) { return true; }
        return false;
    }
    const BulkData & mesh;
    OnlyVisitOnce ovo;
};

struct OnlyVisitGhostsOnce
{
    OnlyVisitGhostsOnce(BulkData & mesh_in) : mesh(mesh_in), ovo(mesh_in) {}
    bool operator()(Entity entity) {
        if (ovo(entity) && mesh.in_receive_ghost(entity)) { return true; }
        return false;
    }
   BulkData & mesh;
   OnlyVisitOnce ovo;
};

template<class DO_THIS_FOR_ENTITY_IN_CLOSURE>
void VisitClosure(
        const stk::mesh::BulkData & mesh,
        stk::mesh::Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}


template<class DO_THIS_FOR_ENTITY_IN_CLOSURE, typename FORWARD_ITERATOR>
void VisitClosure(
        const stk::mesh::BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitClosureGeneral(mesh,start,finish,do_this,ovo);
}


// cyclomatic complexity = 6
template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, class DESIRED_ENTITY>
void VisitUpwardClosureGeneral(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    if (desired_entity(entity_of_interest)) {
        do_this(entity_of_interest);
        if (mesh.is_valid(entity_of_interest)) {
            EntityRank endRank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());
            EntityRank entity_of_interest_rank = mesh.entity_rank(entity_of_interest);
            EntityVector entities_of_rank_up;
            for (EntityRank rank_up = EntityRank(endRank-1) ; rank_up > entity_of_interest_rank ; --rank_up) {
                size_t num_entities_of_rank_up = mesh.num_connectivity(entity_of_interest,rank_up);
                const Entity * entity_up_it = mesh.begin(entity_of_interest,rank_up);

                for (size_t j=0 ; j<num_entities_of_rank_up ; ++j, ++entity_up_it) {
                    VisitUpwardClosureGeneral(mesh,*entity_up_it,do_this,desired_entity);
                }
            }
        }
    }
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitUpwardClosureGeneral(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    for (FORWARD_ITERATOR entity_iterator = start ; entity_iterator != finish ; ++entity_iterator)
    {
        VisitUpwardClosureGeneral<DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE,DESIRED_ENTITY>(mesh,*entity_iterator,do_this,desired_entity);
    }
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE>
void VisitUpwardClosure(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitUpwardClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE, typename FORWARD_ITERATOR>
void VisitUpwardClosure(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_UPWARD_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitUpwardClosureGeneral(mesh,start,finish,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, typename FORWARD_ITERATOR, class DESIRED_ENTITY>
void VisitAuraClosureGeneral(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    std::set<Entity> entity_set;
    StoreInSet<std::set<Entity> > sis(entity_set);
    VisitClosure(mesh,start,finish,sis);
    VisitUpwardClosure(mesh,entity_set.begin(),entity_set.end(),sis);
    VisitClosureGeneral(mesh,entity_set.begin(),entity_set.end(),do_this,desired_entity);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, class DESIRED_ENTITY>
void VisitAuraClosureGeneral(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this,
        DESIRED_ENTITY & desired_entity)
{
    Entity * start = &entity_of_interest;
    Entity * finish = start+1;
    VisitAuraClosureGeneral(mesh,start,finish,do_this,desired_entity);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE, typename FORWARD_ITERATOR>
void VisitAuraClosure(
        const BulkData & mesh,
        const FORWARD_ITERATOR & start,
        const FORWARD_ITERATOR & finish,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitAuraClosureGeneral(mesh,start,finish,do_this,ovo);
}

template<class DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE>
void VisitAuraClosure(
        const BulkData & mesh,
        Entity entity_of_interest,
        DO_THIS_FOR_ENTITY_IN_AURA_CLOSURE & do_this)
{
    OnlyVisitOnce ovo(mesh);
    VisitAuraClosureGeneral(mesh,entity_of_interest,do_this,ovo);
}

struct StoreInEntityProcMapping {
    StoreInEntityProcMapping(BulkData & mesh_in, EntityProcMapping& epm_in)
    :mesh(mesh_in)
    ,myMapping(epm_in)
    {}

    void operator()(Entity entity) {
      myMapping.addEntityProc(entity, proc);
    }

    BulkData & mesh;
    EntityProcMapping& myMapping;
    int proc;
};

struct StoreInEntityProcSet {
    StoreInEntityProcSet(
            BulkData & mesh_in,
            std::set<stk::mesh::EntityProc, stk::mesh::EntityLess> & set_in)
    :mesh(mesh_in)
    ,myset(set_in)
    ,alreadyGhostedToProc(mesh_in.get_size_of_entity_index_space(), -1) { }

    void operator()(Entity entity) {
      if (proc != alreadyGhostedToProc[entity.local_offset()]) {
        alreadyGhostedToProc[entity.local_offset()] = proc;
        myset.insert(stk::mesh::EntityProc(entity,proc));
      }
    }

    BulkData & mesh;
    std::set<stk::mesh::EntityProc , stk::mesh::EntityLess> & myset;
    int proc;
    std::vector<int> alreadyGhostedToProc;
};

struct OnlyGhosts  {
    OnlyGhosts(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity) {
      if (mesh.is_valid(entity)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          const bool isSharedWithProc = mesh.in_shared(entity, proc);
          return !isSharedWithProc;
        }
      }
      return false;
    }
    BulkData & mesh;
    int proc;
};

struct OnlyGhostsEPM  {
    OnlyGhostsEPM(BulkData & mesh_in, const EntityProcMapping& epm_in, const EntityProcMapping& entityShr)
    : mesh(mesh_in), myMapping(epm_in), entitySharing(entityShr) {}
    bool operator()(Entity entity) {
      if (!myMapping.find(entity, proc)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          const bool isSharedWithProc = entitySharing.find(entity, proc);
          return !isSharedWithProc;
        }
      }
      return false;
    }
    BulkData & mesh;
    const EntityProcMapping& myMapping;
    const EntityProcMapping& entitySharing;
    int proc;
};

struct OnlyNewGhosts  {
    OnlyNewGhosts(const BulkData & mesh_in, const Ghosting& ghosting_in) : mesh(mesh_in), ghosting(ghosting_in) {}
    bool operator()(Entity entity) {
      if (mesh.is_valid(entity)) {
        if (proc != mesh.parallel_owner_rank(entity)) {
          if (!mesh.in_ghost(ghosting, entity, proc)) {
            const bool isSharedWithProc = mesh.in_shared(entity, proc);
            return !isSharedWithProc;
          }
        }
      }
      return false;
    }
    const BulkData& mesh;
    const Ghosting& ghosting;
    int proc;
};

struct OnlyRecvGhosts {
  OnlyRecvGhosts(const BulkData& mesh_in, const Ghosting& ghost, const std::vector<bool>& status)
  : mesh(mesh_in), ghosting(ghost), ghostStatus(status) {}
  bool operator()(Entity entity) {
    return mesh.is_valid(entity) && mesh.in_receive_ghost(ghosting, entity) && !ghostStatus[entity.local_offset()];
  }
  const BulkData& mesh;
  const Ghosting& ghosting;
  const std::vector<bool>& ghostStatus;
};

struct VecPushBack {
  VecPushBack(std::vector<Entity>& rcvGhosts, std::vector<bool>& status)
  : recvGhosts(rcvGhosts), ghostStatus(status) {}
  void operator()(Entity entity) {
    recvGhosts.push_back(entity);
    ghostStatus[entity.local_offset()] = true;
  }
  std::vector<Entity>& recvGhosts;
  std::vector<bool>& ghostStatus;
};

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_base_impl_Visitors_hpp

