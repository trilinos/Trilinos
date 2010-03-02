/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_SetOwners_hpp
#define stk_mesh_SetOwners_hpp

#include <set>
#include <map>

#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

typedef std::less<unsigned> LowestRankSharingProcOwns;
typedef std::greater<unsigned> HighestRankSharingProcOwns;

/** Sets the owner for shared entities according to the template parameter OwnershipRule.
 * OwnershipRule is used as a the comparison operator in a std::set.
 * The default behavior of stk::mesh is to give ownership to the highest-rank sharing proc.
*/
template<class OwnershipRule>
void set_owners(BulkData& mesh_bulk_data)
{
  const std::vector<EntityProc>& shared_entities = mesh_bulk_data.shared_entities();

  typedef std::map<Entity*,std::set<unsigned,OwnershipRule> > entity_set_map;

  entity_set_map entities;

  unsigned local_proc = mesh_bulk_data.parallel_rank();

  for(size_t i=0; i<shared_entities.size(); ++i) {
    Entity* entity = shared_entities[i].first;
    const unsigned& proc = shared_entities[i].second;

    entities[entity].insert(proc);
  }

  std::vector<EntityProc> entity_new_owners;

  typename entity_set_map::iterator
    iter = entities.begin(), iend = entities.end();

  for(; iter!=iend; ++iter) {
    Entity* const entity = iter->first;
    const unsigned& new_owner_proc = *(iter->second.begin());

    if (entity->owner_rank() == local_proc) {
      entity_new_owners.push_back(std::make_pair( entity, new_owner_proc ) );
    }
  }

  mesh_bulk_data.modification_begin();

  mesh_bulk_data.change_entity_owner( entity_new_owners );

  mesh_bulk_data.modification_end();
}

}//namespace mesh
}//namespace stk

#endif // stk_mesh_SetOwner_hpp

