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

namespace stk_classic {
namespace mesh {

typedef std::less<unsigned> LowestRankSharingProcOwns;
typedef std::greater<unsigned> HighestRankSharingProcOwns;

/** Sets the owner for shared entities according to the template parameter OwnershipRule.
 * OwnershipRule is used as a the comparison operator in a std::set.
 * The default behavior of stk_classic::mesh is to give ownership to the highest-rank sharing proc.
*/
template<class OwnershipRule>
void set_owners(BulkData& mesh_bulk_data)
{
  typedef std::set<unsigned,OwnershipRule> ProcSet ;

  const unsigned local_proc = mesh_bulk_data.parallel_rank();

  std::vector<EntityProc> entity_new_owners;

  const std::vector<Entity*>& entity_comm = mesh_bulk_data.entity_comm();

  for ( size_t i=0; i<entity_comm.size(); ++i) {
    Entity * const entity = entity_comm[i] ;

    const PairIterEntityComm sharing = entity->sharing();

    if ( ! sharing.empty() && entity->owner_rank() == local_proc ) {
      ProcSet proc_set ;

      proc_set.insert( local_proc );

      for ( size_t j = 0 ; j < sharing.size() ; ++j ) {
        proc_set.insert( sharing[j].proc );
      }

      const unsigned new_owner_proc = *proc_set.begin();

      entity_new_owners.push_back(std::make_pair( entity, new_owner_proc ) );
    }
  }

  mesh_bulk_data.modification_begin();

  mesh_bulk_data.change_entity_owner( entity_new_owners );

  mesh_bulk_data.modification_end();
}

}//namespace mesh
}//namespace stk_classic

#endif // stk_mesh_SetOwner_hpp

