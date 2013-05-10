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

  const std::vector<EntityCommListInfo>& entity_comm = mesh.comm_list();

  for ( size_t i=0; i<entity_comm.size(); ++i) {
    Entity const entity = entity_comm[i].entity;;

    const PairIterEntityComm sharing = mesh.entity_comm_sharing(entity_comm[i].key);

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
