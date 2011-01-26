/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_rebalance_utils/RebalanceUtils.hpp>
#include <stk_mesh/base/GetEntities.hpp>

//----------------------------------------------------------------------

bool stk::rebalance::verify_dependent_ownership( const stk::mesh::EntityRank & parent_rank,
                                                 stk::mesh::EntityVector & entities, 
                                                 stk::mesh::DefaultFEM & fem )
{
  bool is_with_elem = true;
  for( size_t i = 0; i < entities.size(); ++i )
  {
    is_with_elem = false;

    stk::mesh::Entity * entity = entities[i];
    unsigned owner_proc = entity->owner_rank();
    const stk::mesh::PairIterRelation rel = entity->relations( parent_rank );
    const unsigned num_elems = rel.size();

    for ( unsigned j = 0 ; j < num_elems ; ++j ) 
    {
      stk::mesh::Entity & elem = * rel[j].entity();
      if( owner_proc == elem.owner_rank() )
      {
        is_with_elem = true;
        break;
      }
    }
    if( !is_with_elem )
      return false;
  }

  return is_with_elem;
}
