/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

namespace stk {
namespace mesh {

void boundary_analysis(const BulkData& bulk_data,
                       const std::vector< Entity *> & entities_closure,
                       unsigned closure_rank,
                       EntitySideVector& boundary)
{
  Part& locally_used_part = bulk_data.mesh_meta_data().locally_used_part();

  // find an iterator that points to the last item in the closure that is of a 
  // lower-order than the closure_rank
  std::vector<Entity*>::const_iterator itr = std::lower_bound(entities_closure.begin(),
                                                              entities_closure.end(),
                                                              EntityKey(closure_rank, 0),
                                                              EntityLess());

  // iterate over all the entities in the closure up to the iterator we computed above
  for ( ; itr != entities_closure.end() && (*itr)->entity_type() == closure_rank; ++itr) {
    // some temporaries for clarity
    std::vector<std::pair<Entity*, unsigned> > adjacent_entities;
    Entity& curr_entity = **itr;
    const CellTopologyData* celltopology = get_cell_topology(curr_entity);
    if (celltopology == NULL) {
      continue;
    }

    unsigned subcell_rank = closure_rank - 1;
    PairIterRelation relations = curr_entity.relations(Node);

    // iterate over the subcells of the current entity
    for (unsigned nitr = 0; nitr < celltopology->subcell_count[subcell_rank]; ++nitr) {
      // find the entities (same rank as subcell) adjacent to this subcell
      unsigned subcell_identifier = nitr;
      get_adjacent_entities(curr_entity,
                            closure_rank - 1,
                            subcell_identifier,
                            adjacent_entities);
      
      // try to figure out if we want to keep ((curr_entity, subcell_identifier), 
      //                                       (adjacent_entity, [k]))
      // if so, push into boundary
      
      // it is a keeper if adjacent entities[k] is not in the entities closure
      // AND if either curr_entity OR adjacent entities[k] is a member of the 
      //               bulk_data.locally_used

      if (adjacent_entities.empty()) {
        EntitySide keeper;
        keeper.first.first = &curr_entity;
        keeper.first.second = subcell_identifier;
        keeper.second.first = NULL;
        keeper.second.second = 0;
        boundary.push_back(keeper);
        continue;
      }

      // iterate over adjacent entities (our neighbors)
      for (std::vector<std::pair<Entity*, unsigned> >::const_iterator adj_itr = adjacent_entities.begin(); adj_itr != adjacent_entities.end(); ++adj_itr) {
        // grab a reference to this neighbor for clarity
        const Entity& neighbor = *(adj_itr->first);

        // see if this neighbor is in the closure, if so, not a keeper
        std::vector<Entity*>::const_iterator search_itr = std::lower_bound(entities_closure.begin(), entities_closure.end(), neighbor, EntityLess());
        EntityEqual eq;
        if (search_itr != entities_closure.end() && eq(**search_itr, neighbor)) {
          continue;
        }

        // if neighbor or curr_entity is locally-used, add it to keeper
        if (has_superset(neighbor.bucket(), locally_used_part) ||
            has_superset(curr_entity.bucket(), locally_used_part)) {
          EntitySide keeper;
          keeper.first.first = &curr_entity;
          keeper.first.second = subcell_identifier;
          keeper.second = *adj_itr;

          boundary.push_back(keeper);
        }
      }
    }
  }
}

}
}
