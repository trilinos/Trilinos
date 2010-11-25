/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <set>
#include <algorithm>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

namespace stk {
namespace mesh {

void boundary_analysis(const BulkData& bulk_data,
                       const std::vector< Entity *> & entities_closure,
                       EntityRank closure_rank,
                       EntitySideVector& boundary)
{
  const Selector locally_used = bulk_data.mesh_meta_data().locally_owned_part()
                              | bulk_data.mesh_meta_data().globally_shared_part();

  // find an iterator that points to the last item in the closure that is of a
  // lower-order than the closure_rank
  std::vector<Entity*>::const_iterator itr =
    std::lower_bound(entities_closure.begin(),
                     entities_closure.end(),
                     EntityKey(closure_rank, 0),
                     EntityLess());

  // iterate over all the entities in the closure up to the iterator we computed above
  for ( ; itr != entities_closure.end() && (*itr)->entity_rank() == closure_rank; ++itr) {
    // some temporaries for clarity
    std::vector<EntitySideComponent > adjacent_entities;
    Entity& curr_entity = **itr;
    const CellTopologyData* celltopology = fem::get_cell_topology(curr_entity).getTopologyData();
    if (celltopology == NULL) {
      continue;
    }

    unsigned subcell_rank = closure_rank - 1;
    PairIterRelation relations = curr_entity.relations(fem::NODE_RANK);

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
        keeper.inside.entity = &curr_entity;
        keeper.inside.side_ordinal = subcell_identifier;
        keeper.outside.entity = NULL;
        keeper.outside.side_ordinal = 0;
        boundary.push_back(keeper);
        continue;
      }

      // iterate over adjacent entities (our neighbors)
      for (std::vector<EntitySideComponent >::const_iterator
           adj_itr = adjacent_entities.begin();
           adj_itr != adjacent_entities.end(); ++adj_itr) {
        // grab a reference to this neighbor for clarity
        const Entity& neighbor = *(adj_itr->entity);

        // see if this neighbor is in the closure, if so, not a keeper
        bool neighbor_is_in_closure =
          std::binary_search(entities_closure.begin(),
                             entities_closure.end(),
                             neighbor,
                             EntityLess());

        if (neighbor_is_in_closure) {
          continue;
        }

        // if neighbor or curr_entity is locally-used, add it to keeper
        if ( locally_used( neighbor.bucket()) || locally_used( curr_entity.bucket() ) ) {
          EntitySide keeper;
          keeper.inside.entity = &curr_entity;
          keeper.inside.side_ordinal = subcell_identifier;
          keeper.outside = *adj_itr;

          boundary.push_back(keeper);
        }
      }
    }
  }
}

void get_adjacent_entities( const Entity & entity ,
                            unsigned subcell_rank ,
                            unsigned subcell_identifier ,
                            std::vector< EntitySideComponent> & adjacent_entities )
{
  adjacent_entities.clear();

  // get cell topology
  const CellTopologyData* celltopology = fem::get_cell_topology(entity).getTopologyData();
  if (celltopology == NULL) {
    return;
  }

  // valid ranks fall within the dimension of the cell topology
  bool bad_rank = subcell_rank >= celltopology->dimension;

  // local id should be < number of entities of the desired type
  // (if you have 4 edges, their ids should be 0-3)
  bool bad_id = false;
  if (!bad_rank) {
    bad_id = subcell_identifier >= celltopology->subcell_count[subcell_rank];
  }

  if (bad_rank || bad_id) {
    std::ostringstream msg;
    //parallel consisent throw
    if (bad_rank) {
      msg << "get_adjacent_entities( const Entity& entity, unsigned subcell_rank, ... ) subcell_rank is >= celltopology dimension\n";
    }
    else if (bad_id) {
      msg << "get_adjacent_entities( const Entity& entity, unsigned subcell_rank, unsigned subcell_identifier, ... ) subcell_identifier is >= subcell count\n";
    }

    throw std::runtime_error(msg.str());
  }

  // For the potentially common subcell, get it's nodes and num_nodes
  const unsigned* side_node_local_ids =
    celltopology->subcell[subcell_rank][subcell_identifier].node;

  const CellTopologyData * side_topology =
    celltopology->subcell[subcell_rank][subcell_identifier].topology;
  int num_nodes_in_side = side_topology->node_count;

  // Get all the nodal relationships for this entity. We are guaranteed
  // that, if we make it this far, the entity is guaranteed to have
  // some relationship to nodes (we know it is a higher-order entity
  // than Node).

  // Get the node entities for the nodes that make up the side. We put these
  // in in reverse order so that it has the correct orientation with respect
  // the potential adjacent entities we are evaluating. The relations are
  // right-hand-rule ordered for the owning entity, but we need something
  // that's compatible w/ the adjacent entities.
  std::vector<Entity*> side_node_entities;
  side_node_entities.reserve(num_nodes_in_side);
  PairIterRelation irel = entity.relations(fem::NODE_RANK);
  for (int itr = num_nodes_in_side; itr > 0; ) {
    --itr;
    side_node_entities.push_back(irel[side_node_local_ids[itr]].entity());
  }

  // Get the node entities for the nodes that make up the entity
  std::vector<Entity*> entity_nodes;
  entity_nodes.reserve(irel.size());
  for ( ; !irel.empty(); ++irel ) {
    entity_nodes.push_back(irel->entity());
  }
  std::sort(entity_nodes.begin(), entity_nodes.end());

  // Given the nodes related to the side, find all entities
  // with the same rank that have a relation to all of these nodes
  std::vector<Entity*> elements;
  elements.reserve(8); //big enough that resizing should be rare
  get_entities_through_relations(side_node_entities,
                                 entity.entity_rank(),
                                 elements);

  // Make sure to remove the all superimposed entities from the list
  unsigned num_nodes_in_orig_entity = entity_nodes.size();
  std::vector<Entity*> current_nodes;
  current_nodes.resize(num_nodes_in_orig_entity);
  std::vector<Entity*>::iterator itr = elements.begin();
  while ( itr != elements.end() ) {
    Entity * current_entity = *itr;
    PairIterRelation relations = current_entity->relations(fem::NODE_RANK);

    if (current_entity == &entity) {
      // We do not want to be adjacent to ourself
      itr = elements.erase(itr);
    }
    else if (relations.size() != num_nodes_in_orig_entity) {
      // current_entity has a different number of nodes than entity, they
      // cannot be superimposed
      ++itr;
    }
    else {
      for (unsigned i = 0; relations.first != relations.second;
           ++relations.first, ++i ) {
        current_nodes[i] = relations.first->entity();
      }
      std::sort(current_nodes.begin(), current_nodes.end());

      bool entities_are_superimposed = entity_nodes == current_nodes;
      if (entities_are_superimposed) {
        itr = elements.erase(itr);
      }
      else {
        ++itr;
      }
    }
  }

  // Add the local ids, from the POV of the adj entitiy, to the return value
  for (std::vector<Entity*>::const_iterator eitr = elements.begin();
       eitr != elements.end(); ++eitr) {
    int local_side_num = element_local_side_id(**eitr,
                                               side_topology,
                                               side_node_entities);
    if ( local_side_num != -1) {
      adjacent_entities.push_back(EntitySideComponent(*eitr, local_side_num));
    }
  }
}

}
}
