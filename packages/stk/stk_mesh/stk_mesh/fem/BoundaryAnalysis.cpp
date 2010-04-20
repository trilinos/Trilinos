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
  for ( ; itr != entities_closure.end() && (*itr)->entity_rank() == closure_rank; ++itr) {
    // some temporaries for clarity
    std::vector<EntitySideComponent > adjacent_entities;
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
        keeper.inside.entity = &curr_entity;
        keeper.inside.side_ordinal = subcell_identifier;
        keeper.outside.entity = NULL;
        keeper.outside.side_ordinal = 0;
        boundary.push_back(keeper);
        continue;
      }

      // iterate over adjacent entities (our neighbors)
      for (std::vector<EntitySideComponent >::const_iterator adj_itr = adjacent_entities.begin(); adj_itr != adjacent_entities.end(); ++adj_itr) {
        // grab a reference to this neighbor for clarity
        const Entity& neighbor = *(adj_itr->entity);

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
  const CellTopologyData* celltopology = get_cell_topology(entity);
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
      msg << "stk::mesh::get_adjacent_entities( const Entity& entity, unsigned subcell_rank, ... ) subcell_rank is >= celltopology dimension\n";
    }
    else if (bad_id) {
      msg << "stk::mesh::get_adjacent_entities( const Entity& entity, unsigned subcell_rank, unsigned subcell_identifier, ... ) subcell_identifier is >= subcell count\n";
    }

    throw std::runtime_error(msg.str());
  }

  // For the potentially common subcell, get it's nodes and num_nodes
  const unsigned* nodes = celltopology->subcell[subcell_rank][subcell_identifier].node;
  unsigned num_nodes = celltopology->subcell[subcell_rank][subcell_identifier].topology->node_count;

  // Get all the nodal relationships for this entity. We are guaranteed
  // that, if we make it this far, the entity is guaranteed to have
  // some relationship to nodes (we know it is a higher-order entity
  // than Node).
  PairIterRelation relations = entity.relations(Node);

  // Get the node entities that are related to entity
  std::vector<Entity*> node_entities;
  for (unsigned itr = 0; itr < num_nodes; ++itr) {
    node_entities.push_back(relations[nodes[itr]].entity());
  }

  // Given the nodes related to the original entity, find all entities
  // of similar rank that have some relation to one or more of these nodes
  std::vector<Entity*> elements;
  get_entities_through_relations(node_entities, entity.entity_rank(), elements);

  // Make sure to remove the original entity from the list
  bool found = false;
  for (std::vector<Entity*>::iterator itr = elements.begin();
       itr != elements.end(); ++itr) {
    if (*itr == &entity) {
      elements.erase(itr);
      found = true;
      break;
    }
  }
  // The original entity should be related to the nodes of its subcells
  if (! found) {
    throw std::logic_error( "stk::mesh::get_adjacent_entities");
  }

  // Add the local ids, from the POV of the adj entitiy, to the return value
  for (std::vector<Entity*>::const_iterator itr = elements.begin();
       itr != elements.end(); ++itr) {
    unsigned local_side_num = element_local_side_id(**itr, node_entities);
    adjacent_entities.push_back(EntitySideComponent(*itr, local_side_num));
  }
}

}
}
