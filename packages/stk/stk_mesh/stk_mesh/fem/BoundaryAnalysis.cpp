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
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Part.hpp>

namespace stk {
namespace mesh {

namespace {

void filter_superimposed_entities(const Entity & entity, EntityVector & entities)
{
  // Get the node entities for the nodes that make up the entity, we'll
  // use this to check for superimposed entities
  PairIterRelation irel = entity.relations(fem::NODE_RANK);
  EntityVector entity_nodes;
  entity_nodes.reserve(irel.size());
  for ( ; !irel.empty(); ++irel ) {
    entity_nodes.push_back(irel->entity());
  }
  std::sort(entity_nodes.begin(), entity_nodes.end());

  // Make sure to remove the all superimposed entities from the list
  unsigned num_nodes_in_orig_entity = entity_nodes.size();
  EntityVector current_nodes;
  current_nodes.resize(num_nodes_in_orig_entity);
  EntityVector::iterator itr = entities.begin();
  while ( itr != entities.end() ) {
    Entity * current_entity = *itr;
    PairIterRelation relations = current_entity->relations(fem::NODE_RANK);

    if (current_entity == &entity) {
      // Superimposed with self by definition
      itr = entities.erase(itr);
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
        itr = entities.erase(itr);
      }
      else {
        ++itr;
      }
    }
  }
}

} // unnamed namespace

void boundary_analysis(const BulkData& bulk_data,
                       const EntityVector & entities_closure,
                       EntityRank closure_rank,
                       EntitySideVector& boundary)
{
  const Selector locally_used = bulk_data.mesh_meta_data().locally_owned_part()
                              | bulk_data.mesh_meta_data().globally_shared_part();

  // find an iterator that points to the last item in the closure that is of a
  // lower-order than the closure_rank
  EntityVector::const_iterator itr =
    std::lower_bound(entities_closure.begin(),
                     entities_closure.end(),
                     EntityKey(closure_rank, 0),
                     EntityLess());

  // iterate over all the entities in the closure up to the iterator we computed above
  for ( ; itr != entities_closure.end() && (*itr)->entity_rank() == closure_rank; ++itr) {
    // some temporaries for clarity
    std::vector<EntitySideComponent > adjacent_entities;
    Entity& curr_entity = **itr;
    const CellTopologyData* celltopology = fem::get_cell_topology(curr_entity).getCellTopologyData();
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

const CellTopologyData * get_subcell_nodes(const Entity & entity ,
                                           EntityRank subcell_rank ,
                                           unsigned subcell_identifier ,
                                           EntityVector & subcell_nodes,
                                           bool use_reverse_polarity)
{
  // get cell topology
  const CellTopologyData* celltopology = fem::get_cell_topology(entity).getCellTopologyData();
  if (celltopology == NULL) {
    return NULL;
  }

  // valid ranks fall within the dimension of the cell topology
  bool bad_rank = subcell_rank >= celltopology->dimension;

  // local id should be < number of entities of the desired type
  // (if you have 4 edges, their ids should be 0-3)
  bool bad_id = false;
  if (!bad_rank) {
    bad_id = subcell_identifier >= celltopology->subcell_count[subcell_rank];
  }

  ThrowInvalidArgMsgIf( bad_rank, "subcell_rank is >= celltopology dimension\n");
  ThrowInvalidArgMsgIf( bad_id,   "subcell_id is >= subcell_count\n");

  // For the subcell, get it's nodes and num_nodes

  const unsigned* subcell_node_local_ids =
    celltopology->subcell[subcell_rank][subcell_identifier].node;

  const CellTopologyData * subcell_topology =
    celltopology->subcell[subcell_rank][subcell_identifier].topology;
  int num_nodes_in_subcell = subcell_topology->node_count;

  // Get all the nodal relationships for this entity. We are guaranteed
  // that, if we make it this far, the entity is guaranteed to have
  // some relationship to nodes (we know it is a higher-order entity
  // than Node).
  PairIterRelation irel = entity.relations(fem::NODE_RANK);

  // Get the node entities for the nodes that make up the side. We put these
  // in in reverse order so that it has the correct orientation with respect
  // the potential adjacent entities we are evaluating. The relations are
  // right-hand-rule ordered for the owning entity, but we need something
  // that's compatible w/ the adjacent entities. However, if check_both_polarities
  // is defined, we need to do both.

  subcell_nodes.reserve(num_nodes_in_subcell);
  if (use_reverse_polarity) {
    for (int itr = num_nodes_in_subcell - 1; itr >= 0; --itr) {
      subcell_nodes.push_back(irel[subcell_node_local_ids[itr]].entity());
    }
  }
  else {
    for (int itr = 0; itr < num_nodes_in_subcell; ++itr ) {
      subcell_nodes.push_back(irel[subcell_node_local_ids[itr]].entity());
    }
  }

  return subcell_topology;
}

void get_adjacent_entities( const Entity & entity ,
                            EntityRank subcell_rank ,
                            unsigned subcell_identifier ,
                            std::vector< EntitySideComponent> & adjacent_entities,
                            bool use_reverse_polarity,
                            EntityRank adjacent_entities_rank)
{
  adjacent_entities.clear();

  // Get nodes that make up the subcell we're looking at
  EntityVector subcell_nodes;
  const CellTopologyData * subcell_topology = get_subcell_nodes(entity,
                                                                subcell_rank,
                                                                subcell_identifier,
                                                                subcell_nodes,
                                                                use_reverse_polarity);

  // Given the nodes related to the subcell, find all entities
  // with the same rank that have a relation to all of these nodes
  EntityVector potentially_adjacent_entities;
  EntityRank entity_rank_to_get = adjacent_entities_rank == fem::INVALID_RANK ?
                                  entity.entity_rank() :
                                  adjacent_entities_rank;
  get_entities_through_relations(subcell_nodes,
                                 entity_rank_to_get,
                                 potentially_adjacent_entities);

  // We don't want to include entities that are superimposed with
  // the input entity
  filter_superimposed_entities(entity, potentially_adjacent_entities);

  // Add the local ids, from the POV of the adj entitiy, to the return value
  for (EntityVector::const_iterator eitr = potentially_adjacent_entities.begin();
       eitr != potentially_adjacent_entities.end(); ++eitr) {
    int local_subcell_num = get_entity_subcell_id(**eitr,
                                                  subcell_rank,
                                                  subcell_topology,
                                                  subcell_nodes);
    if ( local_subcell_num != -1) {
      adjacent_entities.push_back(EntitySideComponent(*eitr, local_subcell_num));
    }
  }
}

}
}
