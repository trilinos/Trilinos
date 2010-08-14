
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>


namespace stk {
namespace mesh {

namespace {

// Note that we will return the related entites in a sorted by pointer
// address.
void get_related_entities( unsigned type, const Entity& entity,
                           std::vector<Entity*>& related_entities )
{
  PairIterRelation irel = entity.relations(type);
  for ( ; !irel.empty(); ++irel ) {
    related_entities.push_back(irel->entity());
  }
  std::sort(related_entities.begin(), related_entities.end());
}

}

void skin_mesh( BulkData & mesh, unsigned mesh_rank, Part * skin_part) {
  if (mesh.synchronized_state() ==  BulkData::MODIFIABLE) {
    throw std::runtime_error("stk::mesh::skin_mesh is not SYNCHRONIZED");
  }

  if (0 == mesh_rank) {
    return;
  }

  std::vector<Entity *> entities, entities_closure;

  // select owned
  Selector owned = mesh.mesh_meta_data().locally_owned_part();
  get_selected_entities( owned,
                         mesh.buckets(fem_entity_rank(mesh_rank)),
                         entities);

  // compute owned closure
  find_closure( mesh, entities, entities_closure );

  // compute boundary
  EntitySideVector boundary;
  boundary_analysis( mesh, entities_closure, mesh_rank, boundary);

  // find the sides that need to be created to make a skin. The vector below
  // maps a side index to a vector of entities that will be sharing that side
  std::vector<std::vector<EntitySideComponent> > skin;

  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
       itr != boundary.end(); ++itr) {
    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const unsigned side_ordinal = inside.side_ordinal;
    Entity& inside_entity = *(inside.entity);

    // If this process owns the inside and the outside entity does not exist,
    // we need to apply a skin. This means we need to ensure that the side
    // entity exists at this boundary.
    if ( inside_entity.owner_rank() == mesh.parallel_rank() &&
         outside.entity == NULL ) {
      // search through existing sides
      PairIterRelation existing_sides = inside_entity.relations(mesh_rank -1);
      for (; existing_sides.first != existing_sides.second &&
             existing_sides.first->identifier() != side_ordinal ;
             ++existing_sides.first);

      // the side we need for the skin does not exist on the inside entity,
      // we may need to create it
      if (existing_sides.first == existing_sides.second) {
        // Get the nodes for the inside entity
        std::vector<Entity*> inside_entity_nodes;
        get_related_entities(Node, inside_entity, inside_entity_nodes);

        // We need to be very careful not to create multiple skins on
        // the same boundary. This is easy to do when you have multiple shells
        // on the boundary that are superimposed on each other. To avoid
        // this problem, we need to see if any of the Entities we have already
        // registered with the skin are coincident with this Entity.

        bool need_to_create_additional_side = true;
        for ( size_t i = 0; i < skin.size(); ++i) {
          for ( size_t j = 0; j < skin[i].size(); ++j ) {
            // Get the nodes for the skinned entity
            Entity & skinned_entity = *(skin[i][j].entity);
            std::vector<Entity*> skinned_entity_nodes;
            get_related_entities(Node, skinned_entity, skinned_entity_nodes);

            if ( &skinned_entity == &inside_entity) {
              continue;
            }

            // Check if either entity is coicident with the other
            bool entities_are_superimposed =
              skinned_entity_nodes == inside_entity_nodes;

            // If superimposed, we do not want to create an additional side
            // entity but we need to make sure that the side entity that is
            // already being created gets attached to this entity.
            if (entities_are_superimposed) {
              skin[i].push_back(inside);
              need_to_create_additional_side = false;
              break;
            }
          }
        }

        if (need_to_create_additional_side) {
          std::vector<EntitySideComponent> new_side_vec;
          new_side_vec.push_back(inside);
          skin.push_back(new_side_vec);
        }
      }
    }
  }

  mesh.modification_begin();

  // formulate request ids for the new sides
  std::vector<size_t> requests(mesh.mesh_meta_data().entity_rank_count(), 0);
  requests[mesh_rank -1] = skin.size();

  // create the new sides
  std::vector<Entity *> requested_entities;
  mesh.generate_new_entities(requests, requested_entities);

  // Attach the sides (skin) to the appropriate entities
  for ( size_t i = 0; i < skin.size(); ++i ) {
    Entity & side = *(requested_entities[i]);
    for ( size_t j = 0; j < skin[i].size(); ++j ) {
      Entity & entity              = *(skin[i][j].entity);
      const unsigned side_ordinal  = skin[i][j].side_ordinal;

      if (j == 0) {
        declare_element_side(entity, side, side_ordinal, skin_part);
      }
      else {
        // declare relation between entity and side
        mesh.declare_relation(entity, side, side_ordinal);
      }
    }
  }

  mesh.modification_end();
}

}
}
