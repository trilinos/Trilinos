
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

void skin_mesh( BulkData & mesh, unsigned mesh_rank, Part * part) {

  if (mesh.synchronized_state() ==  BulkData::MODIFIABLE) {
    throw std::runtime_error("stk::mesh::skin_mesh is not SYNCHRONIZED");
  }


  if ( 0 == mesh_rank) {
    return;
  }

  std::vector<Entity *> entities, entities_closure;

  Selector owned = mesh.mesh_meta_data().locally_owned_part();

  get_selected_entities( owned, mesh.buckets(fem_entity_type(mesh_rank)), entities) ; // select owned
  find_closure( mesh, entities, entities_closure);

  EntitySideVector boundary;
  boundary_analysis( mesh, entities_closure, mesh_rank, boundary);

  std::vector<EntitySideComponent> skin;
  //find the sides that need to be created to make a skin
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
       itr != boundary.end(); ++itr) {

    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const unsigned side_ordinal = inside.side_ordinal;

    if ( inside.entity->owner_rank() == mesh.parallel_rank() &&
        outside.entity == NULL ) {

      PairIterRelation existing_sides = inside.entity->relations(mesh_rank -1);

      for (; existing_sides.first != existing_sides.second &&
             existing_sides.first->identifier() != side_ordinal ;
             ++existing_sides.first);
      //new side
      if (existing_sides.first == existing_sides.second) {
       skin.push_back(inside);
      }
    }
  }

  mesh.modification_begin();

  std::vector<size_t> requests(mesh.mesh_meta_data().entity_type_count(), 0);
  std::vector<stk::mesh::Entity *> requested_entities;

  //request ids for the new sides
  requests[mesh_rank -1] = skin.size();

  mesh.generate_new_entities(requests, requested_entities);


  //create the skin
  for ( size_t i = 0; i < skin.size(); ++i) {
    stk::mesh::Entity & entity = *(skin[i].entity);
    const unsigned side_ordinal  = skin[i].side_ordinal;
    stk::mesh::Entity & side   = * (requested_entities[i]);

    stk::mesh::declare_element_side(entity, side, side_ordinal, part);
  }

  mesh.modification_end();
}

}
}

