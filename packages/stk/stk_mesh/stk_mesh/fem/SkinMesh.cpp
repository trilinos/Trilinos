
#include <map>

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

#include <stk_mesh/fem/SkinMesh.hpp>

namespace stk {
namespace mesh {

namespace {

typedef std::multimap< std::vector<EntityId>, EntitySideComponent>  BoundaryMap;
typedef std::pair< std::vector<EntityId>, EntitySideComponent>  BoundaryValue;

// \TODO Without the sorting, this would be a good general utility for 'fem'
void get_elem_side_nodes( const Entity & elem,
                          unsigned side_ordinal,
                          std::vector<EntityId> & nodes
                        )
{
  const CellTopologyData * const elem_top = TopologicalMetaData::get_cell_topology( elem );

  if (elem_top == NULL) {
    std::ostringstream msg ;
    msg << "skin_mesh" ;
    msg << "( Element[" << elem.identifier() << "] has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }

  const CellTopologyData * const side_top = elem_top->side[ side_ordinal ].topology;

  if (side_top == NULL) {
    std::ostringstream msg ;
    msg << "skin_mesh" ;
    msg << "( Element[" << elem.identifier() << "]" ;
    msg << " , side_ordinal = " << side_ordinal << " ) FAILED: " ;
    msg << " Side has no defined topology" ;
    throw std::runtime_error( msg.str() );
  }

  const unsigned * const side_node_map = elem_top->side[ side_ordinal ].node ;

  PairIterRelation relations = elem.relations( BaseEntityRank );

  // Find positive polarity permutation that starts with lowest entity id:
  // We're using this as the unique key in a map for a side.
  const int num_permutations = side_top->permutation_count;
  int lowest_entity_id_permutation = 0;
  for (int p = 0; p < num_permutations; ++p) {

    if (side_top->permutation[p].polarity ==
        CELL_PERMUTATION_POLARITY_POSITIVE) {
     
      const unsigned * const pot_lowest_perm_node =
        side_top->permutation[p].node ;

      const unsigned * const curr_lowest_perm_node = 
        side_top->permutation[lowest_entity_id_permutation].node;

      unsigned first_node_index = 0;

      unsigned current_lowest_entity_id = 
        relations[side_node_map[curr_lowest_perm_node[first_node_index]]].entity()->identifier();

      unsigned potential_lowest_entity_id = 
        relations[side_node_map[pot_lowest_perm_node[first_node_index]]].entity()->identifier();

      if ( potential_lowest_entity_id < current_lowest_entity_id ) {
        lowest_entity_id_permutation = p;
      }
    }
  }
  const unsigned * const perm_node = 
    side_top->permutation[lowest_entity_id_permutation].node;

  nodes.reserve(side_top->node_count);
  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    unsigned node_id = relations[side_node_map[perm_node[i]]].entity()->identifier(); 
    nodes.push_back(node_id);
  }

}

}

void skin_mesh( BulkData & mesh, EntityRank element_rank, Part * skin_part) {
  if (mesh.synchronized_state() ==  BulkData::MODIFIABLE) {
    throw std::runtime_error("stk::mesh::skin_mesh is not SYNCHRONIZED");
  }

  EntityVector owned_elements;

  // select owned
  Selector owned = mesh.mesh_meta_data().locally_owned_part();
  get_selected_entities( owned,
                         mesh.buckets(element_rank),
                         owned_elements);

  reskin_mesh(mesh, element_rank, owned_elements, skin_part);
}

void reskin_mesh( BulkData & mesh, EntityRank element_rank, EntityVector & owned_elements, Part * skin_part) {
  if (mesh.synchronized_state() ==  BulkData::MODIFIABLE) {
    throw std::runtime_error("stk::mesh::skin_mesh is not SYNCHRONIZED");
  }

  EntityVector elements_closure;

  // compute owned closure
  find_closure( mesh, owned_elements, elements_closure );

  // compute boundary
  EntitySideVector boundary;
  boundary_analysis( mesh, elements_closure, element_rank, boundary);

  // find the sides that need to be created to make a skin. The map below
  // maps sorted EntityVectors of nodes to EntitySideComponents
  BoundaryMap skin;

  int num_new_sides = 0;
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
      PairIterRelation existing_sides = inside_entity.relations(element_rank -1);
      for (; existing_sides.first != existing_sides.second &&
             existing_sides.first->identifier() != side_ordinal ;
             ++existing_sides.first);

      // the side we need for the skin does not exist on the inside entity,
      // we may need to create it
      if (existing_sides.first == existing_sides.second) {
        // Get the nodes for the inside entity
        std::vector<EntityId> inside_key;
        get_elem_side_nodes(inside_entity, side_ordinal, inside_key);

        //if a side for these nodes is not present in the map already increment the
        //number of side we need to create
        if (skin.count(inside_key) == 0) {
          ++num_new_sides;
        }

        skin.insert( BoundaryValue(inside_key, inside));

      }
    }
  }

  mesh.modification_begin();

  // formulate request ids for the new sides
  std::vector<size_t> requests(mesh.mesh_meta_data().entity_rank_count(), 0);
  requests[element_rank -1] = num_new_sides;

  // create the new sides
  EntityVector requested_sides;
  mesh.generate_new_entities(requests, requested_sides);

  std::vector<EntityId> previous_key;
  size_t current_side = 0;
  for ( BoundaryMap::const_iterator i = skin.begin(); i!= skin.end(); ++i) {
    const std::vector<EntityId> & key = i->first;
    Entity & entity = *(i->second.entity);
    const unsigned side_ordinal = i->second.side_ordinal;

    if (key != previous_key) {
      Entity & side = *(requested_sides[current_side]);
      declare_element_side(entity, side, side_ordinal, skin_part);

      previous_key = key;
      ++current_side;
    }
    else {
      Entity & side = *(requested_sides[current_side-1]);
      mesh.declare_relation(entity, side, side_ordinal);
    }
  }

  mesh.modification_end();
}

}
}
