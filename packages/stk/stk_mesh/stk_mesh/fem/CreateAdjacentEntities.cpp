
#include <map>
#include <set>
#include <algorithm>

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Relation.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>
#include <stk_mesh/fem/CellTopology.hpp>

#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {
namespace mesh {

namespace {
bool is_degenerate( const fem::CellTopology & topo) {
  return topo.getSideCount() <= 3;
}

bool relation_exist( const Entity & entity, EntityRank subcell_rank, RelationIdentifier subcell_id ) {

  bool found = false;
  PairIterRelation relations = entity.relations(subcell_rank);

  for (; !relations.empty(); ++relations) {
    if (relations->identifier() == subcell_id) {
      found = true;
      break;
    }
  }

  return found;

}

int get_entity_subcell_id( const Entity & entity ,
    const EntityRank subcell_rank,
    const CellTopologyData * subcell_topology,
    const std::vector<Entity*>& subcell_nodes )
{
  const int INVALID_SIDE = -1;

  unsigned num_nodes = subcell_topology->node_count;

  if (num_nodes != subcell_nodes.size()) {
    return INVALID_SIDE;
  }

  const CellTopologyData* entity_topology = fem::get_cell_topology(entity).getCellTopologyData();

  if (entity_topology == NULL) {
    return INVALID_SIDE;
  }

  // get nodal relations for entity
  PairIterRelation relations = entity.relations(fem::NODE_RANK);

  const int num_permutations = subcell_topology->permutation_count;


  // Iterate over the subcells of entity...
  for (unsigned local_subcell_ordinal = 0;
      local_subcell_ordinal < entity_topology->subcell_count[subcell_rank];
      ++local_subcell_ordinal) {

    // get topological data for this subcell
    const CellTopologyData* curr_subcell_topology =
      entity_topology->subcell[subcell_rank][local_subcell_ordinal].topology;

    // If topologies are not the same, there is no way the subcells are the same
    if (subcell_topology == curr_subcell_topology) {

      const unsigned* const subcell_node_map = entity_topology->subcell[subcell_rank][local_subcell_ordinal].node;

      // Taking all positive permutations into account, check if this subcell
      // has the same nodes as the subcell_nodes argument. Note that this
      // implementation preserves the node-order so that we can take
      // entity-orientation into account.
      for (int p = 0; p < num_permutations; ++p) {

        if (curr_subcell_topology->permutation[p].polarity ==
            CELL_PERMUTATION_POLARITY_POSITIVE) {

          const unsigned * const perm_node =
            curr_subcell_topology->permutation[p].node ;

          bool all_match = true;
          for (unsigned j = 0 ; j < num_nodes; ++j ) {
            if (subcell_nodes[j] !=
                relations[subcell_node_map[perm_node[j]]].entity()) {
              all_match = false;
              break;
            }
          }

          // all nodes were the same, we have a match
          if ( all_match ) {
            return local_subcell_ordinal ;
          }
        }
      }
    }
  }

  return INVALID_SIDE;
}

const CellTopologyData * get_subcell_nodes(
    const Entity & entity ,
    unsigned subcell_rank ,
    unsigned subcell_identifier ,
    EntityVector & subcell_nodes,
    bool use_reverse_polarity = true
    )
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

  ThrowRequireMsg( ! bad_rank, "subcell_rank is >= celltopology dimension\n");
  ThrowRequireMsg( ! bad_id, "subcell_id is >= subcell_count\n");


  // For the potentially common subcell, get it's nodes and num_nodes
  const unsigned* subcell_node_local_ids =
    celltopology->subcell[subcell_rank][subcell_identifier].node;

  const CellTopologyData * subcell_topology =
    celltopology->subcell[subcell_rank][subcell_identifier].topology;
  int num_nodes_in_subcell = subcell_topology->node_count;

  // Get all the nodal relationships for this entity. We are guaranteed
  // that, if we make it this far, the entity is guaranteed to have
  // some relationship to nodes (we know it is a higher-order entity
  // than Node).

  subcell_nodes.reserve(num_nodes_in_subcell);
  PairIterRelation irel = entity.relations(fem::NODE_RANK);

  if (use_reverse_polarity) {
    for (int itr = num_nodes_in_subcell; itr > 0; ) {
      --itr;
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

void get_adjacent_entities(
    const CellTopologyData * subcell_topology,
    const EntityRank subcell_rank,
    const EntityVector & subcell_nodes,
    const EntityRank adjacent_entities_rank,
    std::vector< EntitySideComponent> & adjacent_entities
    )
{

  // Given the nodes related to the subcell, find all entities
  // with the given rank that have a relation to all of these nodes
  EntityVector entities;
  get_entities_through_relations(subcell_nodes,
      adjacent_entities_rank,
      entities);


  // Add the local ids, from the POV of the adj entitiy, to the return value
  for (EntityVector::const_iterator eitr = entities.begin();
      eitr != entities.end(); ++eitr) {
    int local_subcell_num = get_entity_subcell_id(
        **eitr,
        subcell_rank,
        subcell_topology,
        subcell_nodes);
    if ( local_subcell_num != -1) {
      adjacent_entities.push_back(EntitySideComponent(*eitr, local_subcell_num));
    }
  }
}

} // un-named namespace

void create_adjacent_entities( BulkData & mesh, PartVector & arg_add_parts)
{
  ThrowErrorMsgIf(mesh.synchronized_state() == BulkData::MODIFIABLE,
                  "stk::mesh::skin_mesh is not SYNCHRONIZED");

  const size_t num_ranks = mesh.mesh_meta_data().entity_rank_count();

  fem::FEMInterface & fem_interface = fem::get_fem_interface(mesh);
  const EntityRank element_rank = fem::element_rank(fem_interface);
  const EntityRank side_rank = fem::side_rank(fem_interface);
  const EntityRank edge_rank = fem::edge_rank(fem_interface);
  const EntityRank node_rank = fem::edge_rank(fem_interface);

  Selector select_owned = mesh.mesh_meta_data().locally_owned_part();

  {

    BucketVector element_buckets;

    get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);

    std::vector<size_t> entities_to_request(num_ranks, 0);

    for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
      for (BucketVector::iterator bitr = element_buckets.begin();
          bitr != element_buckets.end();
          ++bitr)
      {
        Bucket & b = **bitr;
        const fem::CellTopology topo = fem::get_cell_topology(b);

        if ( !is_degenerate(topo) ) { // don't loop over shell elements

          for (size_t i = 0; i<b.size(); ++i) {

            Entity & elem = b[i];


            PairIterRelation node_relations = elem.relations(node_rank);


            PairIterRelation subcell_relations = elem.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( ! relation_exist( elem, subcell_rank, subcell_id) ) { //


                EntityVector subcell_nodes, reverse_subcell_nodes;

                const CellTopologyData * subcell_topology =
                  get_subcell_nodes(
                      elem,
                      subcell_rank,
                      subcell_id,
                      reverse_subcell_nodes,
                      true // reverse polarity
                      );

                get_subcell_nodes(
                    elem,
                    subcell_rank,
                    subcell_id,
                    subcell_nodes,
                    false // positive polarity
                    );

                EntitySideComponentVector adjacent_elements;

                get_adjacent_entities(
                    subcell_topology,
                    subcell_rank,
                    subcell_nodes,
                    element_rank,
                    adjacent_elements
                    );

                get_adjacent_entities(
                    subcell_topology,
                    subcell_rank,
                    reverse_subcell_nodes,
                    element_rank,
                    adjacent_elements
                    );

                bool current_elem_has_lowest_id = true;
                //does this process own the element with the lowest id?

                for (EntitySideComponentVector::iterator adjacent_itr = adjacent_elements.begin();
                    adjacent_itr != adjacent_elements.end();
                    ++adjacent_itr)
                {
                  if (adjacent_itr->entity->identifier() < elem.identifier()) {
                    current_elem_has_lowest_id = false;
                    break;
                  }
                }

                // This process owns the lowest element so
                // needs to generate a request to create
                // the subcell
                if (current_elem_has_lowest_id) {
                  entities_to_request[subcell_rank]++;
                }
              }
            }
          }
        }
      }
    }

    EntityVector requested_entities_flat_vector;

    mesh.modification_begin();
    mesh.generate_new_entities(entities_to_request, requested_entities_flat_vector);

    std::vector< EntityVector > requested_entities(num_ranks);

    //shape the requested_entities vector
    EntityVector::iterator b_itr = requested_entities_flat_vector.begin();

    for (size_t i=0; i<num_ranks; ++i) {
      EntityVector & temp = requested_entities[i];
      temp.insert(temp.begin(), b_itr, b_itr + entities_to_request[i]);
      b_itr += entities_to_request[i];
    }

    ThrowRequire(b_itr == requested_entities_flat_vector.end());

    std::vector<size_t> entities_used(num_ranks, 0);

    for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
      //add the relationship to the correct entities
      for (BucketVector::iterator bitr = element_buckets.begin();
          bitr != element_buckets.end();
          ++bitr)
      {
        Bucket & b = **bitr;
        const fem::CellTopology topo = fem::get_cell_topology(b);

        if ( !is_degenerate(topo) ) { // don't loop over shell elements

          for (size_t i = 0; i<b.size(); ++i) {
            Entity & elem = b[i];

            PairIterRelation node_relations = elem.relations(node_rank);


            PairIterRelation subcell_relations = elem.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( ! relation_exist(elem, subcell_rank, subcell_id) ) { //

                EntityVector subcell_nodes, reverse_subcell_nodes;

                const CellTopologyData * subcell_topology =
                  get_subcell_nodes(
                      elem,
                      subcell_rank,
                      subcell_id,
                      reverse_subcell_nodes,
                      true // reverse polarity
                      );

                get_subcell_nodes(
                    elem,
                    subcell_rank,
                    subcell_id,
                    subcell_nodes,
                    false // positive polarity
                    );

                EntitySideComponentVector adjacent_entities;

                get_adjacent_entities(
                    subcell_topology,
                    subcell_rank,
                    subcell_nodes,
                    element_rank,
                    adjacent_entities
                    );

                get_adjacent_entities(
                    subcell_topology,
                    subcell_rank,
                    reverse_subcell_nodes,
                    element_rank,
                    adjacent_entities
                    );

                bool current_elem_has_lowest_id = true;
                //does this process own the element with the lowest id?

                for (EntitySideComponentVector::iterator adjacent_itr = adjacent_entities.begin();
                    adjacent_itr != adjacent_entities.end();
                    ++adjacent_itr)
                {
                  if (adjacent_itr->entity->identifier() < elem.identifier()) {
                    current_elem_has_lowest_id = false;
                    break;
                  }
                }

                // This process owns the lowest element so
                // needs to generate a request to create
                // the subcell
                if (current_elem_has_lowest_id) {
                  Entity & subcell = * requested_entities[subcell_rank][entities_used[subcell_rank]++];


                  //declare the node relations for this subcell
                  for (size_t i = 0; i<subcell_nodes.size(); ++i) {
                    Entity & node = *subcell_nodes[i];
                    mesh.declare_relation( subcell, node, i);
                  }

                  mesh.declare_relation( elem, subcell, subcell_id);


                  PartVector add_parts, empty_remove_parts;
                  add_parts.push_back( & fem::get_part( mesh.mesh_meta_data(), topo.getTopology(subcell_rank,subcell_id)));

                  mesh.change_entity_parts(subcell, add_parts, empty_remove_parts);


                }
              }
            }
          }
        }
      }
    }

    mesh.modification_end();




    for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
      //add the relationship to the correct entities

      for (EntityRank entity_rank = element_rank; entity_rank > edge_rank; --entity_rank) {

        if (entity_rank <= subcell_rank) {
          continue;
        }

        mesh.modification_begin();
        BucketVector entity_buckets;

        Selector select_owned_or_shared = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();

        get_buckets(select_owned_or_shared, mesh.buckets(entity_rank),entity_buckets);

        for (BucketVector::iterator bitr = entity_buckets.begin();
            bitr != entity_buckets.end();
            ++bitr)
        {
          Bucket & b = **bitr;
          const fem::CellTopology topo = fem::get_cell_topology(b);

          {

            for (size_t i = 0; i<b.size(); ++i) {
              Entity & entity = b[i];

              PairIterRelation node_relations = entity.relations(node_rank);

              PairIterRelation subcell_relations = entity.relations(subcell_rank);

              const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

              for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

                if ( !relation_exist(entity, subcell_rank, subcell_id) ) { //

                  EntityVector subcell_nodes, reverse_subcell_nodes;

                  const CellTopologyData * subcell_topology =
                    get_subcell_nodes(
                        entity,
                        subcell_rank,
                        subcell_id,
                        reverse_subcell_nodes,
                        true // reverse polarity
                        );

                  get_subcell_nodes(
                      entity,
                      subcell_rank,
                      subcell_id,
                      subcell_nodes,
                      false // positive polarity
                      );

                  EntitySideComponentVector adjacent_entity;

                  get_adjacent_entities(
                      subcell_topology,
                      subcell_rank,
                      subcell_nodes,
                      subcell_rank,
                      adjacent_entity
                      );

                  get_adjacent_entities(
                      subcell_topology,
                      subcell_rank,
                      reverse_subcell_nodes,
                      subcell_rank,
                      adjacent_entity
                      );

                  ThrowRequire( adjacent_entity.size() == 1);

                  mesh.declare_relation( entity, *adjacent_entity[0].entity, subcell_id);

                }
              }
            }
          }
        }
        mesh.modification_end();
      }
    }
  }
}

}
}
