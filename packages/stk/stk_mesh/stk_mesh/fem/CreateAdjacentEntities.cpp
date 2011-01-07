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

// Check if 3d element topology topo is degenerate
bool is_degenerate( const fem::CellTopology & topo)
{
  return topo.getSideCount() <= 3;
}

// Check if entity has a specific relation to an entity of subcell_rank
bool relation_exist( const Entity & entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
{
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

void get_entities_with_given_subcell(
  const CellTopologyData * subcell_topology,
  const EntityRank subcell_rank,
  const EntityVector & subcell_nodes,
  const EntityRank entities_rank,
  std::vector< EntitySideComponent> & entities_with_subcell
                                     )
{
  // Get all entities that have relations to all the subcell nodes
  EntityVector entities;
  get_entities_through_relations(subcell_nodes,
                                 entities_rank,
                                 entities);

  // For all such entities, add id info for the subcell if the subcell
  // nodes compose a valid subcell of the entity
  for (EntityVector::const_iterator eitr = entities.begin();
       eitr != entities.end(); ++eitr) {
    int local_subcell_num = get_entity_subcell_id(
      **eitr,
      subcell_rank,
      subcell_topology,
      subcell_nodes);
    if ( local_subcell_num != -1) {
      entities_with_subcell.push_back(EntitySideComponent(*eitr, local_subcell_num));
    }
  }
}

// Returns all the entities that have a subcell matching the specified
// subcell.
void get_entities_with_given_subcell_for_both_polarities(
  const Entity & entity,
  const EntityRank subcell_rank,
  const unsigned subcell_id,
  std::vector< EntitySideComponent> & adjacent_entities,
  const EntityRank adjacent_entities_rank)
{
  adjacent_entities.clear();
  bool reverse_polarity = true;
  do {
    EntityVector subcell_nodes;
    const CellTopologyData * subcell_topology = get_subcell_nodes(entity,
                                                                  subcell_rank,
                                                                  subcell_id,
                                                                  subcell_nodes,
                                                                  reverse_polarity);

    get_entities_with_given_subcell(subcell_topology,
                                    subcell_rank,
                                    subcell_nodes,
                                    adjacent_entities_rank,
                                    adjacent_entities);

    reverse_polarity = !reverse_polarity;
  } while (!reverse_polarity);
}

void fill_out_elements(BulkData & mesh)
{
  const size_t num_ranks = mesh.mesh_meta_data().entity_rank_count();

  fem::FEMInterface & fem_interface = fem::get_fem_interface(mesh);
  const EntityRank element_rank = fem::element_rank(fem_interface);
  const EntityRank side_rank = fem::side_rank(fem_interface);
  const EntityRank edge_rank = fem::edge_rank(fem_interface);

  Selector select_owned = mesh.mesh_meta_data().locally_owned_part();

  BucketVector element_buckets;
  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);

  // only used when actually populating relations, keeps track of what
  // index to use when pulling things out of requested_entities
  std::vector<size_t> entities_used(num_ranks, 0);

  // used to count the number of entities we need to create to fill
  // out the elements
  std::vector<size_t> entities_to_request(num_ranks, 0);

  // A vector of EntityVector. This will serve as a map from
  // entity rank to a vector of newly requested entities of
  // that rank.
  std::vector<EntityVector> requested_entities(num_ranks);

  // The "filling out" will occur in two phases, the first for counting
  // the number of new entities to create and the second to add relations
  // to the new entities.
  bool count_only = true;
  do {
    for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
      for (BucketVector::iterator bitr = element_buckets.begin();
           bitr != element_buckets.end();
           ++bitr) {
        Bucket & bucket = **bitr;
        const fem::CellTopology topo = fem::get_cell_topology(bucket);

        if ( !is_degenerate(topo) ) { // don't loop over shell elements

          for (size_t i = 0; i < bucket.size(); ++i) {
            Entity & elem = bucket[i];

            PairIterRelation subcell_relations = elem.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( ! relation_exist( elem, subcell_rank, subcell_id) ) {
                EntitySideComponentVector adjacent_elements;
                get_entities_with_given_subcell_for_both_polarities(elem,
                                                                    subcell_rank,
                                                                    subcell_id,
                                                                    adjacent_elements,
                                                                    elem.entity_rank());

                bool current_elem_has_lowest_id = true;

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
                  if (count_only) {
                    entities_to_request[subcell_rank]++;
                  }
                  else {
                    Entity & subcell = * requested_entities[subcell_rank][entities_used[subcell_rank]++];

                    // Get nodes for this subcell in positive-polarity order
                    EntityVector subcell_nodes;
                    get_subcell_nodes(elem,
                                      subcell_rank,
                                      subcell_id,
                                      subcell_nodes,
                                      false); // positive polarity

                    // Declare the node relations for this subcell
                    for (size_t i = 0; i < subcell_nodes.size(); ++i) {
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
    }

    if (count_only) {
      // At the end of the counting phase, start a new modification cycle
      // and generate the necessary entities. The new entities will be inserted
      // into requested_entities.

      EntityVector requested_entities_flat_vector;

      mesh.modification_begin();
      mesh.generate_new_entities(entities_to_request, requested_entities_flat_vector);

      //shape the requested_entities vector
      EntityVector::iterator b_itr = requested_entities_flat_vector.begin();

      for (size_t i=0; i<num_ranks; ++i) {
        EntityVector & temp = requested_entities[i];
        temp.insert(temp.begin(), b_itr, b_itr + entities_to_request[i]);
        b_itr += entities_to_request[i];
      }

      ThrowRequire(b_itr == requested_entities_flat_vector.end());
    }
    else {
      // At the end of the relation-population phase, end the modification cycle.
      mesh.modification_end();
    }

    count_only = !count_only;
  } while (!count_only);
}

} // un-named namespace

void create_adjacent_entities( BulkData & mesh, PartVector & arg_add_parts)
{
  ThrowErrorMsgIf(mesh.synchronized_state() == BulkData::MODIFIABLE,
                  "stk::mesh::skin_mesh is not SYNCHRONIZED");

  fem::FEMInterface & fem_interface = fem::get_fem_interface(mesh);
  const EntityRank element_rank = fem::element_rank(fem_interface);
  const EntityRank side_rank = fem::side_rank(fem_interface);
  const EntityRank edge_rank = fem::edge_rank(fem_interface);

  // For all elements, find all of the sides and edges that need to be
  // created to "fill-out" each element during the first pass. During
  // the second pass, actually create them and set up the relation. In
  // order to avoid double-creation of shared sides/edges, we must look
  // all elements sharing a side/edge and make sure it's only counted/created
  // once. We do this by assigning responsibility for the side/edge to
  // the element with the lowest id.

  fill_out_elements(mesh);

  // Add the relationship to the correct entities. So far, we've added the
  // edges/side to the sharing element w/ lowest id, but all other elements
  // that contain that edge/side still need to have the relationship set up.
  // We do that below...

  for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {

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

            PairIterRelation subcell_relations = entity.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( !relation_exist(entity, subcell_rank, subcell_id) ) {
                EntitySideComponentVector adjacent_entities;
                get_entities_with_given_subcell_for_both_polarities(entity,
                                                                    subcell_rank,
                                                                    subcell_id,
                                                                    adjacent_entities,
                                                                    subcell_rank);

                ThrowRequire( adjacent_entities.size() == 1);

                mesh.declare_relation( entity, *adjacent_entities[0].entity, subcell_id);
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
