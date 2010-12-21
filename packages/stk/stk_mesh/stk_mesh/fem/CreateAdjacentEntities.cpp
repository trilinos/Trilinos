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

void get_adjacent_entities_both_polarities(
  const Entity & entity,
  const EntityRank subcell_rank,
  const unsigned subcell_id,
  std::vector< EntitySideComponent> & adjacent_entities,
  const EntityRank adjacent_entities_rank = fem::INVALID_RANK)
{
  adjacent_entities.clear();
  EntitySideComponentVector temp;
  bool reverse_polarity = false;
  do {
    get_adjacent_entities(entity,
                          subcell_rank,
                          subcell_id,
                          temp,
                          reverse_polarity,
                          adjacent_entities_rank);
    adjacent_entities.reserve(adjacent_entities.size() + temp.size());
    adjacent_entities.insert(adjacent_entities.end(), temp.begin(), temp.end());

    reverse_polarity = !reverse_polarity;
  } while (reverse_polarity);
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

  BucketVector element_buckets;
  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);

  std::vector<size_t> entities_to_request(num_ranks, 0);

  // For all elements, find all of the sides and edges that need to be
  // created to "fill-out" each element during the first pass. During
  // the second pass, actually create them and set up the relation. In
  // order to avoid double-creation of shared sides/edges, we must look
  // all elements sharing a side/edge and make sure it's only counted/created
  // once. We do this by assigning responsibility for the side/edge to
  // the element with the lowest id.

  for ( int fill_request = 1; fill_request >= 0; --fill_request ) {
    for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
      for (BucketVector::iterator bitr = element_buckets.begin();
           bitr != element_buckets.end();
           ++bitr) {
        Bucket & bucket = **bitr;
        const fem::CellTopology topo = fem::get_cell_topology(bucket);

        if ( !is_degenerate(topo) ) { // don't loop over shell elements

          for (size_t i = 0; i < bucket.size(); ++i) {
            Entity & elem = bucket[i];

            PairIterRelation node_relations = elem.relations(node_rank);

            PairIterRelation subcell_relations = elem.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( ! relation_exist( elem, subcell_rank, subcell_id) ) {
                EntitySideComponentVector adjacent_elements;
                get_adjacent_entities_both_polarities(elem,
                                                      subcell_rank,
                                                      subcell_id,
                                                      adjacent_elements);

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

            if ( ! relation_exist(elem, subcell_rank, subcell_id) ) {
              EntitySideComponentVector adjacent_entities;
              get_adjacent_entities_both_polarities(elem,
                                                    subcell_rank,
                                                    subcell_id,
                                                    adjacent_entities);

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

  mesh.modification_end();

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

            PairIterRelation node_relations = entity.relations(node_rank);

            PairIterRelation subcell_relations = entity.relations(subcell_rank);

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( !relation_exist(entity, subcell_rank, subcell_id) ) {
                EntitySideComponentVector adjacent_entities;
                get_adjacent_entities_both_polarities(entity,
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
