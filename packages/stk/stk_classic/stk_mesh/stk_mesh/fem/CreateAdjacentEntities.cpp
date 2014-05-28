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

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/CellTopology.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {
namespace mesh {

namespace {


struct EntitySubcellComponent {
  public:
    EntitySubcellComponent()
      : entity(NULL)
      , subcell_rank(0)
      , subcell_id(0)
  {}

    EntitySubcellComponent(
        Entity     * arg_entity,
        EntityRank   arg_subcell_rank,
        unsigned     arg_subcell_id
        )
      : entity(arg_entity)
      , subcell_rank(arg_subcell_rank)
      , subcell_id(arg_subcell_id)
  {}

    Entity     * entity;
    EntityRank subcell_rank;
    unsigned   subcell_id;
};



void get_entities_with_given_subcell(
  const CellTopologyData * subcell_topology,
  const EntityRank subcell_rank,
  const EntityVector & subcell_nodes,
  const EntityRank entities_rank,
  std::vector< EntitySubcellComponent> & entities_with_subcell
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
    int local_subcell_num = fem::get_entity_subcell_id(
      **eitr,
      subcell_rank,
      subcell_topology,
      subcell_nodes);
    if ( local_subcell_num != -1) {
      entities_with_subcell.push_back(EntitySubcellComponent(*eitr, subcell_rank, local_subcell_num));
    }
  }
}



// Check if 3d element topology topo is degenerate
bool is_degenerate( const fem::CellTopology & topo)
{
  return topo.getSideCount() < 3;
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



void internal_count_entities_to_create( BulkData & mesh, std::vector<size_t> & entities_to_request) {

  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(mesh);
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank = fem_meta.side_rank();
  const EntityRank edge_rank = fem_meta.edge_rank();

  Selector select_owned = fem::FEMMetaData::get(mesh).locally_owned_part();


  BucketVector element_buckets;

  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);


  for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {
    for (BucketVector::iterator bitr = element_buckets.begin();
        bitr != element_buckets.end();
        ++bitr)
    {
      Bucket & b = **bitr;
      const fem::CellTopology topo = fem::get_cell_topology(b);

      ThrowErrorMsgIf( is_degenerate(topo),
          "stk::mesh::create_adjacent_entities(...) does not yet support degenerate topologies (i.e. shells and beams)");


      if ( !is_degenerate(topo) ) { // don't loop over shell elements

        for (size_t i = 0; i<b.size(); ++i) {

          Entity & elem = b[i];

          const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

          for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

            if ( ! relation_exist( elem, subcell_rank, subcell_id) ) { //


              EntityVector subcell_nodes;

              const CellTopologyData * subcell_topology =
                fem::get_subcell_nodes(
                    elem,
                    subcell_rank,
                    subcell_id,
                    subcell_nodes
                    );

              std::vector<EntitySubcellComponent> adjacent_elements;

              get_entities_with_given_subcell(
                  subcell_topology,
                  subcell_rank,
                  subcell_nodes,
                  element_rank,
                  adjacent_elements
                  );

              std::reverse( subcell_nodes.begin(), subcell_nodes.end());

              get_entities_with_given_subcell(
                  subcell_topology,
                  subcell_rank,
                  subcell_nodes,
                  element_rank,
                  adjacent_elements
                  );

              bool current_elem_has_lowest_id = true;
              //does this process own the element with the lowest id?

              for (std::vector<EntitySubcellComponent>::iterator adjacent_itr = adjacent_elements.begin();
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

void request_entities(
   BulkData & mesh,
   std::vector<size_t> & entities_to_request,
   std::vector< EntityVector > & requested_entities)
{
  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(mesh);
  const size_t num_ranks = fem_meta.entity_rank_count();

  requested_entities.clear();
  requested_entities.resize(num_ranks);

  EntityVector requested_entities_flat_vector;
  mesh.generate_new_entities(entities_to_request, requested_entities_flat_vector);

  EntityVector::iterator b_itr = requested_entities_flat_vector.begin();

  for (size_t i=0; i<num_ranks; ++i) {
    EntityVector & temp = requested_entities[i];
    temp.insert(temp.begin(), b_itr, b_itr + entities_to_request[i]);
    b_itr += entities_to_request[i];
  }

  ThrowRequire(b_itr == requested_entities_flat_vector.end());

}

void internal_create_adjacent_entities( BulkData & mesh, const PartVector & arg_add_parts, std::vector<size_t> & entities_to_request) {

  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(mesh);
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank = fem_meta.side_rank();
  const EntityRank edge_rank = fem_meta.edge_rank();

  const size_t num_ranks = fem_meta.entity_rank_count();

  Selector select_owned = fem_meta.locally_owned_part();

  BucketVector element_buckets;

  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);


  mesh.modification_begin();


  std::vector< EntityVector > requested_entities;

  request_entities(
      mesh,
      entities_to_request,
      requested_entities
      );

  std::vector<size_t> entities_used(num_ranks, 0);

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

          const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

          for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

            if ( ! relation_exist( elem, subcell_rank, subcell_id) ) { //


              EntityVector subcell_nodes;

              const CellTopologyData * subcell_topology =
                fem::get_subcell_nodes(
                    elem,
                    subcell_rank,
                    subcell_id,
                    subcell_nodes
                    );

              std::vector<EntitySubcellComponent> adjacent_elements;

              get_entities_with_given_subcell(
                  subcell_topology,
                  subcell_rank,
                  subcell_nodes,
                  element_rank,
                  adjacent_elements
                  );

              std::reverse( subcell_nodes.begin(), subcell_nodes.end());

              get_entities_with_given_subcell(
                  subcell_topology,
                  subcell_rank,
                  subcell_nodes,
                  element_rank,
                  adjacent_elements
                  );

              bool current_elem_has_lowest_id = true;
              //does this process own the element with the lowest id?

              for (std::vector<EntitySubcellComponent>::iterator adjacent_itr = adjacent_elements.begin();
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
                Entity & subcell = * requested_entities[subcell_rank][entities_used[subcell_rank]++];


                //declare the node relations for this subcell
                for (size_t n = 0; n<subcell_nodes.size(); ++n) {
                  Entity & node = *subcell_nodes[n];
                  mesh.declare_relation( subcell, node, n);
                }

                mesh.declare_relation( elem, subcell, subcell_id);


                PartVector empty_remove_parts;

                PartVector add_parts = arg_add_parts;
                add_parts.push_back( & fem_meta.get_cell_topology_root_part(topo.getCellTopologyData(subcell_rank,subcell_id)));

                mesh.change_entity_parts(subcell, add_parts, empty_remove_parts);

              }
            }
          }
        }
      }
    }
  }

  mesh.modification_end();
}

void complete_connectivity( BulkData & mesh ) {

  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(mesh);
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank = fem_meta.side_rank();
  const EntityRank edge_rank = fem_meta.edge_rank();

  Selector select_owned_or_shared = fem_meta.locally_owned_part() | fem_meta.globally_shared_part();

  BucketVector element_buckets;

  // Add the relationship to the correct entities. So far, we've added the
  // edges/side to the sharing element w/ lowest id, but all other elements
  // that contain that edge/side still need to have the relationship set up.
  // We do that below...

  for ( EntityRank subcell_rank = side_rank; subcell_rank >= edge_rank; --subcell_rank) {

    mesh.modification_begin();
    for (EntityRank entity_rank = element_rank; entity_rank > subcell_rank; --entity_rank) {


      BucketVector entity_buckets;


      get_buckets(select_owned_or_shared, mesh.buckets(entity_rank),entity_buckets);

      for (BucketVector::iterator bitr = entity_buckets.begin();
          bitr != entity_buckets.end();
          ++bitr)
      {
        Bucket & b = **bitr;
        const fem::CellTopology topo = fem::get_cell_topology(b);

        ThrowErrorMsgIf( is_degenerate(topo),
          "stk::mesh::create_adjacent_entities(...) does not yet support degenerate topologies (i.e. shells and beams)");

        {
          for (size_t i = 0; i<b.size(); ++i) {
            Entity & entity = b[i];

            const unsigned subcell_count = topo.getSubcellCount(subcell_rank);

            for (size_t subcell_id = 0; subcell_id < subcell_count; ++subcell_id ) {

              if ( !relation_exist(entity, subcell_rank, subcell_id) ) {

                EntityVector subcell_nodes;

                const CellTopologyData * subcell_topology =
                  fem::get_subcell_nodes(
                      entity,
                      subcell_rank,
                      subcell_id,
                      subcell_nodes
                      );

                std::vector<EntitySubcellComponent> adjacent_entities;

                // add polarity information to newly created relations
                // polarity information is required to correctly attached
                // degenerate elements to the correct faces and edges

                get_entities_with_given_subcell(
                    subcell_topology,
                    subcell_rank,
                    subcell_nodes,
                    subcell_rank,
                    adjacent_entities
                    );

                std::reverse( subcell_nodes.begin(), subcell_nodes.end());

                get_entities_with_given_subcell(
                    subcell_topology,
                    subcell_rank,
                    subcell_nodes,
                    subcell_rank,
                    adjacent_entities
                    );


                if ( !adjacent_entities.empty()) {

                  mesh.declare_relation( entity, *adjacent_entities[0].entity, subcell_id);
                }
              }
            }
          }
        }
      }
    }
    mesh.modification_end();
  }

}

} // un-named namespace

void create_adjacent_entities( BulkData & mesh, PartVector & arg_add_parts)
{
  ThrowErrorMsgIf(mesh.synchronized_state() == BulkData::MODIFIABLE,
                  "stk::mesh::skin_mesh is not SYNCHRONIZED");

  // to handle degenerate topologies we anticipate the following order of operations
  //
  // complete_connectivity
  // count degenerate entities to create
  // create degenerate entities
  // complete_connectivity
  // count non degenerate entities to create
  // create non degenerate entities
  // complete_connectivity
  //
  // to complete the connectivity (with degenerate elements) we require that
  // polarity information to be stored on each relation


  complete_connectivity(mesh);


  fem::FEMMetaData & fem_meta = fem::FEMMetaData::get(mesh);
  const size_t num_ranks = fem_meta.entity_rank_count();
  std::vector<size_t> entities_to_request(num_ranks, 0);

  internal_count_entities_to_create( mesh, entities_to_request);

  internal_create_adjacent_entities( mesh, arg_add_parts, entities_to_request);


  complete_connectivity(mesh);


}

}
}
