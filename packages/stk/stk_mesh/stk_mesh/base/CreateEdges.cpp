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
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CellTopology.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {
namespace mesh {

namespace {

struct EntitySubcellComponent {
  public:
    EntitySubcellComponent()
      : entity()
      , subcell_rank(0)
      , subcell_id(0)
  {}

    EntitySubcellComponent(
        Entity arg_entity,
        EntityRank   arg_subcell_rank,
        unsigned     arg_subcell_id
        )
      : entity(arg_entity)
      , subcell_rank(arg_subcell_rank)
      , subcell_id(arg_subcell_id)
  {}

    Entity entity;
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
    int local_subcell_num = get_entity_subcell_id(
      *eitr,
      subcell_rank,
      subcell_topology,
      subcell_nodes);
    if ( local_subcell_num != -1) {
      entities_with_subcell.push_back(EntitySubcellComponent(*eitr, subcell_rank, local_subcell_num));
    }
  }
}

// Check if entity has a specific relation to an entity of subcell_rank
bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
{
  bool found = false;
  PairIterRelation relations = entity.relations(subcell_rank);

  for (; !relations.empty(); ++relations) {
    if (relations->relation_ordinal() == subcell_id) {
      found = true;
      break;
    }
  }

  return found;
}

size_t internal_count_edges_to_create( BulkData & mesh)
{
  size_t edges_to_request = 0;

  const EntityRank element_rank = MetaData::ELEMENT_RANK;
  const EntityRank edge_rank = MetaData::EDGE_RANK;

  Selector select_owned = MetaData::get(mesh).locally_owned_part();

  BucketVector element_buckets;

  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);

  for (BucketVector::iterator bitr = element_buckets.begin();
      bitr != element_buckets.end();
      ++bitr)
  {
    Bucket & b = **bitr;
    const CellTopology topo = get_cell_topology(b);


    for (size_t i = 0; i<b.size(); ++i) {

      Entity elem = b[i];

      const unsigned num_edges = topo.getSubcellCount(edge_rank);

      for (size_t edge_id = 0; edge_id < num_edges; ++edge_id ) {

        if ( ! relation_exist( elem, edge_rank, edge_id) ) { //


          EntityVector edge_nodes;

          const CellTopologyData * edge_topology =
            get_subcell_nodes(
                elem,
                edge_rank,
                edge_id,
                edge_nodes
                );

          std::vector<EntitySubcellComponent> adjacent_elements;

          get_entities_with_given_subcell(
              edge_topology,
              edge_rank,
              edge_nodes,
              element_rank,
              adjacent_elements
              );

          std::reverse( edge_nodes.begin(), edge_nodes.end());

          get_entities_with_given_subcell(
              edge_topology,
              edge_rank,
              edge_nodes,
              element_rank,
              adjacent_elements
              );

          bool current_elem_has_lowest_id = true;
          //does this process own the element with the lowest id?

          for (std::vector<EntitySubcellComponent>::iterator adjacent_itr = adjacent_elements.begin();
              adjacent_itr != adjacent_elements.end();
              ++adjacent_itr)
          {
            if (adjacent_itr->entity.identifier() < elem.identifier()) {
              current_elem_has_lowest_id = false;
              break;
            }
          }

          // This process owns the lowest element so
          // needs to generate a request to create
          // the subcell
          if (current_elem_has_lowest_id) {
            ++edges_to_request;
          }
        }
      }
    }
  }
  return edges_to_request;
}

void request_edges(
   BulkData & mesh,
   size_t edges_to_request,
   EntityVector & requested_edges)
{
  MetaData & fem_meta = MetaData::get(mesh);
  const size_t num_ranks = fem_meta.entity_rank_count();
  std::vector<size_t> entities_to_request(num_ranks, 0);

  entities_to_request[MetaData::EDGE_RANK] = edges_to_request;

  requested_edges.clear();
  requested_edges.resize(num_ranks);

  EntityVector requested_entities_flat_vector;
  mesh.generate_new_entities(entities_to_request, requested_edges);
}

void internal_create_edges( BulkData & mesh, const PartVector & arg_add_parts, size_t edges_to_request)
{
  MetaData & fem_meta = MetaData::get(mesh);
  const EntityRank element_rank = MetaData::ELEMENT_RANK;
  const EntityRank edge_rank = MetaData::EDGE_RANK;


  Selector select_owned = fem_meta.locally_owned_part();

  BucketVector element_buckets;

  get_buckets( select_owned, mesh.buckets(element_rank), element_buckets);


  mesh.modification_begin();


  EntityVector requested_edges;

  request_edges(
      mesh,
      edges_to_request,
      requested_edges
      );

  size_t current_edge = 0;

    for (BucketVector::iterator bitr = element_buckets.begin();
        bitr != element_buckets.end();
        ++bitr)
    {
      Bucket & b = **bitr;
      const CellTopology topo = get_cell_topology(b);


        for (size_t i = 0; i<b.size(); ++i) {

          Entity elem = b[i];

          const unsigned num_edges = topo.getSubcellCount(edge_rank);

          for (size_t edge_id = 0; edge_id < num_edges; ++edge_id ) {

            if ( ! relation_exist( elem, edge_rank, edge_id) ) { //

              EntityVector edge_nodes;

              const CellTopologyData * edge_topology =
                get_subcell_nodes(
                    elem,
                    edge_rank,
                    edge_id,
                    edge_nodes
                    );

              std::vector<EntitySubcellComponent> adjacent_elements;

              get_entities_with_given_subcell(
                  edge_topology,
                  edge_rank,
                  edge_nodes,
                  element_rank,
                  adjacent_elements
                  );

              std::reverse( edge_nodes.begin(), edge_nodes.end());

              get_entities_with_given_subcell(
                  edge_topology,
                  edge_rank,
                  edge_nodes,
                  element_rank,
                  adjacent_elements
                  );

              bool current_elem_has_lowest_id = true;
              //does this process own the element with the lowest id?

              for (std::vector<EntitySubcellComponent>::iterator adjacent_itr = adjacent_elements.begin();
                  adjacent_itr != adjacent_elements.end();
                  ++adjacent_itr)
              {
                if (adjacent_itr->entity.identifier() < elem.identifier()) {
                  current_elem_has_lowest_id = false;
                  break;
                }
              }

              // This process owns the lowest element so
              // needs to generate a request to create
              // the edge
              if (current_elem_has_lowest_id) {
                Entity edge = requested_edges[current_edge++];

                //declare the node relations for this edge
                for (size_t n = 0; n<edge_nodes.size(); ++n) {
                  Entity node = edge_nodes[n];
                  mesh.declare_relation( edge, node, n);
                }

                mesh.declare_relation( elem, edge, edge_id);

                PartVector empty_remove_parts;

                PartVector add_parts = arg_add_parts;
                add_parts.push_back( & fem_meta.get_cell_topology_root_part(topo.getCellTopologyData(edge_rank,edge_id)));

                mesh.change_entity_parts(edge, add_parts, empty_remove_parts);
              }
            }
          }
        }
    }

  mesh.modification_end();
}

void connect_edges( BulkData & mesh )
{
  MetaData & fem_meta = MetaData::get(mesh);
  const EntityRank element_rank = MetaData::ELEMENT_RANK;
  const EntityRank edge_rank = MetaData::EDGE_RANK;

  Selector select_owned_or_shared = fem_meta.locally_owned_part() | fem_meta.globally_shared_part();

  BucketVector element_buckets;

  // Add the relationship to the correct entities. So far, we've added the
  // edges/side to the sharing element w/ lowest id, but all other elements
  // that contain that edge/side still need to have the relationship set up.
  // We do that below...

  mesh.modification_begin();

  for (EntityRank entity_rank = element_rank; entity_rank > edge_rank; --entity_rank) {

    BucketVector entity_buckets;

    get_buckets(select_owned_or_shared, mesh.buckets(entity_rank),entity_buckets);

    for (BucketVector::iterator bitr = entity_buckets.begin();
        bitr != entity_buckets.end();
        ++bitr)
    {
      Bucket & b = **bitr;
      const CellTopology topo = get_cell_topology(b);

      for (size_t i = 0; i<b.size(); ++i) {
        Entity entity = b[i];

        const unsigned num_edges = topo.getSubcellCount(edge_rank);

        for (size_t edge_id = 0; edge_id < num_edges; ++edge_id ) {

          if ( !relation_exist(entity, edge_rank, edge_id) ) {

            EntityVector subcell_nodes;

            const CellTopologyData * edge_topology =
              get_subcell_nodes(
                  entity,
                  edge_rank,
                  edge_id,
                  subcell_nodes
                  );

            std::vector<EntitySubcellComponent> adjacent_entities;

            // add polarity information to newly created relations
            // polarity information is required to correctly attached
            // degenerate elements to the correct faces and edges

            get_entities_with_given_subcell(
                edge_topology,
                edge_rank,
                subcell_nodes,
                edge_rank,
                adjacent_entities
                );

            std::reverse( subcell_nodes.begin(), subcell_nodes.end());

            get_entities_with_given_subcell(
                edge_topology,
                edge_rank,
                subcell_nodes,
                edge_rank,
                adjacent_entities
                );


            if ( !adjacent_entities.empty()) {
              mesh.declare_relation( entity, adjacent_entities[0].entity, edge_id);
            }
          }
        }
      }
    }
  }
  mesh.modification_end();
}

} // un-named namespace

void create_edges( BulkData & mesh, const PartVector & arg_add_parts)
{
  ThrowErrorMsgIf(mesh.synchronized_state() == BulkData::MODIFIABLE,
                  "stk::mesh::skin_mesh is not SYNCHRONIZED");

  // to handle degenerate topologies we anticipate the following order of operations
  //
  // connect_edges
  // count degenerate entities to create
  // create degenerate entities
  // connect_edges
  // count non degenerate entities to create
  // create non degenerate entities
  // connect_edges
  //
  // to complete the connectivity (with degenerate elements) we require that
  // polarity information to be stored on each relation

  //connect_edges(mesh);

  size_t edges_to_request = internal_count_edges_to_create( mesh );

  internal_create_edges( mesh, arg_add_parts, edges_to_request);

  connect_edges(mesh);
}

}
}
