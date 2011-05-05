#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <Shards_CellTopologyTraits.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <sstream>
#include <stdexcept>

namespace stk {
namespace mesh {
namespace fem {

namespace {

void verify_declare_element_side(
    const BulkData & mesh,
    const Entity & elem,
    const unsigned local_side_id
    )
{
  const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();

  const CellTopologyData * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  ThrowErrorMsgIf( &mesh != & BulkData::get(elem),
    "For elem " << print_entity_key(elem) <<
    ", Bulkdata for 'elem' and mesh are different");

  ThrowErrorMsgIf( elem_top && local_side_id >= elem_top->side_count,
    "For elem " << print_entity_key(elem) << ", local_side_id " << local_side_id << ", " <<
    "local_side_id exceeds " << elem_top->name << ".side_count = " << elem_top->side_count );

  ThrowErrorMsgIf( side_top == NULL,
    "For elem " << print_entity_key(elem) << ", local_side_id " << local_side_id << ", " <<
    "No element topology found");
}

} // unnamed namespace

Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const EntityId elem_id ,
                          const EntityId node_id[] )
{
  FEMMetaData & fem_meta = FEMMetaData::get(mesh);
  const CellTopologyData * const top = fem_meta.get_cell_topology( part ).getCellTopologyData();

  ThrowErrorMsgIf(top == NULL,
                  "Part " << part.name() << " does not have a local topology");

  PartVector empty ;
  PartVector add( 1 ); add[0] = & part ;

  const EntityRank entity_rank = fem_meta.element_rank();

  Entity & elem = mesh.declare_entity( entity_rank, elem_id, add );

  const EntityRank node_rank = fem_meta.node_rank();

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    //declare node if it doesn't already exist
    Entity * node = mesh.get_entity( node_rank , node_id[i]);
    if ( NULL == node) {
      node = & mesh.declare_entity( node_rank , node_id[i], empty );
    }

    mesh.declare_relation( elem , *node , i );
  }
  return elem ;
}

Entity & declare_element_side(
  Entity & elem ,
  Entity & side,
  const unsigned local_side_id ,
  Part * part )
{
  BulkData & mesh = BulkData::get(side);

  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << elem.identifier() << "] has no defined topology" );

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  ThrowErrorMsgIf( side_top == NULL,
      "Element[" << elem.identifier() << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  PartVector add_parts ;

  if ( part ) { add_parts.push_back( part ); }

  mesh.change_entity_parts(side, add_parts);

  mesh.declare_relation( elem , side , local_side_id );

  PairIterRelation rel = elem.relations( FEMMetaData::NODE_RANK );

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    Entity & node = * rel[ side_node_map[i] ].entity();
    mesh.declare_relation( side , node , i );
  }

  return side ;
}

Entity & declare_element_side(
  BulkData & mesh ,
  const stk::mesh::EntityId global_side_id ,
  Entity & elem ,
  const unsigned local_side_id ,
  Part * part )
{
  verify_declare_element_side(mesh, elem, local_side_id);

  const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();

  ThrowErrorMsgIf( elem_top == NULL,
      "Element[" << elem.identifier() << "] has no defined topology");

  const CellTopologyData * const side_top = elem_top->side[ local_side_id ].topology;

  ThrowErrorMsgIf( side_top == NULL,
      "Element[" << elem.identifier() << "], local_side_id = " <<
      local_side_id << ", side has no defined topology" );

  PartVector empty_parts ;

  Entity & side = mesh.declare_entity( side_top->dimension , global_side_id, empty_parts );
  return declare_element_side( elem, side, local_side_id, part);
}



const CellTopologyData * get_subcell_nodes(const Entity & entity ,
                                           EntityRank subcell_rank ,
                                           unsigned subcell_identifier ,
                                           EntityVector & subcell_nodes)
{
  subcell_nodes.clear();

  // get cell topology
  const CellTopologyData* celltopology = get_cell_topology(entity).getCellTopologyData();

  //error checking
  {
    //no celltopology defined
    if (celltopology == NULL) {
      return NULL;
    }

    // valid ranks fall within the dimension of the cell topology
    const bool bad_rank = subcell_rank >= celltopology->dimension;
    ThrowInvalidArgMsgIf( bad_rank, "subcell_rank is >= celltopology dimension\n");

    // subcell_identifier must be less than the subcell count
    const bool bad_id = subcell_identifier >= celltopology->subcell_count[subcell_rank];
    ThrowInvalidArgMsgIf( bad_id,   "subcell_id is >= subcell_count\n");
  }

  // Get the cell topology of the subcell
  const CellTopologyData * subcell_topology =
    celltopology->subcell[subcell_rank][subcell_identifier].topology;

  const int num_nodes_in_subcell = subcell_topology->node_count;

  // For the subcell, get it's local nodes ids
  const unsigned* subcell_node_local_ids =
    celltopology->subcell[subcell_rank][subcell_identifier].node;

  FEMMetaData & fem_meta = FEMMetaData::get(entity);
  const EntityRank node_rank = fem_meta.node_rank();
  PairIterRelation node_relations = entity.relations(node_rank);

  subcell_nodes.reserve(num_nodes_in_subcell);

  for (int i = 0; i < num_nodes_in_subcell; ++i ) {
    subcell_nodes.push_back( node_relations[subcell_node_local_ids[i]].entity() );
  }

  return subcell_topology;
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

  // get topology of elem
  const CellTopologyData* entity_topology = get_cell_topology(entity).getCellTopologyData();
  if (entity_topology == NULL) {
    return INVALID_SIDE;
  }

  // get nodal relations for entity
  FEMMetaData & fem_meta = FEMMetaData::get(entity);
  const EntityRank node_rank = fem_meta.node_rank();
  PairIterRelation relations = entity.relations(node_rank);

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

bool comm_mesh_counts( BulkData & M ,
                       std::vector<size_t> & counts ,
                       bool local_flag )
{
  const size_t zero = 0 ;

  // Count locally owned entities

  const FEMMetaData & S = FEMMetaData::get(M);
  const unsigned entity_rank_count = S.entity_rank_count();
  const size_t   comm_count        = entity_rank_count + 1 ;

  std::vector<size_t> local(  comm_count , zero );
  std::vector<size_t> global( comm_count , zero );

  ParallelMachine comm = M.parallel();
  Part & owns = S.locally_owned_part();

  for ( unsigned i = 0 ; i < entity_rank_count ; ++i ) {
    const std::vector<Bucket*> & ks = M.buckets( i );

    std::vector<Bucket*>::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( has_superset( **ik , owns ) ) {
        local[i] += (*ik)->size();
      }
    }
  }

  local[ entity_rank_count ] = local_flag ;

  stk::all_reduce_sum( comm , & local[0] , & global[0] , comm_count );

  counts.assign( global.begin() , global.begin() + entity_rank_count );

  return 0 < global[ entity_rank_count ] ;
}

bool element_side_polarity( const Entity & elem ,
                            const Entity & side , int local_side_id )
{
  // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
  FEMMetaData &fem_meta = FEMMetaData::get(elem);
  const bool is_side = side.entity_rank() != fem_meta.edge_rank();
  const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();

  const unsigned side_count = ! elem_top ? 0 : (
                                is_side ? elem_top->side_count
                                        : elem_top->edge_count );

  ThrowErrorMsgIf( elem_top == NULL,
                   "For Element[" << elem.identifier() << "], element has no defined topology");

  ThrowErrorMsgIf( local_side_id < 0 || static_cast<int>(side_count) <= local_side_id,
    "For Element[" << elem.identifier() << "], " <<
    "side: " << print_entity_key(side) << ", " <<
    "local_side_id = " << local_side_id <<
    " ; unsupported local_side_id");

  const CellTopologyData * const side_top =
    is_side ? elem_top->side[ local_side_id ].topology
            : elem_top->edge[ local_side_id ].topology ;

  const unsigned * const side_map =
    is_side ? elem_top->side[ local_side_id ].node
            : elem_top->edge[ local_side_id ].node ;

  const PairIterRelation elem_nodes = elem.relations( FEMMetaData::NODE_RANK );
  const PairIterRelation side_nodes = side.relations( FEMMetaData::NODE_RANK );

  bool good = true ;
  for ( unsigned j = 0 ; good && j < side_top->node_count ; ++j ) {
    good = side_nodes[j].entity() == elem_nodes[ side_map[j] ].entity();
  }
  return good ;
}

}
}
}
