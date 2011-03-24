#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <Shards_CellTopologyTraits.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

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
  const CellTopologyData * const elem_top = get_cell_topology_new( elem ).getCellTopologyData();

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
  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);
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

  const CellTopologyData * const elem_top = get_cell_topology_new( elem ).getCellTopologyData();

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

  const CellTopologyData * const elem_top = get_cell_topology_new( elem ).getCellTopologyData();
  
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

CellTopology get_cell_topology_new( const Bucket & bucket)
{
  const BulkData   &  bulk_data = BulkData::get(bucket);
  const FEMMetaData & fem_meta_data = FEMMetaData::get(bulk_data);
  const PartVector & all_parts = fem_meta_data.get_parts();

  CellTopology cell_topology;

  const std::pair< const unsigned *, const unsigned * > supersets = bucket.superset_part_ordinals();

  if (supersets.first != supersets.second) {
    const Part *first_found_part = NULL;

    for ( const unsigned * it = supersets.first ; it != supersets.second ; ++it ) {

      const Part & part = * all_parts[*it] ;

      if ( part.primary_entity_rank() == bucket.entity_rank() ) {

        CellTopology top = fem_meta_data.get_cell_topology( part );

        if ( ! cell_topology.getCellTopologyData() ) {
          cell_topology = top ;

          if (!first_found_part)
            first_found_part = &part;
        }
        else {
          ThrowErrorMsgIf( top.getCellTopologyData() && top != cell_topology,
            "Cell topology is ambiguously defined. It is defined as " << cell_topology.getName() <<
            " on part " << first_found_part->name() << " and as " << top.getName() << " on its superset part " << part.name() );
        }
      }
    }
  }

  return cell_topology ;  
}

CellTopology get_cell_topology_new( const Entity & entity)
{
  return get_cell_topology_new(entity.bucket());
}

}
}
}
