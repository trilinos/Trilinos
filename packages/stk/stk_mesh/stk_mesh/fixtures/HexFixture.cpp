/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

HexFixture::~HexFixture()
{}


HexFixture::HexFixture(stk::ParallelMachine pm, unsigned nx, unsigned ny, unsigned nz)
  :
    NX(nx)
  , NY(ny)
  , NZ(nz)
  , meta_data( fem_entity_rank_names() )
  , bulk_data( meta_data , pm )
  , hex_part( meta_data.declare_part("hex_part", Element) )
  , coord_field( meta_data.declare_field<CoordFieldType>("Coordinates") )
  , coord_gather_field(
        meta_data.declare_field<CoordGatherFieldType>("GatherCoordinates")
      )
{
  typedef shards::Hexahedron<8> Hex8 ;
  enum { SpatialDim = 3 };
  enum { NodesPerElem = Hex8::node_count };

  // Set topology of the element block part
  set_cell_topology<shards::Hexahedron<8> >(hex_part);

  //put coord-field on all nodes:
  put_field(
      coord_field,
      Node,
      meta_data.universal_part(),
      SpatialDim
      );

  //put coord-gather-field on all elements:
  put_field(
      coord_gather_field,
      Element,
      meta_data.universal_part(),
      NodesPerElem
      );

  // Field relation so coord-gather-field on elements points
  // to coord-field of the element's nodes
  meta_data.declare_field_relation( coord_gather_field, element_node_stencil<Hex8>, coord_field);

}

void HexFixture::generate_mesh() {
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();
  const unsigned num_elems = NX * NY * NZ ;

  const EntityId beg_elem = 1 + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = 1 + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void HexFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor) {

  //sort and unique the input elements
  std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
  std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

  std::sort( ib, ie);
  ib = std::unique( ib, ie);
  element_ids_on_this_processor.erase(ib, ie);


  bulk_data.modification_begin();

  {

    for ( unsigned iz = 0 ; iz < NZ ; ++iz ) {
    for ( unsigned iy = 0 ; iy < NY ; ++iy ) {
    for ( unsigned ix = 0 ; ix < NX ; ++ix ) {
      EntityId i_elem = elem_id( ix, iy, iz);

      std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
      std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
      ib = lower_bound( ib, ie, i_elem);

      if ( ib != ie && *ib ==  i_elem ) {

        stk::mesh::EntityId elem_node[8] ;

        elem_node[0] = node_id( ix   , iy   , iz   );
        elem_node[1] = node_id( ix+1 , iy   , iz   );
        elem_node[2] = node_id( ix+1 , iy   , iz+1 );
        elem_node[3] = node_id( ix   , iy   , iz+1 );
        elem_node[4] = node_id( ix   , iy+1 , iz   );
        elem_node[5] = node_id( ix+1 , iy+1 , iz   );
        elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
        elem_node[7] = node_id( ix   , iy+1 , iz+1 );

        stk::mesh::declare_element( bulk_data, hex_part, elem_id( ix , iy , iz ) , elem_node);
      }
    }
  }
  }
  }

  bulk_data.modification_end();

  for ( unsigned iz = 0 ; iz <= NZ ; ++iz ) {
  for ( unsigned iy = 0 ; iy <= NY ; ++iy ) {
  for ( unsigned ix = 0 ; ix <= NX ; ++ix ) {
    stk::mesh::Entity * const entity =
      bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy,iz) );

    if ( entity ) {
      Scalar * data = stk::mesh::field_data( coord_field , *entity );

      data[0] = ix ;
      data[1] = iy ;
      data[2] = - iz ;
    }
  }
  }
  }

}

} // fixtures
} // mesh
} // stk
