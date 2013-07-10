/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

  HexFixture::HexFixture(   stk::ParallelMachine pm
              , unsigned nx
              , unsigned ny
              , unsigned nz
              , ConnectivityMap * connectivity_map
            )
  : m_spatial_dimension(3),
    m_nx(nx),
    m_ny(ny),
    m_nz(nz),
    m_fem_meta( m_spatial_dimension ),
    m_bulk_data(  m_fem_meta
                , pm
#ifdef SIERRA_MIGRATION
                , false
#endif
                , connectivity_map
               ),
    m_hex_part( declare_part<shards::Hexahedron<8> >(m_fem_meta, "hex_part") ),
    m_node_part( m_fem_meta.declare_part("node_part", MetaData::NODE_RANK ) ),
    m_coord_field( m_fem_meta.declare_field<CoordFieldType>("Coordinates") )
{
  typedef shards::Hexahedron<8> Hex8 ;

  //put coord-field on all nodes:
  put_field(
    m_coord_field,
    MetaData::NODE_RANK,
    m_fem_meta.universal_part(),
    m_spatial_dimension);

}

void HexFixture::generate_mesh()
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny * m_nz ;

  const EntityId beg_elem = 1 + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = 1 + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void HexFixture::node_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const
{
  entity_id -= 1;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id % (m_ny+1);
  entity_id /= (m_ny+1);

  z = entity_id;
}

void HexFixture::elem_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const
{
  entity_id -= 1;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id % m_ny;
  entity_id /= m_ny;

  z = entity_id;
}

void HexFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor)
{
  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  m_bulk_data.modification_begin();

  {
    // Declare the elements that belong on this process

    PartVector add_parts;
    add_parts.push_back(&m_node_part);

    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0, iz = 0;
      elem_x_y_z(entity_id, ix, iy, iz);

      stk::mesh::EntityId elem_node[8] ;

      elem_node[0] = node_id( ix   , iy   , iz   );
      elem_node[1] = node_id( ix+1 , iy   , iz   );
      elem_node[2] = node_id( ix+1 , iy+1 , iz   );
      elem_node[3] = node_id( ix   , iy+1 , iz   );
      elem_node[4] = node_id( ix   , iy   , iz+1 );
      elem_node[5] = node_id( ix+1 , iy   , iz+1 );
      elem_node[6] = node_id( ix+1 , iy+1 , iz+1 );
      elem_node[7] = node_id( ix   , iy+1 , iz+1 );

      stk::mesh::declare_element( m_bulk_data, m_hex_part, elem_id( ix , iy , iz ) , elem_node);

      for (unsigned i = 0; i<8; ++i) {
        stk::mesh::Entity const node = m_bulk_data.get_entity( MetaData::NODE_RANK , elem_node[i] );
        m_bulk_data.change_entity_parts(node, add_parts);

        ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        unsigned nx = 0, ny = 0, nz = 0;
        node_x_y_z(elem_node[i], nx, ny, nz);

        Scalar * data = m_bulk_data.field_data( m_coord_field , node );

        data[0] = nx ;
        data[1] = ny ;
        data[2] = nz ;
      }
    }
  }

  m_bulk_data.modification_end();
}

} // fixtures
} // mesh
} // stk
