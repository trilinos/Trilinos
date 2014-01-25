/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/fixtures/QuadFixture.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                          unsigned nx , unsigned ny,
                          const std::vector<std::string>& rank_names )
  : m_spatial_dimension(2),
    m_meta( m_spatial_dimension, rank_names ),
    m_bulk_data( m_meta, pm ),
    m_quad_part( declare_part<shards::Quadrilateral<4> >(m_meta, "quad_part" ) ),
    m_coord_field( m_meta.declare_field<CoordFieldType>(stk::topology::NODE_RANK, "Coordinates") ),
    m_nx( nx ),
    m_ny( ny )
{
  typedef shards::Quadrilateral<4> Quad4 ;

  //put coord-field on all nodes:
  put_field(
      m_coord_field,
      m_meta.universal_part(),
      m_spatial_dimension
      );
}

void QuadFixture::node_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= 1;

  x = entity_id % (m_nx+1);
  entity_id /= (m_nx+1);

  y = entity_id;
}

void QuadFixture::elem_x_y( EntityId entity_id, unsigned &x , unsigned &y ) const
{
  entity_id -= 1;

  x = entity_id % m_nx;
  entity_id /= m_nx;

  y = entity_id;
}


void QuadFixture::generate_mesh()
{
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = m_bulk_data.parallel_size();
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned num_elems = m_nx * m_ny;

  const EntityId beg_elem = 1 + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = 1 + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void QuadFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor)
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

    std::vector<EntityId>::const_iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::const_iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_x_y(entity_id, ix, iy);

      stk::mesh::EntityId elem_nodes[4] ;

      elem_nodes[0] = node_id( ix   , iy );
      elem_nodes[1] = node_id( ix+1 , iy );
      elem_nodes[2] = node_id( ix+1 , iy+1 );
      elem_nodes[3] = node_id( ix   , iy+1 );

      stk::mesh::declare_element( m_bulk_data, m_quad_part, elem_id( ix , iy ) , elem_nodes);
      for (unsigned i = 0; i<4; ++i) {
        stk::mesh::Entity const node = m_bulk_data.get_entity( stk::topology::NODE_RANK , elem_nodes[i] );

        ThrowRequireMsg( m_bulk_data.is_valid(node),
          "This process should know about the nodes that make up its element");

        // Compute and assign coordinates to the node
        unsigned nx = 0, ny = 0;
        node_x_y(elem_nodes[i], nx, ny);

        Scalar * data = stk::mesh::field_data( m_coord_field , node );

        data[0] = nx ;
        data[1] = ny ;
      }
    }
  }

  m_bulk_data.modification_end();
}


} // fixtures
} // mesh
} // stk
