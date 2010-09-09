/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>

#include <stk_mesh/fixtures/QuadFixture.hpp>

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

QuadFixture::~QuadFixture()
{}

QuadFixture::QuadFixture( stk::ParallelMachine pm ,
                        unsigned nx , unsigned ny )
  : spatial_dimension(2),
    meta_data( TopologicalMetaData::entity_rank_names(spatial_dimension) ),
    bulk_data( meta_data , pm ),
    top_data( meta_data, spatial_dimension ),
    quad_part( top_data.declare_part<shards::Quadrilateral<4> >("quad_part" ) ),
    coord_field( meta_data.declare_field<CoordFieldType>("Coordinates") ),
    coord_gather_field( meta_data.declare_field<CoordGatherFieldType>("GatherCoordinates") ),
    NX( nx ),
    NY( ny )
{
  typedef shards::Quadrilateral<4> Quad4 ;
  enum { SpatialDim = 2 };
  enum { NodesPerElem = Quad4::node_count };

  //put coord-field on all nodes:
  put_field(
      coord_field,
      top_data.node_rank,
      meta_data.universal_part(),
      SpatialDim
      );

  //put coord-gather-field on all elements:
  put_field(
      coord_gather_field,
      top_data.element_rank,
      meta_data.universal_part(),
      NodesPerElem
      );

  // Field relation so coord-gather-field on elements points
  // to coord-field of the element's nodes
  meta_data.declare_field_relation( coord_gather_field, element_node_stencil<Quad4>, coord_field);

}

void QuadFixture::generate_mesh() {
  std::vector<EntityId> element_ids_on_this_processor;

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();
  const unsigned num_elems = NX * NY;

  const EntityId beg_elem = 1 + ( num_elems * p_rank ) / p_size ;
  const EntityId end_elem = 1 + ( num_elems * ( p_rank + 1 ) ) / p_size ;

  for ( EntityId i = beg_elem; i != end_elem; ++i) {
    element_ids_on_this_processor.push_back(i);
  }

  generate_mesh(element_ids_on_this_processor);
}

void QuadFixture::generate_mesh(std::vector<EntityId> & element_ids_on_this_processor) {

  {
    //sort and unique the input elements
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();

    std::sort( ib, ie);
    ib = std::unique( ib, ie);
    element_ids_on_this_processor.erase(ib, ie);
  }

  bulk_data.modification_begin();

  {
    std::vector<EntityId>::iterator ib = element_ids_on_this_processor.begin();
    const std::vector<EntityId>::iterator ie = element_ids_on_this_processor.end();
    for (; ib != ie; ++ib) {
      EntityId entity_id = *ib;
      unsigned ix = 0, iy = 0;
      elem_ix_iy(entity_id, ix, iy);

      stk::mesh::EntityId elem_node[4] ;

      elem_node[0] = node_id( ix   , iy );
      elem_node[1] = node_id( ix+1 , iy );
      elem_node[2] = node_id( ix+1 , iy+1 );
      elem_node[3] = node_id( ix   , iy+1 );

      stk::mesh::declare_element( bulk_data, quad_part, elem_id( ix , iy ) , elem_node);
      for (unsigned i = 0; i<4; ++i) {
        stk::mesh::Entity * const node =
          bulk_data.get_entity( top_data.node_rank , elem_node[i] );

        if ( node != NULL) {

          unsigned nx = 0, ny = 0;
          node_ix_iy(elem_node[i], nx, ny);

          Scalar * data = stk::mesh::field_data( coord_field , *node );

          data[0] = nx ;
          data[1] = ny ;
        }
      }
    }
  }

  bulk_data.modification_end();

}


} // fixtures
} // mesh
} // stk
