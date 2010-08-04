/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

class HexFixture
{
public:

  typedef double Scalar ;
  typedef Field<Scalar, Cartesian>     CoordFieldType;
  typedef Field<Scalar*,ElementNode>   CoordGatherFieldType;

  HexFixture(stk::ParallelMachine pm, unsigned nx, unsigned ny, unsigned nz);

  ~HexFixture();

  const unsigned NX;
  const unsigned NY;
  const unsigned NZ;

  MetaData  meta_data;
  BulkData  bulk_data;

  Part    & hex_part;

  CoordFieldType       & coord_field ;
  CoordGatherFieldType & coord_gather_field ;

  EntityId node_id( unsigned ix , unsigned iy , unsigned iz ) const  {
    return 1 + ix + ( NX + 1 ) * ( iy + ( NY + 1 ) * iz );
  }

  EntityId elem_id( unsigned ix , unsigned iy , unsigned iz ) const  {
    return 1 + ix + NX * ( iy + NY * iz );
  }

  Entity * node( unsigned ix , unsigned iy , unsigned iz ) const {
    return bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy,iz) );
  }

  Entity * elem( unsigned ix , unsigned iy , unsigned iz ) const {
    return bulk_data.get_entity( stk::mesh::Element , elem_id(ix,iy,iz) );
  }

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor );

  void generate_mesh();
private:

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator = (const HexFixture &);
};


} // fixtures
} // mesh
} // stk
#endif
