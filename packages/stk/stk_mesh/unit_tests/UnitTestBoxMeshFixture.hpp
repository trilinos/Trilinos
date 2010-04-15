/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_UNITTEST_BOX_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_BOX_MESH_FIXTURE_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/DataTraits.hpp>

class BoxMeshFixture
{
public:
  typedef int Scalar ;
  typedef stk::mesh::Field<Scalar, stk::mesh::Cartesian>     CoordFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::ElementNode>   CoordGatherFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::QuadratureTag> QuadFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::BasisTag>      BasisFieldType;

  BoxMeshFixture(stk::ParallelMachine pm);

  void fill_mesh();

  ~BoxMeshFixture();

  stk::mesh::MetaData    m_meta_data ;
  stk::mesh::BulkData    m_bulk_data ;
  stk::mesh::Part      & m_elem_block ;
  CoordFieldType       & m_coord_field ;
  CoordGatherFieldType & m_coord_gather_field ;
  QuadFieldType        & m_quad_field ;
  BasisFieldType       & m_basis_field ;

  //  Box of 8 elements:
  //  Up to 8 processes, so each process owns at least one element.
  //  Then the automatic ghosting will fully populate the nodes and elements.

  stk::mesh::Entity * m_nodes[3][3][3] ;
  stk::mesh::Entity * m_elems[2][2][2] ;
  stk::mesh::EntityId m_node_id[3][3][3] ;
  Scalar              m_node_coord[3][3][3][3] ;

private:
  BoxMeshFixture( const BoxMeshFixture & );
  BoxMeshFixture & operator = ( const BoxMeshFixture & );
};

#endif

