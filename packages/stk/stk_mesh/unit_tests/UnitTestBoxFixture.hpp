/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_UNITTEST_BOX_FIXTURE_HPP
#define STK_MESH_UNITTEST_BOX_FIXTURE_HPP

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

class BoxFixture
{
public:
  typedef int Scalar ;
  typedef stk::mesh::Field<Scalar, stk::mesh::Cartesian> CoordFieldType;

  BoxFixture( stk::ParallelMachine pm,
              unsigned elem_nx , unsigned elem_ny , unsigned elem_nz );

  ~BoxFixture();

  stk::mesh::MetaData m_meta_data ;
  stk::mesh::BulkData m_bulk_data ;
  stk::mesh::Part   & m_elem_block ;
  CoordFieldType    & m_coord_field ;
  const unsigned      m_elem_nx ;
  const unsigned      m_elem_ny ;
  const unsigned      m_elem_nz ;

  stk::mesh::EntityId node_id( unsigned ix , unsigned iy , unsigned iz ) const
    { return 1 + ix + ( m_elem_nx + 1 ) * ( iy + ( m_elem_ny + 1 ) * iz ); }

  stk::mesh::EntityId elem_id( unsigned ix , unsigned iy , unsigned iz ) const
    { return 1 + ix + m_elem_nx * ( iy + m_elem_ny * iz ); }

  stk::mesh::Entity * node( unsigned ix , unsigned iy , unsigned iz ) const
    { return m_bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy,iz) ); }

  stk::mesh::Entity * elem( unsigned ix , unsigned iy , unsigned iz ) const
    { return m_bulk_data.get_entity( stk::mesh::Element , elem_id(ix,iy,iz) ); }

private:

  BoxFixture( const BoxFixture & );
  BoxFixture & operator = ( const BoxFixture & );
};

class QuadFixture
{
public:
  typedef int Scalar ;
  typedef stk::mesh::Field<Scalar, stk::mesh::Cartesian> CoordFieldType;

  QuadFixture( stk::ParallelMachine pm, unsigned elem_nx , unsigned elem_ny );

  ~QuadFixture();

  stk::mesh::MetaData m_meta_data ;
  stk::mesh::BulkData m_bulk_data ;
  stk::mesh::Part   & m_elem_block ;
  CoordFieldType    & m_coord_field ;
  const unsigned      m_elem_nx ;
  const unsigned      m_elem_ny ;

  stk::mesh::EntityId node_id( unsigned ix , unsigned iy ) const
    { return 1 + ix + ( m_elem_nx + 1 ) * iy ; }

  stk::mesh::EntityId elem_id( unsigned ix , unsigned iy ) const
    { return 1 + ix + m_elem_nx * iy ; }

  stk::mesh::Entity * node( unsigned ix , unsigned iy ) const
    { return m_bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy) ); }

  stk::mesh::Entity * elem( unsigned ix , unsigned iy ) const
    { return m_bulk_data.get_entity( stk::mesh::Element , elem_id(ix,iy) ); }

private:

  QuadFixture( const QuadFixture & );
  QuadFixture & operator = ( const QuadFixture & );
};

#endif

