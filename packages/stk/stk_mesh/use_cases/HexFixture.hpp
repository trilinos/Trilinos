/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_UNITTEST_HEX_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_HEX_MESH_FIXTURE_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/DataTraits.hpp>

template < unsigned NX, unsigned NY, unsigned NZ>
class HexFixture
{
public:

  typedef double Scalar ;
  typedef stk::mesh::Field<Scalar, stk::mesh::Cartesian>     CoordFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::ElementNode>   CoordGatherFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::QuadratureTag> QuadFieldType;
  typedef stk::mesh::Field<Scalar*,stk::mesh::BasisTag>      BasisFieldType;


  HexFixture(stk::ParallelMachine pm);

  ~HexFixture() {};

  stk::mesh::MetaData  m_meta_data;
  stk::mesh::BulkData  m_bulk_data;
  stk::mesh::Part    & m_hex_part;
  stk::mesh::Part    & m_skin_part;

  CoordFieldType       & m_coord_field ;
  CoordGatherFieldType & m_coord_gather_field ;
  QuadFieldType        & m_quad_field ;
  BasisFieldType       & m_basis_field ;



  stk::mesh::EntityId m_elems_id[NZ][NY][NX] ;
  stk::mesh::EntityId m_node_id[NZ+1][NY+1][NX+1] ;
  Scalar              m_node_coord[NZ+1][NY+1][NX+1][3] ;

private:
  void generate_mesh();

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator = (const HexFixture &);
};

typedef HexFixture<3,3,3> HexFixture3x3x3;


template < unsigned NX, unsigned NY, unsigned NZ>
HexFixture<NX,NY,NZ>::HexFixture(stk::ParallelMachine pm)
  : m_meta_data( stk::mesh::fem_entity_rank_names() )
  , m_bulk_data( m_meta_data , pm )
  , m_hex_part( m_meta_data.declare_part("hex_part", stk::mesh::Element) )
  , m_skin_part( m_meta_data.declare_part("skin_part"))
  , m_coord_field( m_meta_data.declare_field<CoordFieldType>("Coordinates") )
  , m_coord_gather_field(
        m_meta_data.declare_field<CoordGatherFieldType>("GatherCoordinates")
      )
  , m_quad_field( m_meta_data.declare_field<QuadFieldType>("Quad") )
  , m_basis_field( m_meta_data.declare_field<BasisFieldType>("Basis") )
{
  typedef shards::Hexahedron<8> Hex8 ;
  enum { SpatialDim = 3 };
  enum { NodesPerElem = Hex8::node_count };

  // Set topology of the element block part
  stk::mesh::set_cell_topology<shards::Hexahedron<8> >(m_hex_part);

  //put coord-field on all nodes:
  stk::mesh::put_field(
      m_coord_field,
      stk::mesh::Node,
      m_meta_data.universal_part(),
      SpatialDim
      );

  //put coord-gather-field on all elements:
  stk::mesh::put_field(
      m_coord_gather_field,
      stk::mesh::Element,
      m_meta_data.universal_part(),
      NodesPerElem
      );

  // Field relation so coord-gather-field on elements points
  // to coord-field of the element's nodes
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_stencil<shards::Hexahedron<8> >, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_stencil<void>, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, stk::mesh::element_node_lock_stencil<void>, m_coord_field);

  m_meta_data.commit();

  generate_mesh();
}

template < unsigned NX, unsigned NY, unsigned NZ>
void HexFixture<NX,NY,NZ>::generate_mesh() {

  enum { num_elems = NX*NY*NZ };
  const unsigned p_rank = m_bulk_data.parallel_rank();
  const unsigned p_size = m_bulk_data.parallel_size();


  //assign node id's and coordinates
  for ( unsigned iz = 0 ; iz < NZ+1 ; ++iz ) {
  for ( unsigned iy = 0 ; iy < NY+1 ; ++iy ) {
  for ( unsigned ix = 0 ; ix < NX+1 ; ++ix ) {
    m_node_id[iz][iy][ix] = 1 + ix + (NX+1) * ( iy + (NY+1) * iz );
    m_node_coord[iz][iy][ix][0] = ix -1.5 ;
    m_node_coord[iz][iy][ix][1] = iy -1.5 ;
    m_node_coord[iz][iy][ix][2] = iz -1.5 ;
  }
  }
  }

  for ( unsigned iz = 0 ; iz < NZ ; ++iz ) {
  for ( unsigned iy = 0 ; iy < NY; ++iy ) {
  for ( unsigned ix = 0 ; ix < NX ; ++ix ) {
    m_elems_id[iz][iy][ix] = 1 + ix + NX * ( iy + NY * iz );
  }
  }
  }

  const unsigned beg_elem = ( num_elems * p_rank ) / p_size ;
  const unsigned end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;

  m_bulk_data.modification_begin();

  unsigned elem = 0 ;
  stk::mesh::EntityId face_id = 0;

  for ( unsigned iz = 0 ; iz < NZ ; ++iz ) {
  for ( unsigned iy = 0 ; iy < NY ; ++iy ) {
  for ( unsigned ix = 0 ; ix < NX ; ++ix , ++elem ) {
    if ( beg_elem <= elem && elem < end_elem ) {
      stk::mesh::EntityId elem_id = m_elems_id[iz][iy][ix];
      stk::mesh::EntityId elem_node[8] ;
      Scalar * node_coord[8] ;
      face_id = elem_id;

      elem_node[0] = m_node_id[iz  ][iy  ][ix  ] ;
      elem_node[1] = m_node_id[iz  ][iy  ][ix+1] ;
      elem_node[2] = m_node_id[iz+1][iy  ][ix+1] ;
      elem_node[3] = m_node_id[iz+1][iy  ][ix  ] ;
      elem_node[4] = m_node_id[iz  ][iy+1][ix  ] ;
      elem_node[5] = m_node_id[iz  ][iy+1][ix+1] ;
      elem_node[6] = m_node_id[iz+1][iy+1][ix+1] ;
      elem_node[7] = m_node_id[iz+1][iy+1][ix  ] ;

      node_coord[0] = m_node_coord[iz  ][iy  ][ix  ] ;
      node_coord[1] = m_node_coord[iz  ][iy  ][ix+1] ;
      node_coord[2] = m_node_coord[iz+1][iy  ][ix+1] ;
      node_coord[3] = m_node_coord[iz+1][iy  ][ix  ] ;
      node_coord[4] = m_node_coord[iz  ][iy+1][ix  ] ;
      node_coord[5] = m_node_coord[iz  ][iy+1][ix+1] ;
      node_coord[6] = m_node_coord[iz+1][iy+1][ix+1] ;
      node_coord[7] = m_node_coord[iz+1][iy+1][ix  ] ;

      stk::mesh::Entity& element = stk::mesh::declare_element( m_bulk_data, m_hex_part, elem_id, elem_node);

      stk::mesh::PairIterRelation rel = element.relations(stk::mesh::Node);

      for(unsigned i=0; i< 8 ; ++i ) {
        stk::mesh::Entity& node = *rel[i].entity();
        Scalar* data = stk::mesh::field_data( m_coord_field, node);
        data[0] = node_coord[i][0];
        data[1] = node_coord[i][1];
        data[2] = node_coord[i][2];
      }
    }
  }
  }
  }

  m_bulk_data.modification_end();

  skin_mesh(m_bulk_data, stk::mesh::Element, &m_skin_part);

}
#endif
