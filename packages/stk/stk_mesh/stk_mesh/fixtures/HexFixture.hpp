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

namespace stk {
namespace mesh {
namespace fixtures {

template < unsigned NX, unsigned NY, unsigned NZ>
class HexFixture
{
public:

  typedef double Scalar ;
  typedef Field<Scalar, Cartesian>     CoordFieldType;
  typedef Field<Scalar*,ElementNode>   CoordGatherFieldType;
  typedef Field<Scalar*,QuadratureTag> QuadFieldType;
  typedef Field<Scalar*,BasisTag>      BasisFieldType;


  HexFixture(stk::ParallelMachine pm, const std::vector<unsigned> * element_ids_this_processor = NULL);

  ~HexFixture() {};

  MetaData  m_meta_data;
  BulkData  m_bulk_data;
  Part    & m_hex_part;
  Part    & m_skin_part;

  CoordFieldType       & m_coord_field ;
  CoordGatherFieldType & m_coord_gather_field ;
  QuadFieldType        & m_quad_field ;
  BasisFieldType       & m_basis_field ;



  EntityId m_elems_id[NZ][NY][NX] ;
  EntityId m_node_id[NZ+1][NY+1][NX+1] ;
  Scalar              m_node_coord[NZ+1][NY+1][NX+1][3] ;

private:
  void generate_mesh(const std::vector<unsigned> * element_ids_this_processor);

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator = (const HexFixture &);
};

typedef HexFixture<3,3,3> HexFixture3x3x3;


template < unsigned NX, unsigned NY, unsigned NZ>
HexFixture<NX,NY,NZ>::HexFixture(stk::ParallelMachine pm, const std::vector<unsigned> * element_ids_this_processor)
  : m_meta_data( fem_entity_rank_names() )
  , m_bulk_data( m_meta_data , pm )
  , m_hex_part( m_meta_data.declare_part("hex_part", Element) )
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
  set_cell_topology<shards::Hexahedron<8> >(m_hex_part);

  //put coord-field on all nodes:
  put_field(
      m_coord_field,
      Node,
      m_meta_data.universal_part(),
      SpatialDim
      );

  //put coord-gather-field on all elements:
  put_field(
      m_coord_gather_field,
      Element,
      m_meta_data.universal_part(),
      NodesPerElem
      );

  // Field relation so coord-gather-field on elements points
  // to coord-field of the element's nodes
  m_meta_data.declare_field_relation( m_coord_gather_field, element_node_stencil<shards::Hexahedron<8> >, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, element_node_stencil<void>, m_coord_field);
  m_meta_data.declare_field_relation( m_coord_gather_field, element_node_lock_stencil<void>, m_coord_field);

  m_meta_data.commit();

  generate_mesh(element_ids_this_processor);
}

template < unsigned NX, unsigned NY, unsigned NZ>
void HexFixture<NX,NY,NZ>::generate_mesh(const std::vector<unsigned> * element_ids_this_processor) {

  enum { num_elems = NX*NY*NZ };


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

  std::map<unsigned,std::vector<unsigned> > elems_id_map;
  for ( unsigned iz = 0 ; iz < NZ ; ++iz ) {
  for ( unsigned iy = 0 ; iy < NY; ++iy ) {
  for ( unsigned ix = 0 ; ix < NX ; ++ix ) {
    const unsigned id = 1 + ix + NX * ( iy + NY * iz );
    // DEBUG:
    //std::cout << "id = " << id << " (" << ix << "," << iy << "," << iz << ")" << std::endl;;
    m_elems_id[iz][iy][ix] = id;
    std::vector<unsigned> id_index;
    id_index.push_back(ix);
    id_index.push_back(iy);
    id_index.push_back(iz);
    elems_id_map[id] = id_index;
  }
  }
  }


  std::vector<unsigned> local_element_ids_this_processor;
  if (element_ids_this_processor != NULL) {
    local_element_ids_this_processor = *element_ids_this_processor;
  }
  else {
    const unsigned p_rank = m_bulk_data.parallel_rank();
    const unsigned p_size = m_bulk_data.parallel_size();
    const unsigned beg_elem = ( num_elems * p_rank ) / p_size ;
    const unsigned end_elem = ( num_elems * ( p_rank + 1 ) ) / p_size ;
    for (unsigned elem_index = beg_elem+1 ; elem_index <= end_elem ; ++elem_index) {
      local_element_ids_this_processor.push_back(elem_index);
    }
  }
  // DEBUG:
  //for (unsigned i=0 ; i<local_element_ids_this_processor.size() ; ++i) {
  //  std::cout << "local_element_ids_this_processor["<<i<<"] = " << local_element_ids_this_processor[i] << std::endl;
  //}

  m_bulk_data.modification_begin();

  const std::vector<unsigned>::iterator elem_it_begin = local_element_ids_this_processor.begin();
  std::vector<unsigned>::iterator elem_it = elem_it_begin;
  const std::vector<unsigned>::iterator elem_it_end = local_element_ids_this_processor.end();

  for ( ; elem_it != elem_it_end ; elem_it++ ) {
    EntityId elem_id = *elem_it;
    // figure out what ix,iy,iz is for this element_id
    const unsigned ix = elems_id_map[elem_id][0];
    const unsigned iy = elems_id_map[elem_id][1];
    const unsigned iz = elems_id_map[elem_id][2];

    EntityId elem_node[8] ;
    Scalar * node_coord[8] ;

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

    Entity& element = declare_element( m_bulk_data, m_hex_part, elem_id, elem_node);

    PairIterRelation rel = element.relations(Node);

    for(unsigned i=0; i< 8 ; ++i ) {
      Entity& node = *rel[i].entity();
      Scalar* data = field_data( m_coord_field, node);
      data[0] = node_coord[i][0];
      data[1] = node_coord[i][1];
      data[2] = node_coord[i][2];
    }

  }

  m_bulk_data.modification_end();

  skin_mesh(m_bulk_data, Element, &m_skin_part);

}

} // fixtures
} // mesh
} // stk
#endif
