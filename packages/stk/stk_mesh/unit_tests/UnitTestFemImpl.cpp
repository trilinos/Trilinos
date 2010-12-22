/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>

namespace {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * Really just a meta-data, fem-object pair. They share the same
 * spatial dimension and the meta-data is associated with the fem-object.
 */
class TestFixture {
public:
  typedef stk::mesh::MetaData           BaseMetaData ;
  typedef stk::mesh::DefaultFEM         TopoMetaData ;

  BaseMetaData m_meta_data ;
  TopoMetaData m_fem ;

  TestFixture(unsigned spatial_dimension)
    : m_meta_data(stk::mesh::fem::entity_rank_names(spatial_dimension)),
      m_fem(m_meta_data, spatial_dimension)
  {}

  stk::mesh::EntityRank get_entity_rank(const CellTopologyData * top) const {
    return m_fem.get_entity_rank(top);
  }
};

STKUNIT_UNIT_TEST(UnitTestDefaultFEM, entity_rank)
{
  TestFixture test1D(1u);
  TestFixture test2D(2u);
  TestFixture test3D(3u);
  {
    stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(test1D.m_meta_data);

    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::node_rank(fem), 0u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::edge_rank(fem), stk::mesh::fem::INVALID_RANK);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::face_rank(fem), stk::mesh::fem::INVALID_RANK);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::side_rank(fem), 0u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::element_rank(fem), 1u);
//    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::patch_rank(fem), 2u);
  }

  {
    stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(test2D.m_meta_data);

    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::NODE_RANK, 0u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::edge_rank(fem), 1u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::face_rank(fem), stk::mesh::fem::INVALID_RANK);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::side_rank(fem), 1u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::element_rank(fem), 2u);
//    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::patch_rank(fem), 3u);
  }

  {
    stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(test3D.m_meta_data);

    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::node_rank(fem), 0u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::edge_rank(fem), 1u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::side_rank(fem), 2u);
    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::element_rank(fem), 3u);
//    STKUNIT_EXPECT_EQUAL(stk::mesh::fem::patch_rank(fem), 4u);
  }
}


STKUNIT_UNIT_TEST(UnitTestDefaultFEM, cellTopology)
{
  const CellTopologyData * node  = shards::getCellTopologyData< shards::Node >();

  const CellTopologyData * line2 = shards::getCellTopologyData< shards::Line<2> >();
  const CellTopologyData * line3 = shards::getCellTopologyData< shards::Line<3> >();

  const CellTopologyData * tri3 = shards::getCellTopologyData< shards::Triangle<3> >();
  const CellTopologyData * tri6 = shards::getCellTopologyData< shards::Triangle<6> >();
  const CellTopologyData * tri4 = shards::getCellTopologyData< shards::Triangle<4> >();

  const CellTopologyData * quad4 = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  const CellTopologyData * quad8 = shards::getCellTopologyData< shards::Quadrilateral<8> >();
  const CellTopologyData * quad9 = shards::getCellTopologyData< shards::Quadrilateral<9> >();

  const CellTopologyData * tet4  = shards::getCellTopologyData< shards::Tetrahedron<4> >();
  const CellTopologyData * tet10 = shards::getCellTopologyData< shards::Tetrahedron<10> >();
  const CellTopologyData * tet8  = shards::getCellTopologyData< shards::Tetrahedron<8> >();

  const CellTopologyData * pyr5  = shards::getCellTopologyData< shards::Pyramid<5> >();
  const CellTopologyData * pyr13 = shards::getCellTopologyData< shards::Pyramid<13> >();
  const CellTopologyData * pyr14 = shards::getCellTopologyData< shards::Pyramid<14> >();

  const CellTopologyData * wedge6  = shards::getCellTopologyData< shards::Wedge<6> >();
  const CellTopologyData * wedge15 = shards::getCellTopologyData< shards::Wedge<15> >();
  const CellTopologyData * wedge18 = shards::getCellTopologyData< shards::Wedge<18> >();

  const CellTopologyData * hex8  = shards::getCellTopologyData< shards::Hexahedron<8> >();
  const CellTopologyData * hex20 = shards::getCellTopologyData< shards::Hexahedron<20> >();
  const CellTopologyData * hex27 = shards::getCellTopologyData< shards::Hexahedron<27> >();


  const CellTopologyData * particle = shards::getCellTopologyData< shards::Particle >();

  const CellTopologyData * beam2 = shards::getCellTopologyData< shards::Beam<2> >();
  const CellTopologyData * beam3 = shards::getCellTopologyData< shards::Beam<3> >();

  const CellTopologyData * shellLine2 = shards::getCellTopologyData< shards::ShellLine<2> >();
  const CellTopologyData * shellLine3 = shards::getCellTopologyData< shards::ShellLine<3> >();

  const CellTopologyData * shellTri3 = shards::getCellTopologyData< shards::ShellTriangle<3> >();
  const CellTopologyData * shellTri6 = shards::getCellTopologyData< shards::ShellTriangle<6> >();

  const CellTopologyData * shellQuad4 = shards::getCellTopologyData< shards::ShellQuadrilateral<4> >();
  const CellTopologyData * shellQuad8 = shards::getCellTopologyData< shards::ShellQuadrilateral<8> >();
  const CellTopologyData * shellQuad9 = shards::getCellTopologyData< shards::ShellQuadrilateral<9> >();


  for (unsigned spatial_dimension = 1 ; spatial_dimension < 4 ; ++spatial_dimension) {
    TestFixture test(spatial_dimension);

    STKUNIT_EXPECT_EQUAL(0u, test.get_entity_rank(node));
    STKUNIT_EXPECT_EQUAL(1u, test.get_entity_rank(line2));
    STKUNIT_EXPECT_EQUAL(1u, test.get_entity_rank(line3));
    STKUNIT_EXPECT_EQUAL(spatial_dimension, test.get_entity_rank(particle));
    if (1 < spatial_dimension) {
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(tri3));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(tri6));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(tri4));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(quad4));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(quad8));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(quad9));

      STKUNIT_EXPECT_EQUAL(spatial_dimension, test.get_entity_rank(beam2));
      STKUNIT_EXPECT_EQUAL(spatial_dimension, test.get_entity_rank(beam3));
    }
    else {
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tri3));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tri4));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tri6));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(quad4));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(quad8));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(quad9));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(beam2));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(beam3));
    }

    if (2 == spatial_dimension) {
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(shellLine2));
      STKUNIT_EXPECT_EQUAL(2u, test.get_entity_rank(shellLine3));
    }
    else {
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellLine2));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellLine3));
    }

    if (2 < spatial_dimension) {
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(tet4));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(tet10));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(tet8));

      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(pyr5));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(pyr13));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(pyr14));

      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(wedge6));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(wedge15));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(wedge18));

      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(hex8));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(hex20));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(hex27));

      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(shellTri3));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(shellTri6));

      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(shellQuad4));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(shellQuad8));
      STKUNIT_EXPECT_EQUAL(3u, test.get_entity_rank(shellQuad9));
    }
    else {
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tet4));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tet10));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(tet8));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(pyr5));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(pyr13));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(pyr14));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(wedge6));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(wedge15));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(wedge18));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(hex8));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(hex20));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(hex27));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellTri3));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellTri6));

      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellQuad4));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellQuad8));
      STKUNIT_EXPECT_EQUAL(stk::mesh::fem::INVALID_RANK, test.get_entity_rank(shellQuad9));
    }
  }
}


STKUNIT_UNIT_TEST(UnitTestDefaultFEM, part_supersets)
{
  const unsigned spatial_dimension = 3 ;

  TestFixture test(spatial_dimension);

  const CellTopologyData * const top_hex8 = shards::getCellTopologyData< shards::Hexahedron<8> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  stk::mesh::Part & hex8 = stk::mesh::declare_part(test.m_meta_data, top_hex8->name, top_hex8);
  stk::mesh::Part & tet4 = stk::mesh::declare_part(test.m_meta_data, top_tet4->name, top_tet4);

  // stk::mesh::fem::set_cell_topology(hex8, get_entity_rank(top_hex8), top_hex8);

  stk::mesh::Part & part_hex_1 = declare_part(test.m_meta_data,  "block_1", spatial_dimension);
  stk::mesh::Part & part_hex_2 = declare_part(test.m_meta_data,  "block_2", spatial_dimension);
  stk::mesh::Part & part_hex_3 = declare_part(test.m_meta_data,  "block_3", spatial_dimension);
  stk::mesh::Part & part_tet_1 = declare_part(test.m_meta_data,  "block_4", spatial_dimension);
  stk::mesh::Part & part_tet_2 = declare_part(test.m_meta_data,  "block_5", spatial_dimension);

  test.m_meta_data.declare_part_subset(hex8, part_hex_1);
  test.m_meta_data.declare_part_subset(hex8, part_hex_2);
  test.m_meta_data.declare_part_subset(hex8, part_hex_3);

  test.m_meta_data.declare_part_subset(tet4, part_tet_1);
  test.m_meta_data.declare_part_subset(tet4, part_tet_2);

  STKUNIT_EXPECT_EQUAL(top_hex8, stk::mesh::fem::get_cell_topology(hex8).getCellTopologyData());
  STKUNIT_EXPECT_EQUAL(top_hex8, stk::mesh::fem::get_cell_topology(part_hex_1).getCellTopologyData());
  STKUNIT_EXPECT_EQUAL(top_hex8, stk::mesh::fem::get_cell_topology(part_hex_2).getCellTopologyData());
  STKUNIT_EXPECT_EQUAL(top_hex8, stk::mesh::fem::get_cell_topology(part_hex_3).getCellTopologyData());

  STKUNIT_EXPECT_EQUAL(top_tet4, stk::mesh::fem::get_cell_topology(tet4).getCellTopologyData());
  STKUNIT_EXPECT_EQUAL(top_tet4, stk::mesh::fem::get_cell_topology(part_tet_1).getCellTopologyData());
  STKUNIT_EXPECT_EQUAL(top_tet4, stk::mesh::fem::get_cell_topology(part_tet_2).getCellTopologyData());

}


STKUNIT_UNIT_TEST(UnitTestDefaultFEM, expected_throws)
{
  const unsigned spatial_dimension = 2 ;

  TestFixture test(spatial_dimension);

  const CellTopologyData * const top_quad4 = shards::getCellTopologyData< shards::Quadrilateral<4> >();
  const CellTopologyData * const top_tri3 = shards::getCellTopologyData< shards::Triangle<3> >();
  const CellTopologyData * const top_tet4 = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  STKUNIT_ASSERT_THROW(test.m_fem.register_cell_topology(top_quad4, 0), std::runtime_error);
  STKUNIT_ASSERT_THROW(test.m_fem.register_cell_topology(top_quad4, 1), std::runtime_error);
  STKUNIT_ASSERT_THROW(test.m_fem.register_cell_topology(top_quad4, 3), std::runtime_error);

  STKUNIT_ASSERT_THROW(test.m_fem.register_cell_topology(top_tri3, 3), std::runtime_error);
  STKUNIT_ASSERT_THROW(test.m_fem.register_cell_topology(top_tet4, 3), std::runtime_error);

  STKUNIT_ASSERT_THROW(stk::mesh::fem::set_spatial_dimension(test.m_meta_data, 3), std::runtime_error);
}

} // namespace
