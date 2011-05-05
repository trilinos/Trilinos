/*------------------------------------------------------------------------*/
/*         _        Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Transaction.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_mesh/baseImpl/BucketImpl.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <Shards_BasicTopologies.hpp>

using stk::ParallelMachine;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;
using stk::mesh::EntitySideComponent;
using stk::mesh::PairIterRelation;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::fem::FEMMetaData;

class TopologyHelpersTestingFixture
{
 public:
  TopologyHelpersTestingFixture(ParallelMachine pm);
  ~TopologyHelpersTestingFixture() {}

  const int spatial_dimension;
  FEMMetaData meta;
  BulkData bulk;
  const EntityRank element_rank;
  const EntityRank side_rank;
  Part & generic_element_part;
  Part & element_tet_part;
  Part & element_wedge_part;
  Part & generic_face_part;
  Part & another_generic_face_part;
  Part & face_quad_part;
  Part & another_generic_element_part;

  EntityId nextEntityId()
  { return psize*(++entity_id)+prank; }

  Entity & create_entity( EntityRank rank, Part& part_membership)
  {
    PartVector part_intersection;
    part_intersection.push_back ( &part_membership );
    return bulk.declare_entity(rank, nextEntityId(), part_intersection);
  }

 private:
  EntityId entity_id;
  const int psize;
  const int prank;
};

TopologyHelpersTestingFixture::TopologyHelpersTestingFixture(ParallelMachine pm)
  : spatial_dimension( 3 )
  , meta( spatial_dimension )
  , bulk( FEMMetaData::get_meta_data(meta), pm, 100 )
  , element_rank( meta.element_rank())
  , side_rank( meta.side_rank())
  , generic_element_part( meta.declare_part("another part", element_rank ) )
  , element_tet_part( stk::mesh::fem::declare_part<shards::Tetrahedron<4> >( meta, "block_left_1" ) )
  , element_wedge_part( stk::mesh::fem::declare_part<shards::Wedge<15> >(meta, "block_left_2" ) )
  , generic_face_part( stk::mesh::fem::declare_part<shards::Quadrilateral<4> >(meta, "A_1" ) )
  , another_generic_face_part( meta.declare_part("A_2", side_rank ) )
  , face_quad_part( meta.declare_part("A_3", side_rank ) )
  , another_generic_element_part( meta.declare_part("B_3", element_rank ) )
  , entity_id(0u)
  , psize(bulk.parallel_size())
  , prank(bulk.parallel_rank())
{
  meta.commit();
}

namespace {

const EntityRank NODE_RANK = FEMMetaData::NODE_RANK;

STKUNIT_UNIT_TEST( testTopologyHelpers, get_cell_topology_based_on_part)
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();
  Entity & elem1  = fix.create_entity( fix.side_rank, fix.generic_face_part );

  PartVector tmp(1);
  tmp[0] = & fix.face_quad_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  STKUNIT_ASSERT_EQUAL( stk::mesh::fem::get_cell_topology(elem1).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  fix.bulk.change_entity_parts ( elem1 , tmp );
  STKUNIT_ASSERT_EQUAL( stk::mesh::fem::get_cell_topology(elem1).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  tmp[0] = & fix.another_generic_face_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  STKUNIT_ASSERT_EQUAL( stk::mesh::fem::get_cell_topology(elem1).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  STKUNIT_ASSERT_NE( stk::mesh::fem::get_cell_topology( elem1).getCellTopologyData() , shards::getCellTopologyData< shards::Wedge<15> >() );

  fix.bulk.modification_end();
}

STKUNIT_UNIT_TEST( testTopologyHelpers, get_cell_topology_multiple_topologies )
{
  // Coverage for get_cell_topology in TopologyHelpers.cpp; (FAILED WITH MULTIPLE LOCAL TOPOLOGIES)
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();
  Entity & elem  = fix.create_entity( fix.element_rank, fix.generic_element_part );
  PartVector add_parts;
  add_parts.push_back( &fix.element_tet_part );
  add_parts.push_back( &fix.element_wedge_part );
  fix.bulk.change_entity_parts( elem, add_parts );
  fix.bulk.modification_end();
  STKUNIT_ASSERT_THROW( stk::mesh::fem::get_cell_topology( elem ).getCellTopologyData(), std::runtime_error );
}

// No longer in the public API
// STKUNIT_UNIT_TEST( testTopologyHelpers, get_adjacent_entities_trivial )
// {
//   // Element, elem2, has NULL topology
//   TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
//
//   if ( 1 == fix.bulk.parallel_size() ) {
//
//     fix.bulk.modification_begin();
//     Entity & elem2  = fix.create_entity( fix.element_rank, fix.generic_element_part );
//     fix.bulk.modification_end();
//
//     std::vector<EntitySideComponent> adjacent_entities;
//     const EntityRank subcell_rank = fix.element_rank;
//     const EntityId subcell_identifier = 1;
//     get_adjacent_entities( elem2 , subcell_rank, subcell_identifier, adjacent_entities);
//     STKUNIT_ASSERT_TRUE( true );
//   }
// }
//
// STKUNIT_UNIT_TEST( testTopologyHelpers, get_adjacent_entities_invalid )
// {
//   TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
//   fix.bulk.modification_begin();
//   Entity & elem3  = fix.create_entity( fix.element_rank , fix.generic_element_part );
//
//   PartVector add_parts;
//   add_parts.push_back( & fix.element_tet_part );
//   fix.bulk.change_entity_parts ( elem3 , add_parts );
//   fix.bulk.modification_end();
//   std::vector<EntitySideComponent> adjacent_entities2;
//   {
//     const EntityRank invalid_subcell_rank = 4;
//     const EntityId valid_subcell_identifier = 0;
//     STKUNIT_ASSERT_THROW(
//         get_adjacent_entities( elem3 , invalid_subcell_rank, valid_subcell_identifier, adjacent_entities2),
//         std::invalid_argument
//         );
//   }
//   {
//     const EntityRank valid_subcell_rank = 1;
//     const EntityId invalid_subcell_identifier = 8;
//     STKUNIT_ASSERT_THROW(
//       get_adjacent_entities( elem3 , valid_subcell_rank, invalid_subcell_identifier, adjacent_entities2),
//       std::invalid_argument
//       );
//   }
// }

STKUNIT_UNIT_TEST( testTopologyHelpers, declare_element_side_no_topology )
{
  // Coverage for declare_element_side - TopologyHelpers.cpp - "Cannot discern element topology"
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();
  Entity & elem4  = fix.create_entity( fix.element_rank , fix.generic_element_part );
  STKUNIT_ASSERT_THROW(
    stk::mesh::fem::declare_element_side( fix.bulk, fix.element_rank, elem4, fix.nextEntityId(), &fix.element_wedge_part ),
    std::runtime_error
      );
  fix.bulk.modification_end();


  {
    EntityId elem_node[4];
    elem_node[0] = 1;
    elem_node[1] = 2;
    elem_node[2] = 3;
    elem_node[3] = 4;
    fix.bulk.modification_begin();
    // Cannot declare an element without a topology defined
    STKUNIT_ASSERT_THROW(
      stk::mesh::fem::declare_element(fix.bulk, fix.generic_element_part, fix.nextEntityId(), elem_node),
        std::runtime_error
        );
    fix.bulk.modification_end();
  }
}

STKUNIT_UNIT_TEST( testTopologyHelpers, declare_element_side_wrong_bulk_data)
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "BulkData for 'elem' and 'side' are different"
  TopologyHelpersTestingFixture fix1(MPI_COMM_WORLD);

  fix1.bulk.modification_begin();

  TopologyHelpersTestingFixture fix2(MPI_COMM_WORLD);
  fix2.bulk.modification_begin();
  Entity & elem4_2  = fix2.create_entity( fix2.element_rank , fix2.generic_element_part );
  fix2.bulk.modification_end();

  STKUNIT_ASSERT_THROW(
    stk::mesh::fem::declare_element_side( fix1.bulk, fix1.element_rank, elem4_2, fix1.nextEntityId(), &fix1.element_wedge_part),
    std::runtime_error
      );
    fix1.bulk.modification_end();
}

STKUNIT_UNIT_TEST( testTopologyHelpers, declare_element_side_no_topology_2 )
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "No element topology found and cell side id exceeds..."
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;
  Entity & element  = stk::mesh::fem::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node);
  const CellTopologyData * const elem_top = stk::mesh::fem::get_cell_topology( element ).getCellTopologyData();
  const EntityId nSideCount = elem_top->side_count + 10 ;
  STKUNIT_ASSERT_THROW(
    stk::mesh::fem::declare_element_side( fix.bulk, fix.nextEntityId(), element, nSideCount, &fix.element_tet_part ),
    std::runtime_error
      );
  fix.bulk.modification_end();
}

STKUNIT_UNIT_TEST( testTopologyHelpers, declare_element_side_full )
{
  // Go all way the through declare_element_side - use new element
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  Entity& element = stk::mesh::fem::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );

  const EntityId zero_side_count = 0;
  Entity& face2 = stk::mesh::fem::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count);
  fix.bulk.modification_end();

  PairIterRelation rel2 = face2.relations(NODE_RANK);

  STKUNIT_ASSERT_TRUE( true );
}

STKUNIT_UNIT_TEST( testTopologyHelpers, element_side_polarity_valid )
{
  // Coverage of element_side_polarity in TopologyHelpers.cpp 168-181 and 200-215
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  fix.bulk.modification_begin();
  Entity & element = stk::mesh::fem::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  const EntityId zero_side_count = 0;
  Entity& face2 = stk::mesh::fem::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count);
  fix.bulk.modification_end();

  const int local_side_id = 0;
  STKUNIT_ASSERT_TRUE( stk::mesh::fem::element_side_polarity( element, face2, local_side_id) );

}

STKUNIT_UNIT_TEST( testTopologyHelpers, element_side_polarity_invalid_1 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  // Coverage of element_side_polarity in TopologyHelpers.cpp
  {
    fix.bulk.modification_begin();
    Entity & element = stk::mesh::fem::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
    const EntityId zero_side_count = 0;
    Entity& face = stk::mesh::fem::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count);
    fix.bulk.modification_end();

    const int invalid_local_side_id = -1;
    // Hits "Unsuported local_side_id" error condition:
    STKUNIT_ASSERT_THROW(
        stk::mesh::fem::element_side_polarity( element, face, invalid_local_side_id),
        std::runtime_error
        );
  }
}

STKUNIT_UNIT_TEST( testTopologyHelpers, element_side_polarity_invalid_2 )
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  // Coverage of element_side_polarity in TopologyHelpers.cpp - NULL = elem_top
  fix.bulk.modification_begin();

  PartVector part_intersection;
  part_intersection.push_back ( &fix.generic_element_part);
  Entity & element = fix.bulk.declare_entity(fix.element_rank, fix.nextEntityId(), part_intersection);
  STKUNIT_ASSERT_TRUE( stk::mesh::fem::get_cell_topology( element ).getCellTopologyData() == NULL );

  Entity & element_with_top = stk::mesh::fem::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  STKUNIT_ASSERT_TRUE( stk::mesh::fem::get_cell_topology( element_with_top ).getCellTopologyData() != NULL );

  const EntityId zero_side_count = 0;
  Entity& face_with_top = stk::mesh::fem::declare_element_side( fix.bulk, fix.nextEntityId(), element_with_top, zero_side_count);

  fix.bulk.modification_end();

  const int valid_local_side_id = 0;
  // Hits "Element has no defined topology" error condition:
  STKUNIT_ASSERT_TRUE( stk::mesh::fem::get_cell_topology( element ).getCellTopologyData() == NULL );
  STKUNIT_ASSERT_THROW(
      stk::mesh::fem::element_side_polarity( element, face_with_top, valid_local_side_id),
      std::runtime_error
      );

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}
