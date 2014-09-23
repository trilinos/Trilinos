/*------------------------------------------------------------------------*/
/*         _        Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for NULL
#include <Shards_BasicTopologies.hpp>   // for getCellTopologyData, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element, etc
#include <stk_mesh/base/MetaData.hpp>   // for get_cell_topology, MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include "Shards_CellTopologyData.h"    // for CellTopologyData
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct EntitySideComponent; } }

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

class TopologyHelpersTestingFixture
{
 public:
  TopologyHelpersTestingFixture(ParallelMachine pm);
  ~TopologyHelpersTestingFixture() {}

  const int spatial_dimension;
  MetaData meta;
  BulkData bulk;
  const EntityRank element_rank;
  const EntityRank side_rank;
  Part & generic_element_part;
  Part & element_tet_part;
  Part & element_wedge_part;
  Part & generic_face_part;
  Part & another_generic_face_part;
  Part & face_tri_part;
  Part & face_quad_part;
  Part & another_generic_element_part;

  EntityId nextEntityId()
  { return psize*(++entity_id)+prank; }

  Entity create_entity( EntityRank rank, Part& part_membership)
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
  , bulk( meta, pm, 100 )
  , element_rank( stk::topology::ELEMENT_RANK )
  , side_rank( meta.side_rank())
  , generic_element_part( meta.declare_part("another part", element_rank ) )
  , element_tet_part( meta.declare_part_with_topology( "block_left_1", stk::topology::TET_4 ) )
  , element_wedge_part( meta.declare_part_with_topology( "block_left_2", stk::topology::WEDGE_15 ) )
  , generic_face_part( meta.declare_part_with_topology( "A_1", stk::topology::QUAD_4 ) )
  , another_generic_face_part( meta.declare_part("A_2", side_rank ) )
  , face_tri_part( meta.declare_part_with_topology("A_3_0", stk::topology::TRI_3))
  , face_quad_part( meta.declare_part("A_3", side_rank ) )
  , another_generic_element_part( meta.declare_part("B_3", element_rank ) )
  , entity_id(0u)
  , psize(bulk.parallel_size())
  , prank(bulk.parallel_rank())
{
  meta.commit();
}

namespace {

TEST( testTopologyHelpers, get_cell_topology_based_on_part)
{
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();
  Entity elem1  = fix.create_entity( fix.side_rank, fix.generic_face_part );

  std::vector<Entity> elem_node(4);
  for (int i = 0; i < 4; ++i) {
    elem_node[i] = fix.bulk.declare_entity(stk::topology::NODE_RANK, 100 + i);
    fix.bulk.declare_relation(elem1, elem_node[i], i);
  }

  PartVector tmp(1);
  tmp[0] = & fix.face_quad_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( stk::mesh::get_cell_topology(fix.bulk.bucket(elem1)).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( stk::mesh::get_cell_topology(fix.bulk.bucket(elem1)).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  tmp[0] = & fix.another_generic_face_part;
  fix.bulk.change_entity_parts ( elem1 , tmp );
  ASSERT_EQ( stk::mesh::get_cell_topology(fix.bulk.bucket(elem1)).getCellTopologyData(), shards::getCellTopologyData< shards::Quadrilateral<4> >() );
  ASSERT_NE( stk::mesh::get_cell_topology(fix.bulk.bucket(elem1)).getCellTopologyData() , shards::getCellTopologyData< shards::Wedge<15> >() );

  fix.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_no_topology )
{
  // Coverage for declare_element_side - TopologyHelpers.cpp - "Cannot discern element topology"
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();
  Entity elem4  = fix.create_entity( fix.element_rank , fix.generic_element_part );
  ASSERT_THROW(
    stk::mesh::declare_element_side( fix.bulk, fix.element_rank, elem4, fix.nextEntityId(), &fix.element_wedge_part ),
    std::runtime_error
      );
  //fix.bulk.modification_end();


  {
    EntityId elem_node[4];
    elem_node[0] = 1;
    elem_node[1] = 2;
    elem_node[2] = 3;
    elem_node[3] = 4;
    //fix.bulk.modification_begin();
    // Cannot declare an element without a topology defined
    ASSERT_THROW(
      stk::mesh::declare_element(fix.bulk, fix.generic_element_part, fix.nextEntityId(), elem_node),
        std::runtime_error
        );
  }
}

TEST( testTopologyHelpers, declare_element_side_wrong_bulk_data)
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "BulkData for 'elem' and 'side' are different"
  TopologyHelpersTestingFixture fix1(MPI_COMM_WORLD);

  fix1.bulk.modification_begin();

  TopologyHelpersTestingFixture fix2(MPI_COMM_WORLD);
  fix2.bulk.modification_begin();
  fix2.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_no_topology_2 )
{
  // Coverage for verify_declare_element_side - in TopologyHelpers.cpp - "No element topology found and cell side id exceeds..."
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;
  Entity element  = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node);
  const CellTopologyData * const elem_top = stk::mesh::get_cell_topology( fix.bulk.bucket(element) ).getCellTopologyData();
  const EntityId nSideCount = elem_top->side_count + 10 ;
  ASSERT_THROW(
    stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, nSideCount, &fix.element_tet_part ),
    std::runtime_error
      );
  fix.bulk.modification_end();
}

TEST( testTopologyHelpers, declare_element_side_full )
{
  // Go all way the through declare_element_side - use new element
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);

  fix.bulk.modification_begin();

  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );

  const EntityId zero_side_count = 0;
  Entity face2 = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                  &fix.face_tri_part);
  fix.bulk.modification_end();

  stk::mesh::Entity const *rel2_nodes = fix.bulk.begin_nodes(face2);
  ASSERT_TRUE(rel2_nodes != 0);

  ASSERT_TRUE( true );  // This test is to check compilation.
}

TEST( testTopologyHelpers, element_side_polarity_valid )
{
  // Coverage of element_side_polarity in TopologyHelpers.cpp 168-181 and 200-215
  TopologyHelpersTestingFixture fix(MPI_COMM_WORLD);
  EntityId elem_node[4];
  elem_node[0] = 1;
  elem_node[1] = 2;
  elem_node[2] = 3;
  elem_node[3] = 4;

  fix.bulk.modification_begin();
  Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  const EntityId zero_side_count = 0;
  Entity face2 = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                  &fix.face_tri_part);
  fix.bulk.modification_end();

  const int local_side_id = 0;
  ASSERT_TRUE( fix.bulk.element_side_polarity( element, face2, local_side_id) );

}

TEST( testTopologyHelpers, element_side_polarity_invalid_1 )
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
    Entity element = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
    const EntityId zero_side_count = 0;
    Entity face = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element, zero_side_count,
                                                   &fix.face_tri_part);
    fix.bulk.modification_end();

    const unsigned invalid_local_side_id = static_cast<unsigned>(-1);
    // Hits "Unsuported local_side_id" error condition:
    ASSERT_THROW(
        fix.bulk.element_side_polarity( element, face, invalid_local_side_id),
        std::runtime_error
        );
  }
}

TEST( testTopologyHelpers, element_side_polarity_invalid_2 )
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
  Entity element = fix.bulk.declare_entity(fix.element_rank, fix.nextEntityId(), part_intersection);
  ASSERT_TRUE( stk::mesh::get_cell_topology(fix.bulk.bucket(element) ).getCellTopologyData() == NULL );

  Entity element_with_top = stk::mesh::declare_element(fix.bulk, fix.element_tet_part, fix.nextEntityId(), elem_node );
  ASSERT_TRUE( stk::mesh::get_cell_topology( fix.bulk.bucket(element_with_top) ).getCellTopologyData() != NULL );

  const EntityId zero_side_count = 0;
  Entity face_with_top = stk::mesh::declare_element_side( fix.bulk, fix.nextEntityId(), element_with_top, zero_side_count,
                                                          &fix.face_tri_part);
  const int valid_local_side_id = 0;
  ASSERT_THROW(
      fix.bulk.element_side_polarity( element, face_with_top, valid_local_side_id),
      std::runtime_error
      );

  //modification_end is not called due to difference in expected behavior for release and debug builds - debug should throw, release should not
  //difference occurs within check_for_connected_nodes method
  //ASSERT_THROW(fix.bulk.modification_end(), std::logic_error);

  // Hits "Element has no defined topology" error condition:
  //ASSERT_TRUE( stk::mesh::get_cell_topology( fix.bulk.bucket(element) ).getCellTopologyData() == NULL );



}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}
