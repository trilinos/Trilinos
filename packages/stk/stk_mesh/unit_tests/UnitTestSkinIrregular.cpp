/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <iostream>

using stk::mesh::EntityId;
using stk::mesh::EntityRank;

//---------------------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestSkin, SkinPocket)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.declare_part( "hex_part", hex_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  //create and skin a 2 hex-element mesh with a pocket
  //in a normal mesh 6 and 13 would be the same node
  //
  //    8-------7-------12
  //   /|      /|\     /|
  //  / |     / | \   / |
  // 5-------6  | 13-11 |
  // |  4----|--3/---|--10
  // | /     | //    | /
  // |/      |//     |/
  // 1-------2-------9
  //

  fem_meta.commit();

  bulk_data.modification_begin();

  // declare left element on first process
  if (p_rank == 0)
  {
    EntityId element_id = 1;
    EntityId node_ids[8] = { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::fem::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  // declare right element on last process
  if (p_rank == p_size -1)
  {
    EntityId element_id = 2;
    EntityId node_ids[8] = { 2, 9, 10, 3, 13, 11, 12, 7};

    stk::mesh::fem::declare_element( bulk_data, hex_part, element_id, node_ids);
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //each element should have 6 faces attached to it
  for (EntityId element_id = 1; element_id < 3; ++element_id) {
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, element_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 6);
    }
  }
}

STKUNIT_UNIT_TEST( UnitTestSkin, SkinTwoStackedShells)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part( "shell_part", shell_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  fem_meta.commit();

  //create and skin stacked shells
  // 4-------3
  // |       |
  // |       |
  // |       |
  // 1-------2
  //
  // create following 8 shells
  // shells 1 defined with right hand rule
  // shells 2 defined with left hand rule
  //
  // shell_id:  node_list
  // 1: (1,2,3,4)
  // 2: (1,2,3,4)

  EntityId node_ids[4] = { 1, 2, 3, 4};

  bulk_data.modification_begin();

  {
    bool create_shell_on_proc = static_cast<int>(p_rank) == 0;
    if (create_shell_on_proc) {
      EntityId element_id = static_cast<EntityId>(1);
      stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
    }
  }

  {
    bool create_shell_on_proc = static_cast<int>(p_rank) == (std::max(0,static_cast<int>(p_size)-1));
    if (create_shell_on_proc) {
      EntityId element_id = static_cast<EntityId>(2);
      stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
    }
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //count number of sides in mesh
  {
    stk::mesh::Selector select_sides = meta_data.locally_owned_part()  ;
    const std::vector<stk::mesh::Bucket*>& side_buckets = bulk_data.buckets(side_rank);
    int num_sides = stk::mesh::count_selected_entities( select_sides, side_buckets);


    stk::all_reduce(MPI_COMM_WORLD, stk::ReduceSum<1>(&num_sides));

    // Verify that the correct 2 sides are present.

    STKUNIT_ASSERT_EQUAL( num_sides, 2 );
  }
}

//---------------------------------------------------------------------------------------
STKUNIT_UNIT_TEST( UnitTestSkin, SkinStackedShells)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part( "shell_part", shell_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  fem_meta.commit();

  //create and skin stacked shells
  // 4-------3
  // |       |
  // |       |
  // |       |
  // 1-------2
  //
  // create following 8 shells
  // shells 1-4 defined with right hand rule
  // shells 5-6 defined with left hand rule
  //
  // shell_id:  node_list
  // 1: (1,2,3,4)
  // 2: (2,3,4,1)
  // 3: (3,4,1,2)
  // 4: (4,1,2,3)
  // 5: (4,3,2,1)
  // 6: (3,2,1,4)
  // 7: (2,1,4,3)
  // 8: (1,4,3,2)

  EntityId node_ids[8] = { 1, 2, 3, 4, 1, 2, 3, 4};
  EntityId reverse_node_ids[8] = { 4, 3, 2, 1, 4, 3, 2, 1};

  bulk_data.modification_begin();

  //create shells
  for ( int i = 0; i<4; i++ ) {

    bool create_shell_on_proc = static_cast<int>(p_rank) == (std::max(0,static_cast<int>(p_size)-1-i));
    if (create_shell_on_proc) {
      EntityId element_id = static_cast<EntityId>(i+1);
      stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids+i);
    }

    bool create_reverse_shell_on_proc = static_cast<int>(p_rank) == (std::max(0,static_cast<int>(p_size)-1-4-i));
    if (create_reverse_shell_on_proc) {
      EntityId reverse_element_id = static_cast<EntityId>(i+1+4);
      stk::mesh::fem::declare_element( bulk_data, shell_part, reverse_element_id, reverse_node_ids+i);
    }
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //count number of sides in mesh
  {
    stk::mesh::Selector select_sides = meta_data.locally_owned_part()  ;
    const std::vector<stk::mesh::Bucket*>& side_buckets = bulk_data.buckets( side_rank);
    int num_sides = stk::mesh::count_selected_entities( select_sides, side_buckets);


    stk::all_reduce(MPI_COMM_WORLD, stk::ReduceSum<1>(&num_sides));

    // Verify that the correct 2 sides are present.

    STKUNIT_ASSERT_EQUAL( num_sides, 2 );
  }

  //check that faces are attached to correct sides
  {
    EntityId face_1_id = 0; //invalid face id
    EntityId face_2_id = 0; //invalid face id
    for (EntityId shell_id = 1; shell_id < 5; ++shell_id) {
      stk::mesh::Entity * shell = bulk_data.get_entity( element_rank, shell_id);
      if ( shell != NULL) {
        stk::mesh::PairIterRelation shell_side_relations = shell->relations(side_rank);

        STKUNIT_ASSERT_TRUE( shell_side_relations.size() == 2);

        // verify that only one side has been created
        // and that all stacked shells reference this side
        if (face_1_id == 0) {
          face_1_id = shell_side_relations->entity()->identifier();
        }
        else {
          STKUNIT_EXPECT_TRUE( face_1_id == shell_side_relations->entity()->identifier());
        }

        //check that the side is one the correct local side of the shell
        STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 0);

        ++shell_side_relations;

        if (face_2_id == 0) {
          face_2_id = shell_side_relations->entity()->identifier();
        }
        else {
          STKUNIT_EXPECT_TRUE( face_2_id == shell_side_relations->entity()->identifier());
        }

        //check that the side is one the correct local side of the shell
        STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 1);

      }
    }

    for (EntityId shell_id = 5; shell_id < 9; ++shell_id) {
      stk::mesh::Entity * shell = bulk_data.get_entity( element_rank, shell_id);
      if ( shell != NULL) {
        stk::mesh::PairIterRelation shell_side_relations = shell->relations(side_rank);

        STKUNIT_ASSERT_TRUE( shell_side_relations.size() == 2);

        // verify that only one side has been created
        // and that all stacked shells reference this side
        if (face_2_id == 0) {
          face_2_id = shell_side_relations->entity()->identifier();
        }
        else {
          STKUNIT_EXPECT_TRUE( face_2_id == shell_side_relations->entity()->identifier());
        }

        //check that the side is one the correct local side of the shell
        STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 0);

        ++shell_side_relations;

        if (face_1_id == 0) {
          face_1_id = shell_side_relations->entity()->identifier();
        }
        else {
          STKUNIT_EXPECT_TRUE( face_1_id == shell_side_relations->entity()->identifier());
        }

        //check that the side is one the correct local side of the shell
        STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 1);

      }
    }
  }
}

//---------------------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestSkin, SkinShellOnHex)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.declare_part( "hex_part", hex_top );
  stk::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part( "shell_part", shell_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  //create and skin a hex element mesh with a shell on the first side of the hex
  // Using a shell defined by the nodes (1, 2, 6, 5) produces an orientated shell
  //
  // Therefore the shell will be skinned with 1 side and the hex will have 5.
  //
  //    8-------7
  //   /|      /|
  //  / |     / |
  // 5=======6  |
  // || 4----||-3
  // ||/     ||/
  // |/      |/
  // 1=======2
  //

  fem_meta.commit();

  bulk_data.modification_begin();

  // declare hex element on first process
  if (p_rank == 0)
  {
    EntityId element_id = 1;
    EntityId node_ids[8] = { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::fem::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  // declare shell element on last process
  if (p_rank == p_size -1)
  {
    EntityId element_id = 2;
    EntityId node_ids[8] = { 1, 2, 6, 5};

    stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //check hex
  {
    EntityId hex_id = 1;
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, hex_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 5);
      for (; !element_side_relations.empty(); ++element_side_relations) {
        unsigned local_side_id = element_side_relations->identifier();
        bool correct_side_skinned = local_side_id > 0 && local_side_id < 6;
        STKUNIT_EXPECT_TRUE (correct_side_skinned);
        std::cout << "Hex local side id: " << local_side_id << std::endl;
      }
    }
  }

  //check shell
  {
    EntityId shell_id = 2;
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, shell_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 1);
      for (; !element_side_relations.empty(); ++element_side_relations) {
        unsigned local_side_id = element_side_relations->identifier();
        bool correct_side_skinned = local_side_id == 0;
        STKUNIT_EXPECT_TRUE (correct_side_skinned);
        std::cout << "Shell local side id: " << local_side_id << std::endl;
      }
    }
  }
}

//---------------------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestSkin, SkinInvertedShellOnHex)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.declare_part( "hex_part", hex_top );
  stk::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part( "shell_part", shell_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  //create and skin a hex element mesh with an inverted shell
  // Using a shell defined by the nodes (1, 2, 5, 6) produces an inverted shell
  // with no valid orientation
  //
  // Therefore the shell will be skinned with 2 sides and the hex will have 6.
  //
  //    8-------7
  //   /|      /|
  //  / |     / |
  // 5=======6  |
  // || 4----||-3
  // ||/     ||/
  // |/      |/
  // 1=======2
  //

  fem_meta.commit();

  bulk_data.modification_begin();

  // declare hex element on first process
  if (p_rank == 0)
  {
    EntityId element_id = 1;
    EntityId node_ids[8] = { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::fem::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  // declare shell element on last process
  if (p_rank == p_size -1)
  {
    EntityId element_id = 2;
    EntityId node_ids[8] = { 1, 2, 5, 6};

    stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //check hex
  {
    EntityId hex_id = 1;
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, hex_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 6);
      for (; !element_side_relations.empty(); ++element_side_relations) {
        unsigned local_side_id = element_side_relations->identifier();
        bool correct_side_skinned = local_side_id < 6;
        STKUNIT_EXPECT_TRUE (correct_side_skinned);
        std::cout << "Hex local side id: " << local_side_id << std::endl;
      }
    }
  }

  //check shell
  {
    EntityId shell_id = 2;
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, shell_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 2);
      for (; !element_side_relations.empty(); ++element_side_relations) {
        unsigned local_side_id = element_side_relations->identifier();
        bool correct_side_skinned = local_side_id < 2;
        STKUNIT_EXPECT_TRUE (correct_side_skinned);
        std::cout << "Shell local side id: " << element_side_relations->identifier() << std::endl;
      }
    }
  }
}

//---------------------------------------------------------------------------------------

STKUNIT_UNIT_TEST( UnitTestSkin, SkinStackedShellOnHex)
{
  enum { SpatialDim = 3 };

  stk::ParallelMachine pm = MPI_COMM_WORLD ;
  const unsigned p_rank = stk::parallel_machine_rank( pm );
  const unsigned p_size = stk::parallel_machine_size( pm );

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim));
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data , pm );
  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::Part & hex_part = fem_meta.declare_part( "hex_part", hex_top );
  stk::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk::mesh::Part & shell_part = fem_meta.declare_part( "shell_part", shell_top );
  const EntityRank element_rank = fem_meta.element_rank();
  const EntityRank side_rank    = fem_meta.side_rank();

  //create and skin a hex element mesh with 3 shells on the first side of the hex
  // Using shells defined by the nodes (1, 2, 6, 5), (6, 5, 1, 2), and (1, 5, 6, 2)
  // produces an orientated shells.
  //
  // Therefore the shells will all have a relation to the same side.
  //
  //    8-------7
  //   /|      /|
  //  / |     / |
  // 5=======6  |
  // || 4----||-3
  // ||/     ||/
  // |/      |/
  // 1=======2
  //

  fem_meta.commit();

  bulk_data.modification_begin();

  bool create_hex_this_proc = (p_rank == 0);
  bool create_shell_1_this_proc = static_cast<int>(p_rank) == (std::max(0,static_cast<int>(p_size)-3));
  bool create_shell_2_this_proc = static_cast<int>(p_rank) == (std::max(0,static_cast<int>(p_size)-2));
  bool create_shell_3_this_proc = (p_rank == p_size -1);

  if (create_hex_this_proc)
  {
    EntityId element_id = 1;
    EntityId node_ids[8] = { 1, 2, 3, 4, 5, 6, 7, 8};

    stk::mesh::fem::declare_element( bulk_data, hex_part, element_id, node_ids);

  }

  if (create_shell_1_this_proc)
  {
    EntityId element_id = 2;
    EntityId node_ids[8] = { 1, 2, 6, 5};

    stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
  }

  if (create_shell_2_this_proc)
  {
    EntityId element_id = 3;
    EntityId node_ids[8] = { 6, 5, 1, 2};

    stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
  }

  if (create_shell_3_this_proc)
  {
    EntityId element_id = 4;
    EntityId node_ids[8] = { 1, 5, 6, 2};

    stk::mesh::fem::declare_element( bulk_data, shell_part, element_id, node_ids);
  }

  bulk_data.modification_end();

  //skin the mesh
  stk::mesh::skin_mesh(bulk_data, element_rank);

  //check hex
  {
    EntityId hex_id = 1;
    stk::mesh::Entity * element = bulk_data.get_entity( element_rank, hex_id);
    if ( element != NULL) {
      stk::mesh::PairIterRelation element_side_relations = element->relations(side_rank);
      STKUNIT_EXPECT_TRUE( element_side_relations.size() == 5);
      for (; !element_side_relations.empty(); ++element_side_relations) {
        unsigned local_side_id = element_side_relations->identifier();
        bool correct_side_skinned = local_side_id > 0 && local_side_id < 6;
        STKUNIT_EXPECT_TRUE (correct_side_skinned);
        std::cout << "Hex local side id: " << local_side_id << std::endl;
      }
    }
  }

  //check shells
  {
    EntityId face_id = 0; //invalid face id
    for (EntityId shell_id = 2; shell_id < 5; ++shell_id) {
      stk::mesh::Entity * shell = bulk_data.get_entity( element_rank, shell_id);
      if ( shell != NULL) {
        stk::mesh::PairIterRelation shell_side_relations = shell->relations(side_rank);

        STKUNIT_EXPECT_TRUE( shell_side_relations.size() == 1);

        // verify that only one side has been created
        // and that all stacked shells reference this side
        if (face_id == 0) {
          face_id = shell_side_relations->entity()->identifier();
        }
        else {
          STKUNIT_EXPECT_TRUE( face_id == shell_side_relations->entity()->identifier());
          std::cout << "Shell: " << shell_id
                    << "\tFace_id: " << face_id
                    << "\tFace_id: " << shell_side_relations->entity()->identifier()
                    << std::endl;
        }

        //check that the side is one the correct local side of the shell
        //shells 1 and 2 follow the right hand rule so the side should be on
        //local_side_id 0.
        //shell 3 should have it's side on local_side_1d 1
        if (shell_id != 4) {
          STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 0);
        }
        else {
          STKUNIT_EXPECT_TRUE( shell_side_relations->identifier() == 1);
        }
      }
    }
  }
}
