/*

 */


#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FEMHelpers.hpp"  // for EntityKey
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>



TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown1ChownNoNodes )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/0---4/0---5/0
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
  else {
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown1ChownUncommonNodesFlip)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //                                           then flip
  //   1/0---4/0---5/0      1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |        |     |     |          |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |        |     |     |          |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/1        2/1---3/0---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part, element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part, element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(elem1) && bulk.bucket_ptr(elem1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  } else {
    EXPECT_TRUE(bulk.bucket_ptr(elem2) && bulk.bucket_ptr(elem2)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }

  //okay now flip
  //    1/0---4/0---5/1          1/1---4/0---5/0
  //     |     |     |            |     |     |
  //     | 1/0 | 2/1 |       =>   | 1/1 | 2/0 |
  //     |     |     |            |     |     |
  //    2/0---3/0---6/1          2/1---3/0---6/0
  stk::mesh::EntityProcVec entity_procs_flip;
  if (p_rank == 0) {
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 2), 1));
  } else {
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::ELEM_RANK, 2), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 0));
  }
  bulk.change_entity_owner(entity_procs_flip);

  //now do checks for flip
  elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(elem1));
  EXPECT_TRUE(bulk.is_valid(elem2));
  EXPECT_TRUE(bulk.is_valid(node1));
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));
  if (p_rank == 0) {
      EXPECT_TRUE(bulk.bucket_ptr(elem2) && bulk.bucket_ptr(elem2)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(elem1) && !bulk.bucket_ptr(elem1)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
      EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
      EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->shared());
    } else {
      EXPECT_TRUE(bulk.bucket_ptr(elem2) && !bulk.bucket_ptr(elem2)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(elem1) && bulk.bucket_ptr(elem1)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
      EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
      EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
      EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->shared());
    }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown1ChownItsNodes)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/1---6/1

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 3), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
  else {
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown2ChownNoNodes )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/0---4/0---5/0
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/1 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/0---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[0], 1));
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node1));
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
  else {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown2ChownCommonNodes)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/0---4/1---5/0
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/1 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/1---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[0], 1));
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 3), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node1));
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
  else {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->shared());
  }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown2ChownUncommonNodes)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/1---4/0---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/1 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/1---3/0---6/1

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[0], 1));
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 2), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node1));
  EXPECT_TRUE(bulk.is_valid(node3));
  EXPECT_TRUE(bulk.is_valid(node5));

  if (p_rank == 0) {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
  else {
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
}



TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown2ChownAllNodesOf1)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/0---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/1 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/0---3/1---6/1

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[0], 1));
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 3), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.is_valid(node1));
  EXPECT_TRUE(bulk.is_valid(node3));

  if (p_rank == 0) {
    EXPECT_FALSE(bulk.is_valid(node5));
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && ! bulk.bucket_ptr(node3)->shared());
  }
  else {
    EXPECT_TRUE(bulk.is_valid(node5));
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
}


TEST( UnitTestOfBulkData, testChangeEntityOwnerShared_2EltsChown2ChownAllNodes)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::mesh::BulkData bulk( meta, pm);

  //   id/owner_proc 
  //
  //   1/0---4/0---5/0      1/1---4/1---5/1
  //    |     |     |        |     |     |
  //    | 1/0 | 2/0 |   =>   | 1/1 | 2/1 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/0      2/1---3/1---6/1

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  // stk::mesh::Part &node_part = meta.declare_part_with_topology("node_part",stk::topology::NODE);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[0], 1));
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 2), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 3), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
  stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);

  if (p_rank == 0) {
    EXPECT_FALSE(bulk.is_valid(node1));
    EXPECT_FALSE(bulk.is_valid(node3));
    EXPECT_FALSE(bulk.is_valid(node5));
  }
  else {
    EXPECT_TRUE(bulk.is_valid(node1));
    EXPECT_TRUE(bulk.is_valid(node3));
    EXPECT_TRUE(bulk.is_valid(node5));
    EXPECT_TRUE(bulk.bucket_ptr(node1) && bulk.bucket_ptr(node1)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node1) && !bulk.bucket_ptr(node1)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && bulk.bucket_ptr(node3)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node3) && !bulk.bucket_ptr(node3)->shared());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && bulk.bucket_ptr(node5)->owned());
    EXPECT_TRUE(bulk.bucket_ptr(node5) && !bulk.bucket_ptr(node5)->shared());
  }
}
