/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


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

void setupKeyholeMesh2D_case1(stk::mesh::BulkData& bulk)
{
//
//   proc 0      proc 1
//            |
//            | block_2
//            | 10---9
// block_1    | | 3  |
//    4----3  | 3----8 
//    | 1  |  |
//    1----2  | 2----7
//            | | 2  |
//            | 5----6
//            |
//
//shared nodes 2 and 3 should be members of block_1 and block_2 on both procs
//

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  meta.commit();

  bulk.modification_begin();

  const unsigned num_nodes = 10;
  std::vector<stk::mesh::Entity> nodes(num_nodes+1, stk::mesh::Entity());

  for(unsigned i=1; i<=num_nodes; ++i) {
    if (bulk.parallel_rank() == 0) {
      if (i > 4) break;
    }
    if (bulk.parallel_rank() == 1) {
      if (i==1 || i==4) continue;
    }
    nodes[i] = bulk.declare_entity(stk::topology::NODE_RANK, static_cast<stk::mesh::EntityId>(i));
  }

  stk::mesh::EntityId id = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity elem1 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_1);
    bulk.declare_relation(elem1, nodes[1], 0);
    bulk.declare_relation(elem1, nodes[2], 1);
    bulk.declare_relation(elem1, nodes[3], 2);
    bulk.declare_relation(elem1, nodes[4], 3);
  }

  if (bulk.parallel_rank() == 1) {
    id = 2;
    stk::mesh::Entity elem2 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_2);
    bulk.declare_relation(elem2, nodes[5], 0);
    bulk.declare_relation(elem2, nodes[6], 1);
    bulk.declare_relation(elem2, nodes[7], 2);
    bulk.declare_relation(elem2, nodes[2], 3);

    id = 3;
    stk::mesh::Entity elem3 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_2);
    bulk.declare_relation(elem3, nodes[3], 0);
    bulk.declare_relation(elem3, nodes[8], 1);
    bulk.declare_relation(elem3, nodes[9], 2);
    bulk.declare_relation(elem3, nodes[10], 3);
  }

  bulk.modification_end();
}

void setupKeyholeMesh2D_case2(stk::mesh::BulkData& bulk)
{
//
//   proc 0      proc 1
//            |
//            | block_2 block_3
//            |
//            |         12---11
// block_1    |         | 4  |
//    4----3  | 3----6  6----10
//    | 1  |  | |  2 |
//    1----2  | 2----5  5----9
//            |         | 3  |
//            |         7----8
//            |
//
//nodes 5 and 6 are ghosts (aura) on proc 0,
//and should be members of block_2 and block_3 on proc 0
//

  stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::QUAD_4_2D);
  meta.commit();

  bulk.modification_begin();

  const unsigned num_nodes = 12;
  std::vector<stk::mesh::Entity> nodes(num_nodes+1, stk::mesh::Entity());

  for(unsigned i=1; i<=num_nodes; ++i) {
    if (bulk.parallel_rank() == 0) {
      if (i > 4) break;
    }
    if (bulk.parallel_rank() == 1) {
      if (i==1 || i==4) continue;
    }
    nodes[i] = bulk.declare_entity(stk::topology::NODE_RANK, static_cast<stk::mesh::EntityId>(i));
  }

  stk::mesh::EntityId id = 1;
  if (bulk.parallel_rank() == 0) {
    stk::mesh::Entity elem1 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_1);
    bulk.declare_relation(elem1, nodes[1], 0);
    bulk.declare_relation(elem1, nodes[2], 1);
    bulk.declare_relation(elem1, nodes[3], 2);
    bulk.declare_relation(elem1, nodes[4], 3);
  }

  if (bulk.parallel_rank() == 1) {
    id = 2;
    stk::mesh::Entity elem2 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_2);
    bulk.declare_relation(elem2, nodes[2], 0);
    bulk.declare_relation(elem2, nodes[5], 1);
    bulk.declare_relation(elem2, nodes[6], 2);
    bulk.declare_relation(elem2, nodes[3], 3);

    id = 3;
    stk::mesh::Entity elem3 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_3);
    bulk.declare_relation(elem3, nodes[7], 0);
    bulk.declare_relation(elem3, nodes[8], 1);
    bulk.declare_relation(elem3, nodes[9], 2);
    bulk.declare_relation(elem3, nodes[5], 3);

    id = 4;
    stk::mesh::Entity elem4 = bulk.declare_entity(stk::topology::ELEM_RANK, id, block_3);
    bulk.declare_relation(elem4, nodes[6], 0);
    bulk.declare_relation(elem4, nodes[10], 1);
    bulk.declare_relation(elem4, nodes[11], 2);
    bulk.declare_relation(elem4, nodes[12], 3);
  }

  bulk.modification_end();
}

TEST(UnitTestKeyhole, NodeParts_case1)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case1(bulk);

  stk::mesh::Part& shared = meta.globally_shared_part();
  const stk::mesh::BucketVector& shared_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, shared);
  stk::mesh::PartVector blocks(2);
  blocks[0] = meta.get_part("block_1");
  blocks[1] = meta.get_part("block_2");
  unsigned num_shared_nodes = 0;
  for(size_t i=0; i<shared_node_buckets.size(); ++i) {
    num_shared_nodes += shared_node_buckets[i]->size();
    const stk::mesh::Bucket& bucket = *shared_node_buckets[i];
    std::ostringstream oss;
    oss<<"proc "<<bulk.parallel_rank()<<", shared node ids: ";
    for(size_t j=0; j<bucket.size(); ++j) oss <<bulk.identifier(bucket[j])<<" ";
    std::cerr<<oss.str()<<std::endl;
    bool in_both_blocks = bucket.member_all(blocks);
    EXPECT_TRUE(in_both_blocks);
  }

  const unsigned expected_num_shared_nodes = 2;
  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);
}

TEST(UnitTestKeyhole, NodeParts_case2)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 2) {
    return;
  }

  const unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setupKeyholeMesh2D_case2(bulk);

  if (bulk.parallel_rank() == 0) {
    stk::mesh::Part& aura = meta.aura_part();
    const stk::mesh::BucketVector& aura_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, aura);
    stk::mesh::PartVector blocks(2);
    blocks[0] = meta.get_part("block_2");
    blocks[1] = meta.get_part("block_3");
    unsigned num_aura_nodes = 0;
    for(size_t i=0; i<aura_node_buckets.size(); ++i) {
      num_aura_nodes += aura_node_buckets[i]->size();
      const stk::mesh::Bucket& bucket = *aura_node_buckets[i];
      std::cerr<<"proc 0, aura node ids: ";
      for(size_t j=0; j<bucket.size(); ++j) std::cerr<<bulk.identifier(bucket[j])<<" ";
      std::cerr<<std::endl;
      bool in_both_blocks = bucket.member_all(blocks);
      EXPECT_TRUE(in_both_blocks);
    }
  
    const unsigned expected_num_aura_nodes = 2;
    EXPECT_EQ(expected_num_aura_nodes, num_aura_nodes);
  }
}

