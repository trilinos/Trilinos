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
#include "mpi.h"                        // for MPI_COMM_WORLD, MPI_Barrier, etc
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

void testNodesAreSelected(stk::mesh::BulkData &stkMeshBulkData,
                     const stk::mesh::EntityVector &nodes,
                     stk::mesh::Selector part1Selector)
{
    stk::mesh::EntityVector selectedNodes;
    stk::mesh::get_selected_entities(part1Selector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), selectedNodes);
    EXPECT_EQ(selectedNodes.size(), nodes.size());
    for(size_t i = 0; i < nodes.size(); i++)
    {
        EXPECT_EQ(nodes[i], selectedNodes[i]);
    }
}

void testAddingNodesToPart(stk::mesh::BulkData &stkMeshBulkData,
                           const stk::mesh::EntityVector &nodes,
                           stk::mesh::Part &nodePart,
                           stk::mesh::Selector part1Selector,
                           stk::mesh::Field<double> &nodeField1)
{
    stkMeshBulkData.modification_begin();
    stk::mesh::PartVector addParts(1, &nodePart);
    for(size_t i = 0; i < nodes.size(); ++i)
    {
        if(stkMeshBulkData.parallel_owner_rank(nodes[i]) == stkMeshBulkData.parallel_rank())
        {
            stkMeshBulkData.change_entity_parts(nodes[i], addParts);
        }
    }
    stkMeshBulkData.modification_end();

    for(size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_TRUE(stkMeshBulkData.bucket(nodes[i]).member(nodePart));
        double* fieldPtr = stk::mesh::field_data(nodeField1, nodes[i]);
        EXPECT_EQ(*static_cast<double*>(nodeField1.get_initial_value()), *fieldPtr);
    }

    testNodesAreSelected(stkMeshBulkData, nodes, part1Selector);
}

TEST(UnitTestParts, CreateAfterCommit)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  const std::string generatedMeshSpecification = "generated:1x1x4";
  stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();

  stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
  stk::mesh::Part& nodePart1 = stkMeshMetaData.declare_part("nodePart1");
  stk::mesh::Field<double>& nodeField1 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodeField1");
  const double initialValue = 3.14;
  stk::mesh::put_field(nodeField1, nodePart1, &initialValue);
  stkMeshIoBroker.populate_bulk_data();

  stk::mesh::BulkData& stkMeshBulkData = stkMeshIoBroker.bulk_data();

  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);

  stk::mesh::Selector part1Selector = nodePart1;
  testAddingNodesToPart(stkMeshBulkData, nodes, nodePart1, part1Selector, nodeField1);

  EXPECT_TRUE(stkMeshMetaData.is_commit());

  stk::mesh::Part& new_part = stkMeshMetaData.declare_part("new_part");

  testAddingNodesToPart(stkMeshBulkData, nodes, new_part, part1Selector, nodeField1);

  stk::mesh::Selector newPartSelector = new_part;
  testNodesAreSelected(stkMeshBulkData, nodes, newPartSelector);
}

