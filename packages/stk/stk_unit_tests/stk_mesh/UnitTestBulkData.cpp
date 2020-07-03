// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//



#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <algorithm>                    // for sort
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

#include "Setup8Quad4ProcMesh.hpp"
#include "UnitTestCEO2Elem.hpp"
#include "UnitTestCEO3Elem.hpp"
#include "UnitTestCEO4ElemEdge.hpp"
#include "UnitTestCEO4ElemRotate.hpp"
#include "UnitTestCEO8Elem.hpp"
#include "UnitTestCEOCommonUtils.hpp"
#include "UnitTestRingFixture.hpp"  // for test_shift_ring
#include "stk_io/StkMeshIoBroker.hpp"
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
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/MeshCommImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp"  // for BoxFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/RingFixture.hpp"  // for RingFixture
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include <gtest/gtest.h>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/FieldBLAS.hpp>  // for stk::mesh::field_fill
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/FaceTestingUtils.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

//====================
extern int gl_argc;
extern char** gl_argv;

namespace
{

void donate_one_element(stk::unit_test_util::BulkDataTester & mesh)
{
    const int p_rank = mesh.parallel_rank();

    Selector select_owned(mesh.mesh_meta_data().locally_owned_part());

    std::vector<size_t> before_count;
    std::vector<size_t> after_count;

    count_entities(select_owned, mesh, before_count);

    // Change owner of an element on a process boundary
    // from P0 to P1, and then recount to confirm ownership change

    std::vector<EntityProc> change;

    // A shared node:
    EntityKey node_key;
    Entity elem = Entity();

    for(stk::mesh::EntityCommListInfoVector::const_iterator i = mesh.my_internal_comm_list().begin(); i != mesh.my_internal_comm_list().end(); ++i)
    {
        if(mesh.in_shared(i->key) && i->key.rank() == stk::topology::NODE_RANK)
        {
            node_key = i->key;
            break;
        }
    }

    ASSERT_TRUE( node_key.is_valid());

    Entity node = mesh.get_entity(node_key);
    ASSERT_TRUE( mesh.is_valid(node));

    Entity const *node_elems_i = mesh.begin_elements(node);
    Entity const *node_elems_e = mesh.end_elements(node);
    for(; (node_elems_i != node_elems_e) && !mesh.is_valid(elem); ++node_elems_i)
    {
        elem = *node_elems_i;
        if(mesh.parallel_owner_rank(elem) != p_rank)
        {
            elem = Entity();
        }
    }

    ASSERT_TRUE( mesh.is_valid(elem));

    unsigned donated_nodes = 0;

    // Only process #0 donates an element and its owned nodes:
    if(0 == p_rank)
    {
        EntityProc entry;
        entry.first = elem;
        std::vector<int> shared_procs;
        mesh.comm_shared_procs(mesh.entity_key(node),shared_procs);
        entry.second = shared_procs[0];
        change.push_back(entry);

        Entity const *elem_nodes_i = mesh.begin_nodes(elem);
        Entity const *elem_nodes_e = mesh.end_nodes(elem);
        for(; elem_nodes_i != elem_nodes_e; ++elem_nodes_i)
        {
            if(mesh.parallel_owner_rank(*elem_nodes_i) == p_rank)
            {
                entry.first = *elem_nodes_i;
                change.push_back(entry);
                ++donated_nodes;
            }
        }
    }

    mesh.change_entity_owner(change);

    count_entities(select_owned, mesh, after_count);

    if(0 == p_rank)
    {
        ASSERT_EQ( before_count[3] - 1, after_count[3]);
        ASSERT_EQ( before_count[0] - donated_nodes, after_count[0]);
    }
}

void donate_all_shared_nodes(stk::unit_test_util::BulkDataTester & mesh)
{
    const int p_rank = mesh.parallel_rank();

    const MetaData& meta = mesh.mesh_meta_data();
    const Selector select_used = meta.locally_owned_part() | meta.globally_shared_part();

    std::vector<size_t> before_count;
    std::vector<size_t> after_count;

    count_entities(select_used, mesh, before_count);

    // Donate owned shared nodes to first sharing process.

    const stk::mesh::EntityCommListInfoVector & entity_comm = mesh.my_internal_comm_list();

    ASSERT_TRUE( ! entity_comm.empty());

    std::vector<EntityProc> change;

    for(stk::mesh::EntityCommListInfoVector::const_iterator i = entity_comm.begin();
            i != entity_comm.end() && i->key.rank() == stk::topology::NODE_RANK; ++i)
            {
        Entity const node = i->entity;
        std::vector<int> shared_procs;
        mesh.comm_shared_procs(i->key,shared_procs);

        if(mesh.parallel_owner_rank(node) == p_rank && !shared_procs.empty())
        {
            change.push_back(EntityProc(node, shared_procs[0]));
        }
    }

    mesh.change_entity_owner(change);

    count_entities(select_used, mesh, after_count);

    ASSERT_TRUE( 3 <= after_count.size());
    ASSERT_EQ( before_count[0], after_count[0]);
    ASSERT_EQ( before_count[1], after_count[1]);
    ASSERT_EQ( before_count[2], after_count[2]);
    ASSERT_EQ( before_count[3], after_count[3]);
}

//----------------------------------------------------------------------
// Testing for mesh entities without relations

TEST(BulkData, testChangeOwner_nodes)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    enum
    {
        nPerProc = 10
    };
    const int p_rank = stk::parallel_machine_rank(pm);
    const int p_size = stk::parallel_machine_size(pm);
    const unsigned id_total = nPerProc * p_size;
    const unsigned id_begin = nPerProc * p_rank;
    const unsigned id_end = nPerProc * (p_rank + 1);

    const int spatial_dimension = 3;
    MetaData meta(spatial_dimension);
    BulkData bulk(meta, pm, stk::mesh::BulkData::AUTO_AURA);

    const stk::mesh::ConstPartVector no_parts;

    meta.commit();
    bulk.modification_begin();

    // Ids for all entities (all entities have type 0):

    std::vector<EntityId> ids(id_total);

    for(unsigned i = 0; i < id_total; ++i)
    {
        ids[i] = i + 1;
    }

    // Declare just those entities in my range of ids:

    for(unsigned i = id_begin; i < id_end; ++i)
    {
        bulk.declare_node(ids[i], no_parts);
    }

    ASSERT_TRUE( bulk.modification_end());

    // Verify that I only have entities in my range:

    for(unsigned i = 0; i < id_total; ++i)
    {
        Entity e = bulk.get_entity(stk::topology::NODE_RANK, ids[i]);
        if(id_begin <= i && i < id_end)
        {
            ASSERT_TRUE( bulk.is_valid(e));
        }
        else
        {
            ASSERT_TRUE( !bulk.is_valid(e));
        }
    }

    // Test change owner no-op first:

    std::vector<EntityProc> change;

    bulk.change_entity_owner(change);

    for(unsigned i = 0; i < id_total; ++i)
    {
        Entity e = bulk.get_entity(stk::topology::NODE_RANK, ids[i]);
        if(id_begin <= i && i < id_end)
        {
            ASSERT_TRUE( bulk.is_valid(e));
        }
        else
        {
            ASSERT_TRUE( !bulk.is_valid(e));
        }
    }

    // Can only test changing owner in parallel.

    if(1 < p_size)
    {
        // Give my last two ids to the next process
        // Get the previous process' last two ids

        const int p_give = (p_rank + 1) % p_size;
        const unsigned id_give = id_end - 2;
        const unsigned id_get = (id_begin + id_total - 2) % id_total;

        ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_give] )));
        ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_give+1] )));
        ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_get] )));
        ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_get+1] )));

        change.resize(2);
        change[0].first = bulk.get_entity(stk::topology::NODE_RANK, ids[id_give]);
        change[0].second = p_give;
        change[1].first = bulk.get_entity(stk::topology::NODE_RANK, ids[id_give + 1]);
        change[1].second = p_give;

        bulk.change_entity_owner(change);

        ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_get] )));
        ASSERT_TRUE( bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_get+1] )));

        // Entities given away are destroyed until the next modification cycle
        {
            Entity const e0 = bulk.get_entity(stk::topology::NODE_RANK, ids[id_give]);
            Entity const e1 = bulk.get_entity(stk::topology::NODE_RANK, ids[id_give + 1]);
            ASSERT_TRUE( !bulk.is_valid(e0));
            ASSERT_TRUE( !bulk.is_valid(e1));
        }

        ASSERT_TRUE( bulk.modification_begin());
        ASSERT_TRUE( bulk.modification_end());

        ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_give] )));
        ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK, ids[id_give+1] )));
    }
}

//----------------------------------------------------------------------
// Testing for creating existing mesh entities.

TEST(BulkData, testCreateMoreExistingOwnershipIsKept)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  enum { nPerProc = 10 };

  const int p_size = stk::parallel_machine_size( pm );
  const int p_rank = stk::parallel_machine_rank( pm );

  if (p_size != 2) {
    return;
  }


  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );

  const int spatial_dimension = 3;
  MetaData meta( spatial_dimension );

  const stk::mesh::ConstPartVector no_parts ;
  Part* edge_part = &meta.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  Part* tri_part  = &meta.declare_part_with_topology("tri_part", stk::topology::TRIANGLE_3);
  Part* shell_part  = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_TRI_3);

  meta.commit();

  BulkData bulk( meta , pm , stk::mesh::BulkData::AUTO_AURA );

  bulk.modification_begin();

  // Ids for all entities (all entities have type 0):

  std::vector<EntityId> ids( id_total );

  for ( unsigned i = 0 ; i < id_total ; ++i ) { ids[i] = i + 1; }

  // Declare just those entities in my range of ids:
  //    P0: node ids  1,...,10
  //    P1: node ids 11,...,20
  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    bulk.declare_node(ids[i], no_parts );
  }

  ASSERT_TRUE( bulk.modification_end() );

  // Quick check that the ranks have different sets.
  const unsigned id_get  = ( id_begin + id_total - 2 ) % id_total ;
  ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get] )) );
  ASSERT_TRUE( !bulk.is_valid(bulk.get_entity( stk::topology::NODE_RANK , ids[id_get+1] )) );

  {
    ASSERT_TRUE( bulk.modification_begin() );

    if (0 == p_rank) {
      Entity n10 = bulk.get_entity( stk::topology::NODE_RANK , 10 );
      Entity n11 = bulk.declare_node(11, no_parts);
      bulk.add_node_sharing(n10, 1);
      bulk.add_node_sharing(n11, 1);
    }
    else if ( 1 == p_rank ) {
      Entity n10 = bulk.declare_node(10, no_parts);
      Entity n11 = bulk.get_entity( stk::topology::NODE_RANK , 11 );
      bulk.add_node_sharing(n10, 0);
      bulk.add_node_sharing(n11, 0);
    }

    ASSERT_TRUE(bulk.modification_end());

    {
      Entity n10 = bulk.get_entity( stk::topology::NODE_RANK , 10 );
      Entity n11 = bulk.get_entity( stk::topology::NODE_RANK , 11);
      ASSERT_TRUE( bulk.in_shared(bulk.entity_key(n10)));
      ASSERT_TRUE( bulk.in_shared(bulk.entity_key(n11)));
      ASSERT_TRUE( 0 == bulk.parallel_owner_rank(n10) );
      ASSERT_TRUE( 1 == bulk.parallel_owner_rank(n11) );
    }
  }

  {
    ASSERT_TRUE(bulk.modification_begin());

    if ( 0 == p_rank ) {

      Entity n8 = bulk.get_entity( stk::topology::NODE_RANK , 8 );
      Entity n9 = bulk.get_entity( stk::topology::NODE_RANK , 9 );
      Entity n20 = bulk.declare_node(20, no_parts);
      bulk.add_node_sharing(n9, 1);
      bulk.add_node_sharing(n20, 1);

      Entity e_9_20 = bulk.declare_edge(1, stk::mesh::ConstPartVector{edge_part});
      bulk.declare_relation( e_9_20 , n9 , 0 );
      bulk.declare_relation( e_9_20 , n20 , 1 );

      Entity tri_shell = bulk.declare_element(1, stk::mesh::ConstPartVector{shell_part});
      bulk.declare_relation( tri_shell , n8 , 0 );
      bulk.declare_relation( tri_shell , n9 , 1 );
      bulk.declare_relation( tri_shell , n20 , 2 );
      bulk.declare_relation( tri_shell , e_9_20, 1 );
      bulk.declare_element_side( tri_shell , 0 , stk::mesh::ConstPartVector{tri_part} );
    }
    else if (1 == p_rank) {

      Entity n18 = bulk.get_entity( stk::topology::NODE_RANK , 18 );
      Entity n9 = bulk.declare_node(9, no_parts);
      Entity n20 = bulk.get_entity( stk::topology::NODE_RANK , 20 );
      bulk.add_node_sharing(n9, 0);
      bulk.add_node_sharing(n20, 0);

      Entity e_9_20 = bulk.declare_edge(1, stk::mesh::ConstPartVector{edge_part});
      bulk.declare_relation( e_9_20 , n9 , 0 );
      bulk.declare_relation( e_9_20 , n20 , 1 );

      Entity tri_shell = bulk.declare_element(11, stk::mesh::ConstPartVector{shell_part});
      bulk.declare_relation( tri_shell , n18 , 0 );
      bulk.declare_relation( tri_shell , n9 , 1 );
      bulk.declare_relation( tri_shell , n20 , 2 );
      bulk.declare_relation( tri_shell , e_9_20, 1 );
      bulk.declare_element_side( tri_shell, 0, stk::mesh::ConstPartVector{tri_part} );
    }

    ASSERT_TRUE(bulk.modification_end());

    {
      Entity n9 = bulk.get_entity( stk::topology::NODE_RANK , 9 );
      Entity n20 = bulk.get_entity( stk::topology::NODE_RANK , 20 );
      Entity e_9_20 = bulk.get_entity( stk::topology::EDGE_RANK , 1);
      ASSERT_TRUE(bulk.in_shared(bulk.entity_key(n9)));
      ASSERT_TRUE(bulk.in_shared(bulk.entity_key(n20)));
      ASSERT_TRUE(bulk.in_shared(bulk.entity_key(e_9_20)));
      ASSERT_TRUE( 0 == bulk.parallel_owner_rank(n9) );
      ASSERT_TRUE( 1 == bulk.parallel_owner_rank(n20) );
      ASSERT_TRUE( 0 == bulk.parallel_owner_rank(e_9_20) );
    }
  }

  {
    // Trip an error condition: an edge or face can only be shared
    // if there is a downward relation to it from a locally owned entity.
    ASSERT_TRUE(bulk.modification_begin());
    if ( 0 == p_rank ) {

      Entity n6 = bulk.get_entity( stk::topology::NODE_RANK , 6 );
      Entity n17 = bulk.declare_node(17, no_parts);
      bulk.add_node_sharing(n6, 1);
      bulk.add_node_sharing(n17, 1);

      Entity e_6_17 = bulk.declare_edge(2, stk::mesh::ConstPartVector{edge_part});
      bulk.declare_relation( e_6_17 , n6 , 0 );
      bulk.declare_relation( e_6_17 , n17 , 1 );
    }
    else if (1 == p_rank) {

      Entity n6 = bulk.declare_node(6, no_parts );
      Entity n17 = bulk.get_entity( stk::topology::NODE_RANK , 17 );
      bulk.add_node_sharing(n6, 0);
      bulk.add_node_sharing(n17, 0);

      Entity e_6_17 = bulk.declare_edge(2, stk::mesh::ConstPartVector{edge_part});
      bulk.declare_relation( e_6_17 , n6 , 0 );
      bulk.declare_relation( e_6_17 , n17 , 1 );
    }
    ASSERT_THROW(bulk.modification_end(), std::runtime_error);
  }
}

TEST(BulkData, inducedPartsOnFacesWorks)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);

    if(numProcs == 2)
    {
        const int spatial_dimension = 3;
        MetaData meta(spatial_dimension);
        //Part& node_part = meta.declare_part("node_part", stk::topology::NODE_RANK);
        Part& face_part = meta.declare_part_with_topology("quad", stk::topology::QUAD_4);
        Part& shell_part = meta.declare_part_with_topology("shell", stk::topology::SHELL_QUAD_4);
        Part& messyPart = meta.declare_part_with_topology("messyPart", stk::topology::QUAD_4);
        meta.commit();

        BulkData bulkData(meta, comm);
        bulkData.modification_begin();
        stk::mesh::Entity face;
        if (bulkData.parallel_rank() == 0)
        {
            stk::mesh::Entity node1 = bulkData.declare_node(1);
            bulkData.add_node_sharing(node1, 1);
        }
        else
        {
            stk::mesh::Entity node1 = bulkData.declare_node(1);
            stk::mesh::Entity node2 = bulkData.declare_node(2);
            stk::mesh::Entity node3 = bulkData.declare_node(3);
            stk::mesh::Entity node4 = bulkData.declare_node(4);

            stk::mesh::Entity shell = bulkData.declare_element(1, stk::mesh::ConstPartVector{&shell_part});
            bulkData.declare_relation(shell, node1, 0);
            bulkData.declare_relation(shell, node2, 1);
            bulkData.declare_relation(shell, node3, 2);
            bulkData.declare_relation(shell, node4, 3);
            bulkData.add_node_sharing(node1, 0);
            face = bulkData.declare_element_side(shell, 0, stk::mesh::ConstPartVector{&face_part});
        }

        bulkData.modification_end();

        stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, 1);
        EXPECT_TRUE(bulkData.bucket(node1).shared());

        bulkData.modification_begin();

        if ( bulkData.parallel_rank() == 1 )
        {
            stk::mesh::ConstPartVector addParts;
            addParts.push_back(&messyPart);
            bulkData.change_entity_parts(face, addParts, stk::mesh::ConstPartVector());
        }

        bulkData.modification_end();

        EXPECT_TRUE(bulkData.bucket(node1).member(messyPart));

        bulkData.modification_begin();

        if ( bulkData.parallel_rank() == 1 )
        {
            stk::mesh::ConstPartVector rmParts;
            rmParts.push_back(&messyPart);
            bulkData.change_entity_parts(face, stk::mesh::ConstPartVector(), rmParts);
        }

        bulkData.modification_end();

        EXPECT_TRUE(!bulkData.bucket(node1).member(messyPart));
    }
}

TEST(BulkData, inducedPartsOnFacesThrowsTicket12896)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);

    if(numProcs == 2)
    {
        const int spatial_dimension = 3;
        MetaData meta(spatial_dimension);
        Part& shell_part = meta.declare_part_with_topology("shell", stk::topology::SHELL_QUAD_4);
        Part& face_part = meta.declare_part_with_topology("quad", stk::topology::QUAD_4);
        Part& messyPart = meta.declare_part_with_topology("messyPart", stk::topology::QUAD_4);
        meta.commit();

        BulkData bulkData(meta, comm);
        bulkData.modification_begin();
        stk::mesh::Entity face;
        if (bulkData.parallel_rank() == 0)
        {
            stk::mesh::Entity node1 = bulkData.declare_node(1);
            bulkData.add_node_sharing(node1, 1);
        }
        else
        {
            stk::mesh::Entity node1 = bulkData.declare_node(1);
            stk::mesh::Entity node2 = bulkData.declare_node(2);
            stk::mesh::Entity node3 = bulkData.declare_node(3);
            stk::mesh::Entity node4 = bulkData.declare_node(4);

            stk::mesh::Entity shell = bulkData.declare_element(1, stk::mesh::ConstPartVector{&shell_part});
            bulkData.declare_relation(shell, node1, 0);
            bulkData.declare_relation(shell, node2, 1);
            bulkData.declare_relation(shell, node3, 2);
            bulkData.declare_relation(shell, node4, 3);
            bulkData.add_node_sharing(node1, 0);
            face = bulkData.declare_element_side(shell, 0, stk::mesh::ConstPartVector{&face_part});
        }

        bulkData.modification_end();

        stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, 1);
        EXPECT_TRUE(bulkData.bucket(node1).shared());

        bulkData.modification_begin();

        if ( bulkData.parallel_rank() == 1 )
        {
            stk::mesh::ConstPartVector addParts;
            addParts.push_back(&messyPart);
            bulkData.change_entity_parts(face, addParts, stk::mesh::ConstPartVector());
        }

        if ( bulkData.parallel_rank() == 1 )
        {
            stk::mesh::ConstPartVector rmParts;
            rmParts.push_back(&messyPart);
            bulkData.change_entity_parts(face, stk::mesh::ConstPartVector(), rmParts);
        }

        EXPECT_NO_THROW(bulkData.modification_end());

        EXPECT_TRUE(!bulkData.bucket(node1).member(messyPart));
    }
}

//----------------------------------------------------------------------
TEST(BulkData, testBulkDataRankBeginEnd)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    const int p_size = stk::parallel_machine_size(pm);
    if(p_size != 1)
    {
        return;
    }

    const size_t spatial_dim = 3;
    MetaData meta(spatial_dim, stk::mesh::entity_rank_names());
    BulkData bulk(meta, pm);
    bulk.modification_begin();
    stk::mesh::const_entity_iterator iter = bulk.begin_entities(stk::topology::NODE_RANK);
    stk::mesh::const_entity_iterator end = bulk.end_entities(stk::topology::NODE_RANK);

    ASSERT_TRUE(iter == end);

    EntityId node_id = 1;
    bulk.declare_node(node_id);

    iter = bulk.begin_entities(stk::topology::NODE_RANK);
    end = bulk.end_entities(stk::topology::NODE_RANK);

    //insist that there is 1 node:
    ASSERT_TRUE(iter != end);
    ASSERT_TRUE(std::distance(iter,end) == 1u);

    //now declare an edge...
    EntityId edge_id = 1;
    bulk.declare_edge(edge_id);

    iter = bulk.begin_entities(stk::topology::NODE_RANK);
    end = bulk.end_entities(stk::topology::NODE_RANK);

    //insist that there is still 1 node:
    ASSERT_TRUE(iter != end);
    ASSERT_TRUE(std::distance(iter,end) == 1u);

    iter = bulk.begin_entities(stk::topology::EDGE_RANK);
    end = bulk.end_entities(stk::topology::EDGE_RANK);

    //insist that there is 1 edge:
    ASSERT_TRUE(iter != end);
    ASSERT_TRUE(std::distance(iter,end) == 1u);

    node_id = 2;
    bulk.declare_node(node_id);

    iter = bulk.begin_entities(stk::topology::NODE_RANK);
    end = bulk.end_entities(stk::topology::NODE_RANK);

    //insist that there are 2 nodes:
    ASSERT_TRUE(iter != end);
    ASSERT_TRUE(std::distance(iter,end) == 2u);

    iter = bulk.begin_entities(stk::topology::EDGE_RANK);
    end = bulk.end_entities(stk::topology::EDGE_RANK);

    //insist that there is still 1 edge:
    ASSERT_TRUE(iter != end);
    ASSERT_TRUE(std::distance(iter,end) == 1u);

    iter = bulk.begin_entities(stk::topology::FACE_RANK);
    end = bulk.end_entities(stk::topology::FACE_RANK);

    //insist that there are no faces:
    ASSERT_TRUE(iter == end);
}
//----------------------------------------------------------------------

TEST(BulkData, testChangeOwner_ring)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    enum
    {
        nPerProc = 10
    };
    const int p_rank = stk::parallel_machine_rank(pm);
    const int p_size = stk::parallel_machine_size(pm);
    const unsigned nLocalNode = nPerProc + (1 < p_size ? 1 : 0);
    const unsigned nLocalElement = nPerProc;

    std::vector<size_t> local_count;

    //------------------------------
    {
        RingFixture ring_mesh(pm, nPerProc, false /* no element parts */, stk::mesh::BulkData::NO_AUTO_AURA);
        BulkData & bulk = ring_mesh.m_bulk_data;
        ring_mesh.m_meta_data.commit();

        bulk.modification_begin();
        ring_mesh.generate_mesh();
        ASSERT_TRUE(bulk.modification_end());

        ring_mesh.fixup_node_ownership();

        const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() | ring_mesh.m_meta_data.globally_shared_part();
        const Selector select_all = ring_mesh.m_meta_data.universal_part();

        count_entities(select_used, ring_mesh.m_bulk_data, local_count);
        ASSERT_EQ( local_count[stk::topology::NODE_RANK], nLocalNode);
        ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK], nLocalElement);

        count_entities(select_all, ring_mesh.m_bulk_data, local_count);
        ASSERT_EQ( local_count[stk::topology::NODE_RANK], nLocalNode);
        ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK], nLocalElement);

        if(1 < p_size)
        {
            // Shift ring by two nodes and elements.

            stk::unit_test::test_shift_ring(ring_mesh);

            count_entities(select_used, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement);

            count_entities(select_all, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement);
        }
    }

    //------------------------------
    // Test shift starting with ghosting but not regenerated ghosting.
    {
        RingFixture ring_mesh(pm, nPerProc, false /* no element parts */, stk::mesh::BulkData::NO_AUTO_AURA);
        BulkData& bulk = ring_mesh.m_bulk_data;
        ring_mesh.m_meta_data.commit();

        bulk.modification_begin();
        ring_mesh.generate_mesh();
        ASSERT_TRUE(bulk.modification_end());

        ring_mesh.fixup_node_ownership();

        const Selector select_owned(ring_mesh.m_meta_data.locally_owned_part());
        const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() | ring_mesh.m_meta_data.globally_shared_part();
        const Selector select_all(ring_mesh.m_meta_data.universal_part());

        count_entities(select_used, ring_mesh.m_bulk_data, local_count);
        ASSERT_EQ( local_count[stk::topology::NODE_RANK], nLocalNode);
        ASSERT_EQ( local_count[stk::topology::ELEMENT_RANK], nLocalElement);

        if(1 < p_size)
        {
            stk::unit_test::test_shift_ring(ring_mesh);

            count_entities(select_owned, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nPerProc);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nPerProc);

            count_entities(select_used, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement);

            // All of my ghosts were disrupted and therefore deleted:
            count_entities(select_all, ring_mesh.m_bulk_data, local_count);
            ASSERT_EQ( nLocalElement, local_count[stk::topology::ELEMENT_RANK]);
            ASSERT_EQ( nLocalNode, local_count[stk::topology::NODE_RANK]);
        }
    }
    //------------------------------
    // Test shift starting with ghosting and regenerating ghosting.
    {
        RingFixture ring_mesh(pm, nPerProc, false /* no element parts */);
        BulkData& bulk = ring_mesh.m_bulk_data;
        ring_mesh.m_meta_data.commit();

        bulk.modification_begin();
        ring_mesh.generate_mesh();
        ASSERT_TRUE(bulk.modification_end());

        ring_mesh.fixup_node_ownership();

        const Selector select_owned(ring_mesh.m_meta_data.locally_owned_part());
        const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() | ring_mesh.m_meta_data.globally_shared_part();
        const Selector select_all(ring_mesh.m_meta_data.universal_part());

        count_entities(select_used, ring_mesh.m_bulk_data, local_count);
        ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode);
        ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement);

        count_entities(select_all, ring_mesh.m_bulk_data, local_count);
        const unsigned n_extra = 1 < p_size ? 2 : 0;
        ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode + n_extra);
        ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement + n_extra);

        if(1 < p_size)
        {
            stk::unit_test::test_shift_ring(ring_mesh);

            count_entities(select_owned, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nPerProc);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nPerProc);

            count_entities(select_used, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement);

            // All of my ghosts were regenerated:
            count_entities(select_all, ring_mesh.m_bulk_data, local_count);
            ASSERT_TRUE( local_count[stk::topology::NODE_RANK] == nLocalNode + n_extra);
            ASSERT_TRUE( local_count[stk::topology::ELEMENT_RANK] == nLocalElement + n_extra);
        }
    }
    //------------------------------
    // Test bad owner change catching:
    if(1 < p_size)
    {
        RingFixture ring_mesh(pm, nPerProc, false /* no element parts */);
        BulkData& bulk = ring_mesh.m_bulk_data;
        ring_mesh.m_meta_data.commit();

        bulk.modification_begin();
        ring_mesh.generate_mesh();
        ASSERT_TRUE(bulk.modification_end());

        ring_mesh.fixup_node_ownership();

        std::vector<EntityProc> change;

        if(0 == p_rank)
        {
            change.resize(4);
            // Error to change to bad owner:
            change[0].first = ring_mesh.m_bulk_data.get_entity(stk::topology::NODE_RANK, ring_mesh.m_node_ids[1]);
            change[0].second = p_size;
            // Error to change a ghost:
            for(stk::mesh::EntityCommListInfoVector::const_iterator ec = ring_mesh.m_bulk_data.my_internal_comm_list().begin();
                    ec != ring_mesh.m_bulk_data.my_internal_comm_list().end(); ++ec)
                    {
                if(bulk.in_receive_ghost(ec->key))
                {
                    change[1].first = ec->entity;
                    break;
                }
            }
            change[1].second = p_rank;
            // Error to change to multiple owners:
            change[2].first = ring_mesh.m_bulk_data.get_entity(stk::topology::NODE_RANK, ring_mesh.m_node_ids[1]);
            change[2].second = (p_rank + 1) % p_size;
            change[3].first = change[2].first;
            change[3].second = (p_rank + 2) % p_size;
        }

        ASSERT_THROW( ring_mesh.m_bulk_data.change_entity_owner( change ), std::runtime_error);
    }
    //------------------------------
    // Test move one element with initial ghosting but not regenerated ghosting:
    // last processor give its shared node to P0
    if(1 < p_size)
    {
        RingFixture ring_mesh(pm, nPerProc, false /* no element parts */, stk::mesh::BulkData::NO_AUTO_AURA);
        BulkData& bulk = ring_mesh.m_bulk_data;
        ring_mesh.m_meta_data.commit();

        bulk.modification_begin();
        ring_mesh.generate_mesh();
        ASSERT_TRUE(bulk.modification_end());

        ring_mesh.fixup_node_ownership();

        const Selector select_owned(ring_mesh.m_meta_data.locally_owned_part());
        const Selector select_used = ring_mesh.m_meta_data.locally_owned_part() | ring_mesh.m_meta_data.globally_shared_part();
        const Selector select_all(ring_mesh.m_meta_data.universal_part());

        std::vector<EntityProc> change;

        if(p_rank + 1 == p_size)
        {
            EntityProc entry;
            entry.first = ring_mesh.m_bulk_data.get_entity(stk::topology::NODE_RANK, ring_mesh.m_node_ids[0]);
            entry.second = 0;
            ASSERT_EQ( p_rank, bulk.parallel_owner_rank(entry.first));
            change.push_back(entry);
        }

        ring_mesh.m_bulk_data.change_entity_owner(change);

        count_entities(select_owned, ring_mesh.m_bulk_data, local_count);
        const unsigned n_node = p_rank == 0 ? nPerProc + 1 : (p_rank + 1 == p_size ? nPerProc - 1 : nPerProc );

        ASSERT_EQ( n_node, local_count[stk::topology::NODE_RANK]);
        ASSERT_EQ( static_cast<unsigned>(nPerProc), local_count[stk::topology::ELEMENT_RANK]);

        count_entities(select_used, ring_mesh.m_bulk_data, local_count);
        ASSERT_EQ( nLocalNode, local_count[stk::topology::NODE_RANK]);
        ASSERT_EQ( nLocalElement, local_count[stk::topology::ELEMENT_RANK]);
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing for collection of boxes

TEST(BulkData, testChangeOwner_box)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    const int root_box[3][2] = { {0, 4}, {0, 5}, {0, 6}};

    const int p_size = stk::parallel_machine_size(pm);

    const int spatial_dimension = 3;
    MetaData meta(spatial_dimension);

    meta.commit();

    //------------------------------
    {
        BoxFixture fixture(pm, stk::mesh::BulkData::AUTO_AURA, 100);
        fixture.fem_meta().commit();
        stk::unit_test_util::BulkDataTester & bulk = fixture.bulk_data();
        int local_box[3][2] = { {0, 0}, {0, 0}, {0, 0}};

        bulk.modification_begin();
        fixture.generate_boxes(root_box, local_box);
        ASSERT_TRUE(bulk.modification_end());

        if(1 < p_size)
        {
            donate_one_element(bulk);
        }
    }

    if(1 < p_size)
    {
        BoxFixture fixture(pm, stk::mesh::BulkData::AUTO_AURA, 100);
        fixture.fem_meta().commit();
        stk::unit_test_util::BulkDataTester & bulk = fixture.bulk_data();
        int local_box[3][2] = { {0, 0}, {0, 0}, {0, 0}};

        bulk.modification_begin();
        fixture.generate_boxes(root_box, local_box);
        ASSERT_TRUE(bulk.modification_end());

        donate_all_shared_nodes(bulk);
    }
    //------------------------------
    if(1 < p_size)
    {
        BoxFixture fixture(pm, stk::mesh::BulkData::AUTO_AURA, 100);
        fixture.fem_meta().commit();
        stk::unit_test_util::BulkDataTester & bulk = fixture.bulk_data();
        int local_box[3][2] = { {0, 0}, {0, 0}, {0, 0}};

        bulk.modification_begin();
        fixture.generate_boxes(root_box, local_box);
        ASSERT_TRUE(bulk.modification_end());

        donate_one_element(bulk);
    }
    //------------------------------
    // Introduce ghosts:
    if(1 < p_size)
    {
        BoxFixture fixture(pm, stk::mesh::BulkData::NO_AUTO_AURA, 100);
        stk::unit_test_util::BulkDataTester & bulk = fixture.bulk_data();
        MetaData & box_meta = fixture.fem_meta();
        box_meta.commit();
        int local_box[3][2] = { {0, 0}, {0, 0}, {0, 0}};

        bulk.modification_begin();
        fixture.generate_boxes(root_box, local_box);
        ASSERT_TRUE(bulk.modification_end());

        std::vector<size_t> used_count;
        std::vector<size_t> all_count;

        const Selector select_owned(box_meta.locally_owned_part());
        const Selector select_used = box_meta.locally_owned_part() | box_meta.globally_shared_part();
        const Selector select_all(box_meta.universal_part());

        count_entities(select_all, bulk, all_count);
        count_entities(select_used, bulk, used_count);

        ASSERT_EQ( used_count[0], all_count[0]);
        ASSERT_EQ( used_count[3], all_count[3]);

        donate_all_shared_nodes(bulk);

        count_entities(select_all, bulk, all_count);
        count_entities(select_used, bulk, used_count);

        ASSERT_EQ( used_count[0], all_count[0]);
        ASSERT_EQ( used_count[3], all_count[3]);
    }
}

TEST(BulkData, testModifyPropagation)
{
    // Our new modification model makes it so the modified status
    // of an entity is propagated up to higher-ranked entities
    // that have relations to the modified entity. We test this
    // by grabbing a node off of a ring mesh, modifying it, and
    // checking that its element also gets marked as modified.

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    const unsigned nPerProc = 2;
    const int p_size = stk::parallel_machine_size(pm);

    // this test only needs to be run w/ one processor
    if(p_size > 1)
        return;

    // Make a ring_mesh and add an extra part
    RingFixture ring_mesh(pm, nPerProc, false /* don't use element parts */);
    stk::mesh::Part& special_part = ring_mesh.m_meta_data.declare_part("special_node_part", stk::topology::NODE_RANK);
    ring_mesh.m_meta_data.commit();
    BulkData& bulk = ring_mesh.m_bulk_data;

    bulk.modification_begin();
    ring_mesh.generate_mesh();
    ASSERT_TRUE(bulk.modification_end());

    ring_mesh.fixup_node_ownership();

    // grab the first element
    EntityVector elements;
    const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;
    stk::mesh::get_entities(ring_mesh.m_bulk_data, element_rank, elements);
    stk::mesh::Entity element = elements.front();

    // get one of the nodes related to this element
    ASSERT_TRUE(bulk.num_nodes(element) > 0);
    stk::mesh::Entity node = *bulk.begin_nodes(element);
    ASSERT_EQ( bulk.entity_rank(node), stk::topology::NODE_RANK);

    // make a modification to the node by changing its parts
    ring_mesh.m_bulk_data.modification_begin();
    stk::mesh::ConstPartVector parts;
    parts.push_back(&special_part);
    bulk.change_entity_parts(node, parts);

    // check that the node AND it's element are marked as modified
    ASSERT_EQ( bulk.state(node), stk::mesh::Modified);
    ASSERT_EQ( bulk.state(element), stk::mesh::Modified);

    ASSERT_TRUE( bulk.modification_end());
}

TEST(BulkData, testChangeEntityOwnerFromSelfToSelf)
{
    // It should be legal to "change" entity ownership from yourself to yourself.
    //
    // 1---3---5
    // | 1 | 2 |
    // 2---4---6
    //
    // To test this, we use the mesh above, with elem 1 going on rank 0 and
    // elem 2 going on rank 1. Nodes 3,4 are shared. After the mesh is set up
    // we change the ownership of a few nodes to the same proc that already
    // owns them.

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    meta_data.commit();
    BulkData mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    // Bail if we only have one proc
    if(p_size == 1)
    {
        return;
    }

    // Begin modification cycle so we can create the entities and relations
    mesh.modification_begin();

    EntityVector nodes;
    const unsigned nodes_per_elem = 4, nodes_per_side = 2;

    if(p_rank < 2)
    {
        // We're just going to add everything to the universal part
        stk::mesh::ConstPartVector node_parts, elem_parts;
        stk::mesh::Part &node_part = meta_data.get_topology_root_part(stk::topology::NODE);
        node_parts.push_back(&node_part);
        stk::mesh::Part &quad4_part = meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
        elem_parts.push_back(&quad4_part);

        // Create element
        Entity elem = mesh.declare_element(p_rank + 1, elem_parts);

        // Create nodes
        const unsigned starting_node_id = p_rank * nodes_per_side + 1;
        for(unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id)
        {
            nodes.push_back(mesh.declare_node(id, node_parts));
        }

        // Add relations to nodes
        unsigned rel_id = 0;
        for(EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
        {
            mesh.declare_relation(elem, *itr, rel_id);
        }
        if (p_rank == 0)
        {
            Entity shared_node0 = nodes[2];
            Entity shared_node1 = nodes[3];
            mesh.add_node_sharing(shared_node0, 1);
            mesh.add_node_sharing(shared_node1, 1);
        }
        else
        {
            Entity shared_node0 = nodes[0];
            Entity shared_node1 = nodes[1];
            mesh.add_node_sharing(shared_node0, 0);
            mesh.add_node_sharing(shared_node1, 0);
        }
    }

    mesh.modification_end();

    std::vector<EntityProc> change;
    if(p_rank < 2)
    {
        // Change ownership of some nodes to the same proc that owns them

        // Add a non-shared node to change list
        if(p_rank == 0)
        {
            EntityProc entry(nodes.front(), p_rank);
            change.push_back(entry);
        }
        else
        {
            EntityProc entry(nodes.back(), p_rank);
            change.push_back(entry);
        }

        // Add a shared node to change list
        Entity shared_node = nodes[p_rank == 0 ? nodes_per_side : 0];
        EntityId expected_id = 3;
        Part& shared_part = meta_data.globally_shared_part();
        ASSERT_TRUE( has_superset(mesh.bucket(shared_node), shared_part));
        ASSERT_EQ(mesh.identifier(shared_node), expected_id);
        if(mesh.parallel_owner_rank(shared_node) == p_rank)
        {
            EntityProc entry(shared_node, p_rank);
            change.push_back(entry);
        }
    }

    mesh.change_entity_owner(change);
}

TEST(BulkData, test_internal_clean_and_verify_parallel_change_trivial)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);
    std::vector<EntityProc> local_change;
    stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change);
    EXPECT_TRUE(local_change.size() == 0);
}
TEST(BulkData, test_internal_clean_and_verify_parallel_change_sort_unique)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    int myRank = stk::parallel_machine_rank(pm);
    if (numProcs < 2) { return; }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);
    mesh.modification_begin();
    Entity node1, node2, node3, node4;
    if (myRank == 0) {
        node1 = mesh.declare_node(1);
        node2 = mesh.declare_node(2);
        node3 = mesh.declare_node(3);
        node4 = mesh.declare_node(4);
    }
    mesh.modification_end();
    std::vector<EntityProc> local_change;
    if (myRank == 0)
    {
        local_change.push_back(EntityProc(node4,0));
        local_change.push_back(EntityProc(node3,1));
        local_change.push_back(EntityProc(node2,1));
        local_change.push_back(EntityProc(node2,1));
        local_change.push_back(EntityProc(node1,1));
    }
    stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change);
    if (myRank == 0)
    {
        EXPECT_EQ(3u,local_change.size());
        EXPECT_TRUE( EntityProc(node1,1) == local_change[0] );
        EXPECT_TRUE( EntityProc(node2,1) == local_change[1] );
        EXPECT_TRUE( EntityProc(node3,1) == local_change[2] );
    }
}
TEST(BulkData, test_internal_clean_and_verify_parallel_change_bad_null)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);
    std::vector<EntityProc> local_change;
    local_change.push_back(EntityProc(Entity(),1));
    EXPECT_THROW( stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change),
                  std::runtime_error );
}
TEST(BulkData, test_internal_clean_and_verify_parallel_change_not_owner)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    int myRank = stk::parallel_machine_rank(pm);
    if (numProcs < 2) { return; }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
    stk::mesh::BulkData mesh(meta, pm);
    mesh.modification_begin();
    Entity node1, node2, node3, node4, node5, node6;
    if (myRank == 0)
    {
        Entity element = mesh.declare_element(1, stk::mesh::ConstPartVector{&block_1});
        node1 = mesh.declare_node(1);
        node2 = mesh.declare_node(2); // shared
        node3 = mesh.declare_node(3); // shared
        node4 = mesh.declare_node(4);
        mesh.declare_relation(element,node1,0);
        mesh.declare_relation(element,node2,1);
        mesh.declare_relation(element,node3,2);
        mesh.declare_relation(element,node4,3);
        mesh.add_node_sharing(node2,1);
        mesh.add_node_sharing(node3,1);
    }
    if (myRank == 1)
    {
        Entity element = mesh.declare_element(2, stk::mesh::ConstPartVector{&block_1});
        node2 = mesh.declare_node(2); // shared
        node5 = mesh.declare_node(5);
        node6 = mesh.declare_node(6);
        node3 = mesh.declare_node(3); // shared
        mesh.declare_relation(element,node2,0);
        mesh.declare_relation(element,node5,1);
        mesh.declare_relation(element,node6,2);
        mesh.declare_relation(element,node3,3);
        mesh.add_node_sharing(node2,0);
        mesh.add_node_sharing(node3,0);
    }
    mesh.modification_end();

    std::vector<EntityProc> local_change;
    if (myRank == 1)
    {
        EXPECT_FALSE( mesh.bucket(node2).owned() );
        local_change.push_back(EntityProc(node2,1));
    }
    EXPECT_THROW( stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change),
                  std::runtime_error );
}
TEST(BulkData, test_internal_clean_and_verify_parallel_change_invalid_owner)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    int myRank = stk::parallel_machine_rank(pm);
    if (numProcs < 2) { return; }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);
    mesh.modification_begin();
    Entity node1;
    if (myRank == 0)
    {
        node1 = mesh.declare_node(1);
    }
    mesh.modification_end();
    std::vector<EntityProc> local_change;
    if (myRank == 0)
    {
        local_change.push_back(EntityProc(node1,numProcs));
    }
    EXPECT_THROW( stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change),
                  std::runtime_error );
}
TEST(BulkData, test_internal_clean_and_verify_parallel_change_send_to_2_owners)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    int myRank = stk::parallel_machine_rank(pm);
    if (numProcs < 3) { return; }

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);
    mesh.modification_begin();
    Entity node1;
    if (myRank == 0) {
        node1 = mesh.declare_node(1);
    }
    mesh.modification_end();
    std::vector<EntityProc> local_change;
    if (myRank == 0)
    {
        local_change.push_back(EntityProc(node1,1));
        local_change.push_back(EntityProc(node1,2));
    }
    EXPECT_THROW( stk::mesh::impl::internal_clean_and_verify_parallel_change(mesh,local_change),
                  std::runtime_error );
}

TEST(BulkData, test_internal_generate_parallel_change_lists_trivial)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;

    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    std::vector<EntityProc> local_change;
    std::vector<EntityProc> ghosted_change;
    std::vector<EntityProc> shared_change;
    stk::mesh::impl::internal_generate_parallel_change_lists(mesh, local_change,
                                                       shared_change, ghosted_change);
    EXPECT_TRUE(shared_change.empty());
    EXPECT_TRUE(ghosted_change.empty());
}

TEST(BulkData, test_internal_generate_parallel_change_lists_2EltsChown1ChownItsNodes)
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
  stk::mesh::EntityIdVector elem_node_ids[] {
      {1, 2, 3, 4},
      {4, 3, 6, 5}
  };

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  std::vector<stk::mesh::Entity> elems;
  bulk.modification_begin();
  if (p_rank == 0) {
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] ) );
    elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec local_change;
  if (p_rank == 0) {
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 3), 1));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
    local_change.push_back(stk::mesh::EntityProc(elems[1], 1));
  }

  std::vector<EntityProc> shared_change;
  std::vector<EntityProc> ghosted_change;
  stk::mesh::impl::internal_generate_parallel_change_lists(bulk, local_change, shared_change, ghosted_change);

  EXPECT_TRUE(shared_change.empty());
  EXPECT_TRUE(ghosted_change.empty());
}

TEST(BulkData, test_internal_generate_parallel_change_lists_2EltsFlip)
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
  //   1/0---4/0---5/1      1/1---4/1---5/0
  //    |     |     |        |     |     |
  //    | 1/0 | 2/1 |   =>   | 1/1 | 2/0 |
  //    |     |     |        |     |     |
  //   2/0---3/0---6/1      2/1---3/0---6/0

  stk::mesh::EntityId element_ids [2] = {1, 2};
  stk::mesh::EntityIdVector elem_node_ids [] = {
      {1, 2, 3, 4},
      {4, 3, 6, 5}
  };

  stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part",stk::topology::QUAD_4_2D);
  meta.commit();

  // Start with all entities on proc 0
  Entity elem;
  bulk.modification_begin();
  if (p_rank == 0) {
    elem = stk::mesh::declare_element(bulk, elem_part ,element_ids[0], elem_node_ids[0] );
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 3), 1);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 4), 1);
  }
  else if (p_rank == 1) {
    elem = stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] );
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 3), 0);
    bulk.add_node_sharing(bulk.get_entity(stk::topology::NODE_RANK, 4), 0);
  }
  bulk.modification_end();

  stk::mesh::EntityProcVec local_change;
  if (p_rank == 0) {
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 1), 1));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 2), 1));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    local_change.push_back(stk::mesh::EntityProc(elem, 1));
  }
  else if (p_rank == 1) {
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 0));
    local_change.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 0));
    local_change.push_back(stk::mesh::EntityProc(elem, 0));
  }

  std::vector<EntityProc> shared_change;
  std::vector<EntityProc> ghosted_change;
  stk::mesh::impl::internal_generate_parallel_change_lists(bulk, local_change, shared_change, ghosted_change);


  if (p_rank == 0) {
    ASSERT_EQ(0u, shared_change.size());
    ASSERT_EQ(3u, ghosted_change.size());
    EntityProc node5_ghosted = ghosted_change[0];
    EntityProc node6_ghosted = ghosted_change[1];
    EntityProc elt2_ghosted  = ghosted_change[2];
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::NODE_RANK, 5), bulk.entity_key(node5_ghosted.first));
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::NODE_RANK, 6), bulk.entity_key(node6_ghosted.first));
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::ELEMENT_RANK, 2), bulk.entity_key(elt2_ghosted.first));
    EXPECT_EQ(0, node5_ghosted.second);
    EXPECT_EQ(0, node6_ghosted.second);
    EXPECT_EQ(0, elt2_ghosted.second);
  }
  else if (p_rank == 1) {
    ASSERT_EQ(1u, shared_change.size());
    ASSERT_EQ(3u, ghosted_change.size());
    EntityProc node4_shared  = shared_change[0];
    EntityProc node1_ghosted = ghosted_change[0];
    EntityProc node2_ghosted = ghosted_change[1];
    EntityProc elt1_ghosted  = ghosted_change[2];
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::NODE_RANK, 4), bulk.entity_key(node4_shared.first));
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::NODE_RANK, 1), bulk.entity_key(node1_ghosted.first));
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::NODE_RANK, 2), bulk.entity_key(node2_ghosted.first));
    EXPECT_EQ(stk::mesh::EntityKey(stk::topology::ELEMENT_RANK, 1), bulk.entity_key(elt1_ghosted.first));
    EXPECT_EQ(1, node4_shared.second);
    EXPECT_EQ(1, node1_ghosted.second);
    EXPECT_EQ(1, node2_ghosted.second);
    EXPECT_EQ(1, elt1_ghosted.second);
  }
}

TEST(BulkData, testFamilyTreeGhosting)
{
    // A family tree is a higher-rank entity (rank = element_rank() + 1) that
    // has down-relations to elements used, for example, to hold parent/child
    // relations in an adapted mesh.
    //
    // 1---3---5---7
    // | 1 | 2 | 3 | ...
    // 2---4---6---8
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc.
    // After the mesh is set up we add rank-3 (family tree) entities and have them point down to
    // just the single rank-2 elements.  Then we check that they are properly
    // ghosted after modification_end.

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;

    std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
    entity_rank_names.push_back("FAMILY_TREE");

    MetaData meta_data(spatial_dim, entity_rank_names);
    const unsigned nodes_per_elem = 4, nodes_per_side = 2;
    Part &elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    meta_data.commit();
    BulkData mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    Part& owned = meta_data.locally_owned_part();
    Part& shared = meta_data.globally_shared_part();

    //
    // Begin modification cycle so we can create the entities and relations
    //

    mesh.modification_begin();

    EntityVector nodes;
    const EntityRank family_tree_rank = static_cast<EntityRank>(stk::topology::ELEMENT_RANK + 1);
    const EntityId my_family_tree_id = p_rank + 1;

    // We're just going to add everything to the universal part
    stk::mesh::ConstPartVector empty_parts;
    stk::mesh::ConstPartVector elem_parts;
    elem_parts.push_back(&elem_part);

    // Create element
    Entity elem = mesh.declare_element(p_rank + 1, elem_parts);

    // Create nodes
    const unsigned starting_node_id = p_rank * nodes_per_side + 1;
    for(unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id)
    {
        nodes.push_back(mesh.declare_node(id, empty_parts));
    }
    if(p_rank > 0)
    {
        mesh.add_node_sharing(nodes[0], p_rank - 1);
        mesh.add_node_sharing(nodes[1], p_rank - 1);
    }
    if(p_rank < (p_size - 1))
    {
        mesh.add_node_sharing(nodes[2], p_rank + 1);
        mesh.add_node_sharing(nodes[3], p_rank + 1);
    }

    // Add relations to nodes
    unsigned rel_id = 0;
    for(EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
    {
        mesh.declare_relation(elem, *itr, rel_id);
    }

    // Create family tree
    Entity family_tree = mesh.declare_constraint(my_family_tree_id, empty_parts);
    // Add relation to element
    unsigned downward_ordinal = 0; // we only have 1 down relation, it has ordinal 0
    mesh.declare_relation(family_tree, elem, downward_ordinal);

    mesh.modification_end();

    //
    // Test correctness of ghosting: Check that adjacent family-trees are ghosted on this proc
    //

    // Compute and store ids of adjacent family-trees
    std::vector<EntityId> family_tree_ghost_ids;
    if(p_rank > 0)
    {
        family_tree_ghost_ids.push_back(my_family_tree_id - 1);
    }
    if(p_rank < p_size - 1)
    {
        family_tree_ghost_ids.push_back(my_family_tree_id + 1);
    }

    // Check that my_family_tree exists and I own it
    Entity my_family_tree = mesh.get_entity(family_tree_rank, my_family_tree_id);
    ASSERT_TRUE(mesh.is_valid(my_family_tree));
    ASSERT_TRUE( (p_rank) == mesh.parallel_owner_rank(my_family_tree));

    // Check that adjacent family-trees exist and are ghosted
    for(std::vector<EntityId>::const_iterator itr = family_tree_ghost_ids.begin(); itr != family_tree_ghost_ids.end(); ++itr)
    {
        int expected_ghosted_family_tree_id = *itr;

        Entity expected_ghosted_family_tree = mesh.get_entity(family_tree_rank, expected_ghosted_family_tree_id);
        ASSERT_TRUE(mesh.is_valid(expected_ghosted_family_tree));
        ASSERT_TRUE(expected_ghosted_family_tree_id - 1 == mesh.parallel_owner_rank(expected_ghosted_family_tree));

        stk::mesh::Bucket& bucket = mesh.bucket(expected_ghosted_family_tree);
        ASSERT_TRUE(!bucket.member(owned) && !bucket.member(shared));
    }
}

TEST(BulkData, testChangeEntityPartsOfShared)
{
    //
    // This unit-test is designed to test what happens when a shared entity
    // is moved on one processor during the same modification cycle in which
    // it was declared.
    //
    //   p0  p1
    // 1---3---5
    // | 1 | 2 |
    // 2---4---6
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. Node 3 is the node we'll be testing.
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    const EntityRank node_rank = stk::topology::NODE_RANK;

    stk::mesh::Part& extra_node_part = meta_data.declare_part("extra_node_part", node_rank);
    meta_data.commit();

    BulkData mesh(meta_data, pm);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    // Bail unless in parallel
    if(p_size == 1)
    {
        return;
    }

    // Begin modification cycle so we can create the entities and relations
    if(p_rank < 2)
    {
        mesh.modification_begin();

        const unsigned nodes_per_elem = 4, nodes_per_side = 2;
        EntityKey node_key_to_move(node_rank, 3 /*id*/);

        // We're just going to add everything to the universal part
        stk::mesh::ConstPartVector empty_parts, elem_parts;
        stk::mesh::Part &quad4_part = meta_data.get_topology_root_part(stk::topology::QUAD_4_2D);
        elem_parts.push_back(&quad4_part);

        // Create element
        Entity elem = mesh.declare_element(p_rank + 1, elem_parts);

        // Create nodes
        EntityVector nodes;
        const unsigned starting_node_id = p_rank * nodes_per_side + 1;
        for(unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id)
        {
            nodes.push_back(mesh.declare_node(id, empty_parts));
        }

        // Add relations to nodes
        unsigned rel_id = 0;
        for(EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
        {
            mesh.declare_relation(elem, *itr, rel_id);
        }

        // On the processor that does *not* end up as the owner of the node, change its parts
        Entity changing_node = mesh.get_entity(node_key_to_move);
        if(p_rank == 1)
        {
            PartVector add_parts(1, &extra_node_part);
            mesh.change_entity_parts(changing_node, add_parts);
        }
        if (p_rank == 0)
        {
            Entity shared_node0 = nodes[2];
            Entity shared_node1 = nodes[3];
            mesh.add_node_sharing(shared_node0, 1);
            mesh.add_node_sharing(shared_node1, 1);
        }
        else
        {
            Entity shared_node0 = nodes[0];
            Entity shared_node1 = nodes[1];
            mesh.add_node_sharing(shared_node0, 0);
            mesh.add_node_sharing(shared_node1, 0);
        }
        mesh.modification_end();

        // Expect that this is a shared node
        {
            std::vector<int> shared_procs;
            mesh.comm_shared_procs(mesh.entity_key(changing_node),shared_procs);
            EXPECT_FALSE(shared_procs.empty());
        }

        EXPECT_TRUE(mesh.bucket(changing_node).member(extra_node_part));

        mesh.modification_begin();

        // On the processor that owns the node, change its parts
        if(p_rank == 0)
        {
            PartVector add_parts(1, &extra_node_part);
            mesh.change_entity_parts(changing_node, add_parts);
        }

        mesh.modification_end();

        // Expect that the part change *did* have an impact
        EXPECT_TRUE(mesh.bucket(changing_node).member(extra_node_part));
    }
    else
    {
        // On extra procs, do bare minimum
        mesh.modification_begin();
        mesh.modification_end();
        mesh.modification_begin();
        mesh.modification_end();
    }
}

bool is_entity_key_shared(const stk::mesh::BulkData & mesh, stk::mesh::EntityKey key) {
    std::vector<int> shared_procs;
    mesh.comm_shared_procs(key,shared_procs);
    return !shared_procs.empty();
}

void testParallelSideCreation(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
    //
    // This unit-test is designed to test what happens when a shared sides are created on
    // both processors that share the side with different global ids.  Then synchonization
    // is done as a second step.
    //
    // 1---4---5
    // | 1 | 2 |
    // 2---3---6

    // element 1 conn = { 1, 2, 3, 4 }
    // element 2 conn = { 3, 4, 5, 6 }

    // edge is nodes 3-4 (element 1, side id 2, perm 0)
    // edge is nodes 3-4 (element 2, side id 0, perm 0)
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. Node 3 is the node we'll be testing.
    //

    stk::ParallelMachine pm = MPI_COMM_WORLD;

    // Set up meta and bulk data
    const unsigned spatial_dim = 2;
    MetaData meta_data(spatial_dim);
    const EntityRank node_rank = stk::topology::NODE_RANK;
    const EntityRank side_rank = stk::topology::EDGE_RANK;

    stk::mesh::Part& side_part = meta_data.declare_part_with_topology("side_part", stk::topology::LINE_2);
    stk::mesh::Part& elem_part = meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);

    meta_data.commit();

    BulkData mesh(meta_data, pm, autoAuraOption);
    int p_rank = mesh.parallel_rank();
    int p_size = mesh.parallel_size();

    // Bail unless in parallel
    if(p_size == 1)
    {
        return;
    }

    // Begin modification cycle so we can create the entities and relations
    if(p_rank < 2)
    {
        mesh.modification_begin();

        const unsigned nodes_per_elem = 4, nodes_per_side = 2;
        EntityKey node_key_to_move(node_rank, 3 /*id*/);

        // We're just going to add everything to the universal part
        stk::mesh::ConstPartVector empty_parts;

        // Create nodes
        EntityVector nodes;
        const unsigned starting_node_id = p_rank * nodes_per_side + 1;
        for(unsigned id = starting_node_id; id < starting_node_id + nodes_per_elem; ++id)
        {
            nodes.push_back(mesh.declare_node(id, empty_parts));
        }

        // 1-2, 1
        // 2-3, 2
        // 3-4, 3
        // 4-1, 4

        // 3-4, 1
        // 4-5, 2
        // 5-6, 3
        // 6-3, 4

        // Create element
        const EntityId elem_id = p_rank + 1;
        Entity elem = mesh.declare_element(elem_id, stk::mesh::ConstPartVector{&elem_part});

        // Add element relations to nodes
        unsigned elem_rel_id = 0;
        for(EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++elem_rel_id)
        {
            mesh.declare_relation(elem, *itr, elem_rel_id);
        }

        // Add side relations to nodes and element
        EntityVector side_nodes;
        side_nodes.push_back(mesh.get_entity(stk::topology::NODE_RANK, 3));
        side_nodes.push_back(mesh.get_entity(stk::topology::NODE_RANK, 4));
        mesh.add_node_sharing(side_nodes[0], (p_rank == 0 ? 1 : 0));
        mesh.add_node_sharing(side_nodes[1], (p_rank == 0 ? 1 : 0));
        stk::topology elem_top = mesh.bucket(elem).topology();
        unsigned local_side_ordinal = 2;
        if (p_rank == 1)
        {
            local_side_ordinal = 0;
        }

        // Create local version of side on each proc
        Entity side = mesh.declare_element_side(elem, local_side_ordinal, stk::mesh::ConstPartVector{&side_part});

        stk::mesh::Permutation perm1 = mesh.find_permutation(elem_top, &nodes[0], elem_top.side_topology(local_side_ordinal), &side_nodes[0], local_side_ordinal);
        ASSERT_TRUE(perm1 != stk::mesh::Permutation::INVALID_PERMUTATION);
        mesh.modification_end();

        // Expect that the side is not shared, but the nodes of side are shared
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side)));
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side_nodes[0])));
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side_nodes[1])));

        // Now "detect" that there is a duplicate aura side using the side nodes
        EntityVector sides;
        get_entities_through_relations(mesh, side_nodes, side_rank, sides);
        EXPECT_EQ(1u, sides.size());

        mesh.modification_begin();

        // Delete the local side and create new, shared side
        side = sides[0];
        bool destroyrelationship = mesh.destroy_relation(elem, side, local_side_ordinal);
        EXPECT_TRUE(destroyrelationship);
        mesh.modification_end();
        //must call this here to delete ghosts, kills relationship between side and ghost of elem on other proc, allows side to be deleted in next phase
        mesh.modification_begin();
        bool successfully_destroyed = mesh.destroy_entity(side);
        if (p_rank == 0) {
            EXPECT_TRUE(successfully_destroyed);
        }
        else {
            EXPECT_FALSE(successfully_destroyed);
        }
        side = mesh.declare_element_side(elem, local_side_ordinal, stk::mesh::ConstPartVector{&side_part});

        stk::mesh::Permutation perm2 = mesh.find_permutation(elem_top, &nodes[0], elem_top.side_topology(local_side_ordinal), &side_nodes[0], local_side_ordinal);
        ASSERT_TRUE(perm2 != stk::mesh::Permutation::INVALID_PERMUTATION);

        mesh.modification_end();

        // Expect that the side is shared, and nodes of side are shared
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side)));
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side_nodes[0])));
        EXPECT_TRUE(is_entity_key_shared(mesh,mesh.entity_key(side_nodes[1])));

        // Check that there is only a single side using the side nodes
        get_entities_through_relations(mesh, side_nodes, side_rank, sides);
        EXPECT_EQ(1u, sides.size());
    }
    else
    {
        // On extra procs, do bare minimum, collective calls must be made on each proc
        mesh.modification_begin();
        mesh.modification_end();
        mesh.modification_begin();
        mesh.modification_end();
        mesh.modification_begin();
        mesh.modification_end();
    }
}

TEST(BulkData, testParallelSideCreationWithAura)
{
    testParallelSideCreation(stk::mesh::BulkData::AUTO_AURA);
}

TEST(BulkData, testParallelSideCreationWithoutAura)
{
    testParallelSideCreation(stk::mesh::BulkData::NO_AUTO_AURA);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of field_data_footprint(.)
TEST(BulkData, test_total_field_data_footprint )
{
    // Test 3x1x1 HexFixture structure
    const unsigned NX = 3;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
    hf.m_meta.commit();
    hf.generate_mesh();

    const stk::mesh::BulkData &mesh = hf.m_bulk_data;

    // Call function we're testing
    size_t field_data_footprint = mesh.total_field_data_footprint(stk::topology::NODE_RANK);

    // Alternative computation explicitly gathers buckets.
    size_t node_fields_footprint = 0;
    const stk::mesh::BucketVector &node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    for(size_t i = 0; i < node_buckets.size(); ++i)
    {
        node_fields_footprint += node_buckets[i]->capacity() * field_bytes_per_entity(hf.m_coord_field, *node_buckets[i]);
    }

    EXPECT_EQ(node_fields_footprint, field_data_footprint);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of get_buckets and get_entities functions
TEST(BulkData, test_get_entities )
{
    // Test 3x4x4 HexFixture structure
    const unsigned NX = 3;
    const unsigned NY = 4;
    const unsigned NZ = 40;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);

    hf.m_meta.commit();
    hf.generate_mesh();
    const stk::mesh::BulkData &mesh = hf.m_bulk_data;

    Selector select_owned(mesh.mesh_meta_data().locally_owned_part());
    const stk::mesh::BucketVector &bucket_ptrs = mesh.get_buckets(stk::topology::NODE_RANK, select_owned);
    stk::mesh::EntityVector entities;
    mesh.get_entities(stk::topology::NODE_RANK, select_owned, entities);
    //
    //  Confirm that the number of entities exracted by either bucket or entity access is identical
    //
    int numBucketEntities = 0;
    std::map<stk::mesh::Entity, int> entityMap;
    for(unsigned int ibucket = 0; ibucket < bucket_ptrs.size(); ++ibucket)
    {
        numBucketEntities += bucket_ptrs[ibucket]->size();
        for(unsigned iobj = 0; iobj < bucket_ptrs[ibucket]->size(); ++iobj)
        {
            entityMap[(*bucket_ptrs[ibucket])[iobj]] = 1;
        }
    }
    int numEntities = entities.size();
    EXPECT_EQ(numBucketEntities, numEntities);
    //
    //  Confirm that the actual contents of the entity lists are identical
    //
    for(unsigned int iobj = 0; iobj < entities.size(); ++iobj)
    {
        ASSERT_TRUE(entityMap.find(entities[iobj]) != entityMap.end());
    }
    //
    //  confirm the total number of entities is the expected (41*5*4) = 820, the total number of unique nodes in the mesh
    //
    int globalNumEntities = numEntities;
    stk::all_reduce(MPI_COMM_WORLD, stk::ReduceSum<1>(&globalNumEntities));

    EXPECT_EQ(globalNumEntities, 820);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing of communicate_field_data

typedef stk::mesh::Field<int> PressureFieldType;

static void test_sync_1(stk::mesh::BulkData& eMesh, PressureFieldType& pressure_field, bool sync_shared, bool sync_aura)
{
    unsigned p_rank = eMesh.parallel_rank();
    unsigned p_size = eMesh.parallel_size();
    static_cast<void>(p_size);

    const stk::mesh::BucketVector & buckets = eMesh.buckets(stk::topology::NODE_RANK);

    enum Type
    {
        Owned, Shared, Ghost
    };

    for(stk::mesh::BucketVector::const_iterator k = buckets.begin(); k != buckets.end(); ++k)
    {
        {
            stk::mesh::Bucket & bucket = **k;

            const unsigned num_elements_in_bucket = bucket.size();

            for(unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
            {
                stk::mesh::Entity entity = bucket[iEntity];
                int * const p = stk::mesh::field_data(pressure_field, entity);
                stk::mesh::EntityId id = eMesh.identifier(entity);

                if(bucket.owned())
                {
                    p[0] = (p_rank + 1) * 100 + id;
                }
                else if(bucket.shared())
                {
                    p[0] = -((eMesh.parallel_owner_rank(entity) + 1) * 100 + id);
                }
                else
                {
                    p[0] = ((p_rank + 1) * 1000 + id);
                }

            }
        }
    }

    {
        std::vector<const stk::mesh::FieldBase *> fields;
        fields.push_back(&pressure_field);

        // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
        if(sync_aura)
            stk::mesh::communicate_field_data(eMesh.aura_ghosting(), fields);

        // the shared part (just the shared boundary)
        if(sync_shared)
            stk::mesh::copy_owned_to_shared(eMesh, fields);
    }

    for(stk::mesh::BucketVector::const_iterator k = buckets.begin(); k != buckets.end(); ++k)
    {
        {
            stk::mesh::Bucket & bucket = **k;

            const unsigned num_elements_in_bucket = bucket.size();

            for(unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
            {
                stk::mesh::Entity entity = bucket[iEntity];
                stk::mesh::EntityId id = eMesh.identifier(entity);
                int * const p = stk::mesh::field_data(pressure_field, entity);
                double p_e = (p_rank + 1) * 100 + id;
                if(bucket.owned())
                {
                    ASSERT_EQ(p[0], p_e);
                }
                else if(bucket.shared())
                {
                    p_e = ((eMesh.parallel_owner_rank(entity) + 1) * 100 + id);
                    if(sync_shared)
                    {
                        ASSERT_EQ(p[0], p_e);
                    }
                }
                else
                {
                    p_e = ((eMesh.parallel_owner_rank(entity) + 1) * 100 + id);
                    if(sync_aura)
                    {
                        ASSERT_EQ(p[0], p_e);
                    }
                }
            }
        }
    }
}

TEST(BulkData, testFieldComm)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    // run this with exercise_field_sync_bug = true, and 3 <= nprocs <= 4 to show the possible bug
    bool exercise_field_sync_bug = true;

    const unsigned p_size = stk::parallel_machine_size(pm);

    const int spatial_dimension = 3;
    MetaData meta(spatial_dimension);

    meta.commit();

    //------------------------------
    // 3d seems to be fine...
    if(p_size <= 4)
    {
        const int root_box[3][2] = { {0, 2}, {0, 2}, {0, 1}};  // emulate 2d box

        BoxFixture fixture(pm, stk::mesh::BulkData::AUTO_AURA, 100);
        PressureFieldType& p_field = fixture.fem_meta().declare_field<PressureFieldType>(stk::topology::NODE_RANK, "p");
        stk::mesh::put_field_on_mesh(p_field, fixture.fem_meta().universal_part(),
                                     (stk::mesh::FieldTraits<PressureFieldType>::data_type*) nullptr);
        fixture.fem_meta().commit();
        BulkData & bulk = fixture.bulk_data();
        int local_box[3][2] = { {0, 0}, {0, 0}, {0, 0}};

        bulk.modification_begin();
        fixture.generate_boxes(root_box, local_box);
        bulk.modification_end();

        {
            bool shared_aura = false;
            bool shared = false;
            test_sync_1(bulk, p_field, shared, shared_aura);
            test_sync_1(bulk, p_field, false, true);
            if(exercise_field_sync_bug || p_size <= 2)
            {
                test_sync_1(bulk, p_field, true, false);
                test_sync_1(bulk, p_field, true, true);
            }
        }
    }

    //------------------------------
    // 2d, not so much
    if(p_size <= 4)
    {
        stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
        PressureFieldType& p_field = fixture.m_meta.declare_field<PressureFieldType>(stk::topology::NODE_RANK, "p");
        stk::mesh::put_field_on_mesh(p_field, fixture.m_meta.universal_part(),
                                     (stk::mesh::FieldTraits<PressureFieldType>::data_type*) nullptr);
        fixture.m_meta.commit();
        fixture.generate_mesh();
        stk::mesh::BulkData & bulk = fixture.m_bulk_data;

        {
            bool shared_aura = false;
            bool shared = false;
            test_sync_1(bulk, p_field, shared, shared_aura);
            test_sync_1(bulk, p_field, false, true);
            if(exercise_field_sync_bug || p_size <= 2)
            {
                test_sync_1(bulk, p_field, true, false);
                test_sync_1(bulk, p_field, true, true);
            }
        }
    }
}

// testing comm lists and custom ghosting

TEST(BulkData, testCommList)
{
    /**
     * This is a boiled-down version of a stk_adapt situation that is failing
     *   where a custom ghosted node is later shared because it becomes part
     *   of an element (this occurs during hanging-node refinement).  It is
     *   a definite edge case, but exposed an assumption in the comm_mesh_verify_parallel_consistency()
     *   function, that a comm list can't have both a shared and custom ghosted node.
     * This test currently just verifies the issue exists.  When comm_mesh_verify_parallel_consistency()
     *   is fixed, this test should be modified to enforce no failure under debug mode.
     *
     *  Mesh
     *    7---8---9  P0 owns nodes 1,2,4,5; P, elem 1
     *    | 3 | 4 |  P1 : 3,6, elem 2
     *    4---5---6  P2 : 7,8, elem 3
     *    | 1 | 2 |  P3 : 9,   elem 4
     *    1---2---3
     *
     *  node 5 ghosted to proc 3, node 9 ghosted to proc 0
     *
     *  Comm list looks like this after the ghosting operations (obtained from print_comm_list()):
     *    legend: (ghost_id, proc)
     *    P3: NODE[9] owner(3)     (1,0) (1,1) (1,2) (2,0)
     *    P0: NODE[9] owner(3)     (1,3) (2,3)
     *
     *  elem 1 modified to replace relation elem[1] -> node[1] with node[9] (previously ghosted)
     *
     *  This last step induces a comm list like this
     *
     *    P3: NODE[9] owner(3) mod (0,0) (1,1) (1,2) (2,0)
     *    P0: NODE[9] owner(3) mod (0,3) (2,3)
     *
     *  Note the repetition of proc 3 in both the shared (ghost_id=0) and ghosted (id=2) list
     *  for P0 (and repetition of proc 0 in P3's list).  This causes the
     *  comm_mesh_verify_parallel_consistency() test to fail, although
     *  this should be a valid mesh, e.g., for periodic b.c.'s we might want this situation.
     *
     */

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    MPI_Barrier( pm);

    const unsigned p_size = stk::parallel_machine_size(pm);
    const unsigned p_rank = stk::parallel_machine_rank(pm);

    //------------------------------
    if(p_size != 4)
        return;

    //------------------------------
    // test begin/end pair
    {
        stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
        fixture.m_meta.commit();
        fixture.generate_mesh();
        stk::mesh::BulkData & bulk = fixture.m_bulk_data;
        bulk.modification_begin();
        bulk.modification_end();
    }

    //------------------------------
    // test begin/end pair with mesh mods
    {
        stk::mesh::fixtures::QuadFixture fixture(pm, 2 /*nx*/, 2 /*ny*/);
        fixture.m_meta.commit();
        fixture.generate_mesh();
        stk::mesh::BulkData & bulk = fixture.m_bulk_data;
        bulk.modification_begin();
        bulk.modification_end();

        bulk.modification_begin();

        // add some custom ghosting
        stk::mesh::Ghosting & ghosting = bulk.create_ghosting(std::string("new_nodes"));
        std::vector<stk::mesh::EntityKey> receive;
        ghosting.receive_list(receive);

        std::vector<stk::mesh::EntityProc> nodes_to_ghost;

        if(p_rank == 0)
        {
            stk::mesh::Entity node_5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
            EXPECT_TRUE(bulk.is_valid(node_5));
            nodes_to_ghost.push_back(stk::mesh::EntityProc(node_5, 3));
        }

        if(p_rank == 3)
        {
            stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, 9);
            EXPECT_TRUE(bulk.is_valid(node));
            nodes_to_ghost.push_back(stk::mesh::EntityProc(node, 0));
        }

        bulk.change_ghosting(ghosting, nodes_to_ghost, receive);

        bulk.modification_end();

        bulk.modification_begin();

        stk::mesh::Entity node9 = bulk.get_entity(stk::topology::NODE_RANK, 9);
        EXPECT_TRUE(bulk.is_valid(node9));
        if(p_rank == 0)
        {
            stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
            EXPECT_TRUE(bulk.is_valid(node1));
            stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, 1);
            EXPECT_TRUE(bulk.is_valid(elem));
            bulk.destroy_relation(elem, node1, 0);
            bulk.declare_relation(elem, node9, 0);
            bulk.add_node_sharing(node9, 3);
        }

        if (p_rank == 3)
        {
            bulk.add_node_sharing(node9, 0);
        }

        bool failed = false;
        try
        {
            bulk.modification_end();
        }
        catch(const std::exception & X)
        {
            std::cout<<X.what()<<std::endl;
            failed = true;
        }
        EXPECT_FALSE(failed);

    }

}

std::string printGhostData(stk::mesh::BulkData & bulkData, stk::mesh::Entity entity)
{
    std::ostringstream oss;
    std::vector<stk::mesh::EntityGhostData> egd;
    stk::mesh::impl::get_ghost_data(bulkData, entity, egd);
    for(size_t z = 0; z < egd.size(); ++z)
    {
        oss << "P" << bulkData.parallel_rank() << ":  " << egd[z] << std::endl;
    }
    return oss.str();
}

std::string printGhostDataByRank(stk::mesh::BulkData & bulkData, stk::topology::rank_t rank)
{
    const stk::mesh::BucketVector & buckets = bulkData.buckets(rank);
    std::ostringstream oss;
    oss << "P" << bulkData.parallel_rank() << ":  rank=" << rank << std::endl;
    for(size_t k = 0; k < buckets.size(); ++k)
    {
        const stk::mesh::Bucket::iterator begin = buckets[k]->begin();
        const stk::mesh::Bucket::iterator end = buckets[k]->end();
        for(stk::mesh::Bucket::iterator it = begin; it != end; ++it)
        {
            oss << printGhostData(bulkData, *it);
        }
    }
    return oss.str();
}

TEST(BulkData, EntityGhostData)
{
    std::string gold_result = "(Entity_lid=0, direction=SEND, processor=128, ghosting level=LOCALLY_OWNED)";
    stk::mesh::EntityGhostData data;
    data.direction = stk::mesh::EntityGhostData::SEND;
    data.ghostingLevel = stk::mesh::EntityGhostData::LOCALLY_OWNED;
    data.processor = 128;
    std::ostringstream oss;
    oss << data;
    EXPECT_EQ( gold_result, oss.str());
}

TEST(BulkData, get_ghost_data)
{
    using std::string;
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);

    if(psize == 3)
    { // Skip unless we're on 3 processors
        stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
        const string generatedMeshSpecification = "generated:1x1x3";
        stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

        if(stkMeshBulkData.parallel_rank() == 0)
        {
            std::ostringstream oss;
            for(stk::topology::rank_t rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; ++rank)
            {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            string gold_result =
                    string("P0:  rank=NODE_RANK\n")
                    + string("P0:  (Entity_gid=1, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=1, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=2, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=2, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=3, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=3, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=4, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=4, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=9, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=10, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=11, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=12, rank=0, direction=RECEIVE, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=5, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=5, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n")
                    + string("P0:  (Entity_gid=5, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=6, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=6, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n")
                    + string("P0:  (Entity_gid=6, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=7, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=7, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n")
                    + string("P0:  (Entity_gid=7, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=8, rank=0, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=8, rank=0, direction=SEND, processor=1, ghosting level=SHARED)\n")
                    + string("P0:  (Entity_gid=8, rank=0, direction=SEND, processor=2, ghosting level=AURA)\n")
                    + string("P0:  rank=EDGE_RANK\n")
                    + string("P0:  rank=FACE_RANK\n")
                    + string("P0:  rank=ELEMENT_RANK\n")
                    + string("P0:  (Entity_gid=1, rank=3, direction=NONE, processor=0, ghosting level=LOCALLY_OWNED)\n")
                    + string("P0:  (Entity_gid=1, rank=3, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P0:  (Entity_gid=2, rank=3, direction=RECEIVE, processor=1, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str());
        }
        else if(stkMeshBulkData.parallel_rank() == 1)
        {
            std::ostringstream oss;
            for(stk::topology::rank_t rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; ++rank)
            {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            std::string gold_result =
                    string("P1:  rank=NODE_RANK\n")
                    + string("P1:  (Entity_gid=5, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=6, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=7, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=8, rank=0, direction=RECEIVE, processor=0, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=1, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=2, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=3, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=4, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=13, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=14, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=15, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=16, rank=0, direction=RECEIVE, processor=2, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=9, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n")
                    + string("P1:  (Entity_gid=9, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=9, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=10, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n")
                    + string("P1:  (Entity_gid=10, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=10, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=11, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n")
                    + string("P1:  (Entity_gid=11, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=11, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=12, rank=0, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n")
                    + string("P1:  (Entity_gid=12, rank=0, direction=SEND, processor=2, ghosting level=SHARED)\n")
                    + string("P1:  (Entity_gid=12, rank=0, direction=SEND, processor=0, ghosting level=AURA)\n")
                    + string("P1:  rank=EDGE_RANK\n")
                    + string("P1:  rank=FACE_RANK\n")
                    + string("P1:  rank=ELEMENT_RANK\n")
                    + string("P1:  (Entity_gid=2, rank=3, direction=NONE, processor=1, ghosting level=LOCALLY_OWNED)\n")
                    + string("P1:  (Entity_gid=2, rank=3, direction=SEND, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=2, rank=3, direction=SEND, processor=2, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=1, rank=3, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P1:  (Entity_gid=3, rank=3, direction=RECEIVE, processor=2, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str());
        }
        else
        { // if (stkMeshBulkData.parallel_rank() == 2)
            std::ostringstream oss;
            for(stk::topology::rank_t rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; ++rank)
            {
                oss << printGhostDataByRank(stkMeshBulkData, rank);
            }
            std::string gold_result =
                    string("P2:  rank=NODE_RANK\n")
                    + string("P2:  (Entity_gid=13, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n")
                    + string("P2:  (Entity_gid=13, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=14, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n")
                    + string("P2:  (Entity_gid=14, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=15, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n")
                    + string("P2:  (Entity_gid=15, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=16, rank=0, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n")
                    + string("P2:  (Entity_gid=16, rank=0, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=9, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n")
                    + string("P2:  (Entity_gid=10, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n")
                    + string("P2:  (Entity_gid=11, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n")
                    + string("P2:  (Entity_gid=12, rank=0, direction=RECEIVE, processor=1, ghosting level=SHARED)\n")
                    + string("P2:  (Entity_gid=5, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=6, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=7, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=8, rank=0, direction=RECEIVE, processor=0, ghosting level=AURA)\n")
                    + string("P2:  rank=EDGE_RANK\n")
                    + string("P2:  rank=FACE_RANK\n")
                    + string("P2:  rank=ELEMENT_RANK\n")
                    + string("P2:  (Entity_gid=3, rank=3, direction=NONE, processor=2, ghosting level=LOCALLY_OWNED)\n")
                    + string("P2:  (Entity_gid=3, rank=3, direction=SEND, processor=1, ghosting level=AURA)\n")
                    + string("P2:  (Entity_gid=2, rank=3, direction=RECEIVE, processor=1, ghosting level=AURA)\n");
            EXPECT_EQ( gold_result, oss.str());
        }
    }
}

// Add to documentation tests for modification_end FAQ:
// 1.  How does change owner work?  Can I change the owner of a shared entity when I'm not the owner?
//     A:  Only the owner can give ownership to another processor
// 2.  Can I delete a shared and not-locally owned entity just on one processor?
// 3.  Can I delete a shared and locally owned entity just on one processor?
// 4.  Can I change parts on a shared entity differently on different sharing processors and will modification_end figure out the right parts?
//     A:  Only the owner can change parts on the entity.
// 5.  Can I add relations to shared entities differently on different sharing processors?
// 6.  Is there a difference in any of these between "creation" and "modification" cycles of the modification?
//     For example, can I create shared entities with different parts and it works?  vs.  Modifying shared entities to produce different parts?
TEST(DocTestBulkData, onlyTheOwnerCanChangeEntityParts)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs != 2)
    {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Part & myPart = stkMeshMetaData.declare_part("new_part");
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    stk::mesh::Selector allSharedSelector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> allSharedNodes;
    stk::mesh::get_selected_entities(allSharedSelector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), allSharedNodes);
    stkMeshBulkData.modification_begin();
    stk::mesh::PartVector addParts;
    addParts.push_back(&myPart);
    const int myRank = stk::parallel_machine_rank(communicator);
    for(size_t i = 0; i < allSharedNodes.size(); ++i)
    {
        if(stkMeshBulkData.parallel_owner_rank(allSharedNodes[i]) == myRank)
        {
            stkMeshBulkData.change_entity_parts(allSharedNodes[i], addParts);
        }
        else
        {
            EXPECT_THROW(stkMeshBulkData.change_entity_parts(allSharedNodes[i],addParts), std::logic_error);
        }
    }EXPECT_NO_THROW( stkMeshBulkData.modification_end());

    // Verify parts are correct on all processors.
    stk::mesh::Selector shared_selector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> shared_nodes;
    stk::mesh::get_selected_entities(shared_selector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), shared_nodes);
    for(size_t i = 0; i < shared_nodes.size(); ++i)
    {
        EXPECT_TRUE(stkMeshBulkData.bucket(shared_nodes[i]).member(myPart));
    }
}

TEST(BulkData, onlyKeepTheOwnersParts)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs != 2)
    {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA", stk::topology::NODE_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB", stk::topology::NODE_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::EntityKey node0Key;
    const stk::mesh::EntityId node_id = 2000;
    stkMeshBulkData.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;
    if(myRank == 0)
    {
        stk::mesh::Entity element0 = stkMeshBulkData.declare_element(element_id_0, elem_parts);
        stk::mesh::Entity node0 = stkMeshBulkData.declare_node(node_id, stk::mesh::ConstPartVector{&partA});
        int sharingProc = 1;
        stkMeshBulkData.add_node_sharing(node0, sharingProc);
        node0Key = stkMeshBulkData.entity_key(node0);
        for(stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
        {
            stkMeshBulkData.declare_relation(element0, node0, node_rel_id);
        }
    }
    else  // myRank == 1
    {
        stk::mesh::Entity element1 = stkMeshBulkData.declare_element(element_id_1, elem_parts);
        stk::mesh::Entity node0 = stkMeshBulkData.declare_node(node_id, stk::mesh::ConstPartVector{&partB});
        int sharingProc = 0;
        stkMeshBulkData.add_node_sharing(node0, sharingProc);
        node0Key = stkMeshBulkData.entity_key(node0);
        for(stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
        {
            stkMeshBulkData.declare_relation(element1, node0, node_rel_id);
        }
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Entity node0 = stkMeshBulkData.get_entity(node0Key);
    stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(node0);

    const int nodeOwner = 0;
    EXPECT_EQ( nodeOwner, stkMeshBulkData.parallel_owner_rank(node0));
    EXPECT_TRUE( nodeBucket.member(stkMeshMetaData.globally_shared_part()));

    EXPECT_TRUE( nodeBucket.member(partA));
    EXPECT_TRUE( nodeBucket.member(partB));
    EXPECT_TRUE( stk::mesh::has_superset(nodeBucket,partA));
    EXPECT_TRUE( stk::mesh::has_superset(nodeBucket,partB));

    stk::mesh::Entity element0 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_0);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_1);
    const int element0Owner = 0;
    const int element1Owner = 1;
    EXPECT_EQ( element0Owner, stkMeshBulkData.parallel_owner_rank(element0));
    EXPECT_EQ( element1Owner, stkMeshBulkData.parallel_owner_rank(element1));
}

TEST(BulkData, newSharedNodeGetMergedPartsFromElements)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs != 2)
    {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA", stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB", stk::topology::ELEMENT_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::EntityKey node0Key;
    const stk::mesh::EntityId node_id = 2000;
    stkMeshBulkData.modification_begin();
    if(myRank == 0)
    {
        const stk::mesh::EntityId element_id = 1000;
        stk::mesh::Entity element0 = stkMeshBulkData.declare_element(element_id, stk::mesh::ConstPartVector{&partA});
        stk::mesh::Entity node0 = stkMeshBulkData.declare_node(node_id);
        int sharingProc = 1;
        stkMeshBulkData.add_node_sharing(node0, sharingProc);
        node0Key = stkMeshBulkData.entity_key(node0);
        stkMeshBulkData.change_entity_parts(element0, elem_parts, empty_parts);
        for(stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
        {
            stkMeshBulkData.declare_relation(element0, node0, node_rel_id);
        }
    }
    else  // myRank == 1
    {
        const stk::mesh::EntityId element_id = 1001;
        stk::mesh::Entity element1 = stkMeshBulkData.declare_element(element_id, stk::mesh::ConstPartVector{&partB});
        stk::mesh::Entity node0 = stkMeshBulkData.declare_node(node_id);
        int sharingProc = 0;
        stkMeshBulkData.add_node_sharing(node0, sharingProc);
        node0Key = stkMeshBulkData.entity_key(node0);
        stkMeshBulkData.change_entity_parts(element1, elem_parts, empty_parts);
        for(stk::mesh::RelationIdentifier node_rel_id = 0; node_rel_id < elem_topo.num_nodes(); ++node_rel_id)
        {
            stkMeshBulkData.declare_relation(element1, node0, node_rel_id);
        }
    }
    stkMeshBulkData.modification_end();

    stk::mesh::Entity node0 = stkMeshBulkData.get_entity(node0Key);
    stk::mesh::Bucket & node0Bucket = stkMeshBulkData.bucket(node0);

    const int node0Owner = 0;
    EXPECT_EQ( node0Owner, stkMeshBulkData.parallel_owner_rank(node0));
    EXPECT_TRUE( node0Bucket.member(stkMeshMetaData.globally_shared_part()));

    EXPECT_TRUE( node0Bucket.member(partA));
    EXPECT_TRUE( node0Bucket.member(partB));
    EXPECT_TRUE( stk::mesh::has_superset(node0Bucket,partA));
    EXPECT_TRUE( stk::mesh::has_superset(node0Bucket,partB));
}

TEST(BulkData, mayCreateRelationsToNodesDifferently)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs != 2)
    {
        return;
    }

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    const std::string generatedMeshSpecification = "generated:1x1x2";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    const int myRank = stk::parallel_machine_rank(communicator);
    stk::mesh::Part & partA = stkMeshMetaData.declare_part("PartA", stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = stkMeshMetaData.declare_part("PartB", stk::topology::ELEMENT_RANK);
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (stkMeshMetaData.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = stkMeshMetaData.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::Selector allSharedSelector = stkMeshMetaData.globally_shared_part();
    std::vector<stk::mesh::Entity> allSharedNodes;
    stk::mesh::get_selected_entities(allSharedSelector, stkMeshBulkData.buckets(stk::topology::NODE_RANK), allSharedNodes);

    std::sort(allSharedNodes.begin(), allSharedNodes.end(), stk::mesh::EntityLess(stkMeshBulkData));
    stk::mesh::Entity sharedNode0 = allSharedNodes[0];
    stk::mesh::Entity sharedNode1 = allSharedNodes[1];

    stkMeshBulkData.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;
    stk::mesh::Entity filler_node = stkMeshBulkData.declare_node(123456);
    if(myRank == 0)
    {
        int sharingProc = 1;
        stkMeshBulkData.add_node_sharing(filler_node, sharingProc);
        stk::mesh::Entity element0 = stkMeshBulkData.declare_element(element_id_0, stk::mesh::ConstPartVector{&partA});
        stkMeshBulkData.change_entity_parts(element0, elem_parts, empty_parts);
        stk::mesh::RelationIdentifier node_rel_id = 0;
        stkMeshBulkData.declare_relation(element0, sharedNode0, node_rel_id);
        for(stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
        {
            stkMeshBulkData.declare_relation(element0, filler_node, filler_rel_id);
        }
    }
    else  // myRank == 1
    {
        int sharingProc = 0;
        stkMeshBulkData.add_node_sharing(filler_node, sharingProc);
        stk::mesh::Entity element1 = stkMeshBulkData.declare_element(element_id_1, stk::mesh::ConstPartVector{&partB});
        stkMeshBulkData.change_entity_parts(element1, elem_parts, empty_parts);
        const stk::mesh::RelationIdentifier node_rel_id = 0;
        stkMeshBulkData.declare_relation(element1, sharedNode0, node_rel_id);
        for(stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
        {
            stkMeshBulkData.declare_relation(element1, filler_node, filler_rel_id);
        }
    }EXPECT_NO_THROW( stkMeshBulkData.modification_end());
    {
        stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(sharedNode0);
        EXPECT_TRUE( nodeBucket.member(partA));
        EXPECT_TRUE( nodeBucket.member(partB));
        EXPECT_EQ( 4u, stkMeshBulkData.num_elements(sharedNode0));
    }
    stk::mesh::Entity element0 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_0);
    stk::mesh::Entity element1 = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, element_id_1);
    stkMeshBulkData.modification_begin();
    if(myRank == 0)
    {
        stk::mesh::RelationIdentifier node_rel_id = 1;
        stkMeshBulkData.declare_relation(element0, sharedNode1, node_rel_id);
    }
    else // myRank == 1
    {
        stk::mesh::RelationIdentifier node_rel_id = 1;
        stkMeshBulkData.declare_relation(element1, sharedNode1, node_rel_id);

    }EXPECT_NO_THROW( stkMeshBulkData.modification_end());

    {
        stk::mesh::Bucket & nodeBucket = stkMeshBulkData.bucket(sharedNode1);
        EXPECT_TRUE( nodeBucket.member(partA));
        EXPECT_TRUE( nodeBucket.member(partB));
        EXPECT_EQ( 4u, stkMeshBulkData.num_elements(sharedNode1));
    }
}

stk::mesh::PartVector setupFixture(stk::io::StkMeshIoBroker & io)
{
    const std::string generatedMeshSpecification = "generated:1x1x2";
    io.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    io.create_input_mesh();
    stk::mesh::Part & partA = io.meta_data().declare_part("PartA", stk::topology::ELEMENT_RANK);
    stk::mesh::Part & partB = io.meta_data().declare_part("PartB", stk::topology::ELEMENT_RANK);
    io.populate_bulk_data();
    stk::mesh::PartVector pv;
    pv.push_back(&partA);
    pv.push_back(&partB);
    return pv;
}

stk::mesh::EntityVector getSortedNodes(stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), nodes);
    std::sort(nodes.begin(), nodes.end(), stk::mesh::EntityLess(bulk));
    return nodes;
}

bool no_aura_declare_relation(stk::mesh::BulkData &bulk,
                              stk::mesh::Entity e_from,
                              stk::mesh::Entity e_to,
                              const stk::mesh::RelationIdentifier local_id)
{
    stk::mesh::Bucket const *b_from = bulk.bucket_ptr(e_from);
    stk::mesh::Bucket const *b_to = bulk.bucket_ptr(e_to);
    const bool ownership_ok = (b_from && b_to && !b_from->in_aura() && b_to->in_aura());

    if(ownership_ok)
    {
        bulk.declare_relation(e_from, e_to, local_id);
        return true;
    }
    else
    {
        return false;
    }
}

TEST(DocTestBulkData, inducedPartMembershipIgnoredForNonOwnedHigherRankedEntities)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs != 2)
    {
        return;
    }
    const int myRank = stk::parallel_machine_rank(communicator);

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    stk::mesh::PartVector pv = setupFixture(stkMeshIoBroker);
    stk::mesh::Part & partA = *pv[0];
    stk::mesh::Part & partB = *pv[1];

    stk::mesh::MetaData &meta = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &bulk = stkMeshIoBroker.bulk_data();

    stk::topology elem_topo = (meta.spatial_dimension() == 2 ? stk::topology::QUAD_4_2D : stk::topology::HEX_8);
    stk::mesh::Part &elem_part = meta.get_topology_root_part(elem_topo);
    stk::mesh::PartVector elem_parts, empty_parts;
    elem_parts.push_back(&elem_part);

    stk::mesh::EntityVector ev = getSortedNodes(bulk, stk::mesh::Selector(meta.globally_shared_part()));
    ASSERT_TRUE( ev.size() >= 2);
    stk::mesh::Entity sharedNodeA(ev[0]);
    stk::mesh::Entity sharedNodeB(ev[1]);

    stk::mesh::Bucket & nodeABucket = bulk.bucket(sharedNodeA);
    EXPECT_FALSE( nodeABucket.member(partA));
    EXPECT_FALSE( nodeABucket.member(partB));
    EXPECT_EQ( 2u, bulk.num_elements(sharedNodeA));
    EXPECT_EQ( 2u, bulk.num_elements(sharedNodeB));

    // First we create off processor elements
    bulk.modification_begin();
    const stk::mesh::EntityId element_id_0 = 1000;
    const stk::mesh::EntityId element_id_1 = 1001;

    const stk::mesh::EntityId filler_node_id = 123456;
    stk::mesh::Entity filler_node = bulk.declare_node(filler_node_id);

    bulk.add_node_sharing(filler_node, (myRank == 1 ? 0 : 1));

    if(myRank == 0)
    {
        const stk::mesh::RelationIdentifier node_rel_id = 0;
        stk::mesh::Entity element = bulk.declare_element(element_id_0, stk::mesh::ConstPartVector{&partA});
        bulk.change_entity_parts(element, elem_parts, empty_parts);
        bulk.declare_relation(element, sharedNodeA, node_rel_id);
        for(stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
        {
            bulk.declare_relation(element, filler_node, filler_rel_id);
        }
    }
    else  // myRank == 1
    {
        const stk::mesh::RelationIdentifier node_rel_id = 0;
        stk::mesh::Entity element = bulk.declare_element(element_id_1, stk::mesh::ConstPartVector{&partB});
        bulk.change_entity_parts(element, elem_parts, empty_parts);
        bulk.declare_relation(element, sharedNodeA, node_rel_id);
        for(stk::mesh::RelationIdentifier filler_rel_id = node_rel_id + 1; filler_rel_id < elem_topo.num_nodes(); ++filler_rel_id)
        {
            bulk.declare_relation(element, filler_node, filler_rel_id);
        }
    }

    EXPECT_NO_THROW( bulk.modification_end());
    {
        stk::mesh::Bucket & nodeABkt = bulk.bucket(sharedNodeA);
        EXPECT_TRUE( nodeABkt.member(partA));
        EXPECT_TRUE( nodeABkt.member(partB));
        EXPECT_EQ( 4u, bulk.num_elements(sharedNodeA));
    }

    stk::mesh::Entity element0 = bulk.get_entity(stk::topology::ELEMENT_RANK, element_id_0);
    stk::mesh::Entity element1 = bulk.get_entity(stk::topology::ELEMENT_RANK, element_id_1);

    stk::mesh::Part &locally_owned_part = meta.locally_owned_part();
    stk::mesh::Bucket &element0_bucket = bulk.bucket(element0);
    stk::mesh::Bucket &element1_bucket = bulk.bucket(element1);

    if(myRank == 0)
    {
        EXPECT_TRUE(element0_bucket.member(locally_owned_part));
        EXPECT_FALSE(element1_bucket.member(locally_owned_part));
    }
    else // myRank == 1
    {
        EXPECT_FALSE(element0_bucket.member(locally_owned_part));
        EXPECT_TRUE(element1_bucket.member(locally_owned_part));
    }

    // Disallowing declaring relations to/from an entity in the aura can make it
    // easier to keep the mesh consistent.
    bulk.modification_begin();
    if(myRank == 0)
    {
        const stk::mesh::RelationIdentifier node_rel_id = 1;
        EXPECT_FALSE(no_aura_declare_relation(bulk, element1, sharedNodeB, node_rel_id));
    }
    else // myRank == 1
    {
        const stk::mesh::RelationIdentifier node_rel_id = 1;
        EXPECT_FALSE(no_aura_declare_relation(bulk, element0, sharedNodeB, node_rel_id));
    }EXPECT_NO_THROW( bulk.modification_end());

    {
        stk::mesh::Bucket & nodeBBucket = bulk.bucket(sharedNodeB);
        EXPECT_FALSE( nodeBBucket.member(partA));
        EXPECT_FALSE( nodeBBucket.member(partB));
        EXPECT_EQ( 2u, bulk.num_elements(sharedNodeB));

        const int element0Owner = 0;
        EXPECT_EQ( element0Owner, bulk.parallel_owner_rank(element0));
        const int element1Owner = 1;
        EXPECT_EQ( element1Owner, bulk.parallel_owner_rank(element1));

        if(myRank == 0)
        {
            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);

            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
            EXPECT_TRUE( filler_node == bulk.begin_nodes(element1)[1]);
        }
        else // myRank == 1
        {
            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);
            EXPECT_TRUE( filler_node == bulk.begin_nodes(element0)[1]);

            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
        }
    }

    // Again we try to create relations from ghost element to shared node B.
    // BulkData::declare_relation(..) will do this, but the owner relations
    // will reset the ghost's relations in modification_end().
    //
    // Thus, more work needs to be done to enforce the relationship reciprocity
    // aspect of mesh consistency.
    bulk.modification_begin();
    if(myRank == 0)
    {
        const stk::mesh::RelationIdentifier node_rel_id = 1;
        bulk.declare_relation(element1, sharedNodeB, node_rel_id);
        EXPECT_TRUE( sharedNodeB == bulk.begin_nodes(element1)[1]);
    }
    else // myRank == 1
    {
        const stk::mesh::RelationIdentifier node_rel_id = 1;
        bulk.declare_relation(element0, sharedNodeB, node_rel_id);
        EXPECT_TRUE( sharedNodeB == bulk.begin_nodes(element0)[1]);
    }EXPECT_NO_THROW( bulk.modification_end());

    {
        stk::mesh::Bucket & nodeBBucket = bulk.bucket(sharedNodeB);
        EXPECT_FALSE( nodeBBucket.member(partA));
        EXPECT_FALSE( nodeBBucket.member(partB));
        EXPECT_EQ( 3u, bulk.num_elements(sharedNodeB));

        if(myRank == 0)
        {
            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);

            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
            EXPECT_TRUE( filler_node == bulk.begin_nodes(element1)[1]);
        }
        else // myRank == 1
        {
            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element0)[0]);
            EXPECT_TRUE( filler_node == bulk.begin_nodes(element0)[1]);

            EXPECT_TRUE( sharedNodeA == bulk.begin_nodes(element1)[0]);
        }
    }
}

TEST(BulkData, ModificationEnd)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);

    if(numProcs == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        stk::unit_test_util::BulkDataTester *stkMeshBulkData = new stk::unit_test_util::BulkDataTester(stkMeshMetaData, communicator);

        std::string exodusFileName = stk::unit_test_util::get_option("-i", "generated:1x1x4");

        // STK IO module will be described in separate chapter.
        // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
        // The order of the following lines in {} are important
        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);

            // Inform STK IO which STK Mesh objects to populate later
            exodusFileReader.set_bulk_data(*stkMeshBulkData);

            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

            // Populate the MetaData which has the descriptions of the Parts and Fields.
            exodusFileReader.create_input_mesh();

            // Populate entities in STK Mesh from Exodus file
            exodusFileReader.populate_bulk_data();
        }

        int elementToMove = 3;
        int nodeToCheck = 9;

        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK, nodeToCheck);
        stk::mesh::EntityKey entityToMoveKey(stk::topology::ELEMENT_RANK, elementToMove);

        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData->my_internal_comm_list().begin(),
                                                                                    stkMeshBulkData->my_internal_comm_list().end(),
                                                                                    nodeEntityKey);

        ASSERT_TRUE(iter != stkMeshBulkData->my_internal_comm_list().end());
        EXPECT_EQ(nodeEntityKey, iter->key);
        EXPECT_TRUE(stkMeshBulkData->is_valid(iter->entity));

        stkMeshBulkData->modification_begin();

        ASSERT_TRUE( stkMeshBulkData->is_valid(stkMeshBulkData->get_entity(entityToMoveKey)));

        if(stkMeshBulkData->parallel_rank() == 1)
        {
            stkMeshBulkData->destroy_entity(stkMeshBulkData->get_entity(entityToMoveKey));
        }

        // Really testing destroy_entity
        stkMeshBulkData->my_delete_shared_entities_which_are_no_longer_in_owned_closure();

        iter = std::lower_bound(stkMeshBulkData->my_internal_comm_list().begin(), stkMeshBulkData->my_internal_comm_list().end(), nodeEntityKey);

        ASSERT_TRUE(iter != stkMeshBulkData->my_internal_comm_list().end());
        EXPECT_EQ(nodeEntityKey, iter->key);

        if(stkMeshBulkData->parallel_rank() == 0)
        {
            EXPECT_TRUE(stkMeshBulkData->is_valid(iter->entity));
        }
        else
        {
            EXPECT_FALSE(stkMeshBulkData->is_valid(iter->entity));
        }

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(*stkMeshBulkData, globalCounts);

        delete stkMeshBulkData;
    }
}

TEST(BulkData, resolve_ownership_of_modified_entities_trivial)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    const int myRank = stk::parallel_machine_rank(communicator);

    if(numProcs != 3)
    {
        return;
    }

    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester mesh(stkMeshMetaData, communicator);
    std::string exodusFileName = stk::unit_test_util::get_option("-i", "generated:1x1x3");
    {
        stk::io::StkMeshIoBroker exodusFileReader(communicator);
        exodusFileReader.set_bulk_data(mesh);
        exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
        exodusFileReader.create_input_mesh();
        exodusFileReader.populate_bulk_data();
    }
    std::vector<Entity> modified_entities;
    if (myRank == 0) {
        modified_entities.push_back(mesh.get_entity(stk::topology::NODE_RANK, 1));
    }
    else if (myRank == 1) {
        modified_entities.push_back(mesh.get_entity(stk::topology::NODE_RANK, 9));
    }
    else {
        modified_entities.push_back(mesh.get_entity(stk::topology::NODE_RANK, 13));
    }
    mesh.my_resolve_ownership_of_modified_entities(modified_entities);
    if (myRank == 0) {
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 1), CEOUtils::STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 9), CEOUtils::STATE_OWNED, 1));
    }
    else if (myRank == 1) {
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 1), CEOUtils::STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 9), CEOUtils::STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 13), CEOUtils::STATE_OWNED, 2));
    }
    else {
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 9), CEOUtils::STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(stk::topology::NODE_RANK, 13), CEOUtils::STATE_OWNED, 2));
    }
}

TEST(BulkData, verify_closure_count_is_correct)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    const int myRank = stk::parallel_machine_rank(communicator);

    if(numProcs == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        stk::unit_test_util::BulkDataTester *stkMeshBulkData = new stk::unit_test_util::BulkDataTester(stkMeshMetaData, communicator);

        std::string exodusFileName = stk::unit_test_util::get_option("-i", "generated:1x1x2");

        // STK IO module will be described in separate chapter.
        // It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
        // The order of the following lines in {} are important
        {
            stk::io::StkMeshIoBroker exodusFileReader(communicator);

            // Inform STK IO which STK Mesh objects to populate later
            exodusFileReader.set_bulk_data(*stkMeshBulkData);

            exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

            // Populate the MetaData which has the descriptions of the Parts and Fields.
            exodusFileReader.create_input_mesh();

            // Populate entities in STK Mesh from Exodus file
            exodusFileReader.populate_bulk_data();
        }

        stk::mesh::EntityKey element_1_key(stk::topology::ELEMENT_RANK, 1);
        stk::mesh::EntityKey element_2_key(stk::topology::ELEMENT_RANK, 2);

        stk::mesh::EntityKey node_1_key(stk::topology::NODE_RANK, 1);
        stk::mesh::EntityKey node_2_key(stk::topology::NODE_RANK, 2);
        stk::mesh::EntityKey node_3_key(stk::topology::NODE_RANK, 3);
        stk::mesh::EntityKey node_4_key(stk::topology::NODE_RANK, 4);
        stk::mesh::EntityKey node_5_key(stk::topology::NODE_RANK, 5);
        stk::mesh::EntityKey node_6_key(stk::topology::NODE_RANK, 6);
        stk::mesh::EntityKey node_7_key(stk::topology::NODE_RANK, 7);
        stk::mesh::EntityKey node_8_key(stk::topology::NODE_RANK, 8);
        stk::mesh::EntityKey node_9_key(stk::topology::NODE_RANK, 9);
        stk::mesh::EntityKey node_10_key(stk::topology::NODE_RANK, 10);
        stk::mesh::EntityKey node_11_key(stk::topology::NODE_RANK, 11);
        stk::mesh::EntityKey node_12_key(stk::topology::NODE_RANK, 12);

        if(myRank == 0)
        {
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)));

            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)));
        }
        else // myRank == 1
        {
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)));

            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)));
        }

        std::vector<EntityProc> change_owner_vector;
        if(myRank == 1)
        {
            const int target_processor = 0;
            stk::mesh::Selector my_entities_selector = stkMeshBulkData->mesh_meta_data().locally_owned_part();
            std::vector<Entity> nodes;
            stk::mesh::get_selected_entities(my_entities_selector, stkMeshBulkData->buckets(stk::topology::NODE_RANK), nodes);
            for(size_t i = 0; i < nodes.size(); ++i)
            {
                change_owner_vector.push_back(EntityProc(nodes[i], target_processor));
            }
        }
        stkMeshBulkData->change_entity_owner(change_owner_vector);

        if(myRank == 0)
        {
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)));

            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)));
            EXPECT_EQ( 2u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)));
        }
        else // myRank == 1
        {
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_1_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(element_2_key)));

            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_1_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_2_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_3_key)));
            EXPECT_EQ( 0u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_4_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_5_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_6_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_7_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_8_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_9_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_10_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_11_key)));
            EXPECT_EQ( 1u, stkMeshBulkData->closure_count(stkMeshBulkData->get_entity(node_12_key)));
        }

        delete stkMeshBulkData;
    }

}

TEST(BulkData, orphaned_node_closure_count_shared_nodes_non_owner_adds_element)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  const int myRank = stk::parallel_machine_rank(communicator);

  if (numProcs != 2) { return; }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta,communicator);

  stk::mesh::Part& element_part = meta.declare_part_with_topology("Beam2Part", stk::topology::BEAM_2);

  bulk.modification_begin();
  stk::mesh::Entity node1 = bulk.declare_node(1);
  stk::mesh::Entity node2 = bulk.declare_node(2);
  int other_proc = 1-myRank;
  bulk.add_node_sharing(node1,other_proc);
  bulk.add_node_sharing(node2,other_proc);
  EXPECT_EQ(bulk.my_orphaned_node_marking()+1,bulk.closure_count(node1));
  EXPECT_EQ(bulk.my_orphaned_node_marking()+1,bulk.closure_count(node2));

  ASSERT_NO_THROW(bulk.modification_end());

  if (myRank == 0)
  {
    EXPECT_EQ(1u,bulk.closure_count(node1));
    EXPECT_EQ(1u,bulk.closure_count(node2));
  }
  else
  {
    EXPECT_EQ(bulk.my_orphaned_node_marking()+0,bulk.closure_count(node1));
    EXPECT_EQ(bulk.my_orphaned_node_marking()+0,bulk.closure_count(node2));
  }

  bulk.modification_begin();
  if (myRank == 1)
  {
    stk::mesh::Entity element = bulk.declare_element(1, stk::mesh::ConstPartVector{&element_part});
    bulk.declare_relation(element,node1,0);
    EXPECT_EQ(1u,bulk.closure_count(node1));
    bulk.declare_relation(element,node2,1);
    EXPECT_EQ(1u,bulk.closure_count(node2));
  }
  EXPECT_NO_THROW(bulk.modification_end());

  EXPECT_EQ(1u,bulk.closure_count(node1));
  EXPECT_EQ(1u,bulk.closure_count(node2));

  bulk.modification_begin();
  if (myRank == 0)
  {
    std::vector<Entity> element_vector;
    bulk.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), element_vector);
    ASSERT_EQ(1u, element_vector.size());
    EXPECT_TRUE(bulk.destroy_entity(element_vector[0]));
    EXPECT_TRUE(bulk.destroy_entity(node1));
    EXPECT_TRUE(bulk.destroy_entity(node2));
    EXPECT_FALSE(bulk.is_valid(node1));
    EXPECT_FALSE(bulk.is_valid(node2));
  }
  EXPECT_NO_THROW(bulk.modification_end());
  if (myRank == 0)
  {
    EXPECT_FALSE(bulk.is_valid(node1));
    EXPECT_FALSE(bulk.is_valid(node2));
  }
  else // myRank == 1
  {
    EXPECT_TRUE(bulk.bucket(node1).owned());
    EXPECT_TRUE(bulk.bucket(node2).owned());
    EXPECT_FALSE(bulk.bucket(node1).shared());
    EXPECT_FALSE(bulk.bucket(node2).shared());
  }
}

TEST(BulkData, orphaned_node_closure_count_shared_nodes_owner_deletes)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  const int myRank = stk::parallel_machine_rank(communicator);

  if (numProcs != 2) { return; }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta,communicator);

  bulk.modification_begin();
  stk::mesh::Entity node1 = bulk.declare_node(1);
  int other_proc = 1-myRank;
  bulk.add_node_sharing(node1,other_proc);

  bulk.modification_end();

  bulk.modification_begin();

  if (myRank == 0)
  {
    EXPECT_TRUE(bulk.destroy_entity(node1));
    EXPECT_FALSE(bulk.is_valid(node1));
  }
  EXPECT_NO_THROW(bulk.modification_end());
  if (myRank == 0)
  {
    EXPECT_FALSE(bulk.is_valid(node1));
  }
  else // myRank == 1
  {
    EXPECT_EQ(1u,bulk.closure_count(node1));
    EXPECT_TRUE(bulk.bucket(node1).owned());
    EXPECT_FALSE(bulk.bucket(node1).shared());
  }
}

TEST(BulkData, orphaned_node_closure_count_shared_nodes_change_entity_owner_3proc)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  const int myRank = stk::parallel_machine_rank(communicator);

  if (numProcs != 3) { return; }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta,communicator);

  bulk.modification_begin();
  stk::mesh::Entity node1;
  if (myRank != 2) {
    node1 = bulk.declare_node(1);
    int other_proc = 1-myRank;
    bulk.add_node_sharing(node1,other_proc);
  }
  bulk.modification_end();
  if (myRank == 0)
  {
    EXPECT_EQ(1u, bulk.closure_count(node1));
  }
  else if (myRank == 1)
  {
    EXPECT_EQ(bulk.my_orphaned_node_marking(), bulk.closure_count(node1));
  }
  else // myRank == 2
  {
    EXPECT_FALSE(bulk.is_valid(node1));
  }

  std::vector<EntityProc> new_owners;
  if (myRank == 0)
  {
    new_owners.push_back(stk::mesh::EntityProc(node1,2));
  }
  EXPECT_NO_THROW(bulk.change_entity_owner(new_owners));

  if (myRank == 0)
  {
    EXPECT_FALSE(bulk.is_valid(node1));
  }
  else if (myRank == 1)
  {
    node1 = bulk.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK,1));
    EXPECT_EQ(bulk.my_orphaned_node_marking()+0,bulk.closure_count(node1));
    EXPECT_FALSE(bulk.bucket(node1).owned());
    EXPECT_TRUE(bulk.bucket(node1).shared());
  }
  else // myRank == 2
  {
    node1 = bulk.get_entity(stk::mesh::EntityKey(stk::topology::NODE_RANK,1));
    EXPECT_EQ(1u,bulk.closure_count(node1));
    EXPECT_TRUE(bulk.bucket(node1).owned());
    EXPECT_TRUE(bulk.bucket(node1).shared());
  }
}

TEST(BulkData, orphaned_node_closure_count_shared_nodes_change_entity_owner_2proc)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  const int myRank = stk::parallel_machine_rank(communicator);

  if (numProcs != 2) { return; }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta,communicator);

  bulk.modification_begin();
  stk::mesh::Entity node1 = bulk.declare_node(1);
  int other_proc = 1-myRank;
  bulk.add_node_sharing(node1,other_proc);
  bulk.modification_end();
  if (myRank == 0)
  {
    EXPECT_EQ(1u,bulk.closure_count(node1));
  }
  else // myRank == 1
  {
    EXPECT_EQ(bulk.my_orphaned_node_marking(),bulk.closure_count(node1));
  }

  std::vector<EntityProc> new_owners;
  if (myRank == 0)
  {
    new_owners.push_back(stk::mesh::EntityProc(node1,1));
  }
  EXPECT_NO_THROW(bulk.change_entity_owner(new_owners));

  if (myRank == 0)
  {
    EXPECT_EQ(0, bulk.closure_count(node1));
    EXPECT_FALSE(bulk.is_valid(node1));
  }
  else if (myRank == 1)
  {
    EXPECT_EQ(1u,bulk.closure_count(node1));
    EXPECT_TRUE(bulk.bucket(node1).owned());
    EXPECT_FALSE(bulk.bucket(node1).shared());
  }
}

TEST(BulkData, orphaned_node_closure_count_shared_nodes_owner_adds_element)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  const int myRank = stk::parallel_machine_rank(communicator);

  if (numProcs != 2) { return; }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::unit_test_util::BulkDataTester bulk(meta,communicator);

  stk::mesh::Part& element_part = meta.declare_part_with_topology("Beam2Part", stk::topology::BEAM_2);

  bulk.modification_begin();
  stk::mesh::Entity node1 = bulk.declare_node(1);
  stk::mesh::Entity node2 = bulk.declare_node(2);
  int other_proc = 1-myRank;
  bulk.add_node_sharing(node1,other_proc);
  bulk.add_node_sharing(node2,other_proc);

  bulk.modification_end();

  bulk.modification_begin();
  if (myRank == 0)
  {
    stk::mesh::Entity element = bulk.declare_element(1, stk::mesh::ConstPartVector{&element_part});
    bulk.declare_relation(element,node1,0);
    EXPECT_EQ(2u,bulk.closure_count(node1));
    bulk.declare_relation(element,node2,1);
    EXPECT_EQ(2u,bulk.closure_count(node2));
  }
  EXPECT_NO_THROW(bulk.modification_end());
  if (myRank == 1)
  {
    EXPECT_EQ(bulk.my_orphaned_node_marking(),bulk.closure_count(node1));
    EXPECT_EQ(bulk.my_orphaned_node_marking(),bulk.closure_count(node2));
  }
}

TEST(BulkData, change_entity_owner_no_aura_check)
{
//   id/owner_proc
//
//   1/0---4/0---5/0      1/0---4/1---5/1
//    |     |     |        |     |     |
//    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
//    |     |     |        |     |     |
//   2/0---3/0---6/0      2/0---3/0---6/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::unit_test_util::BulkDataTester bulk( meta, pm, stk::mesh::BulkData::NO_AUTO_AURA);

  std::vector<stk::mesh::Entity> elems;
  CEOUtils::fillMeshfor2Elem2ProcMoveAndTest(bulk, meta, elems);

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  CEOUtils::checkStatesAfterCEOME_2Elem2ProcMove_no_ghost(bulk);
}

TEST(BulkData, modification_end_and_change_entity_owner_no_aura_check)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::unit_test_util::BulkDataTester mesh( meta, pm, stk::mesh::BulkData::NO_AUTO_AURA);

  CEOUtils::fillMeshfor2Elem2ProcFlipAndTest_no_ghost(mesh, meta);

  stk::mesh::EntityProcVec entity_procs_flip;
  if (p_rank == 0) {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 2), 1));
  } else {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 2), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 5), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 6), 0));
  }
  mesh.change_entity_owner(entity_procs_flip);

  CEOUtils::checkStatesAfterCEOME_2Elem2ProcFlip_no_ghost(mesh);
}

//==============================================================================
//These 8 tests thoroughly test change_entity_owner() and consequently:
//in_receive_ghost()
//is_valid()
//parallel_owner_rank() : assuming entity you pass is valid
//entity_comm_map_shared(): this is called for every possible state
//entity_comm_map_owner(): haven't explicitly tested for invalid processor condition
//aura_ghosting:  trivial - returns 1 - no obvious need to unit test
//declare_element: tested ubiquitously throughout the code - including stk mesh unit tests
//add_node_sharing: tested ubiquitously throughout the code - including specific add node sharing test in UnitTestBulkDataSharing
//connectEntityToEdge: exercised through create_edges() and directly within 2 change_entity_owner acceptance tests - not certain what happens if you pass it garbage
//declare_part_with_topology: exercised ubiquitously and directly throughout the unit tests


TEST(BulkData, change_entity_owner_2Elem2ProcMove)
{
//   id/owner_proc
//
//   1/0---4/0---5/0      1/0---4/1---5/1
//    |     |     |        |     |     |
//    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
//    |     |     |        |     |     |
//   2/0---3/0---6/0      2/0---3/0---6/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::unit_test_util::BulkDataTester bulk( meta, pm);

  std::vector<stk::mesh::Entity> elems;
  CEOUtils::fillMeshfor2Elem2ProcMoveAndTest(bulk, meta, elems);

  stk::mesh::EntityProcVec entity_procs;
  if (p_rank == 0) {
    entity_procs.push_back(stk::mesh::EntityProc(elems[1], 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 4), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 5), 1));
    entity_procs.push_back(stk::mesh::EntityProc(bulk.get_entity(stk::topology::NODE_RANK, 6), 1));
  }
  bulk.change_entity_owner(entity_procs);

  CEOUtils::checkStatesAfterCEOME_2Elem2ProcMove(bulk);
}

TEST(BulkData, change_entity_owner_2Elem2ProcFlip)
{
  //   id/owner_proc
  //
  //   1/0---4/0---5/1        1/1---4/0---5/0
  //    |     |     |          |     |     |
  //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
  //    |     |     |          |     |     |
  //   2/0---3/0---6/1        2/1---3/0---6/0

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_rank = stk::parallel_machine_rank( pm );
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 2) {
    return;
  }

  const int spatial_dimension = 2;
  stk::mesh::MetaData meta( spatial_dimension );
  stk::unit_test_util::BulkDataTester mesh( meta, pm);

  CEOUtils::fillMeshfor2Elem2ProcFlipAndTest(mesh, meta);

  stk::mesh::EntityProcVec entity_procs_flip;
  if (p_rank == 0) {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 1), 1));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 2), 1));
  } else {
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::ELEM_RANK, 2), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 5), 0));
    entity_procs_flip.push_back(stk::mesh::EntityProc(mesh.get_entity(stk::topology::NODE_RANK, 6), 0));
  }
  mesh.change_entity_owner(entity_procs_flip);

  CEOUtils::checkStatesAfterCEOME_2Elem2ProcFlip(mesh);
}

TEST(BulkData, change_entity_owner_3Elem2ProcMoveRight)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  stk::unit_test_util::BulkDataTester mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 2)
  {
    return;
  }

  EntityVector nodes;
  EntityVector elements;

  CEOUtils::fillMeshfor3Elem2ProcMoveRightAndTest(mesh, meta_data, nodes, elements);

  std::vector<EntityProc> change;
  if(p_rank == 0)
  {
    change.push_back(EntityProc(elements[1], 1));
    change.push_back(EntityProc(nodes[4], 1));
    change.push_back(EntityProc(nodes[5], 1));
  }

  mesh.change_entity_owner(change);

  CEOUtils::checkStatesAfterCEOME_3Elem2ProcMoveRight(mesh);
}

TEST(BulkData, change_entity_owner_3Elem2ProcMoveLeft)
{
  //   id/owner_proc
  //
  //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
  //    |     |     |     |           |     |     |     |
  //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
  //    |     |     |     |           |     |     |     |
  //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  stk::unit_test_util::BulkDataTester mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 2)
  {
    return;
  }

  EntityVector nodes;
  EntityVector elements;

  CEOUtils::fillMeshfor3Elem2ProcMoveLeftAndTest(mesh, meta_data, nodes, elements);

  std::vector<EntityProc> change;
  if(p_rank == 1)
  {
    change.push_back(EntityProc(elements[0], 0));
    change.push_back(EntityProc(nodes[2], 0));
    change.push_back(EntityProc(nodes[3], 0));
  }

  mesh.change_entity_owner(change);

  CEOUtils::checkStatesAfterCEOME_3Elem2ProcMoveLeft(mesh);
}

TEST(BulkData, change_entity_owner_4Elem4ProcEdge)
{
  // This unit-test is designed to test the conditions that results that
  // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
  // it will test the changing-of-ownership of a shared edge to a proc that
  // either ghosted it or did not know about it.
  //
  //         id/proc                             id/proc
  //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
  //        |      |     |    ||     |          |      |     |    ||     |
  //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
  //        |      |     |    ||     |          |      |     |    ||     |
  //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
  //  this edge moves to p0 --^
  //  element 3 moves to proc 0.
  //  nodes 7&8 move to proc 0.
  //  proc 2 forgets everything.
  //
  // To test this, we use the mesh above, with each elem going on a separate
  // proc, one elem per proc. We will take the edge shared by the last
  // two (rightmost) elements and change the ownership to proc 0.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  stk::unit_test_util::BulkDataTester mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if(p_size != 4)
  {
    return;
  }

  stk::mesh::EntityKey elem_key_chg_own;
  CEOUtils::fillMeshfor4Elem4ProcEdgeAndTest(mesh, meta_data, elem_key_chg_own);

  std::vector<EntityProc> change;
  if(p_rank == 2)
  {
    // Change ownership of changing elem and all entities in it's closure that
    // we own to proc 0.

    Entity changing_elem = mesh.get_entity(elem_key_chg_own);
    ASSERT_TRUE( mesh.is_valid(changing_elem));
    EntityProc eproc(changing_elem, 0 /*new owner*/);
    change.push_back(eproc);

    const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
    for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
    {
      stk::mesh::Entity const *to_i = mesh.begin(changing_elem, irank);
      stk::mesh::Entity const *to_e = mesh.end(changing_elem, irank);
      for (; to_i != to_e; ++to_i)
      {
        if (mesh.parallel_owner_rank(*to_i) == p_rank)
        {
          EntityProc eproc_new(*to_i, 0 /*new owner*/);
          change.push_back(eproc_new);
        }
      }
    }
  }

  mesh.change_entity_owner(change);

  CEOUtils::checkStatesAfterCEOME_4Elem4ProcEdge(mesh);
}

TEST(BulkData, change_entity_owner_8Elem4ProcMoveTop)
{
  //
  //     id/proc                           id/proc
  //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
  //       |     |     |     |     |         |     |     |     |     |
  //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
  //       |     |     |     |     |         |     |     |     |     |
  //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
  //
  // This test moves ownership of elements 6 and 7 (as well as their locally-owned
  // nodes) to procs 3 and 0, respectively.
  //

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 4)
  {
    return;
  }

  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::unit_test_util::BulkDataTester mesh(meta, pm);

  CEOUtils::fillMeshfor8Elem4ProcMoveTopAndTest(mesh, meta);

  std::vector<stk::mesh::EntityProc> entities_to_move;
  if(mesh.parallel_rank() == 1)
  {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 6);
    int dest_proc = 3;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  if(mesh.parallel_rank() == 2)
  {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 7);
    int dest_proc = 0;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }

  mesh.change_entity_owner(entities_to_move);

  CEOUtils::checkStatesAfterCEOME_8Elem4ProcMoveTop(mesh);
}

TEST(BulkData, change_entity_owner_4Elem4ProcRotate)
{
  //
  //     id/proc                id/proc
  //      7/3---8/2---9/2        7/2---8/1---9/1
  //       |     |     |          |     |     |
  //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
  //       |     |     |          |     |     |
  //      4/0---5/0---6/1  -->   4/3---5/3---6/0
  //       |     |     |          |     |     |
  //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
  //       |     |     |          |     |     |
  //      1/0---2/0---3/1        1/3---2/3---3/0

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 4)
  {
    return;
  }

  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::unit_test_util::BulkDataTester mesh(meta, pm);
  const int p_rank = mesh.parallel_rank();

  CEOUtils::fillMeshfor4Elem4ProcRotateAndTest(mesh, meta);

  std::vector<stk::mesh::EntityProc> entities_to_move;
  if (p_rank == 0) {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
    int dest_proc = 3;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  else if (p_rank == 1) {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
    int dest_proc = 0;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  else if (p_rank == 2) {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
    int dest_proc = 1;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  else if (p_rank == 3) {
    stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 4);
    int dest_proc = 2;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }

  mesh.change_entity_owner(entities_to_move);

  CEOUtils::checkStatesAfterCEOME_4Elem4ProcRotate(mesh, meta);
}

TEST(BulkData, change_entity_owner_3Elem4Proc1Edge3D)
{
  //  ID.proc
  //                    15.2--------16.2                      15.1--------16.1
  //                     /|          /|                        /|          /|
  //                    / |         / |                       / |         / |
  //                  7.2---------8.2 |                     7.1---------8.1 |
  //                   |  |  3.2   |  |                      |  |  3.1   |  |
  //                   |  |        |  |                      |  |        |  |
  //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
  //         /|        | *|        | /|   -->      /|        | *|        | /|
  //        / |        |* |        |/ |           / |        |* |        |/ |
  //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
  //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
  //       |  |        |  |        |  |          |  |        |  |        |  |
  //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
  //       | /         | /         | /           | /         | /         | /
  //       |/          |/          |/            |/          |/          |/
  //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
  //
  //      (*)edge: 1.0                          (*)edge: 1.1

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 4)
  {
    return;
  }

  unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::unit_test_util::BulkDataTester mesh(meta, pm);
  const int p_rank = mesh.parallel_rank();
  CEOUtils::fillMeshfor3Elem4Proc1Edge3DAndTest(mesh, meta);

  std::vector<stk::mesh::EntityProc> entities_to_move;
  if (p_rank == 0) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
    int dest_proc = 3;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);

    elem = mesh.get_entity(stk::topology::EDGE_RANK, 1);
    dest_proc = 1;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
  }
  else if (p_rank == 1) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
    int dest_proc = 0;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }
  else if (p_rank == 2) {
    Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 3);
    int dest_proc = 1;
    entities_to_move.push_back(stk::mesh::EntityProc(elem, dest_proc));
    CEOUtils::add_nodes_to_move(mesh, elem, dest_proc, entities_to_move);
  }

  mesh.change_entity_owner(entities_to_move);

  CEOUtils::checkStatesAfterCEOME_3Elem4Proc1Edge3D(mesh);
}

TEST(BulkData, test_find_ghosted_nodes_that_need_to_be_shared)
{
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& elem_part = meta.declare_part_with_topology("beam2", stk::topology::BEAM_2);
    meta.commit();

    stk::unit_test_util::BulkDataTester bulk(meta, MPI_COMM_WORLD);
    if ( bulk.parallel_size() == 2 )
    {
        bulk.modification_begin();
        stk::mesh::Ghosting &ghosting = bulk.create_ghosting("Ghost Both nodes");
        stk::mesh::Part& ghost_part = bulk.ghosting_part(ghosting);

        std::vector< std::pair<stk::mesh::Entity, int> > ghostingStruct;

        if ( bulk.parallel_rank() == 0 )
        {
            stk::mesh::EntityId nodeId1 = 1;
            stk::mesh::Entity node1 = bulk.declare_node(nodeId1);

            stk::mesh::EntityId nodeId2 = 2;
            stk::mesh::Entity node2 = bulk.declare_node(nodeId2);

            ghostingStruct.push_back(std::make_pair(node1, 1));
            ghostingStruct.push_back(std::make_pair(node2, 1));
        }

        bulk.change_ghosting(ghosting, ghostingStruct);

        stk::mesh::EntityVector ghosted_nodes_that_need_to_be_shared;
        stk::mesh::find_ghosted_nodes_that_need_to_be_shared(bulk, ghosted_nodes_that_need_to_be_shared);

        EXPECT_EQ(0u, ghosted_nodes_that_need_to_be_shared.size());

        ASSERT_NO_THROW(bulk.modification_end());

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(bulk, counts);

        size_t validNumberOfNodes = 2;
        EXPECT_EQ(validNumberOfNodes,counts[stk::topology::NODE_RANK]);

        bulk.modification_begin();

        if ( bulk.parallel_rank() == 1 )
        {
            stk::mesh::EntityId elementId = 1;
            stk::mesh::Entity element1 = bulk.declare_element(elementId, stk::mesh::ConstPartVector{&elem_part});
            stk::mesh::Entity ghostedNode1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
            stk::mesh::Entity ghostedNode2 = bulk.get_entity(stk::topology::NODE_RANK, 2);

            ASSERT_FALSE(bulk.bucket(ghostedNode1).in_aura());
            ASSERT_FALSE(bulk.bucket(ghostedNode2).in_aura());

            ASSERT_TRUE(bulk.bucket(ghostedNode1).member(ghost_part));
            ASSERT_TRUE(bulk.bucket(ghostedNode2).member(ghost_part));

            bulk.declare_relation( element1, ghostedNode1, 0 );
            bulk.declare_relation( element1, ghostedNode2, 1 );
        }

        stk::mesh::find_ghosted_nodes_that_need_to_be_shared(bulk, ghosted_nodes_that_need_to_be_shared);

        if ( bulk.parallel_rank() == 1 )
        {
            EXPECT_EQ(2u, ghosted_nodes_that_need_to_be_shared.size());

            stk::mesh::Entity ghostedNode1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
            EXPECT_EQ(ghostedNode1, ghosted_nodes_that_need_to_be_shared[0]);

            stk::mesh::Entity ghostedNode2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
            EXPECT_EQ(ghostedNode2, ghosted_nodes_that_need_to_be_shared[1]);
        }
        else
        {
            EXPECT_EQ(0u, ghosted_nodes_that_need_to_be_shared.size());
        }

        stk::mesh::fixup_ghosted_to_shared_nodes(bulk);
        ASSERT_NO_THROW(bulk.modification_end());

        stk::mesh::Entity shared_node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
        ASSERT_TRUE(bulk.bucket(shared_node1).shared());

        stk::mesh::Entity shared_node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
        ASSERT_TRUE(bulk.bucket(shared_node2).shared());

        stk::mesh::comm_mesh_counts(bulk, counts);
        EXPECT_EQ(validNumberOfNodes,counts[stk::topology::NODE_RANK]);
    }
}


TEST(BulkData, show_how_one_could_add_a_shared_node)
{
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& elem_part = meta.declare_part_with_topology("triangle", stk::topology::SHELL_TRIANGLE_3);
    meta.commit();

    stk::unit_test_util::BulkDataTester bulk(meta, MPI_COMM_WORLD);

    if ( bulk.parallel_size() == 2 )
    {
        bulk.modification_begin();

        stk::mesh::EntityId nodeId1 = 1;
        bulk.declare_node(nodeId1);

        stk::mesh::EntityId nodeId2 = 2;
        bulk.declare_node(nodeId2);

        if ( bulk.parallel_rank() == 0 )
        {
            stk::mesh::EntityId nodeId3 = 3;
            bulk.declare_node(nodeId3);
        }
        else
        {
            stk::mesh::EntityId nodeId4 = 4;
            bulk.declare_node(nodeId4);
        }

        ASSERT_NO_THROW(bulk.modification_end());

        /// Proc 0 has nodes 1, 2, 3
        /// Proc 1 has nodes 1, 2, 4

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(bulk, counts);

        size_t inValidNumberOfNodes = 6;
        EXPECT_EQ(inValidNumberOfNodes,counts[stk::topology::NODE_RANK]);

        stk::mesh::Entity shared_node1 = bulk.get_entity(stk::topology::NODE_RANK, nodeId1);
        ASSERT_FALSE(bulk.bucket(shared_node1).shared());
        ASSERT_TRUE(bulk.bucket(shared_node1).owned());

        stk::mesh::Entity shared_node2 = bulk.get_entity(stk::topology::NODE_RANK, nodeId2);
        ASSERT_FALSE(bulk.bucket(shared_node2).shared());
        ASSERT_TRUE(bulk.bucket(shared_node2).owned());

        bulk.modification_begin();

        int otherProc = 1-bulk.parallel_rank();
        std::vector<Entity> shared_modified ;

        if ( bulk.parallel_rank() == 0 )
        {
            stk::mesh::EntityId elementId = 1;
            stk::mesh::Entity element = bulk.declare_element(elementId, stk::mesh::ConstPartVector{&elem_part});

            stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
            stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
            stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);

            ASSERT_TRUE(bulk.is_valid(node1));
            ASSERT_TRUE(bulk.is_valid(node2));
            ASSERT_TRUE(bulk.is_valid(node3));
            ASSERT_TRUE(bulk.is_valid(element));

            bulk.declare_relation( element, node1, 0 );
            bulk.declare_relation( element, node2, 1 );
            bulk.declare_relation( element, node3, 2 );

            bulk.add_node_sharing(node1, otherProc);
            bulk.add_node_sharing(node2, otherProc);
            shared_modified.push_back(node1);
            shared_modified.push_back(node2);
        }
        else
        {
            stk::mesh::EntityId elementId = 2;
            stk::mesh::Entity element = bulk.declare_element(elementId, stk::mesh::ConstPartVector{&elem_part});

            stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
            stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2);
            stk::mesh::Entity node4 = bulk.get_entity(stk::topology::NODE_RANK, 4);

            ASSERT_TRUE(bulk.is_valid(node1));
            ASSERT_TRUE(bulk.is_valid(node2));
            ASSERT_TRUE(bulk.is_valid(node4));
            ASSERT_TRUE(bulk.is_valid(element));

            bulk.declare_relation( element, node1, 0 );
            bulk.declare_relation( element, node2, 1 );
            bulk.declare_relation( element, node4, 2 );

            bulk.add_node_sharing(node1, otherProc);
            bulk.add_node_sharing(node2, otherProc);
            shared_modified.push_back(node1);
            shared_modified.push_back(node2);
        }


        for (size_t i=0;i<shared_modified.size();i++)
        {
            std::vector<int> procs;
            bulk.comm_shared_procs(bulk.entity_key(shared_modified[i]), procs);
            int owner = bulk.parallel_rank();
            for (size_t j=0;j<procs.size();j++)
            {
               owner = std::min(owner, procs[j]);
            }
            bulk.my_internal_change_owner_in_comm_data(shared_modified[i], owner);
        }

        ASSERT_NO_THROW(bulk.modification_end());

        stk::mesh::comm_mesh_counts(bulk, counts);
        EXPECT_EQ(4u,counts[stk::topology::NODE_RANK]);
    }
}


void write_mesh(const std::string& filename, stk::mesh::BulkData& mesh)
{
  stk::io::StkMeshIoBroker writer(mesh.parallel());
  writer.set_bulk_data(mesh);
  size_t output_handle = writer.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  writer.write_output_mesh(output_handle);
}

void test_nodes(stk::mesh::BulkData &mesh)
{
    stk::mesh::EntityId node1Id = 1;
    stk::mesh::EntityId node2Id = 2;
    stk::mesh::EntityId node3Id = 3;
    stk::mesh::EntityId node4Id = 4;

    stk::mesh::Entity node1 = mesh.get_entity(stk::topology::NODE_RANK, node1Id);
    stk::mesh::Entity node2 = mesh.get_entity(stk::topology::NODE_RANK, node2Id);
    stk::mesh::Entity node3 = mesh.get_entity(stk::topology::NODE_RANK, node3Id);
    stk::mesh::Entity node4 = mesh.get_entity(stk::topology::NODE_RANK, node4Id);

    stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);

    if ( mesh.parallel_rank() == 0)
    {
        EXPECT_TRUE(mesh.bucket(node1).owned());
        EXPECT_TRUE(mesh.bucket(node2).owned());
        EXPECT_TRUE(mesh.bucket(node3).owned());
        EXPECT_FALSE(mesh.is_valid(node4));

        EXPECT_TRUE(mesh.bucket(node1).shared());
        EXPECT_TRUE(mesh.bucket(node2).shared());
        EXPECT_FALSE(mesh.bucket(node3).shared());

        EXPECT_FALSE(mesh.bucket(node1).in_aura());
        EXPECT_FALSE(mesh.bucket(node2).in_aura());
        EXPECT_FALSE(mesh.bucket(node3).in_aura());

        EXPECT_TRUE(mesh.bucket(element1).owned());
    }
    else if ( mesh.parallel_rank() == 1)
    {
        EXPECT_FALSE(mesh.bucket(node1).owned());
        EXPECT_FALSE(mesh.bucket(node2).owned());
        EXPECT_FALSE(mesh.bucket(node3).owned());
        EXPECT_TRUE(mesh.bucket(node4).owned());

        EXPECT_TRUE(mesh.bucket(node1).shared());
        EXPECT_TRUE(mesh.bucket(node2).shared());
        EXPECT_FALSE(mesh.bucket(node3).shared());
        EXPECT_FALSE(mesh.bucket(node4).shared());

        EXPECT_FALSE(mesh.bucket(node1).in_aura());
        EXPECT_FALSE(mesh.bucket(node2).in_aura());
        EXPECT_TRUE(mesh.bucket(node3).in_aura());
        EXPECT_FALSE(mesh.bucket(node4).in_aura());

        EXPECT_FALSE(mesh.bucket(element1).owned());
        EXPECT_TRUE(mesh.bucket(element1).in_aura());
    }
}

TEST(BulkData, can_we_create_shared_nodes)
{
    int numProcs = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if ( numProcs == 2 )
    {
        std::string filename = "test.exo";

        {
            unsigned spatialDim = 3;
            stk::mesh::MetaData meta(spatialDim);

            stk::mesh::Selector all_nodes = meta.universal_part();
            typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> CoordFieldType;
            CoordFieldType& coordField = meta.declare_field<CoordFieldType>(stk::topology::NODE_RANK, "model_coordinates");
            stk::mesh::put_field_on_mesh(coordField, all_nodes,
                                         (stk::mesh::FieldTraits<CoordFieldType>::data_type*) nullptr);

            stk::mesh::Part& elem_part = meta.declare_part_with_topology("block_1", stk::topology::BEAM_2);
            stk::io::put_io_part_attribute(elem_part);

            typedef stk::mesh::Field<double> ScalarField;
            ScalarField& temperatureField = meta.declare_field<ScalarField>(stk::topology::NODE_RANK, "temperature");

            double initialTemperatureValue = 25.0;
            stk::mesh::put_field_on_entire_mesh_with_initial_value(temperatureField, &initialTemperatureValue);

            meta.commit();

            //////////////////////////////////////////
            /*
                    Proc 0   |     Proc 1
                             |
                    3    1   |    1       4
                    o----o   |    o       o
                      E1     |
                             |
                         2   |    2
                         o   |    o
                             |
             */

            stk::unit_test_util::BulkDataTester bulk(meta, MPI_COMM_WORLD);

            bulk.modification_begin();

            double xCoords[4] = {  0.,  1.,  1.,  1. };
            double yCoords[4] = {  0.,  0.,  1.,  1. };
            double zCoords[4] = {  0.,  0.,  0.,  1. };

            stk::mesh::EntityId nodeId1 = 1;
            stk::mesh::Entity node1 = bulk.declare_node(nodeId1);

            stk::mesh::EntityId nodeId2 = 2;
            stk::mesh::Entity node2 = bulk.declare_node(nodeId2);

            stk::mesh::Entity thirdNode;
            if ( bulk.parallel_rank() == 0 )
            {
                stk::mesh::EntityId nodeId3 = 3;
                thirdNode = bulk.declare_node(nodeId3);
            }
            else
            {
                stk::mesh::EntityId nodeId4 = 4;
                thirdNode = bulk.declare_node(nodeId4);
            }

            int otherProc = 1 - bulk.parallel_rank();
            bulk.add_node_sharing(node1, otherProc);
            bulk.add_node_sharing(node2, otherProc);

            if ( bulk.parallel_rank() == 0 )
            {
                stk::mesh::EntityId elementId = 1 + bulk.parallel_rank();
                stk::mesh::Entity element = bulk.declare_element(elementId, stk::mesh::ConstPartVector{&elem_part});
                bulk.declare_relation(element, node1, 0);
                bulk.declare_relation(element, thirdNode, 1);
            }

            ASSERT_NO_THROW(bulk.modification_end());

            test_nodes(bulk);

            ////////////////////////
            // Made it here. nodes 2 and 3 are shared, but are not connected to any higher
            // ranked entities on processor 1

            // Add coordinates so mesh can be written out to file
            int node_id = bulk.identifier(node1);
            double* coords = stk::mesh::field_data(coordField, node1);
            coords[0] = xCoords[node_id-1];
            coords[1] = yCoords[node_id-1];
            coords[2] = zCoords[node_id-1];

            node_id = bulk.identifier(node2);
            coords = stk::mesh::field_data(coordField, node2);
            coords[0] = xCoords[node_id-1];
            coords[1] = yCoords[node_id-1];
            coords[2] = zCoords[node_id-1];

            node_id = bulk.identifier(thirdNode);
            coords = stk::mesh::field_data(coordField, thirdNode);
            coords[0] = xCoords[node_id-1];
            coords[1] = yCoords[node_id-1];
            coords[2] = zCoords[node_id-1];

            std::vector<size_t> counts;
            stk::mesh::comm_mesh_counts(bulk, counts);
            size_t goldNumberNodes = 4;
            EXPECT_EQ(goldNumberNodes, counts[stk::topology::NODE_RANK]);
            write_mesh(filename, bulk);

            // Add temperature field to make sure parallel field communications work

            double owned_value = 10;
            double shared_value = 1;

            double *tempField = stk::mesh::field_data(temperatureField, node1);
            if ( bulk.bucket(node1).owned() )
            {
                *tempField = owned_value;
            }
            else
            {
                *tempField = shared_value;
            }

            tempField = stk::mesh::field_data(temperatureField, node2);
            if ( bulk.bucket(node2).owned() )
            {
                *tempField = owned_value;
            }
            else
            {
                *tempField = shared_value;
            }

            tempField = stk::mesh::field_data(temperatureField, thirdNode);
            if ( bulk.bucket(thirdNode).owned() )
            {
                *tempField = owned_value;
            }
            else
            {
                *tempField = shared_value;
            }

            std::vector<const stk::mesh::FieldBase*> fields(1, &temperatureField);
            // send data, owned to shared
            stk::mesh::communicate_field_data(bulk.shared_ghosting(), fields); /* IMPORTANT PART TO TEST */

            tempField = stk::mesh::field_data(temperatureField, node1);
            EXPECT_NEAR(owned_value, *tempField, 1e-6);
            tempField = stk::mesh::field_data(temperatureField, node2);
            EXPECT_NEAR(owned_value, *tempField, 1e-6);
            tempField = stk::mesh::field_data(temperatureField, thirdNode);
            EXPECT_NEAR(owned_value, *tempField, 1e-6);

            // node 3 on processor 1 is ghosted, so it still has initial value
            if ( bulk.parallel_rank() == 1 )
            {
                stk::mesh::Entity ghostedNode = bulk.get_entity(stk::topology::NODE_RANK, 3);
                tempField = stk::mesh::field_data(temperatureField, ghostedNode);
                EXPECT_NEAR(initialTemperatureValue, *tempField, 1e-6);
            }

            // send data, owned to ghosted
            stk::mesh::communicate_field_data(bulk.aura_ghosting(), fields); /* IMPORTANT PART TO TEST */

            if ( bulk.parallel_rank() == 1 )
            {
                stk::mesh::Entity ghostedNode = bulk.get_entity(stk::topology::NODE_RANK, 3);
                tempField = stk::mesh::field_data(temperatureField, ghostedNode);
                EXPECT_NEAR(owned_value, *tempField, 1e-6);
            }

            std::vector<const stk::mesh::FieldBase*> fields1(1, &temperatureField);
            stk::mesh::parallel_sum(bulk, fields1); /* IMPORTANT PART TO TEST */

            double summed_value = 2*owned_value;
            tempField = stk::mesh::field_data(temperatureField, node1);
            EXPECT_NEAR(summed_value, *tempField, 1e-6);
            tempField = stk::mesh::field_data(temperatureField, node2);
            EXPECT_NEAR(summed_value, *tempField, 1e-6);

            double min_value = 1;
            if ( bulk.parallel_rank() == 1 )
            {
                tempField = stk::mesh::field_data(temperatureField, node1);
                *tempField = min_value;
            }

            stk::mesh::parallel_min(bulk, fields1); /* IMPORTANT PART TO TEST */
            tempField = stk::mesh::field_data(temperatureField, node1);
            EXPECT_NEAR(min_value, *tempField, 1e-6);

            double max_value = 100;
            if ( bulk.parallel_rank() == 1 )
            {
                tempField = stk::mesh::field_data(temperatureField, node1);
                *tempField = max_value;
            }

            stk::mesh::parallel_max(bulk, fields1); /* IMPORTANT PART TO TEST */
            tempField = stk::mesh::field_data(temperatureField, node1);
            EXPECT_NEAR(max_value, *tempField, 1e-6);

            std::vector<EntityProc> entities_to_change;
            if ( bulk.parallel_rank() == 0 )
            {
                entities_to_change.push_back(std::make_pair(node2, 1));
            }

            bulk.change_entity_owner(entities_to_change);

            if (bulk.parallel_rank() == 0)
            {
                EXPECT_FALSE(bulk.is_valid(node2));
            }
            else // bulk.parallel_rank() == 1
            {
                ASSERT_TRUE(bulk.is_valid(node2));
                EXPECT_FALSE(bulk.bucket(node2).shared());
                EXPECT_EQ(1u,bulk.closure_count(node2));
            }

        }

        // Now read in the mesh and see that it passes same tests

        {
            unsigned spatialDim = 3;
            stk::mesh::MetaData meta(spatialDim);
            stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
            stk::io::StkMeshIoBroker reader(mesh.parallel());
            reader.set_bulk_data(mesh);
            reader.add_mesh_database(filename, stk::io::READ_MESH);
            reader.create_input_mesh();
            reader.populate_bulk_data();

            test_nodes(mesh);

            std::vector<size_t> counts;
            stk::mesh::comm_mesh_counts(mesh, counts);
            EXPECT_EQ(4u,counts[stk::topology::NODE_RANK]);

            std::ostringstream os;
            os << filename << ".2." << mesh.parallel_rank();
            std::string fileToRemove = os.str();
            unlink(fileToRemove.c_str());
        }
    }
}

typedef std::pair< std::vector<const stk::mesh::Entity*>, stk::mesh::Entity* > ChildNodeRequest;

struct edgeToBeRefined
{
    stk::mesh::EntityId m_node1;
    stk::mesh::EntityId m_node2;
    std::vector<int> m_sharing_procs;
    std::vector<std::pair<int, stk::mesh::EntityId> > m_id_proc_pairs_from_all_procs;
    stk::mesh::BulkData *m_stk_mesh;
    bool m_id_procs_pairs_have_been_sorted;
    stk::mesh::Entity **m_node_entity;

    edgeToBeRefined(stk::mesh::BulkData& bulkData, stk::mesh::EntityId nodeA, stk::mesh::EntityId nodeB, stk::mesh::Entity **entity_place_holder=NULL)
    : m_node1(0), m_node2(0), m_sharing_procs(), m_id_proc_pairs_from_all_procs(), m_stk_mesh(&bulkData),
      m_id_procs_pairs_have_been_sorted(false), m_node_entity(entity_place_holder)
    {
        if ( nodeA < nodeB )
        {
            m_node1 = nodeA;
            m_node2 = nodeB;
        }
        else
        {
            m_node1 = nodeB;
            m_node2 = nodeA;
        }
    }

    void add_proc_id_pair(int proc_id, stk::mesh::EntityId id)
    {
        m_id_proc_pairs_from_all_procs.push_back(std::make_pair(proc_id, id));
    }

    void calculate_sharing_procs()
    {
        std::vector<int> sharingProcsA;
        stk::mesh::EntityKey keyA(stk::topology::NODE_RANK, m_node1);
        m_stk_mesh->comm_shared_procs(keyA, sharingProcsA);

        std::vector<int> sharingProcsB;
        stk::mesh::EntityKey keyB(stk::topology::NODE_RANK, m_node2);
        m_stk_mesh->comm_shared_procs(keyB, sharingProcsB);

        std::sort(sharingProcsA.begin(), sharingProcsA.end());
        std::sort(sharingProcsB.begin(), sharingProcsB.end());

        std::set_intersection(sharingProcsA.begin(),sharingProcsA.end(),sharingProcsB.begin(),sharingProcsB.end(),std::back_inserter(m_sharing_procs));
    }

    size_t num_sharing_procs() const
    {
        return m_sharing_procs.size();
    }

    int sharing_proc(int index) const
    {
        return m_sharing_procs[index];
    }

    stk::mesh::EntityId suggested_node_id() const
    {
        ThrowRequireMsg(!m_id_procs_pairs_have_been_sorted, "Invalid use of edge calculation. Contact sierra-help");
        return m_id_proc_pairs_from_all_procs[0].second;
    }

    stk::mesh::EntityId node1() const
    {
        return m_node1;
    }

    stk::mesh::EntityId node2() const
    {
        return m_node2;
    }

    void sort_id_proc_pairs()
    {
        m_id_procs_pairs_have_been_sorted = true;
        std::sort(m_id_proc_pairs_from_all_procs.begin(), m_id_proc_pairs_from_all_procs.end());
    }

    stk::mesh::EntityId get_id_for_edge() const
    {
        ThrowRequireMsg(m_id_procs_pairs_have_been_sorted, "Invalid use of edge calculation. Contact sierra-help");
        return m_id_proc_pairs_from_all_procs[0].second;
    }

    void set_node_entity_for_edge()
    {
        this->sort_id_proc_pairs();
        stk::mesh::EntityId id_for_edge = this->get_id_for_edge();
        *(*m_node_entity) = m_stk_mesh->declare_node(id_for_edge);
        for (size_t i=0;i<m_id_proc_pairs_from_all_procs.size();++i)
        {
            if ( m_id_proc_pairs_from_all_procs[i].first != m_stk_mesh->parallel_rank() )
            {
                m_stk_mesh->add_node_sharing(**m_node_entity, m_id_proc_pairs_from_all_procs[i].first);
            }
        }
    }

    bool operator==(const edgeToBeRefined& otherEdge) const
    {
        if ( this->node1() == otherEdge.node1() &&
                this->node2() == otherEdge.node2() )
        {
            return true;
        }
        return false;
    }
};

struct edge_finder
{
    bool operator()(const edgeToBeRefined& edge1, const edgeToBeRefined& edge2) const
    {
        if ( edge1.node1() != edge2.node1() )
        {
            return edge1.node1() < edge2.node1();
        }
        else
        {
            return edge1.node2() < edge2.node2();
        }
    }
};

/*
TEST(BulkData, test_edge_finder)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    int num_procs = stk::parallel_machine_size(comm);

    std::ostringstream os;
    os << "generated:2x2x" << num_procs;

    std::string filename = os.str();

    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker reader(mesh.parallel());

    reader.set_bulk_data(mesh);
    reader.add_mesh_database("generated:2x2x2", stk::io::READ_MESH);
    reader.create_input_mesh();
    reader.populate_bulk_data();

    edge_finder edge(mesh, 10, 11);


}
*/

void batch_create_child_nodes_new(BulkData & mesh, std::vector< ChildNodeRequest > & child_node_requests)
{
    bool communicate_nodes = true;
    std::vector<bool> communicate_edge(child_node_requests.size(), false);

    unsigned num_nodes_requested = child_node_requests.size();
    std::vector<stk::mesh::EntityId> available_node_ids;
    mesh.generate_new_ids(stk::topology::NODE_RANK, num_nodes_requested, available_node_ids);

    while ( communicate_nodes )
    {
        communicate_nodes = false;

        std::vector<edgeToBeRefined> edgesToBeRefined;

        for (unsigned it_req=0; it_req<child_node_requests.size(); ++it_req)
        {
            ChildNodeRequest & request = child_node_requests[it_req];
            std::vector<const stk::mesh::Entity*> & request_parents = request.first;

            ThrowRequireMsg(request_parents.size() == 2, "Invalid size of request, needed exactly 2 parents, found " << request_parents.size() << ". Contact sierra-help.");
            if (mesh.is_valid(*request_parents[0]) && mesh.is_valid(*request_parents[1]) && communicate_edge[it_req] == false )
            {
                stk::mesh::Entity **request_child = &request.second;

                communicate_edge[it_req] = true;
                communicate_nodes = true;
                stk::mesh::EntityId node1 = mesh.identifier(*request_parents[0]);
                stk::mesh::EntityId node2 = mesh.identifier(*request_parents[1]);

                edgesToBeRefined.push_back(edgeToBeRefined(mesh, node1, node2, request_child));
                edgesToBeRefined.back().add_proc_id_pair(mesh.parallel_rank(), available_node_ids[it_req]);
                edgesToBeRefined.back().calculate_sharing_procs();
            }
        }

        int more_work_to_be_done = communicate_nodes ? 1 : 0;
        int global_result = 0;
        stk::all_reduce_max(mesh.parallel(), &more_work_to_be_done, &global_result, 1);
        if ( global_result == 0 ) break;

        std::sort(edgesToBeRefined.begin(), edgesToBeRefined.end(), edge_finder());

        // By this point, I have list of edges that I need to communicate

        stk::CommSparse comm_spec(mesh.parallel());

        for (int phase=0;phase<2;++phase)
        {
            for (size_t edge_index=0;edge_index<edgesToBeRefined.size();++edge_index)
            {
                for (size_t proc_index=0;proc_index<edgesToBeRefined[edge_index].num_sharing_procs();++proc_index)
                {
                    int other_proc = edgesToBeRefined[edge_index].sharing_proc(proc_index);
                    stk::mesh::EntityId node_1 = edgesToBeRefined[edge_index].node1();
                    stk::mesh::EntityId node_2 = edgesToBeRefined[edge_index].node2();
                    stk::mesh::EntityId this_procs_suggested_id = edgesToBeRefined[edge_index].suggested_node_id();
                    comm_spec.send_buffer(other_proc).pack<stk::mesh::EntityId>(node_1).pack<stk::mesh::EntityId>(node_2);
                    comm_spec.send_buffer(other_proc).pack<stk::mesh::EntityId>(this_procs_suggested_id);
                }
            }

            if ( phase == 0 )
            {
                comm_spec.allocate_buffers();
            }
            else
            {
                comm_spec.communicate();
            }
        }

        for(int i = 0; i < mesh.parallel_size(); ++i)
        {
            if(i != mesh.parallel_rank())
            {
                while(comm_spec.recv_buffer(i).remaining())
                {
                    stk::mesh::EntityId node1;
                    stk::mesh::EntityId node2;
                    stk::mesh::EntityId suggested_node_id;
                    comm_spec.recv_buffer(i).unpack<stk::mesh::EntityId>(node1);
                    comm_spec.recv_buffer(i).unpack<stk::mesh::EntityId>(node2);
                    comm_spec.recv_buffer(i).unpack<stk::mesh::EntityId>(suggested_node_id);

                    edgeToBeRefined from_other_proc(mesh, node1, node2);
                    std::vector<edgeToBeRefined>::iterator iter = std::lower_bound(edgesToBeRefined.begin(), edgesToBeRefined.end(),
                                   from_other_proc, edge_finder());

                    if ( iter != edgesToBeRefined.end() && *iter == from_other_proc)
                    {
                        iter->add_proc_id_pair(i, suggested_node_id);
                    }
                }
            }
        }

        for (size_t edge_index=0;edge_index<edgesToBeRefined.size();++edge_index)
        {
            edgesToBeRefined[edge_index].set_node_entity_for_edge();
        }
    }

    std::vector<bool>::iterator iter = std::find(communicate_edge.begin(), communicate_edge.end(), false);
    ThrowRequireMsg(iter == communicate_edge.end(), "Invalid edge requests. Contact sierra-help.");
}

void set_coords_on_new_node(stk::mesh::MetaData& meta, stk::mesh::Entity nodeA, stk::mesh::Entity nodeB, stk::mesh::Entity new_node)
{
    stk::mesh::FieldBase const * coord = meta.coordinate_field();

    //typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> CoordFieldType;
    //CoordFieldType& field = meta.declare_field<CoordFieldType>(stk::topology::NODE_RANK, "model_coordinates");

    double x = 0, y = 0, z = 0;

    double *fieldValue = static_cast<double *>(stk::mesh::field_data(*coord, nodeA));
    x += fieldValue[0]*0.5;
    y += fieldValue[1]*0.5;
    z += fieldValue[2]*0.5;

    fieldValue = static_cast<double *>(stk::mesh::field_data(*coord, nodeB));
    x += fieldValue[0]*0.5;
    y += fieldValue[1]*0.5;
    z += fieldValue[2]*0.5;

    fieldValue = static_cast<double *>(stk::mesh::field_data(*coord, new_node));
    fieldValue[0] = x;
    fieldValue[1] = y;
    fieldValue[2] = z;
}

TEST(BulkData, create_vigilante_nodes_along_shared_edge)
{
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& node_part = meta.declare_part_with_topology("nodelist_1", stk::topology::NODE);
    stk::io::put_io_part_attribute(node_part);
//    meta.commit();

    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker reader(mesh.parallel());

    if ( mesh.parallel_size() == 2 )
    {
        reader.set_bulk_data(mesh);
        reader.add_mesh_database("generated:2x2x2", stk::io::READ_MESH);
        reader.create_input_mesh();
        reader.populate_bulk_data();

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(mesh, counts);

        EXPECT_EQ(27u,counts[stk::topology::NODE_RANK]);
        EXPECT_EQ(8u,counts[stk::topology::ELEMENT_RANK]);

        // nodes 16, 17, 18
        // nodes 13, 14, 15
        // nodes 10, 11, 12

        stk::mesh::Entity node_10 = mesh.get_entity(stk::topology::NODE_RANK, 10);
        stk::mesh::Entity node_11 = mesh.get_entity(stk::topology::NODE_RANK, 11);

        ASSERT_TRUE(mesh.bucket(node_10).shared());
        ASSERT_TRUE(mesh.bucket(node_11).shared());

        {
            std::vector<ChildNodeRequest> child_node_requests;

            stk::mesh::Entity node_between_10_and_11= stk::mesh::Entity();
            {
                std::vector<const stk::mesh::Entity*> node_parents;
                node_parents.push_back(&node_10);
                node_parents.push_back(&node_11);
                ChildNodeRequest node_request(node_parents, &node_between_10_and_11);
                child_node_requests.push_back(node_request);
            }

            stk::mesh::Entity node_between_10_and_new_node = stk::mesh::Entity();
            {
                std::vector<const stk::mesh::Entity*> node_parents;
                node_parents.push_back(&node_10);
                node_parents.push_back(&node_between_10_and_11);
                ChildNodeRequest node_request(node_parents, &node_between_10_and_new_node);
                child_node_requests.push_back(node_request);
            }

            stk::mesh::Entity node_between_10_and_second_new_node = stk::mesh::Entity();
            {
                std::vector<const stk::mesh::Entity*> node_parents;
                node_parents.push_back(&node_10);
                node_parents.push_back(&node_between_10_and_new_node);
                ChildNodeRequest node_request(node_parents, &node_between_10_and_second_new_node);
                child_node_requests.push_back(node_request);
            }


            mesh.modification_begin();

            batch_create_child_nodes_new(mesh, child_node_requests);

            stk::mesh::PartVector add_parts_vec;
            stk::mesh::PartVector rem_parts_vec;
            add_parts_vec.push_back(&node_part);
            mesh.change_entity_parts(node_between_10_and_11, add_parts_vec, rem_parts_vec);
            mesh.change_entity_parts(node_between_10_and_new_node, add_parts_vec, rem_parts_vec);
            mesh.change_entity_parts(node_between_10_and_second_new_node, add_parts_vec, rem_parts_vec);


            ASSERT_NO_THROW(mesh.modification_end());

            for (size_t i=0;i<child_node_requests.size();++i)
            {
                ChildNodeRequest & request = child_node_requests[i];
                std::vector<const stk::mesh::Entity*> & request_parents = request.first;
                stk::mesh::Entity *new_node = request.second;
                set_coords_on_new_node(mesh.mesh_meta_data(), *request_parents[0], *request_parents[1], *new_node);
            }
        }

        stk::mesh::comm_mesh_counts(mesh, counts);

        EXPECT_EQ(30u,counts[stk::topology::NODE_RANK]);
        EXPECT_EQ(8u,counts[stk::topology::ELEMENT_RANK]);

        //std::string filename = "test.exo";
        //write_mesh(filename, mesh);
    }
}

TEST(BulkData, renegade_nodes_along_every_edge)
{
    unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& node_part = meta.declare_part_with_topology("nodelist_1", stk::topology::NODE);
    stk::io::put_io_part_attribute(node_part);

    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::io::StkMeshIoBroker reader(mesh.parallel());

    if(mesh.parallel_size() == 2)
    {
        reader.set_bulk_data(mesh);
        reader.add_mesh_database("generated:2x2x2", stk::io::READ_MESH);
        reader.create_input_mesh();
        reader.populate_bulk_data();

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(mesh, counts);

        size_t num_original_nodes = counts[stk::topology::NODE_RANK];

        EXPECT_EQ(27u, num_original_nodes);
        EXPECT_EQ(8u, counts[stk::topology::ELEMENT_RANK]);

        stk::mesh::create_edges(mesh);

        stk::mesh::comm_mesh_counts(mesh, counts);
        size_t num_edges = counts[stk::topology::EDGE_RANK];
        EXPECT_EQ(54u, num_edges);
        std::vector<stk::mesh::Entity> newNodes(3*num_edges, stk::mesh::Entity());
        size_t new_node_counter = 0;

        std::vector<ChildNodeRequest> child_node_requests;

        const stk::mesh::BucketVector &buckets = mesh.buckets(stk::topology::EDGE_RANK);

        for(size_t b = 0; b < buckets.size(); ++b)
        {
            const stk::mesh::Bucket &bucket = *buckets[b];
            if(!bucket.in_aura())
            {
                for (size_t e=0;e<bucket.size();++e)
                {
                    const stk::mesh::Entity* pnodes = mesh.begin_nodes(bucket[e]);

                    int child_index = 0;
                    {
                        std::vector<const stk::mesh::Entity*> node_parents;
                        node_parents.push_back(pnodes);
                        node_parents.push_back(pnodes+1);
                        ChildNodeRequest node_request(node_parents, &newNodes[new_node_counter]);
                        child_index = new_node_counter;
                        ++new_node_counter;
                        child_node_requests.push_back(node_request);
                    }
                    {
                        std::vector<const stk::mesh::Entity*> node_parents;
                        node_parents.push_back(&newNodes[child_index]);
                        node_parents.push_back(pnodes+1);
                        ChildNodeRequest node_request(node_parents, &newNodes[new_node_counter]);
                        ++new_node_counter;
                        child_node_requests.push_back(node_request);
                    }
                    {
                        std::vector<const stk::mesh::Entity*> node_parents;
                        node_parents.push_back(pnodes);
                        node_parents.push_back(&newNodes[child_index]);
                        ChildNodeRequest node_request(node_parents, &newNodes[new_node_counter]);
                        ++new_node_counter;
                        child_node_requests.push_back(node_request);
                    }
                    ASSERT_TRUE(new_node_counter <= 3*num_edges);
                }
            }
        }

        {
            mesh.modification_begin();

            batch_create_child_nodes_new(mesh, child_node_requests);

            stk::mesh::PartVector add_parts_vec;
            stk::mesh::PartVector rem_parts_vec;
            add_parts_vec.push_back(&node_part);
            for (size_t e=0;e<new_node_counter;++e)
            {
                mesh.change_entity_parts(newNodes[e], add_parts_vec, rem_parts_vec);
            }

            ASSERT_NO_THROW(mesh.modification_end());

            for(size_t i = 0; i < child_node_requests.size(); ++i)
            {
                ChildNodeRequest & request = child_node_requests[i];
                std::vector<const stk::mesh::Entity*> & request_parents = request.first;
                stk::mesh::Entity *new_node = request.second;
                set_coords_on_new_node(mesh.mesh_meta_data(), *request_parents[0], *request_parents[1], *new_node);
            }
        }

        stk::mesh::comm_mesh_counts(mesh, counts);

        size_t gold_number_nodes = num_original_nodes + 3*num_edges;
        EXPECT_EQ(gold_number_nodes, counts[stk::topology::NODE_RANK]);
        EXPECT_EQ(8u, counts[stk::topology::ELEMENT_RANK]);

        //std::string filename = "test.exo";
        //write_mesh(filename, mesh);
    }
}


TEST(BulkData, show_API_for_batch_create_child_nodes)
{
    unsigned spatialDim = 2;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::Part& elem_part = meta.declare_part_with_topology("triangle", stk::topology::TRIANGLE_3_2D);

    meta.commit();

    stk::unit_test_util::BulkDataTester bulk(meta, MPI_COMM_WORLD);

    if ( bulk.parallel_size() != 2 ) return;

    bulk.modification_begin();

    stk::mesh::Entity node1 = bulk.declare_node(1);
    stk::mesh::Entity node2 = bulk.declare_node(2);

    int otherProc = 1-bulk.parallel_rank();
    bulk.add_node_sharing(node1, otherProc);
    bulk.add_node_sharing(node2, otherProc);

    if ( bulk.parallel_rank() == 0 )
    {
      bulk.declare_node(3);

      stk::mesh::EntityIdVector connected_nodes {1, 2, 3 };
      stk::mesh::declare_element(bulk, elem_part, 1, connected_nodes);
    }
    else
    {
      bulk.declare_node(4);

      stk::mesh::EntityIdVector connected_nodes {2, 1, 4 };
      stk::mesh::declare_element(bulk, elem_part, 2, connected_nodes);
    }

    ASSERT_NO_THROW(bulk.modification_end());

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);

    EXPECT_EQ(4u,counts[stk::topology::NODE_RANK]);
    EXPECT_EQ(2u,counts[stk::topology::ELEMENT_RANK]);

    ASSERT_TRUE(bulk.bucket(node1).shared());
    ASSERT_TRUE(bulk.bucket(node2).shared());

    // Nodes 3 and 4 accessible on both processors via aura
    stk::mesh::Entity node3 = bulk.get_entity(stk::topology::NODE_RANK, 3);
    stk::mesh::Entity node4 = bulk.get_entity(stk::topology::NODE_RANK, 4);

    ASSERT_FALSE(bulk.bucket(node3).shared());
    ASSERT_FALSE(bulk.bucket(node4).shared());

    std::vector<ChildNodeRequest> child_node_requests;

    stk::mesh::Entity node5 = stk::mesh::Entity(); // "child" of 1 and 2 is 5
    std::vector<const stk::mesh::Entity*> node5_parents;
    node5_parents.push_back(&node1);
    node5_parents.push_back(&node2);
    ChildNodeRequest node5_request(node5_parents, &node5);
    child_node_requests.push_back(node5_request);

    stk::mesh::Entity node6 = stk::mesh::Entity(); // "child" of 1 and 5 is 6
    std::vector<const stk::mesh::Entity*> node6_parents;
    node6_parents.push_back(&node1);
    node6_parents.push_back(&node5);
    ChildNodeRequest node6_request(node6_parents, &node6);
    child_node_requests.push_back(node6_request);

    // It would be nice if child node creation was processed a batch of requests and didn't have to
    // be part of the modification cycle where the elements are attached.  Currently, this would throw
    // an error because of the shared orphan nodes.
    bulk.modification_begin();
    batch_create_child_nodes_new(bulk, child_node_requests);

    if ( bulk.parallel_rank() == 0 )
    {
      stk::mesh::EntityIdVector connected_nodes(3);
      connected_nodes[0] = 1;
      connected_nodes[1] = bulk.identifier(node6);
      connected_nodes[2] = 3;
      stk::mesh::declare_element(bulk, elem_part, 3, connected_nodes);

      connected_nodes[0] = bulk.identifier(node6);
      connected_nodes[1] = bulk.identifier(node5);
      connected_nodes[2] = 3;
      stk::mesh::declare_element(bulk, elem_part, 4, connected_nodes);

      connected_nodes[0] = bulk.identifier(node5);
      connected_nodes[1] = 2;
      connected_nodes[2] = 3;
      stk::mesh::declare_element(bulk, elem_part, 5, connected_nodes);
    }
    else
    {
      stk::mesh::EntityIdVector connected_nodes(3);
      connected_nodes[0] = bulk.identifier(node6);
      connected_nodes[1] = 1;
      connected_nodes[2] = 4;
      stk::mesh::declare_element(bulk, elem_part, 6, connected_nodes);

      connected_nodes[0] = bulk.identifier(node5);
      connected_nodes[1] = bulk.identifier(node6);
      connected_nodes[2] = 4;
      stk::mesh::declare_element(bulk, elem_part, 7, connected_nodes);

      connected_nodes[0] = 2;
      connected_nodes[1] = bulk.identifier(node5);
      connected_nodes[2] = 4;
      stk::mesh::declare_element(bulk, elem_part, 8, connected_nodes);
    }

    ASSERT_NO_THROW(bulk.modification_end());


    ASSERT_TRUE(bulk.is_valid(node5));
    ASSERT_TRUE(bulk.is_valid(node6));
    ASSERT_TRUE(bulk.bucket(node5).shared());
    ASSERT_TRUE(bulk.bucket(node6).shared());

    stk::mesh::comm_mesh_counts(bulk, counts);

    EXPECT_EQ(6u,counts[stk::topology::NODE_RANK]);
    EXPECT_EQ(8u,counts[stk::topology::ELEMENT_RANK]);

    //std::string filename = "test.exo";
    //write_mesh(filename, bulk);
}

void Test_STK_ParallelPartConsistency_ChangeBlock(stk::mesh::BulkData::AutomaticAuraOption autoAuraOption)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size != 2) return;

  // This test will create a two-element mesh (quad4 elements)
  // in 2 blocks in parallel.  Both elements will start out in block_1
  // and then one of the elements will be moved to block_2.

  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData mesh(meta, pm, autoAuraOption);

  //declare 'block_1' which will hold element 1
  stk::mesh::Part& block_1 = meta.declare_part("block_1", stk::topology::ELEMENT_RANK);
  stk::mesh::set_topology(block_1, stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part("block_2", stk::topology::ELEMENT_RANK);
  stk::mesh::set_topology(block_2, stk::topology::QUAD_4_2D);

  stk::mesh::Selector all_nodes = meta.universal_part();

  //declare a field for coordinates
  typedef stk::mesh::Field<double, stk::mesh::Cartesian2d> CoordFieldType;
  CoordFieldType& coordField = meta.declare_field<CoordFieldType>(stk::topology::NODE_RANK, "model_coordinates");
  stk::mesh::put_field_on_mesh(coordField, all_nodes,
                               (stk::mesh::FieldTraits<CoordFieldType>::data_type*) nullptr);
  stk::mesh::Field<double>& oneField = meta.declare_field< stk::mesh::Field<double> >(stk::topology::NODE_RANK, "field_of_one");
  stk::mesh::put_field_on_mesh(oneField, block_1,
                               (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);

  meta.commit();
  mesh.modification_begin();

  const size_t nodesPerElem = 4;
  const size_t numNodes = 6;

  double xCoords[numNodes] = { -1.,  0.,  1.,  1.,  0., -1. };
  double yCoords[numNodes] = {  0.,  0.,  0.,  1.,  1.,  1. };
  int elem_nodes0[] = {0, 1, 4, 5};
  int elem_nodes1[] = {1, 2, 3, 4};
  int * elem_nodes[] = { elem_nodes0, elem_nodes1 };

  //Next create nodes and set up connectivity to use later for creating the element.
  stk::mesh::EntityIdVector connected_nodes(nodesPerElem);
  for(size_t n=0; n<nodesPerElem; ++n) {
    size_t e = parallel_rank;
    stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n]+1;

    stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
    if (!mesh.is_valid(node))
    {
      node = mesh.declare_node(nodeGlobalId);
    }

    connected_nodes[n] = nodeGlobalId;
  }

  stk::mesh::Entity node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  stk::mesh::Entity node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

  if (parallel_rank == 0)
  {
    mesh.add_node_sharing(node2, 1);
    mesh.add_node_sharing(node5, 1);
  }
  else
  {
    mesh.add_node_sharing(node2, 0);
    mesh.add_node_sharing(node5, 0);
  }

  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(mesh, stk::topology::NODE_RANK, nodes);
  for (size_t n=0; n<nodes.size(); ++n)
  {
    stk::mesh::Entity node = nodes[n];
    int node_id = mesh.identifier(node);

    double* coords = stk::mesh::field_data(coordField, node);
    coords[0] = xCoords[node_id-1];
    coords[1] = yCoords[node_id-1];
  }

  //create 1 element per processor
  stk::mesh::EntityId elemId = parallel_rank + 1;
  stk::mesh::Entity element = stk::mesh::declare_element(mesh, block_1, elemId, connected_nodes );

  mesh.modification_end();

  stk::mesh::field_fill(1.0, oneField);

  EXPECT_TRUE(mesh.is_valid(node2));
  EXPECT_TRUE(mesh.is_valid(node5));

  // check that shared nodes are members of block_1
  EXPECT_TRUE(mesh.bucket(node2).member(block_1));
  EXPECT_TRUE(mesh.bucket(node5).member(block_1));

  // check that all nodes of block_1 have the correct value
  std::vector<stk::mesh::Entity> block_1_nodes;
  stk::mesh::get_selected_entities(stk::mesh::Selector(block_1), mesh.buckets( stk::topology::NODE_RANK ), block_1_nodes);
  for(size_t n=0; n<block_1_nodes.size(); ++n)
  {
    double* data_ptr = stk::mesh::field_data(oneField, block_1_nodes[n]);
    EXPECT_TRUE(data_ptr != NULL);
    double value = (NULL == data_ptr) ? 0.0 : *data_ptr;
    EXPECT_DOUBLE_EQ(1.0, value);
  }

  //
  // now switch the element on proc0 to block_2
  //
  mesh.modification_begin();

  if (0 == parallel_rank)
  {
    stk::mesh::PartVector add_parts(1, &block_2);
    stk::mesh::PartVector remove_parts(1, &block_1);

    element = mesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
    mesh.change_entity_parts(element, add_parts, remove_parts);
  }

  mesh.modification_end();

  // check that shared nodes are now members of both blocks
  EXPECT_TRUE(mesh.bucket(node2).member(block_1));
  EXPECT_TRUE(mesh.bucket(node5).member(block_1));
  EXPECT_TRUE(mesh.bucket(node2).member(block_2));
  EXPECT_TRUE(mesh.bucket(node5).member(block_2));

  // check that all nodes of block_1 have the correct value
  stk::mesh::get_selected_entities(stk::mesh::Selector(block_1), mesh.buckets( stk::topology::NODE_RANK ), block_1_nodes);
  for(size_t n=0; n<block_1_nodes.size(); ++n)
  {
    double* data_ptr = stk::mesh::field_data(oneField, block_1_nodes[n]);
    EXPECT_TRUE(data_ptr != NULL);
    double value = (NULL == data_ptr) ? 0.0 : *data_ptr;
    EXPECT_DOUBLE_EQ(1.0, value);
  }
}

TEST(BulkData, STK_ParallelPartConsistency_ChangeBlock_WithAura)
{
    Test_STK_ParallelPartConsistency_ChangeBlock(stk::mesh::BulkData::AUTO_AURA);
}

TEST(BulkData, STK_ParallelPartConsistency_ChangeBlock_WithoutAura)
{
    Test_STK_ParallelPartConsistency_ChangeBlock(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(BulkData, STK_Deimprint)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size != 1) return;

  // This test will create a two-element mesh (quad4 elements)
  // in 2 blocks in parallel.  Element 1 is in block_1, element 2
  // is in block_2.  The part block_2 is a subset of a part that
  // has no rank (and is therefore non-imprintable).  This test
  // examines the behavior of the parts on the shared nodes when
  // element 2 is destroyed.

  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData mesh(meta, pm);

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2_superset_with_rank = meta.declare_part("block_2_superset_with_rank", stk::topology::ELEMENT_RANK);
  stk::mesh::Part& block_2_superset_without_rank = meta.declare_part("block_2_superset_without_rank");

  meta.declare_part_subset(block_2_superset_with_rank, block_2);
  meta.declare_part_subset(block_2_superset_without_rank, block_2);

  stk::mesh::Selector all_nodes = meta.universal_part();

  meta.commit();
  mesh.modification_begin();

  const size_t nodesPerElem = 4;

  stk::mesh::EntityIdVector elem_nodes[] {
      {1, 2, 5, 6},
      {2, 3, 4, 5}
  };

  const size_t numElem = 2;

  // Create nodes
  for (size_t e=0; e<numElem; ++e)
  {
    for(size_t n=0; n<nodesPerElem; ++n)
    {
      stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];

      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
      if (!mesh.is_valid(node))
      {
        node = mesh.declare_node(nodeGlobalId);
      }
    }
  }

  stk::mesh::Entity node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  stk::mesh::Entity node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

  // Create elements
  for (size_t e=0; e<numElem; ++e)
  {
    stk::mesh::Part& block = (e == 0) ? block_1 : block_2;
    stk::mesh::declare_element(mesh, block, e+1, elem_nodes[e]);
  }

  stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);

  mesh.modification_end();

  EXPECT_TRUE(mesh.is_valid(node2));
  EXPECT_TRUE(mesh.is_valid(node5));
  EXPECT_TRUE(mesh.is_valid(element1));
  EXPECT_TRUE(mesh.is_valid(element2));

  // These two flags should definitely be the same. Otherwise, creation and deletion
  // are not symmetric, and parts of lower rank entities will be corrupted by mesh modification.
  // The code might be simpler and performance might be better if both were false.
  const bool stk_induces_unranked_supersets = false;
  const bool stk_deinduces_unranked_supersets = false;

  // check block membership of shared nodes
  EXPECT_TRUE(mesh.bucket(node2).member(block_1));
  EXPECT_TRUE(mesh.bucket(node5).member(block_1));
  EXPECT_TRUE(mesh.bucket(node2).member(block_2));
  EXPECT_TRUE(mesh.bucket(node5).member(block_2));
  EXPECT_TRUE(mesh.bucket(node2).member(block_2_superset_with_rank));
  EXPECT_TRUE(mesh.bucket(node5).member(block_2_superset_with_rank));
  if (stk_induces_unranked_supersets)
  {
    EXPECT_TRUE(mesh.bucket(node2).member(block_2_superset_without_rank));
    EXPECT_TRUE(mesh.bucket(node5).member(block_2_superset_without_rank));
  }
  else
  {
    EXPECT_FALSE(mesh.bucket(node2).member(block_2_superset_without_rank));
    EXPECT_FALSE(mesh.bucket(node5).member(block_2_superset_without_rank));
  }

  //
  // now delete element 2
  //
  mesh.modification_begin();

  mesh.destroy_entity(element2);

  mesh.modification_end();

  // check block membership of shared nodes
  EXPECT_TRUE(mesh.bucket(node2).member(block_1));
  EXPECT_TRUE(mesh.bucket(node5).member(block_1));
  EXPECT_FALSE(mesh.bucket(node2).member(block_2));
  EXPECT_FALSE(mesh.bucket(node5).member(block_2));
  EXPECT_FALSE(mesh.bucket(node2).member(block_2_superset_with_rank));
  EXPECT_FALSE(mesh.bucket(node5).member(block_2_superset_with_rank));
  if (stk_induces_unranked_supersets && !stk_deinduces_unranked_supersets)
  {
    EXPECT_TRUE(mesh.bucket(node2).member(block_2_superset_without_rank));
    EXPECT_TRUE(mesh.bucket(node5).member(block_2_superset_without_rank));
  }
  else
  {
    EXPECT_FALSE(mesh.bucket(node2).member(block_2_superset_without_rank));
    EXPECT_FALSE(mesh.bucket(node5).member(block_2_superset_without_rank));
  }
}

TEST(BulkData, ChangeAuraElementPart)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size != 4) return;

  // This test will create a four element mesh (quad4 elements)
  // in parallel.

  /*  Mesh
   *    7---8---9  P0 owns nodes 1,2,4,5; P, elem 1
   *    | 3 | 4 |  P1 : 3,6, elem 2
   *    4---5---6  P2 : 7,8, elem 3
   *    | 1 | 2 |  P3 : 9,   elem 4
   *    1---2---3
   */

  unsigned spatialDim = 2;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData mesh(meta, pm);

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);

  stk::mesh::Selector all_nodes = meta.universal_part();

  meta.commit();
  mesh.modification_begin();

  const size_t nodesPerElem = 4;

  stk::mesh::EntityIdVector elem_nodes[] {
      {1, 2, 5, 4},
      {2, 3, 6, 5},
      {4, 5, 8, 7},
      {5, 6, 9, 8}
  };

  int node_sharing[9][4] = { {0}, {0,1}, {1}, {0,2}, {0,1,2,3}, {1,3}, {2}, {2,3}, {3} };
  int num_node_sharing[] = {1, 2, 1, 2, 4, 2, 1, 2, 1};

  // Create nodes
  size_t e = parallel_rank;
  {
    for(size_t n=0; n<nodesPerElem; ++n)
    {
      stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];

      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
      if (!mesh.is_valid(node))
      {
        node = mesh.declare_node(nodeGlobalId);
      }
    }
  }

  // Create elements
  stk::mesh::Entity local_elem = stk::mesh::declare_element(mesh, block_1, e+1, elem_nodes[e]);

  // declare node sharing
  for(size_t n=0; n<nodesPerElem; ++n)
  {
    stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];
    stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
    const int node_index = nodeGlobalId-1;
    for (int share=0; share<num_node_sharing[node_index]; ++share)
    {
      const int share_proc = node_sharing[node_index][share];
      if (parallel_rank != share_proc) mesh.add_node_sharing(node, share_proc);
    }
  }

  mesh.modification_end();

  // All nodes appear everywhere via aura
  stk::mesh::Entity node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);
  stk::mesh::Entity node6 = mesh.get_entity(stk::topology::NODE_RANK, 6);
  EXPECT_TRUE(mesh.is_valid(node5));
  EXPECT_TRUE(mesh.is_valid(node6));

  mesh.modification_begin();
  if (parallel_rank == 1 || parallel_rank == 3)
  {
    stk::mesh::PartVector add_parts(1, &block_2);
    stk::mesh::PartVector remove_parts(1, &block_1);
    mesh.change_entity_parts(local_elem, add_parts, remove_parts);
  }
  mesh.modification_end();

  stk::mesh::Entity elem2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
  stk::mesh::Entity elem4 = mesh.get_entity(stk::topology::ELEMENT_RANK, 4);
  EXPECT_TRUE(mesh.is_valid(elem2));
  EXPECT_TRUE(mesh.is_valid(elem4));

  const bool expect_consistent_parts_on_aura_entities = true;

  if (expect_consistent_parts_on_aura_entities)
  {
    EXPECT_FALSE(mesh.bucket(elem2).member(block_1));
    EXPECT_TRUE(mesh.bucket(elem2).member(block_2));
    EXPECT_FALSE(mesh.bucket(elem4).member(block_1));
    EXPECT_TRUE(mesh.bucket(elem4).member(block_2));

    EXPECT_TRUE(mesh.bucket(node5).member(block_1));
    EXPECT_TRUE(mesh.bucket(node5).member(block_2));
    EXPECT_FALSE(mesh.bucket(node6).member(block_1));
    EXPECT_TRUE(mesh.bucket(node6).member(block_2));
  }
}

TEST(BulkData, generate_new_ids)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    std::ostringstream os;
    os << "generated:10x10x" << psize;
    const std::string generatedMeshSpec = os.str();
    stkMeshIoBroker.add_mesh_database(generatedMeshSpec, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, counts);

    std::vector<stk::mesh::EntityId> requestedIds;
    size_t numIdsNeeded = counts[stk::topology::NODE_RANK] / 10 + 1;
    ASSERT_TRUE(numIdsNeeded>0);
    stkMeshBulkData.generate_new_ids(stk::topology::NODE_RANK, numIdsNeeded, requestedIds);

    std::sort(requestedIds.begin(), requestedIds.end());
    std::vector<stk::mesh::EntityId>::iterator iter = std::unique(requestedIds.begin(), requestedIds.end());

    bool ids_better_be_unique_on_this_proc = ( iter == requestedIds.end() );
    ASSERT_TRUE(ids_better_be_unique_on_this_proc);

    stk::CommSparse comm(stkMeshBulkData.parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        for(int i = 0; i < stkMeshBulkData.parallel_size(); ++i)
        {
            if(i != stkMeshBulkData.parallel_rank())
            {
                for(size_t j = 0; j < requestedIds.size(); ++j)
                {
                    comm.send_buffer(i).pack<stk::mesh::EntityId>(requestedIds[j]);
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int i = 0; i < stkMeshBulkData.parallel_size(); ++i)
    {
        if(i != stkMeshBulkData.parallel_rank())
        {
            while(comm.recv_buffer(i).remaining())
            {
                stk::mesh::EntityId entity_id;
                comm.recv_buffer(i).unpack<stk::mesh::EntityId>(entity_id);
                bool is_other_procs_id_on_this_proc = std::binary_search(requestedIds.begin(), requestedIds.end(), entity_id);
                ASSERT_FALSE(is_other_procs_id_on_this_proc);
            }
        }
    }

}

TEST(BulkData, test_generate_new_entities)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);

    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    std::ostringstream os;
    os << "generated:10x10x" << psize;
    const std::string generatedMeshSpec = os.str();
    stkMeshIoBroker.add_mesh_database(generatedMeshSpec, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, counts);

    std::vector<stk::mesh::Entity> requestedEntities;
    size_t numIdsNeeded = counts[stk::topology::NODE_RANK] / 10 + 1;
    ASSERT_TRUE(numIdsNeeded>0);

    std::vector<size_t> requests(stkMeshBulkData.mesh_meta_data().entity_rank_count(),0);
    requests[0] = numIdsNeeded;

    stkMeshBulkData.modification_begin();
    stkMeshBulkData.generate_new_entities(requests, requestedEntities);
    stkMeshBulkData.modification_end();

    stk::CommSparse comm(stkMeshBulkData.parallel());

    for(int phase = 0; phase < 2; ++phase)
    {
        for(int i = 0; i < stkMeshBulkData.parallel_size(); ++i)
        {
            if(i != stkMeshBulkData.parallel_rank())
            {
                for(size_t j = 0; j < requestedEntities.size(); ++j)
                {
                    comm.send_buffer(i).pack<stk::mesh::EntityKey>(stkMeshBulkData.entity_key(requestedEntities[j]));
                }
            }
        }

        if(phase == 0)
        {
            comm.allocate_buffers();
        }
        else
        {
            comm.communicate();
        }
    }

    for(int i = 0; i < stkMeshBulkData.parallel_size(); ++i)
    {
        if(i != stkMeshBulkData.parallel_rank())
        {
            while(comm.recv_buffer(i).remaining())
            {
                stk::mesh::EntityKey key;
                comm.recv_buffer(i).unpack<stk::mesh::EntityKey>(key);
                ASSERT_FALSE(stkMeshBulkData.is_valid(stkMeshBulkData.get_entity(key)));
            }
        }
    }

}

TEST(BulkData, test_destroy_ghosted_entity_then_create_locally_owned_entity_with_same_identifier)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int psize = stk::parallel_machine_size(communicator);

    if ( psize == 2 )
    {
        stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
        std::ostringstream os;
        os << "generated:2x2x2";
        const std::string generatedMeshSpec = os.str();
        stkMeshIoBroker.add_mesh_database(generatedMeshSpec, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

        stk::mesh::EntityKey key(stk::topology::ELEMENT_RANK,3);

        stkMeshBulkData.modification_begin();
        stkMeshBulkData.destroy_entity(stkMeshBulkData.get_entity(key), false);
        stkMeshBulkData.modification_end();

        stk::mesh::Part *part = stkMeshBulkData.mesh_meta_data().get_part("block_1");

        stkMeshBulkData.modification_begin();
        if ( stkMeshBulkData.parallel_rank() == 1 )
        {
            stk::mesh::Entity element = stkMeshBulkData.declare_element(3, stk::mesh::ConstPartVector{part});
            stk::mesh::EntityId nodes[] = { 100, 101, 102, 103, 104, 105, 106, 107 };
            const int num_nodes = 8;
            for (int i=0;i<num_nodes;++i)
            {
                stk::mesh::Entity node = stkMeshBulkData.declare_node(nodes[i]);
                stkMeshBulkData.declare_relation(element, node, i);
            }
        }
        stkMeshBulkData.modification_end();

        stkMeshBulkData.modification_begin();
        if ( stkMeshBulkData.parallel_rank() == 1 )
        {
            EXPECT_NO_THROW(stkMeshBulkData.destroy_entity(stkMeshBulkData.get_entity(key), false));
        }
        stkMeshBulkData.modification_end();

    }
}

TEST(FaceCreation, test_face_creation_2Hexes_2procs)
{
    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numProcs==2)
    {
        stk::mesh::MetaData meta(3);
        stk::unit_test_util::BulkDataFaceSharingTester mesh(meta, MPI_COMM_WORLD);

        const std::string generatedMeshSpec = "generated:1x1x2";
        stk::io::fill_mesh(generatedMeshSpec, mesh);

        int procId = stk::parallel_machine_rank(MPI_COMM_WORLD);

        unsigned elem_id = procId+1;

        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);

        unsigned face_node_ids[] = { 5, 6, 8, 7 };
        size_t num_nodes_on_entity = 4;
        stk::mesh::EntityVector nodes(num_nodes_on_entity);
        if (procId==0)
        {
            for(size_t n = 0; n < num_nodes_on_entity; ++n)
            {
                nodes[n] = mesh.get_entity(stk::topology::NODE_RANK, face_node_ids[n]);
            }
        }
        else
        {
            for(size_t n = 0; n < num_nodes_on_entity; ++n)
            {
                unsigned index = num_nodes_on_entity - n - 1;
                nodes[n] = mesh.get_entity(stk::topology::NODE_RANK, face_node_ids[index]);
            }
        }

        mesh.modification_begin();

        stk::mesh::Entity side = stk::unit_test_util::declare_element_side_with_nodes(mesh, elem, nodes, 1+procId, meta.get_topology_root_part(stk::topology::QUAD_4));

        EXPECT_TRUE(mesh.is_valid(side));

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(mesh, counts);
        EXPECT_EQ(2u, counts[stk::topology::FACE_RANK]);

        std::vector<stk::mesh::shared_entity_type> potentially_shared_sides;
        mesh.my_markEntitiesForResolvingSharingInfoUsingNodes(stk::topology::FACE_RANK, potentially_shared_sides);

        ASSERT_EQ(1u, potentially_shared_sides.size());

        std::sort(potentially_shared_sides.begin(), potentially_shared_sides.end());

        EXPECT_EQ(side, potentially_shared_sides[0].entity);

        std::vector<std::vector<stk::mesh::shared_entity_type> > shared_entities_by_proc(mesh.parallel_size());
        mesh.my_fillSharedEntities(potentially_shared_sides, shared_entities_by_proc);

        int otherProc = 1 - procId;
        EXPECT_TRUE(shared_entities_by_proc[procId].empty());
        EXPECT_EQ(1u, shared_entities_by_proc[otherProc].size());

        stk::CommSparse comm(mesh.parallel());
        stk::mesh::impl::communicate_shared_entity_info(mesh, comm, shared_entities_by_proc);
        mesh.my_unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(comm, potentially_shared_sides);

        EXPECT_TRUE(mesh.my_internal_is_entity_marked(side) == stk::mesh::BulkData::IS_SHARED);

        mesh.change_entity_key_and_update_sharing_info(potentially_shared_sides);

        mesh.my_modification_end_for_entity_creation({stk::topology::FACE_RANK});

        stk::mesh::comm_mesh_counts(mesh, counts);
        EXPECT_EQ(1u, counts[stk::topology::FACE_RANK]);
    }
}

//
TEST(BulkData, test_parallel_entity_sharing)
{

    stk::mesh::Entity entity;
    stk::mesh::EntityKey quad(stk::mesh::EntityKey(stk::topology::ELEM_RANK, 1));

    entity.set_local_offset(1);
    size_t num_nodes_on_entity = 4;
    std::vector<stk::mesh::EntityKey> keys;
    keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 1));
    keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 2));
    keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 3));
    keys.push_back(stk::mesh::EntityKey(stk::topology::NODE_RANK, 4));

    stk::topology topo = stk::topology::QUAD_4;

    stk::mesh::shared_entity_type sentity(quad, entity, topo);
    sentity.nodes.resize(num_nodes_on_entity);
    for(size_t n = 0; n < num_nodes_on_entity; ++n)
    {
        sentity.nodes[n]=keys[n];
    }

    std::vector<stk::mesh::shared_entity_type> shared_entity_map;
    shared_entity_map.push_back(sentity);

    stk::mesh::shared_entity_type entity_from_other_proc(quad, entity, topo);

    entity_from_other_proc.nodes.resize(num_nodes_on_entity);
    for(size_t n = 0; n < num_nodes_on_entity; ++n)
    {
        int index = num_nodes_on_entity - n - 1;
        entity_from_other_proc.nodes[n]=keys[index];
    }

    int matching_index = stk::unit_test_util::does_entity_exist_in_list(shared_entity_map, entity_from_other_proc);
    EXPECT_TRUE(matching_index >= 0);
}

TEST(BulkData, makeElementWithConflictingTopologies)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size(pm);

  if(p_size != 1)
  {
    return;
  }
  const int spatial_dimension = 2;
  stk::mesh::MetaData meta(spatial_dimension);
  stk::mesh::BulkData mesh(meta, pm);

  stk::mesh::EntityId element_ids[1] = {1};
  stk::mesh::EntityIdVector elem_node_ids[] { {1, 2, 3, 4} };

  stk::mesh::Part * quad_part = &meta.declare_part_with_topology("quad_part", stk::topology::QUAD_4_2D);
  stk::mesh::Part * tri_part  = &meta.declare_part_with_topology( "tri_part", stk::topology::TRI_3_2D);
  meta.commit();

  stk::mesh::PartVector parts;
  parts.push_back(quad_part);
  parts.push_back(tri_part);

  mesh.modification_begin();

  EXPECT_THROW(stk::mesh::declare_element(mesh, parts, element_ids[0], elem_node_ids[0]), std::runtime_error);

  mesh.modification_end();
}

TEST( BulkData, AddSharedNodesInTwoSteps)
{


    stk::ParallelMachine pm = MPI_COMM_WORLD;
    unsigned p_size = stk::parallel_machine_size(pm);
    unsigned p_rank = stk::parallel_machine_rank(pm);

    if(p_size != 3u)
    {
        return;
    }

    int nodeId = 1;
    const unsigned spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData mesh(meta, pm);

    meta.commit();


    mesh.modification_begin();
    Entity node;
    if (0 == p_rank || 1 == p_rank) {
        node = mesh.declare_node(nodeId);
    }
    if (0 == p_rank) {
        mesh.add_node_sharing(node, 1);
    }
    if (1 == p_rank) {
        mesh.add_node_sharing(node, 0);
    }
    mesh.modification_end();



    mesh.modification_begin();
    if (2 == p_rank) {
        node = mesh.declare_node(nodeId);
    }
    if (0 == p_rank) {
        mesh.add_node_sharing(node, 2);
    }
    if (1 == p_rank) {
        mesh.add_node_sharing(node, 2);
    }
    if (2 == p_rank) {
        mesh.add_node_sharing(node, 0);
        mesh.add_node_sharing(node, 1);
    }

    //    EXPECT_THROW(mesh.modification_end(), std::logic_error);
    //this only throws on processor 2, but not in a parallel consistent way
    //this is because you apparently can't create the node on a new processor where it didn't exist before

}

TEST(ChangeEntityId, test_throw_on_shared_node)
{
    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numProcs==2)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);

        const std::string generatedMeshSpec = "generated:1x1x2";
        stk::io::fill_mesh(generatedMeshSpec, mesh);

        stk::mesh::Entity sharedNode5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

        EXPECT_TRUE(mesh.bucket(sharedNode5).shared());

        mesh.modification_begin();

        EXPECT_THROW(mesh.change_entity_id(99, sharedNode5), std::logic_error);

        mesh.modification_end();
    }
}

TEST(AmbiguousTopology, hexRedefinedAsShell)
{
    int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (numProcs==1)
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);

        const std::string generatedMeshSpec = "generated:1x1x1";
        stk::io::fill_mesh(generatedMeshSpec, mesh);

        stk::mesh::Part& shellPart = meta.get_topology_root_part(stk::topology::SHELL_QUAD_4);
        stk::mesh::EntityId elemId = 1;
        mesh.modification_begin();
        ASSERT_THROW(mesh.declare_element(elemId, stk::mesh::ConstPartVector{&shellPart}), std::runtime_error);
    }
}

}// empty namespace

