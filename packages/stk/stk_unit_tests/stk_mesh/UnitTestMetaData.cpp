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

#include <gtest/gtest.h>                // for AssertHelper, TEST, etc
#include <stddef.h>                     // for NULL, size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error, logic_error
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <string>                       // for string, operator==
#include <vector>                       // for vector
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/NamedPair.hpp"
#include <stk_unit_test_utils/stk_mesh_fixtures/TestHexFixture.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>

namespace stk { namespace mesh { class Part; } }

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::MeshBuilder;

namespace {

TEST( UnitTestRootTopology, newPartsWithTopologyAfterCommit )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  MetaData uncommitted_metadata(spatial_dimension);
  MetaData committed_metadata(spatial_dimension);

  committed_metadata.commit();

  EXPECT_NO_THROW(committed_metadata.declare_part_with_topology( std::string("a") , stk::topology::TRI_3  ));

  EXPECT_NO_THROW(uncommitted_metadata.declare_part_with_topology( std::string("a") , stk::topology::TRI_3 ));
  uncommitted_metadata.commit();
}

TEST(UnitTestMetaData, declare_ranked_part_without_spatial_dim)
{
  MetaData meta;
  EXPECT_ANY_THROW(meta.declare_part("myPart",stk::topology::NODE_RANK));
}

TEST(UnitTestMetaData, superElemTopoDeclarePartWithTopology)
{
    const int spatial_dimension = 3;
    MetaData meta(spatial_dimension);
    unsigned numNodes = 11;
    stk::topology superTopo = stk::create_superelement_topology(numNodes);
    Part& part = meta.declare_part_with_topology("super-part", superTopo);
    stk::topology partTopo = meta.get_topology(part);
    EXPECT_TRUE(partTopo.is_superelement());
    EXPECT_EQ(numNodes, partTopo.num_nodes());
}

TEST( UnitTestMetaData, testMetaData )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  MetaData metadata_committed(spatial_dimension);
  MetaData metadata_not_committed(spatial_dimension);
  MetaData metadata(spatial_dimension);

  stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;
  Part &pa = metadata.declare_part( std::string("a") , node_rank );
  Part &pb = metadata.declare_part( std::string("b") , node_rank );
  Part &pc = metadata.declare_part( std::string("c") , node_rank );
  Part &pd = metadata.declare_part( std::string("d") , node_rank );
  Part &pe = metadata.declare_part( std::string("e") , node_rank );
  PartVector part_vector;
  metadata_committed.commit();

  //test get_part with part that does not exist
  std::string test_string = "this_part_does_not_exist";
  ASSERT_THROW( metadata_committed.get_part(test_string,"test_throw"),std::runtime_error);

  //test get_part with valid part
  ASSERT_TRUE( metadata.get_part(std::string("a"),"do_not_throw"));



  part_vector.push_back(& pa);
  part_vector.push_back(& pb);
  part_vector.push_back(& pc);
  part_vector.push_back(& pd);

  //Test declare_part_subset
  ASSERT_THROW(  metadata.declare_part_subset( pe, pe), std::runtime_error);

  metadata.commit();
}

TEST( UnitTestMetaData, rankHigherThanDefined )
{
  //Test function entity_rank_name in MetaData.cpp
  const int spatial_dimension = 3;
  const std::vector<std::string> & rank_names = stk::mesh::entity_rank_names();
  MetaData metadata(spatial_dimension, rank_names);

  const std::string& i_name2 =  metadata.entity_rank_name( stk::topology::EDGE_RANK );

  ASSERT_TRUE( i_name2 == rank_names[stk::topology::EDGE_RANK] );

  EntityRank one_rank_higher_than_defined = static_cast<EntityRank>(rank_names.size());

  ASSERT_THROW(
    metadata.entity_rank_name( one_rank_higher_than_defined ),
    std::runtime_error
                        );
}

TEST( UnitTestMetaData, testEntityKeyMapping )
{
  static const size_t spatial_dimension = 3;

  stk::mesh::MetaData meta ( spatial_dimension );
  stk::mesh::Part & part = meta.declare_part("another part");
  stk::mesh::Part & hex_part = meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);

  meta.commit();

  stk::unit_test_util::BulkDataTester bulk ( meta , MPI_COMM_WORLD );
  std::vector<stk::mesh::Part *>  add_part;
  add_part.push_back ( &part );
  std::vector<stk::mesh::Part *> elem_parts;
  elem_parts.push_back( &part );
  elem_parts.push_back( &hex_part );

  int rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
  int size = stk::parallel_machine_size( MPI_COMM_WORLD );
  PartVector tmp(1);

  bulk.modification_begin();

  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::Entity node = stk::mesh::Entity();
  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = size * id_base + rank;
    node = bulk.declare_node(new_id+1, add_part);
    nodes.push_back(node);
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity elem  = bulk.declare_element(new_id+1, elem_parts);

  for (unsigned ord = 0; ord < 8; ++ord)
  {
    bulk.declare_relation(elem, nodes[ord], ord);
  }

  bulk.my_entity_comm_map_clear(bulk.entity_key(elem));

  bulk.my_entity_comm_map_clear_ghosting(bulk.entity_key(elem));

  const stk::mesh::Ghosting & ghost = bulk.aura_ghosting();

  bulk.modification_end();

  ASSERT_FALSE(bulk.my_entity_comm_map_erase(bulk.entity_key(elem), ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  ASSERT_FALSE(bulk.my_entity_comm_map_erase(bulk.entity_key(elem), comm_info));
  ASSERT_TRUE(bulk.my_entity_comm_map_insert(elem, comm_info).second);
}

TEST( UnitTestMetaData, noEntityTypes )
{
  //MetaData constructor fails because there are no entity types:
  std::vector<std::string> wrong_names(1, "foo");
  ASSERT_THROW(
    MetaData metadata(3 /*dim*/, wrong_names),
    std::runtime_error
    );
}
TEST( UnitTestMetaData, declare_part_with_rank )
{
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  metadata.declare_part("foo");
  ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));
  ASSERT_NO_THROW(metadata.declare_part("foo",stk::topology::EDGE_RANK));

  // Should throw because we're trying to change rank
  ASSERT_THROW(metadata.declare_part("foo",stk::topology::FACE_RANK),std::runtime_error);

  // Should not throw since we did not provide rank
  metadata.declare_part("foo");
}

TEST( UnitTestMetaData, declare_attribute_no_delete )
{
  //Coverage of declare_attribute_no_delete in MetaData.hpp
  const int * singleton = NULL;
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  Part &pa = metadata.declare_part( std::string("a") , stk::topology::NODE_RANK );
  metadata.declare_attribute_no_delete( pa, singleton);
  metadata.commit();
}

TEST(UnitTestMetaData, set_mesh_bulk_data )
{
  const int spatial_dimension = 3;
  MeshBuilder builder(MPI_COMM_WORLD);
  std::shared_ptr<MetaData> meta = builder.set_spatial_dimension(spatial_dimension).create_meta_data();

  std::shared_ptr<BulkData> bulk1 = builder.create(meta);
  ASSERT_THROW(builder.create(meta), std::logic_error);

  //But if we first clear the original BulkData, we should be able to
  //add another one with the same MetaData.
  bulk1.reset();

  std::shared_ptr<BulkData> bulk2 = builder.create(meta);
  ASSERT_TRUE(&meta->mesh_bulk_data() == bulk2.get());
}

class TestHexMeta : public stk::mesh::fixtures::TestHexFixture {};

TEST_F(TestHexMeta, superset_of_shared_part)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs == 2)
    {
        stk::mesh::MetaData &meta = get_meta();

        stk::mesh::Part & mySuperPart = meta.declare_part("my_superset_part_shared");
        meta.declare_part_subset(mySuperPart, meta.globally_shared_part());
        mySuperPart.entity_membership_is_parallel_consistent(false);

        stk::mesh::Part & mySuperPartLocal = meta.declare_part("my_superset_part_local");
        meta.declare_part_subset(mySuperPartLocal, meta.locally_owned_part());
        mySuperPartLocal.entity_membership_is_parallel_consistent(false);

        stk::mesh::Part & userpart = meta.declare_part("userpartsubsettest");
        stk::mesh::Part & usersuper = meta.declare_part("usersuperset");
        meta.declare_part_subset(usersuper, userpart);

        setup_mesh(1, 1, 2);
        stk::mesh::BulkData &mesh = get_bulk();

        bool expect_supersets_to_work_with_shared_part = false;

        if (expect_supersets_to_work_with_shared_part) {

            std::cout << "p[" << mesh.parallel_rank() <<"] num nodes stk shared part=" <<
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, meta.globally_shared_part())
                                                       << std::endl;
            std::cout << "p[" << mesh.parallel_rank() << "] num nodes in superset of stk shared part=" <<
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, mySuperPart)
                                                       << std::endl;
            std::cout << "p[" << mesh.parallel_rank() <<"] num nodes stk local part=" <<
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, meta.locally_owned_part())
                                                        << std::endl;
            std::cout << "p[" << mesh.parallel_rank() << "] num nodes in superset of stk local part=" <<
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, mySuperPartLocal)
                                                        << std::endl;

            EXPECT_EQ(
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, meta.globally_shared_part())
                    ,
                    stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, mySuperPart));
        }

        EXPECT_EQ(stk::mesh::count_selected_entities(meta.locally_owned_part(),
                                                     mesh.buckets(stk::topology::NODE_RANK)),
                  stk::mesh::count_selected_entities(mySuperPartLocal,
                                                     mesh.buckets(stk::topology::NODE_RANK)));

        mesh.modification_begin();

        stk::mesh::Entity node7 = mesh.get_entity(stk::topology::NODE_RANK, 7);
        stk::mesh::PartVector addparts(1, &userpart);
        if (mesh.parallel_rank() == 0) {
            mesh.change_entity_parts(node7, addparts);
        }
        mesh.modification_end();

        EXPECT_EQ(1u, stk::mesh::count_selected_entities(userpart,
                                                mesh.buckets(stk::topology::NODE_RANK)));
        EXPECT_EQ(1u, stk::mesh::count_selected_entities(usersuper,
                                                mesh.buckets(stk::topology::NODE_RANK)));

        //now take the subset part off, hope the superset gets taken off too
        mesh.modification_begin();
        stk::mesh::PartVector addnothing;
        stk::mesh::PartVector removeparts(1, &userpart);
        if (mesh.parallel_rank() == 0) {
            mesh.change_entity_parts(node7, addnothing, removeparts);
        }
        mesh.modification_end();

        EXPECT_EQ(0u, stk::mesh::count_selected_entities(userpart,
                                                mesh.buckets(stk::topology::NODE_RANK)));
        bool expect_entities_removed_from_supersets_when_removed_from_all_subsets = false;

        if (expect_entities_removed_from_supersets_when_removed_from_all_subsets) {
            EXPECT_EQ(0u, stk::mesh::count_selected_entities(usersuper,
                                                mesh.buckets(stk::topology::NODE_RANK)));
        }

    }
}

std::shared_ptr<stk::mesh::BulkData> build_mesh(unsigned spatialDim,
                                                stk::ParallelMachine comm)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  return bulk;
}

TEST(UnitTestMetaData, ConsistentSerialDebugCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  meta.declare_part("part_1", stk::topology::NODE_RANK);
  meta.declare_part("part_2", stk::topology::NODE_RANK);

  meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1");
  meta.declare_field<double>(stk::topology::NODE_RANK, "field_2");

  EXPECT_NO_THROW(bulk.modification_begin());
}

TEST(UnitTestMetaData, ConsistentParallelDebugCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  meta.declare_part("part_1", stk::topology::NODE_RANK);
  meta.declare_part("part_2", stk::topology::NODE_RANK);

  meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1");
  meta.declare_field<double>(stk::topology::NODE_RANK, "field_2");

  EXPECT_NO_THROW(bulk.modification_begin());
}

#define STK_EXPECT_THROW_MSG(runit, myProc, msgProc, goldMsg) \
{ \
  bool threw = false; \
  try { \
    runit; \
  } \
  catch(std::exception& e) { \
    threw = true; \
    if (myProc == msgProc) { \
      EXPECT_TRUE(std::string(e.what()).find(goldMsg) != std::string::npos)<<"failed to find '"<<goldMsg<<"' in exception string '"<<e.what()<<"'"; \
    } \
  } \
  EXPECT_TRUE(threw); \
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadPartNameLength)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
  }
  else {
    meta.declare_part("really_long_part_1", stk::topology::NODE_RANK);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Part name (really_long_part_1) does not match Part name (part_1) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadPartNameText)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
  }
  else {
    meta.declare_part("part_2", stk::topology::NODE_RANK);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Part name (part_2) does not match Part name (part_1) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadPartRank)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
  }
  else {
    meta.declare_part("part_1", stk::topology::ELEM_RANK);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Part part_1 rank (ELEMENT_RANK) does not match Part part_1 rank (NODE_RANK) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadPartTopology)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part_with_topology("part_1", stk::topology::HEX_8);
  }
  else {
    meta.declare_part_with_topology("part_1", stk::topology::TET_4);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Part part_1 topology (TETRAHEDRON_4) does not match Part part_1 topology (HEXAHEDRON_8) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadPartSubset)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
    meta.declare_part("part_2", stk::topology::NODE_RANK);
  }
  else {
    stk::mesh::Part & part_1 = meta.declare_part("part_1", stk::topology::NODE_RANK);
    stk::mesh::Part & part_2 = meta.declare_part("part_2", stk::topology::NODE_RANK);
    meta.declare_part_subset(part_1, part_2);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Part part_1 subset ordinals (46 ) does not match Part part_1 subset ordinals () on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadNumberOfParts_RootTooFew)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
  }
  else {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
    meta.declare_part("part_2", stk::topology::NODE_RANK);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Have extra Part (part_2) that does not exist on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadNumberOfParts_RootTooMany)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
    meta.declare_part("part_2", stk::topology::NODE_RANK);
  }
  else {
    meta.declare_part("part_1", stk::topology::NODE_RANK);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Received extra Part (part_2) from root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadFieldNameLength)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  }
  else {
    meta.declare_field<double>(stk::topology::NODE_RANK, "really_long_field_1");
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Field name (really_long_field_1) does not match Field name (field_1) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadFieldNameText)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  }
  else {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_2");
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Field name (field_2) does not match Field name (field_1) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadFieldRank)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  }
  else {
    meta.declare_field<double>(stk::topology::ELEM_RANK, "field_1");
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Field field_1 rank (ELEMENT_RANK) does not match Field field_1 rank (NODE_RANK) on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadFieldNumberOfStates)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1", 1);
  }
  else {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1", 2);
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Field field_1 number of states (2) does not match Field field_1 number of states (1) on root processor\n[p1] Have extra Field (field_1_STKFS_OLD) that does not exist on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadNumberOfFields_RootTooFew)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  }
  else {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_2");
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Have extra Field (field_2) that does not exist on root processor\n");
}

TEST(UnitTestMetaData, InconsistentParallelDebugCheck_BadNumberOfFields_RootTooMany)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) return;

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  MetaData& meta = bulk.mesh_meta_data();

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_2");
  }
  else {
    meta.declare_field<double>(stk::topology::NODE_RANK, "field_1");
  }

  STK_EXPECT_THROW_MSG(bulk.modification_begin(), bulk.parallel_rank(), 1, "[p1] Received extra Field (field_2) from root processor\n");
}

}
//----------------------------------------------------------------------

