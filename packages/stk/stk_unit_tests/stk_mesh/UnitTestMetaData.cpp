// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stddef.h>                     // for NULL, size_t
#include <exception>                    // for exception
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <string>                       // for string, operator==, etc
#include <vector>                       // for vector
#include "Shards_CellTopologyData.h"    // for CellTopologyData
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/NamedPair.hpp"
#include "unit_tests/BulkDataTester.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/GetEntities.hpp>


namespace stk { namespace mesh { class Part; } }






using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using std::cout;
using std::endl;

//----------------------------------------------------------------------

namespace {

TEST( UnitTestTopologyPart, noNewPartsWithTopologyAfterCommit )
{
  //Test functions in MetaData.cpp
  const int spatial_dimension = 3;
  MetaData uncommited_metadata(spatial_dimension);
  MetaData commited_metadata(spatial_dimension);

  commited_metadata.commit();

  EXPECT_THROW(commited_metadata.declare_part_with_topology( std::string("a") , stk::topology::TRI_3  ), std::logic_error);

  EXPECT_NO_THROW(uncommited_metadata.declare_part_with_topology( std::string("a") , stk::topology::TRI_3 ));
  uncommited_metadata.commit();
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

TEST( UnitTestMetaData, testEntityRepository )
{
  static const size_t spatial_dimension = 3;

  //Test Entity repository - covering EntityRepository.cpp/hpp
  stk::mesh::MetaData meta ( spatial_dimension );
  stk::mesh::Part & part = meta.declare_part("another part");
  stk::mesh::Part & hex_part = meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);

  meta.commit();

  stk::mesh::unit_test::BulkDataTester bulk ( meta , MPI_COMM_WORLD );
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
    node = bulk.declare_entity( stk::topology::NODE_RANK , new_id+1 , add_part );
    nodes.push_back(node);
  }

  int new_id = size * (++id_base) + rank;
  stk::mesh::Entity elem  = bulk.declare_entity( stk::topology::ELEMENT_RANK , new_id+1 , elem_parts );

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
  ASSERT_TRUE(bulk.my_entity_comm_map_insert(elem, comm_info));
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
  //MetaData constructor fails because there are no entity types:
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
  const CellTopologyData * singleton = NULL;
  const int spatial_dimension = 3;
  MetaData metadata(spatial_dimension);
  Part &pa = metadata.declare_part( std::string("a") , stk::topology::NODE_RANK );
  metadata.declare_attribute_no_delete( pa, singleton);
  metadata.commit();
}

TEST(UnitTestMetaData, set_mesh_bulk_data )
{
  const int spatial_dimension = 3;
  MetaData meta(spatial_dimension);
  BulkData* bulk1 = new BulkData(meta, MPI_COMM_WORLD);
  ASSERT_THROW(BulkData bulk2(meta, MPI_COMM_WORLD), std::logic_error);

  //But if we first clear the original BulkData, we should be able to
  //add another one with the same MetaData.
  delete bulk1;
  BulkData bulk2(meta, MPI_COMM_WORLD);
  meta.set_mesh_bulk_data(&bulk2);
  ASSERT_TRUE(&meta.mesh_bulk_data() == &bulk2);
}

TEST(UnitTestMetaData, superset_of_shared_part)
{
    stk::ParallelMachine communicator = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(communicator);
    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
        const std::string generatedMeshSpecification = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stk::mesh::MetaData &meta = stkMeshIoBroker.meta_data();

        stk::mesh::Part & mysupername = meta.declare_part("my_superset_part_shared");
        meta.declare_part_subset(mysupername, meta.globally_shared_part());
        mysupername.entity_membership_is_parallel_consistent(false);

        stk::mesh::Part & mysupernamelocal = meta.declare_part("my_superset_part_local");
        meta.declare_part_subset(mysupernamelocal, meta.locally_owned_part());
        mysupernamelocal.entity_membership_is_parallel_consistent(false);

        stk::mesh::Part & userpart = meta.declare_part("userpartsubsettest");
        stk::mesh::Part & usersuper = meta.declare_part("usersuperset");
        meta.declare_part_subset(usersuper, userpart);

        stkMeshIoBroker.populate_bulk_data();
        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        bool expect_supersets_to_work_with_shared_part = false;

        if (expect_supersets_to_work_with_shared_part) {

            std::cout << "p[" << mesh.parallel_rank() <<"] num nodes stk shared part=" <<
                    stk::mesh::count_selected_entities(meta.globally_shared_part(),
                                                       mesh.buckets(stk::topology::NODE_RANK)) << std::endl;
            std::cout << "p[" << mesh.parallel_rank() << "] num nodes in superset of stk shared part=" <<
                    stk::mesh::count_selected_entities(mysupername,
                                                       mesh.buckets(stk::topology::NODE_RANK)) << std::endl;
            std::cout << "p[" << mesh.parallel_rank() <<"] num nodes stk local part=" <<
                    stk::mesh::count_selected_entities(meta.locally_owned_part(),
                                                       mesh.buckets(stk::topology::NODE_RANK)) << std::endl;
            std::cout << "p[" << mesh.parallel_rank() << "] num nodes in superset of stk local part=" <<
                    stk::mesh::count_selected_entities(mysupernamelocal,
                                                       mesh.buckets(stk::topology::NODE_RANK)) << std::endl;

            EXPECT_EQ(
                    stk::mesh::count_selected_entities(meta.globally_shared_part(),
                                                       mesh.buckets(stk::topology::NODE_RANK)),
                                                       stk::mesh::count_selected_entities(mysupername,
                                                                                          mesh.buckets(stk::topology::NODE_RANK)));
        }

        EXPECT_EQ(stk::mesh::count_selected_entities(meta.locally_owned_part(),
                                                     mesh.buckets(stk::topology::NODE_RANK)),
                  stk::mesh::count_selected_entities(mysupernamelocal,
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

}
//----------------------------------------------------------------------





