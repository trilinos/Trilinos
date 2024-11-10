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

#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/TestHexFixture.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/BoxFixture.hpp"  // for BoxFixture
#include <gtest/gtest.h>
#include <sstream>                      // for ostringstream, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/Bucket.hpp>     // for has_superset, Bucket, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field_on_mesh
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_unit_test_utils/FaceTestingUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc

namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class Part; } }

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;
using stk::mesh::PairIterEntityComm;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::BucketIterator;
using stk::mesh::Selector;
using stk::mesh::Field;
using stk::mesh::FieldBase;
using stk::mesh::put_field_on_mesh;
using stk::mesh::BucketVector;

typedef Field<double> ScalarFieldType;

namespace {

TEST(UnitTestingOfBucket, testBucket)
{
  // Unit test the Part functionality in isolation:

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  // Create a mesh for testing buckets...

  // Create dummy names for entity ranks to be given to MetaData
  std::vector<std::string> entity_names(5);
  for ( size_t i = 0 ; i < 5 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank" << i ;
    entity_names[i] = name.str();
  }

  // Create MetaData, BulkData
  unsigned max_bucket_size = 4;
  stk::mesh::fixtures::BoxFixture fixture(pm, stk::mesh::BulkData::AUTO_AURA, max_bucket_size, entity_names);
  MetaData& meta = fixture.fem_meta();
  BulkData& bulk = fixture.bulk_data();
  // Create two scalar fields, temperature and volume. Put temperature
  // on all the nodes and put volume on all the elements.
  unsigned number_of_states = 4;

  ScalarFieldType & temperature = meta.declare_field<double>(stk::topology::NODE_RANK, "temperature", number_of_states);
  ScalarFieldType & volume =  meta.declare_field<double>(stk::topology::ELEMENT_RANK, "volume", number_of_states);
  Part & universal     = meta.universal_part ();
  put_field_on_mesh ( temperature , universal , nullptr);
  put_field_on_mesh ( volume , universal , nullptr);
  meta.commit();

  // Generate the mesh
  int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  ASSERT_TRUE(bulk.modification_end());

  //  First, test for streaming IO;
  {
    std::string gold1;
    gold1 = "Bucket( EntityRank0 : {UNIVERSAL} {OWNS} {FEM_ROOT_CELL_TOPOLOGY_PART_HEXAHEDRON_8} elem_part )";
    Bucket *b1 = bulk.buckets(stk::topology::NODE_RANK)[0];
    std::stringstream  out1_str;
    out1_str << (*b1);
    bool equal = (gold1 == out1_str.str());
    ASSERT_TRUE(equal);
  }

  // Second, update state of bucket until circular cue is filled
  {
    /* Need to set some data in state, rotate look for it, rotate 3 more times
       and look for it again */
    for ( size_t i = 0 ; i != 10 ; ++i )
      bulk.update_field_data_states ();
  }

  // next, check has_superset (...) and membership functions
  {
    PartVector tmp(2) ;
    tmp[0] = & meta.universal_part();
    tmp[1] = & meta.locally_owned_part();
    ASSERT_TRUE ( has_superset ( *bulk.buckets(stk::topology::NODE_RANK)[0] , tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_any ( tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member_all ( tmp ) );
    ASSERT_TRUE ( bulk.buckets(stk::topology::NODE_RANK)[0]->member ( **meta.get_parts().begin() ) );
  }
}

TEST(UnitTestingOfBucket, bucketSortChangeEntityId)
{
  const unsigned spatialDim=3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::Part& part = meta.declare_part_with_topology("node_part", stk::topology::NODE);
  meta.commit();
  stk::mesh::BulkData& bulk = *bulkPtr;
  if (bulk.parallel_size() > 1) {
    return;
  }
  stk::mesh::EntityId nodeID=1;
  bulk.modification_begin();
  bulk.declare_node(nodeID, stk::mesh::ConstPartVector{&part});
  nodeID=3;
  bulk.declare_node(nodeID, stk::mesh::ConstPartVector{&part});
  nodeID=5;
  bulk.declare_node(nodeID, stk::mesh::ConstPartVector{&part});
  bulk.modification_end();

  const stk::mesh::BucketVector& node_buckets_1 = bulk.get_buckets(stk::topology::NODE_RANK, meta.universal_part());
  size_t expected_num_buckets = 1;
  EXPECT_EQ(expected_num_buckets, node_buckets_1.size());
  size_t expected_bucket_size = 3;
  EXPECT_EQ(expected_bucket_size, node_buckets_1[0]->size());

  stk::mesh::Entity node3 = (*node_buckets_1[0])[1];
  stk::mesh::EntityId node3ID = 3;
  EXPECT_EQ(node3ID, bulk.identifier(node3));

  stk::mesh::Entity node5 = (*node_buckets_1[0])[2];

  bulk.modification_begin();
  stk::mesh::EntityId node2ID = 2;
  bulk.change_entity_id(node2ID, node5);
  bulk.modification_end();

  const stk::mesh::BucketVector& node_buckets_2 = bulk.get_buckets(stk::topology::NODE_RANK, meta.universal_part());

  stk::mesh::Entity node2 = (*node_buckets_2[0])[1];
  EXPECT_EQ(node2ID, bulk.identifier(node2));
}

void test_nodes_and_permutation(stk::mesh::BulkData& bulk, stk::mesh::Entity elem, stk::mesh::Entity side, stk::mesh::EntityVector& nodes)
{
  stk::mesh::EntityRank rank = bulk.entity_rank(side);
  Entity const *rels_itr = bulk.begin_nodes(side);
  unsigned num_nodes = bulk.num_nodes(side);

  for(unsigned i=0;i<num_nodes;++i)
  {
    EXPECT_EQ(nodes[i], rels_itr[i]);
  }

  stk::mesh::Permutation const *perms = bulk.begin_permutations(side, stk::topology::ELEM_RANK);
  std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
      stk::mesh::get_ordinal_and_permutation(bulk, elem, rank, nodes);

  stk::mesh::Permutation gold_permutation = ordinalAndPermutation.second;
  ASSERT_TRUE(gold_permutation!=stk::mesh::INVALID_PERMUTATION);

  unsigned sides_element_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
  unsigned num_elems = bulk.num_elements(side);
  const stk::mesh::Entity *elements = bulk.begin_elements(side);
  for (unsigned i=0;i<num_elems;++i)
  {
    if (elements[i]==elem)
    {
      sides_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
      break;
    }
  }

  EXPECT_EQ(gold_permutation, perms[sides_element_offset]);

  stk::mesh::ConnectivityOrdinal elements_side_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

  unsigned num_sides = bulk.num_connectivity(elem, rank);
  const stk::mesh::Entity *sides = bulk.begin(elem, rank);
  for (unsigned i=0;i<num_sides;++i)
  {
    if (sides[i]==side)
    {
      elements_side_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
      break;
    }
  }

  stk::mesh::Permutation const *perms2 = bulk.begin_permutations(elem, rank);
  EXPECT_EQ(gold_permutation, perms2[elements_side_offset]);
}

bool does_rank_have_permutation(stk::mesh::EntityRank rank)
{
  return rank > stk::topology::NODE_RANK && rank < stk::topology::CONSTRAINT_RANK;
}

class BucketHex : public stk::mesh::fixtures::TestHexFixture {};

TEST_F(BucketHex, testing_valid_permutation_on_various_ranks)
{
  const int num_procs_this_test_works_for = 1;
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == num_procs_this_test_works_for)
  {
    const unsigned num_entity_rank_ranks = 5;
    std::vector< std::string > entity_rank_names(num_entity_rank_ranks);
    entity_rank_names[stk::topology::NODE_RANK] = std::string("NODE");
    entity_rank_names[stk::topology::EDGE_RANK] = std::string("EDGE");
    entity_rank_names[stk::topology::FACE_RANK] = std::string("FACE");
    entity_rank_names[stk::topology::ELEM_RANK] = std::string("ELEMENT");
    entity_rank_names[stk::topology::CONSTRAINT_RANK] = std::string("CONSTRAINT");

    setup_mesh(1, 1, 1, entity_rank_names);
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::BulkData& bulk = get_bulk();

    const unsigned num_nodes_on_one_hex = 8;
    stk::mesh::EntityVector nodes(num_nodes_on_one_hex);
    for (size_t i=0;i<nodes.size();++i)
    {
      stk::mesh::EntityId id = i+1;
      nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, id);
    }

    stk::mesh::EntityId id = 1;
    stk::mesh::EntityVector entities(5);
    entities[stk::topology::NODE_RANK] = bulk.get_entity(stk::topology::NODE_RANK, id);
    entities[stk::topology::ELEM_RANK] = bulk.get_entity(stk::topology::ELEM_RANK, id);

    bulk.modification_begin();

    entities[stk::topology::CONSTRAINT_RANK] = bulk.declare_constraint(id);

    enum node { FIRST_NODE=0, SECOND_NODE=1, THIRD_NODE=2, FOURTH_NODE=3};

    const unsigned num_nodes_on_edge = 2;
    stk::mesh::EntityVector edge_nodes(num_nodes_on_edge);
    edge_nodes[FIRST_NODE]  = nodes[FIRST_NODE];
    edge_nodes[SECOND_NODE] = nodes[SECOND_NODE];

    entities[stk::topology::EDGE_RANK] = stk::unit_test_util::declare_element_to_edge_with_nodes(bulk, entities[stk::topology::ELEM_RANK],
        edge_nodes, id, meta.get_topology_root_part(stk::topology::LINE_2));

    const unsigned num_nodes_on_face = 4;
    stk::mesh::EntityVector face_nodes(num_nodes_on_face);
    face_nodes[FIRST_NODE]  = nodes[FIRST_NODE];
    face_nodes[SECOND_NODE] = nodes[SECOND_NODE];
    face_nodes[THIRD_NODE]  = nodes[FOURTH_NODE];
    face_nodes[FOURTH_NODE] = nodes[THIRD_NODE];

    entities[stk::topology::FACE_RANK] = stk::unit_test_util::declare_element_side_with_nodes(bulk, entities[stk::topology::ELEM_RANK],
        face_nodes, id, meta.get_topology_root_part(stk::topology::QUAD_4));

    bulk.modification_end();

    for (stk::mesh::EntityRank irank=stk::topology::BEGIN_RANK; irank<stk::topology::END_RANK; ++irank)
    {
      for (size_t i=0;i<entities.size();++i)
      {
        ASSERT_TRUE(bulk.is_valid(entities[i]));
        stk::mesh::EntityRank target_rank = irank;
        stk::mesh::EntityRank from_rank = bulk.entity_rank(entities[i]);
        if (does_rank_have_permutation(target_rank) && does_rank_have_permutation(from_rank))
        {
          EXPECT_TRUE(bulk.has_permutation(entities[i], irank));
        }
        else
        {
          EXPECT_FALSE(bulk.has_permutation(entities[i], irank));
        }
      }
    }
  }
}

TEST_F(BucketHex, changing_conn_on_bucket_for_face_to_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_mesh(1, 1, 1);
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::BulkData& bulk = get_bulk();

    unsigned face_node_ids[] = { 5, 6, 8, 7 };
    stk::mesh::EntityVector nodes(4);
    for (size_t i=0;i<nodes.size();++i)
    {
      nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, face_node_ids[i]);
    }

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    bulk.modification_begin();
    stk::mesh::Entity side = stk::unit_test_util::declare_element_side_with_nodes(bulk, elem, nodes, 1, meta.get_topology_root_part(stk::topology::QUAD_4));
    bulk.modification_end();

    test_nodes_and_permutation(bulk, elem, side, nodes);

    unsigned new_face_node_ids[] = { 5, 7, 8, 6 };
    stk::mesh::EntityVector new_nodes(4);
    for (size_t i=0;i<new_nodes.size();++i)
    {
      new_nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, new_face_node_ids[i]);
    }

    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
        stk::mesh::get_ordinal_and_permutation(bulk, elem, stk::topology::FACE_RANK, new_nodes);

    stk::mesh::Permutation new_permutation = ordinalAndPermutation.second;

    unsigned faces_element_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    const stk::mesh::ConnectedEntities elements = bulk.get_connected_entities(side, stk::topology::ELEM_RANK);
    for (unsigned i=0;i<elements.size();++i)
    {
      if (elements[i]==elem)
      {
        faces_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
        break;
      }
    }

    unsigned elements_face_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    const stk::mesh::ConnectedEntities faces = bulk.get_connected_entities(elem, stk::topology::FACE_RANK);
    for (unsigned i=0;i<faces.size();++i)
    {
      if (faces[i]==side)
      {
        elements_face_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
        break;
      }
    }

    ASSERT_TRUE(elements_face_offset!=stk::mesh::INVALID_CONNECTIVITY_ORDINAL);

    stk::mesh::unit_test::BucketTester& bucket_side = static_cast<stk::mesh::unit_test::BucketTester&>(bulk.bucket(side));
    bucket_side.my_change_connected_nodes(bulk.bucket_ordinal(side), &new_nodes[0]);
    bucket_side.my_change_existing_permutation_for_connected_element(bulk.bucket_ordinal(side), faces_element_offset, new_permutation);

    stk::mesh::unit_test::BucketTester& bucket_elem = static_cast<stk::mesh::unit_test::BucketTester&>(bulk.bucket(elem));

    bucket_elem.my_change_existing_permutation_for_connected_face(bulk.bucket_ordinal(elem), elements_face_offset, new_permutation);

    test_nodes_and_permutation(bulk, elem, side, new_nodes);
  }
}

TEST_F(BucketHex, changing_conn_on_bucket_for_edge_to_element)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    setup_mesh(1, 1, 1);
    stk::mesh::MetaData& meta = get_meta();
    stk::mesh::BulkData& bulk = get_bulk();

    unsigned edge_node_ids[] = { 5, 6 };
    stk::mesh::EntityVector nodes(2);
    for (size_t i=0;i<nodes.size();++i)
      nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, edge_node_ids[i]);

    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    bulk.modification_begin();
    stk::mesh::Entity edge = stk::unit_test_util::declare_element_to_edge_with_nodes(bulk, elem, nodes, 1, meta.get_topology_root_part(stk::topology::LINE_2));
    bulk.modification_end();

    test_nodes_and_permutation(bulk, elem, edge, nodes);

    unsigned new_face_node_ids[] = { 6, 5 };
    stk::mesh::EntityVector new_nodes(2);
    for (size_t i=0;i<new_nodes.size();++i)
      new_nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, new_face_node_ids[i]);

    std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation =
        stk::mesh::get_ordinal_and_permutation(bulk, elem, stk::topology::EDGE_RANK, new_nodes);

    stk::mesh::Permutation new_permutation = ordinalAndPermutation.second;

    unsigned edges_element_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    unsigned num_elems = bulk.num_elements(edge);
    const stk::mesh::Entity *elements = bulk.begin_elements(edge);
    for (unsigned i=0;i<num_elems;++i)
    {
      if (elements[i]==elem)
      {
        edges_element_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
        break;
      }
    }

    unsigned elements_edge_offset = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
    unsigned num_edges = bulk.num_edges(elem);
    const stk::mesh::Entity *edges = bulk.begin_edges(elem);
    for (unsigned i=0;i<num_edges;++i)
    {
      if (edges[i]==edge)
      {
        elements_edge_offset = static_cast<stk::mesh::ConnectivityOrdinal>(i);
        break;
      }
    }

    ASSERT_TRUE(elements_edge_offset != stk::mesh::INVALID_CONNECTIVITY_ORDINAL);

    stk::mesh::unit_test::BucketTester& bucket_edge = static_cast<stk::mesh::unit_test::BucketTester&>(bulk.bucket(edge));

    bucket_edge.my_change_connected_nodes(bulk.bucket_ordinal(edge), &new_nodes[0]);
    bucket_edge.my_change_existing_permutation_for_connected_element(bulk.bucket_ordinal(edge), edges_element_offset, new_permutation);

    stk::mesh::unit_test::BucketTester& bucket_elem = static_cast<stk::mesh::unit_test::BucketTester&>(bulk.bucket(elem));

    bucket_elem.my_change_existing_permutation_for_connected_edge(bulk.bucket_ordinal(elem), elements_edge_offset, new_permutation);

    test_nodes_and_permutation(bulk, elem, edge, new_nodes);
  }
}

#ifdef STK_USE_DEVICE_MESH
void do_nonmodifying_debug_check(const stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & coordsField)
{
  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(buckets.size(), 1u);

  buckets[0]->check_size_invariant();

  EXPECT_FALSE(buckets[0]->get_ngp_field_bucket_is_modified(coordsField.mesh_meta_data_ordinal()));
}

void do_modifying_entity_creation(stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & coordsField)
{
  unsigned face_node_ids[] = { 5, 6, 8, 7 };
  stk::mesh::EntityVector nodes(4);
  for (size_t i = 0 ; i < nodes.size(); ++i) {
    nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, face_node_ids[i]);
  }

  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  bulk.modification_begin();
  stk::unit_test_util::declare_element_side_with_nodes(bulk, elem, nodes, 1, meta.get_topology_root_part(stk::topology::QUAD_4));
  bulk.modification_end();

  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(buckets.size(), 2u);
  EXPECT_TRUE(buckets[0]->get_ngp_field_bucket_is_modified(coordsField.mesh_meta_data_ordinal()));
  EXPECT_TRUE(buckets[1]->get_ngp_field_bucket_is_modified(coordsField.mesh_meta_data_ordinal()));
}

TEST_F(BucketHex, checkModifiedStatus)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_mesh(1, 1, 1);
  stk::mesh::MetaData& meta = get_meta();
  stk::mesh::BulkData& bulk = get_bulk();

  const stk::mesh::FieldBase & coordsField = *meta.coordinate_field();
  stk::mesh::get_updated_ngp_field<double>(coordsField);

  do_nonmodifying_debug_check(bulk, coordsField);
  do_modifying_entity_creation(bulk, coordsField);
}
#endif

}
