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

#include <stk_util/stk_config.h>
#if defined ( STK_HAS_MPI )
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/CreateEdges.hpp"  // for create_edges
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include "gtest/gtest.h"                // for EXPECT_EQ
#include <algorithm>                    // for max, min
#include <gtest/gtest.h>
#include <iostream>                     // for basic_ostream::operator<<
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldParallel.hpp>  // for parallel_max, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <vector>                       // for vector
namespace stk { namespace mesh { class Part; } }

namespace stk { namespace mesh { class FieldBase; } }

namespace {

using namespace stk::mesh;
using stk::mesh::fixtures::HexFixture;

template <typename T>
T do_operation(Operation Op, T lhs, T rhs)
{
  switch(Op) {
  case Operation::SUM:
    return lhs + rhs;
  case Operation::MIN:
    return std::min(lhs, rhs);
  case Operation::MAX:
    return std::max(lhs, rhs);
  default:
    STK_ThrowRequire(false);
    return 0;
  }
}



template <typename FieldVector>
void do_assemble(Operation Op, BulkData & bulk, FieldVector const& field_vector)
{
  switch(Op) {
  case Operation::SUM:
    parallel_sum(bulk, field_vector);
    break;
  case Operation::MIN:
    parallel_min(bulk, field_vector);
    break;
  case Operation::MAX:
    parallel_max(bulk, field_vector);
    break;
  default:
    STK_ThrowRequire(false);
  }
}

typedef Field<double> ScalarField;
typedef Field<double> CartesianField;
typedef Field<long double> LongDoubleField;
typedef Field<int> IntField;
typedef Field<unsigned long> IdField;

template <Operation Op>
void do_parallel_assemble()
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  int p_rank = stk::parallel_machine_rank(pm);
  int p_size = stk::parallel_machine_size(pm);

  if (p_size == 1) {
    return;
  }

  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = p_size;

  HexFixture fixture(pm, NX, NY, NZ);

  MetaData& meta = fixture.m_meta;
  BulkData& bulk = fixture.m_bulk_data;

  ScalarField& universal_scalar_node_field       = meta.declare_field<double>( stk::topology::NODE_RANK,    "universal_scalar_node_field" );
  ScalarField& non_universal_scalar_node_field   = meta.declare_field<double>( stk::topology::NODE_RANK,    "non_universal_scalar_node_field" );
  CartesianField& universal_cartesian_node_field = meta.declare_field<double>( stk::topology::NODE_RANK, "universal_cartesian_node_field" );
  IntField& universal_scalar_int_node_field      = meta.declare_field<int>( stk::topology::NODE_RANK,       "universal_scalar_int_node_field" );
  ScalarField& universal_scalar_edge_field       = meta.declare_field<double>( stk::topology::EDGE_RANK,    "universal_scalar_edge_field" );
  LongDoubleField& universal_scalar_long_double_node_field = meta.declare_field<long double>( stk::topology::NODE_RANK, "universal_scalar_long_double_node_field" );
  IdField& universal_scalar_id_node_field        = meta.declare_field<unsigned long>( stk::topology::NODE_RANK, "universal_scalar_id_node_field" );

  Part& center_part = meta.declare_part("center_part", stk::topology::NODE_RANK);

  put_field_on_mesh(universal_scalar_node_field, meta.universal_part(), nullptr);
  put_field_on_mesh(universal_cartesian_node_field, meta.universal_part(), 3, nullptr);
  put_field_on_mesh(universal_scalar_int_node_field, meta.universal_part(), nullptr);
  put_field_on_mesh(universal_scalar_edge_field, meta.universal_part(), nullptr);
  put_field_on_mesh(non_universal_scalar_node_field, center_part, nullptr);
  put_field_on_mesh(universal_scalar_long_double_node_field, meta.universal_part(), nullptr);
  put_field_on_mesh(universal_scalar_id_node_field, meta.universal_part(), nullptr);

  // TODO - Need a field that gets shared between more than 2 procs.

  meta.commit();

  fixture.generate_mesh();

  // Generate edges

  create_edges(bulk);

  Selector shared_sel = meta.globally_shared_part();

  std::vector<size_t> counts;
  count_entities( shared_sel, bulk, counts );

  int multiplier = (p_rank == 0 || p_rank == p_size - 1) ? 1 : 2;
  EXPECT_EQ( multiplier * (NX + 1) * (NY + 1), counts[0]); // nodes
  EXPECT_EQ( multiplier * 2 * NX * (NY + 1),   counts[1]); // edges

  // Move center nodes into center part
  Entity center_element = fixture.elem(NX / 2, NY / 2, p_rank);

  bulk.modification_begin();

  Entity const* nodes = bulk.begin_nodes(center_element);
  EXPECT_EQ(bulk.num_nodes(center_element), 8u);
  PartVector add_center_part(1, &center_part);
  for (int i = 0; i < 8; ++i) {
    if (bulk.parallel_owner_rank(nodes[i]) == p_rank) {
      bulk.change_entity_parts(nodes[i], add_center_part);
    }
  }

  bulk.modification_end();

  // Fill field data

  BucketVector const& all_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, shared_sel);
  for (size_t b = 0, be = all_node_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_node_buckets[b];
    for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
      Entity node = bucket[n];
      EntityId node_id = bulk.identifier(node);

      int field_id = universal_scalar_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_node_field, node) = p_rank + field_id + node_id;

      field_id = universal_cartesian_node_field.mesh_meta_data_ordinal();
      double* data = stk::mesh::field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        data[d] = p_rank + (field_id*d) + node_id;
      }

      field_id = universal_scalar_int_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_int_node_field, node) = p_rank + field_id + node_id;


      field_id = non_universal_scalar_node_field.mesh_meta_data_ordinal();
      if (bucket.member(center_part)) {
        *field_data(non_universal_scalar_node_field, node) = p_rank + field_id + node_id;
      }

      field_id = universal_scalar_long_double_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_long_double_node_field, node) = p_rank + field_id + node_id;

      field_id = universal_scalar_id_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_id_node_field, node) = p_rank + field_id + node_id;
    }
  }

  BucketVector const& all_edge_buckets = bulk.get_buckets(stk::topology::EDGE_RANK, shared_sel);
  for (size_t b = 0, be = all_edge_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_edge_buckets[b];
    for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
      Entity edge = bucket[e];
      EntityId edge_id = bulk.identifier(edge);

      int field_id = universal_scalar_edge_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_edge_field, edge) = p_rank + field_id + edge_id;
    }
  }

  std::vector<const FieldBase*> double_field_vector;
  double_field_vector.push_back(&universal_scalar_node_field);
  double_field_vector.push_back(&non_universal_scalar_node_field);
  double_field_vector.push_back(&universal_cartesian_node_field);
  double_field_vector.push_back(&universal_scalar_edge_field);

  do_assemble(Op, bulk, double_field_vector);

  std::vector<const FieldBase*> int_field_vector(1, &universal_scalar_int_node_field);

  do_assemble(Op, bulk, int_field_vector);

  std::vector<const FieldBase*> long_double_field_vector(1, &universal_scalar_long_double_node_field);

  do_assemble(Op, bulk, long_double_field_vector);

  std::vector<const FieldBase*> id_field_vector(1, &universal_scalar_id_node_field);

  do_assemble(Op, bulk, id_field_vector);

  // Check field values

  for (size_t b = 0, be = all_node_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_node_buckets[b];
    for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
      Entity node = bucket[n];
      EntityId node_id = bulk.identifier(node);

      bool is_left_node = bulk.in_shared(bulk.entity_key(node), p_rank - 1);

      int sharing_rank = is_left_node ? p_rank - 1 : p_rank + 1;
      int field_id = universal_scalar_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( *field_data(universal_scalar_node_field, node),
                 (do_operation<double>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );

      field_id = universal_cartesian_node_field.mesh_meta_data_ordinal();
      double* data = stk::mesh::field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        EXPECT_EQ( data[d], (do_operation<double>(Op, p_rank + field_id*d + node_id, sharing_rank + field_id*d + node_id)) );
      }

      field_id = universal_scalar_int_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( *field_data(universal_scalar_int_node_field, node),
                 (do_operation<int>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );

      field_id = non_universal_scalar_node_field.mesh_meta_data_ordinal();
      if (bucket.member(center_part)) {
        EXPECT_EQ( *field_data(non_universal_scalar_node_field, node),
                   (do_operation<double>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );
      }

      field_id = universal_scalar_long_double_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( *field_data(universal_scalar_long_double_node_field, node),
                 (do_operation<long double>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );

      field_id = universal_scalar_id_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( *field_data(universal_scalar_id_node_field, node),
                 (do_operation<unsigned long>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );
    }
  }

  for (size_t b = 0, be = all_edge_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_edge_buckets[b];
    for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
      Entity edge = bucket[e];
      EntityId edge_id = bulk.identifier(edge);

      bool is_left_edge = bulk.in_shared(bulk.entity_key(edge), p_rank - 1);

      int sharing_rank = is_left_edge ? p_rank - 1 : p_rank + 1;

      int field_id = universal_scalar_edge_field.mesh_meta_data_ordinal();
      EXPECT_EQ( *field_data(universal_scalar_edge_field, edge),
                 (do_operation<double>(Op, p_rank + field_id + edge_id, sharing_rank + field_id + edge_id)) );
    }
  }
}

TEST(FieldParallel, parallel_sum)
{
  do_parallel_assemble<Operation::SUM>();
}

TEST(FieldParallel, parallel_min)
{
  do_parallel_assemble<Operation::MIN>();
}

TEST(FieldParallel, parallel_max)
{
  do_parallel_assemble<Operation::MAX>();
}


template <typename FieldVector>
void do_assemble_including_ghosts(Operation Op, BulkData & bulk, FieldVector const& field_vector)
{
  switch(Op) {
  case Operation::SUM:
    parallel_sum_including_ghosts(bulk, field_vector);
    break;
  case Operation::MIN:
    parallel_min_including_ghosts(bulk, field_vector);
    break;
  case Operation::MAX:
    parallel_max_including_ghosts(bulk, field_vector);
    break;
  default:
    STK_ThrowRequire(false);
  }
}

template <Operation Op>
void do_parallel_assemble_including_ghosts()
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  int p_rank = stk::parallel_machine_rank(pm);
  int p_size = stk::parallel_machine_size(pm);

  if (p_size == 1) {
    return;
  }

  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = p_size;

  const bool autoAuraOn = false;
  HexFixture fixture(pm, NX, NY, NZ, autoAuraOn);

  MetaData& meta = fixture.m_meta;
  BulkData& bulk = fixture.m_bulk_data;

  ScalarField& universal_scalar_node_field       = meta.declare_field<double>( stk::topology::NODE_RANK,    "universal_scalar_node_field" );
  ScalarField& non_universal_scalar_node_field   = meta.declare_field<double>( stk::topology::NODE_RANK,    "non_universal_scalar_node_field" );
  CartesianField& universal_cartesian_node_field = meta.declare_field<double>( stk::topology::NODE_RANK, "universal_cartesian_node_field" );
  IntField& universal_scalar_int_node_field      = meta.declare_field<int>( stk::topology::NODE_RANK,       "universal_scalar_int_node_field" );
  IdField& universal_scalar_id_node_field        = meta.declare_field<unsigned long>( stk::topology::NODE_RANK, "universal_scalar_id_node_field" );

  Part& center_part = meta.declare_part("center_part", stk::topology::NODE_RANK);

  put_field_on_mesh(universal_scalar_node_field, meta.universal_part(), nullptr);
  put_field_on_mesh(universal_cartesian_node_field, meta.universal_part(), 3, nullptr);
  put_field_on_mesh(universal_scalar_int_node_field, meta.universal_part(), nullptr);
  put_field_on_mesh(non_universal_scalar_node_field, center_part, nullptr);
  put_field_on_mesh(universal_scalar_id_node_field, meta.universal_part(), nullptr);

  meta.commit();

  fixture.generate_mesh();

  // Ghost non-shared nodes from the "ends" of the mesh to the opposite side
  std::vector<EntityProc> add_send1;
  std::vector<EntityProc> add_send2;
  Selector non_shared_sel = !meta.globally_shared_part();
  if (p_rank == 0 || p_rank == p_size - 1) {
    const int target_proc = (p_rank == 0) ? (p_size - 1) : 0;  // Opposite ends
    BucketVector const& non_shared_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, non_shared_sel);
    for (size_t b = 0, be = non_shared_node_buckets.size(); b < be; ++b) {
      Bucket& bucket = *non_shared_node_buckets[b];
      for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
        Entity node = bucket[n];
        add_send1.push_back(std::make_pair(node, target_proc));
        add_send2.push_back(std::make_pair(node, target_proc));
      }
    }
  }
  bulk.modification_begin();
  stk::mesh::Ghosting & ghosting1 = bulk.create_ghosting("test ghosting 1");
  stk::mesh::Ghosting & ghosting2 = bulk.create_ghosting("test ghosting 2");
  bulk.change_ghosting(ghosting1, add_send1);
  bulk.change_ghosting(ghosting2, add_send2);
  bulk.modification_end();


  Selector shared_sel = meta.globally_shared_part();

  std::vector<size_t> counts;
  count_entities( shared_sel, bulk, counts );

  int multiplier = (p_rank == 0 || p_rank == p_size - 1) ? 1 : 2;
  EXPECT_EQ( multiplier * (NX + 1) * (NY + 1), counts[0]); // nodes

  // Move center nodes into center part
  Entity center_element = fixture.elem(NX / 2, NY / 2, p_rank);

  bulk.modification_begin();

  Entity const* nodes = bulk.begin_nodes(center_element);
  EXPECT_EQ(bulk.num_nodes(center_element), 8u);
  PartVector add_center_part(1, &center_part);
  for (int i = 0; i < 8; ++i) {
    if (bulk.parallel_owner_rank(nodes[i]) == p_rank) {
      bulk.change_entity_parts(nodes[i], add_center_part);
    }
  }

  bulk.modification_end();

  // Fill field data

  BucketVector const& all_node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, meta.universal_part());
  for (size_t b = 0, be = all_node_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_node_buckets[b];
    for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
      Entity node = bucket[n];
      EntityId node_id = bulk.identifier(node);

      int field_id = universal_scalar_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_node_field, node) = p_rank + field_id + node_id;

      field_id = universal_cartesian_node_field.mesh_meta_data_ordinal();
      double* data = stk::mesh::field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        data[d] = p_rank + (field_id*d) + node_id;
      }

      field_id = universal_scalar_int_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_int_node_field, node) = p_rank + field_id + node_id;


      field_id = non_universal_scalar_node_field.mesh_meta_data_ordinal();
      if (bucket.member(center_part)) {
        *field_data(non_universal_scalar_node_field, node) = p_rank + field_id + node_id;
      }

      field_id = universal_scalar_id_node_field.mesh_meta_data_ordinal();
      *field_data(universal_scalar_id_node_field, node) = p_rank + field_id + node_id;
    }
  }

  std::vector<const FieldBase*> double_field_vector;
  double_field_vector.push_back(&universal_scalar_node_field);
  double_field_vector.push_back(&non_universal_scalar_node_field);
  double_field_vector.push_back(&universal_cartesian_node_field);

  do_assemble_including_ghosts(Op, bulk, double_field_vector);

  std::vector<const FieldBase*> int_field_vector(1, &universal_scalar_int_node_field);

  do_assemble_including_ghosts(Op, bulk, int_field_vector);

  std::vector<const FieldBase*> id_field_vector(1, &universal_scalar_id_node_field);

  do_assemble_including_ghosts(Op, bulk, id_field_vector);

  // Check field values

  for (size_t b = 0, be = all_node_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_node_buckets[b];
    for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
      Entity node = bucket[n];
      EntityId node_id = bulk.identifier(node);

      int other_rank = 0;
      if (bulk.in_shared(bulk.entity_key(node))) {
        bool is_left_node = bulk.in_shared(bulk.entity_key(node), p_rank - 1);
        other_rank = is_left_node ? p_rank - 1 : p_rank + 1;
      }
      else {
        other_rank = (p_rank == 0) ? (p_size - 1) : 0;  // Opposite ends
      }

      int field_id = universal_scalar_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( do_operation<double>(Op, p_rank + field_id + node_id, other_rank + field_id + node_id),
                 *field_data(universal_scalar_node_field, node) );

      field_id = universal_cartesian_node_field.mesh_meta_data_ordinal();
      double* data = stk::mesh::field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        EXPECT_EQ( do_operation<double>(Op, p_rank + field_id*d + node_id, other_rank + field_id*d + node_id), data[d] );
      }

      field_id = universal_scalar_int_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( do_operation<int>(Op, p_rank + field_id + node_id, other_rank + field_id + node_id),
                 *field_data(universal_scalar_int_node_field, node) );

      field_id = non_universal_scalar_node_field.mesh_meta_data_ordinal();
      if (bucket.member(center_part)) {
        EXPECT_EQ( do_operation<double>(Op, p_rank + field_id + node_id, other_rank + field_id + node_id),
                   *field_data(non_universal_scalar_node_field, node) );
      }

      field_id = universal_scalar_id_node_field.mesh_meta_data_ordinal();
      EXPECT_EQ( do_operation<unsigned long>(Op, p_rank + field_id + node_id, other_rank + field_id + node_id),
                 *field_data(universal_scalar_id_node_field, node) );
    }
  }
}

TEST(FieldParallel, parallel_sum_including_ghosts)
{
  do_parallel_assemble_including_ghosts<Operation::SUM>();
}

TEST(FieldParallel, parallel_min_including_ghosts)
{
  do_parallel_assemble_including_ghosts<Operation::MIN>();
}

TEST(FieldParallel, parallel_max_including_ghosts)
{
  do_parallel_assemble_including_ghosts<Operation::MAX>();
}

} //namespace <anonymous>
#endif
