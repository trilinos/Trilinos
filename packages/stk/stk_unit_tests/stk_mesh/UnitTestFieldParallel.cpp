// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stk_util/stk_config.h>
#if defined ( STK_HAS_MPI )
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for max, min
#include <iostream>                     // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldParallel.hpp>  // for parallel_max, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for EXPECT_EQ
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian
#include "stk_mesh/base/CreateEdges.hpp"  // for create_edges
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityId, etc
#include "stk_mesh/fixtures/HexFixture.hpp"  // for HexFixture
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequire
namespace stk { namespace mesh { class Part; } }

namespace stk { namespace mesh { class FieldBase; } }

namespace {

using namespace stk::mesh;
using stk::mesh::fixtures::HexFixture;

enum Operation
{
  SUM,
  MIN,
  MAX
};

template <typename T>
T do_operation(Operation Op, T lhs, T rhs)
{
  switch(Op) {
  case SUM:
    return lhs + rhs;
  case MIN:
    return std::min(lhs, rhs);
  case MAX:
    return std::max(lhs, rhs);
  default:
    ThrowRequire(false);
    return 0;
  }
}



template <typename FieldVector>
void do_assemble(Operation Op, BulkData & bulk, FieldVector const& field_vector)
{
  switch(Op) {
  case SUM:
    parallel_sum(bulk, field_vector);
    break;
  case MIN:
    parallel_min(bulk, field_vector);
    break;
  case MAX:
    parallel_max(bulk, field_vector);
    break;
  default:
    ThrowRequire(false);
  }
}

typedef Field<double> ScalarField;
typedef Field<double, Cartesian> CartesianField;
typedef Field<int> IntField;

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

  ScalarField& universal_scalar_node_field       = meta.declare_field< ScalarField >( stk::topology::NODE_RANK,    "universal_scalar_node_field" );
  ScalarField& non_universal_scalar_node_field   = meta.declare_field< ScalarField >( stk::topology::NODE_RANK,    "non_universal_scalar_node_field" );
  CartesianField& universal_cartesian_node_field = meta.declare_field< CartesianField >( stk::topology::NODE_RANK, "universal_cartesian_node_field" );
  IntField& universal_scalar_int_node_field      = meta.declare_field< IntField >( stk::topology::NODE_RANK,       "universal_scalar_int_node_field" );
  ScalarField& universal_scalar_edge_field       = meta.declare_field< ScalarField >( stk::topology::EDGE_RANK,    "universal_scalar_edge_field" );

  Part& center_part = meta.declare_part("center_part", stk::topology::NODE_RANK);

  put_field( universal_scalar_node_field , meta.universal_part() );
  put_field( universal_cartesian_node_field , meta.universal_part() );
  put_field( universal_scalar_int_node_field , meta.universal_part() );
  put_field( universal_scalar_edge_field , meta.universal_part() );
  put_field( non_universal_scalar_node_field, center_part);

  // TODO - Need a field that gets shared between more than 2 procs.

  meta.commit();

  fixture.generate_mesh();

  // Generate edges

  create_edges(bulk);

  Selector shared_sel = meta.globally_shared_part();

  std::vector<unsigned> counts;
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

      *field_data(universal_scalar_node_field, node) = p_rank + 1.0 + node_id;

      double* data = field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        data[d] = p_rank + (2.0*d) + node_id;
      }

      *field_data(universal_scalar_int_node_field, node) = p_rank + 3 + node_id;

      if (bucket.member(center_part)) {
        *field_data(non_universal_scalar_node_field, node) = p_rank + 4.0 + node_id;
      }
    }
  }

  BucketVector const& all_edge_buckets = bulk.get_buckets(stk::topology::EDGE_RANK, shared_sel);
  for (size_t b = 0, be = all_edge_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_edge_buckets[b];
    for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
      Entity edge = bucket[e];
      EntityId edge_id = bulk.identifier(edge);

      *field_data(universal_scalar_edge_field, edge) = p_rank + 5.0 + edge_id;
    }
  }

  std::vector<FieldBase*> double_field_vector;
  double_field_vector.push_back(&universal_scalar_node_field);
  double_field_vector.push_back(&non_universal_scalar_node_field);
  double_field_vector.push_back(&universal_cartesian_node_field);
  double_field_vector.push_back(&universal_scalar_edge_field);

  do_assemble(Op, bulk, double_field_vector);

  std::vector<FieldBase*> int_field_vector(1, &universal_scalar_int_node_field);

  do_assemble(Op, bulk, int_field_vector);

  // Check field values

  for (size_t b = 0, be = all_node_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_node_buckets[b];
    for (size_t n = 0, ne = bucket.size(); n < ne; ++n) {
      Entity node = bucket[n];
      EntityId node_id = bulk.identifier(node);

      bool is_left_node = bulk.in_shared(bulk.entity_key(node), p_rank - 1);

      int sharing_rank = is_left_node ? p_rank - 1 : p_rank + 1;
      int field_id = 1;

      EXPECT_EQ( *field_data(universal_scalar_node_field, node),
                 (do_operation<double>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );
      ++field_id;

      double* data = field_data(universal_cartesian_node_field, node);
      for (int d = 0; d < 3; ++d) {
        EXPECT_EQ( data[d], (do_operation<double>(Op, p_rank + field_id*d + node_id, sharing_rank + field_id*d + node_id)) );
      }
      ++field_id;

      EXPECT_EQ( *field_data(universal_scalar_int_node_field, node),
                 (do_operation<int>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );
      ++field_id;

      if (bucket.member(center_part)) {
        EXPECT_EQ( *field_data(non_universal_scalar_node_field, node),
                   (do_operation<double>(Op, p_rank + field_id + node_id, sharing_rank + field_id + node_id)) );
      }
    }
  }

  for (size_t b = 0, be = all_edge_buckets.size(); b < be; ++b) {
    Bucket& bucket = *all_edge_buckets[b];
    for (size_t e = 0, ee = bucket.size(); e < ee; ++e) {
      Entity edge = bucket[e];
      EntityId edge_id = bulk.identifier(edge);

      bool is_left_edge = bulk.in_shared(bulk.entity_key(edge), p_rank - 1);

      int sharing_rank = is_left_edge ? p_rank - 1 : p_rank + 1;
      int field_id = 5;

      EXPECT_EQ( *field_data(universal_scalar_edge_field, edge),
                 (do_operation<double>(Op, p_rank + field_id + edge_id, sharing_rank + field_id + edge_id)) );
    }
  }
}

TEST(FieldParallel, parallel_sum)
{
  do_parallel_assemble<SUM>();
}

TEST(FieldParallel, parallel_min)
{
  do_parallel_assemble<MIN>();
}

TEST(FieldParallel, parallel_max)
{
  do_parallel_assemble<MAX>();
}

} //namespace <anonymous>
#endif
