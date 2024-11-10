/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>      // for parallel_machine_rank
#include <stk_util/parallel/ParallelComm.hpp>  // for parallel_data_exchange...
#include <stk_mesh/base/Bucket.hpp>            // for Bucket
#include <stk_mesh/base/BulkData.hpp>          // for BulkData, BulkData::NO...
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/FieldBase.hpp>         // for field_data, FieldBase
#include <stk_mesh/base/FieldParallel.hpp>     // for parallel_sum, parallel...
#include <stk_mesh/base/GetEntities.hpp>       // for count_selected_entities
#include <stk_mesh/base/Part.hpp>              // for Part
#include <stk_mesh/base/Selector.hpp>          // for Selector, operator|
#include <stk_topology/topology.hpp>           // for topology, topology::NO...
#include <stk_util/environment/CPUTime.hpp>    // for cpu_time
#include "stk_io/FillMesh.hpp"                 // for fill_mesh
#include "stk_mesh/base/Entity.hpp"            // for Entity
#include "stk_mesh/base/Field.hpp"             // for Field
#include "stk_mesh/base/MetaData.hpp"          // for MetaData, put_field_on...
#include "stk_mesh/base/Types.hpp"             // for BucketVector, EntityPr...
#include "stk_util/environment/Env.hpp"        // for parallel_rank, paralle...
#include "stk_util/environment/perf_util.hpp"  // for get_max_hwm_across_procs
#include <cstddef>                             // for size_t
#include <algorithm>                           // for fill
#include <iostream>                            // for operator<<, basic_ostr...
#include <string>                              // for operator<<, char_traits
#include <vector>                              // for vector
namespace stk { namespace mesh { class Ghosting; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace STKperf {

static const int X_DIM = 500; //num_elems_x
static const int Y_DIM = 500; //num_elems_y

static const int NUM_FIELDS = 40;

static const int NUM_ITERS = 40;

stk::mesh::EntityId node_id( unsigned x , unsigned y , unsigned z, unsigned nx, unsigned ny )
{
  return 1 + x + ( nx + 1 ) * ( y + ( ny + 1 ) * z );
}

void do_stk_test(bool with_ghosts=false)
{
  using namespace stk::mesh;

  typedef Field<double> ScalarField;

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);
  if (parallel_size < 3)
  {
      if (parallel_rank == 0)
      {
          std::cerr << "returning, this test only runs on 3 or more processors." << std::endl;
      }
      return;
  }
  int z_dim = parallel_size*2;

  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.

  std::ostringstream oss;
  oss << "generated:" << X_DIM << "x" << Y_DIM << "x" << z_dim;

  stk::mesh::MeshBuilder builder(pm);
  unsigned spatialDim = 3;
  builder.set_spatial_dimension(spatialDim);
  builder.set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;
  if (parallel_rank == 0)
  {
      std::cerr << "Mesh: " << oss.str() << std::endl;
  }

  std::vector<const FieldBase*> fields(NUM_FIELDS);
  for (int i = 0; i < NUM_FIELDS; ++i) {
    std::ostringstream oss2;
    oss2 << "field_" << i;
    FieldBase* field = &meta.declare_field<double>(stk::topology::NODE_RANK, oss2.str());
    fields[i] = field;
    stk::mesh::put_field_on_mesh(*field, meta.universal_part(), nullptr);
  }

  PartVector hex_topo(1, &meta.declare_part_with_topology("hex_part", stk::topology::HEX_8));

  stk::io::fill_mesh(oss.str(),bulk);

  if (parallel_rank == 0) {
    std::cout << "STK mesh constructed" << std::endl;
  }

  int last_proc = parallel_size - 1;
  int next_to_last_proc = last_proc - 1;

  stk::mesh::Selector shared_sel = meta.globally_shared_part();
  stk::mesh::Selector communicated_nodes = shared_sel;
  unsigned num_ghost_nodes = 0;

  if (with_ghosts)
  {
      bulk.modification_begin();
      stk::mesh::Ghosting& node_ghosting = bulk.create_ghosting("node ghosting");
      bulk.modification_end();

      stk::mesh::EntityProcVec nodes_to_ghost;


      stk::mesh::Selector send_ghost_selector;
      if (parallel_rank == 0)
      {
          send_ghost_selector = meta.locally_owned_part() & !shared_sel;
          const stk::mesh::BucketVector& node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, send_ghost_selector);
          for(const stk::mesh::Bucket* bucket : node_buckets)
          {
              for(stk::mesh::Entity node : *bucket)
              {
                  nodes_to_ghost.push_back(stk::mesh::EntityProc(node, last_proc));
                  nodes_to_ghost.push_back(stk::mesh::EntityProc(node, next_to_last_proc));
              }
          }
      }

      bulk.batch_add_to_ghosting(node_ghosting, nodes_to_ghost);

      num_ghost_nodes = nodes_to_ghost.size()/2; //because each node is going to 2 destination procs
      if (parallel_rank == last_proc)
      {
          num_ghost_nodes = stk::mesh::count_selected_entities(bulk.ghosting_part(node_ghosting), bulk.buckets(stk::topology::NODE_RANK));
          EXPECT_GT(num_ghost_nodes, 0u);
      }

      communicated_nodes = send_ghost_selector | bulk.ghosting_part(node_ghosting) | shared_sel;

      if (parallel_rank != 0)
      {
          num_ghost_nodes = stk::mesh::count_selected_entities(bulk.ghosting_part(node_ghosting), bulk.buckets(stk::topology::NODE_RANK));
      }

      if (parallel_rank == 0) {
        std::cout << "Ghosts added to mesh" << std::endl;
      }
  }

  // populate field data
  stk::mesh::BucketVector const& node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, communicated_nodes);
  for (int b = 0, be = node_buckets.size(); b < be; ++b) {
    stk::mesh::Bucket const& bucket = *node_buckets[b];
    for (int i = 0; i < NUM_FIELDS; ++i) {
      const ScalarField& field = dynamic_cast<const ScalarField&>(*fields[i]);
      double* data = stk::mesh::field_data(field, bucket);
      for (int n = 0, ne = bucket.size(); n < ne; ++n) {
        data[n] = static_cast<double>(i+1);
      }
    }
  }

  MPI_Barrier(pm);

  double start_time = stk::cpu_time();

  for (int t = 0; t < NUM_ITERS; ++t) {
      if (with_ghosts)
      {
          stk::mesh::parallel_sum_including_ghosts(bulk, fields);
      }
      else
      {
          stk::mesh::parallel_sum(bulk, fields);
      }
  }

  double stk_sum_time = stk::cpu_time() - start_time;

  double max_time;
  MPI_Reduce(static_cast<void*>(&stk_sum_time), static_cast<void*>(&max_time), 1, MPI_DOUBLE, MPI_MAX, 0 /*root*/, MPI_COMM_WORLD);

  double power2 = std::pow(2,NUM_ITERS);
  double power3 = std::pow(3,NUM_ITERS);
  const double tolerance = 1.e-8;

  // Sanity check
  size_t num_comm_nodes = 0;
  for (int b = 0, be = node_buckets.size(); b < be; ++b) {
    stk::mesh::Bucket const& bucket = *node_buckets[b];
    const bool isShared = bucket.shared();
    num_comm_nodes += bucket.size();
    for (int f = 0; f < NUM_FIELDS; ++f) {
      const ScalarField& field = dynamic_cast<const ScalarField&>(*fields[f]);
      const double* stk_data = stk::mesh::field_data(field, bucket);
      const double expected_shared_value = static_cast<double>(f+1) * power2;
      const double expected_ghosted_value = static_cast<double>(f+1) * power3;
      const double expected = isShared ? expected_shared_value : expected_ghosted_value;
      for (int n = 0, ne = bucket.size(); n < ne; ++n) {
        const double relativeError = std::abs(stk_data[n] - expected) / expected;
        EXPECT_NEAR(0.0, relativeError, tolerance);
      }
    }
  }

  size_t expected_num_shared_nodes = (X_DIM+1) * (Y_DIM+1);
  if (parallel_rank != 0 && parallel_rank != last_proc)
  {
      expected_num_shared_nodes += (X_DIM+1) * (Y_DIM+1);
  }
  if (parallel_rank == 0 || parallel_rank == last_proc || parallel_rank == next_to_last_proc)
  {
      expected_num_shared_nodes += num_ghost_nodes;
  }
  EXPECT_EQ(expected_num_shared_nodes, num_comm_nodes);

  size_t hwm_max = stk::get_max_hwm_across_procs(pm);

  stk::print_stats_for_performance_compare(std::cout, max_time, hwm_max, NUM_ITERS, pm);
}

TEST(ParallelDataExchange, test_nonsym_known_sizes_from_proc0_to_all_other_procs)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int num_procs = stk::parallel_machine_size(comm);
    if (num_procs >= 2) {
        int this_proc = stk::parallel_machine_rank(comm);

        std::vector<int> sendOffsets(num_procs+1);
        std::vector<int> recvOffsets(num_procs+1);
        std::vector<double> send_data(num_procs-1);
        std::vector<double> recv_data(1);
        if (this_proc == 0) {
            sendOffsets[0] = 0;
            sendOffsets[1] = 0;  // Don't send anything to self
            for (int iproc = 1; iproc < num_procs; ++iproc) {
                sendOffsets[iproc+1] = sendOffsets[iproc]+1;
            }

            std::fill(send_data.begin(), send_data.end(), this_proc+1);
        }
        else {
            std::fill(recvOffsets.begin(), recvOffsets.end(), 1);
            recvOffsets[0] = 0;  // Only receive from p0
        }

        stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), send_data.data(),
                                                         recvOffsets.data(), recv_data.data(), comm );

        if (this_proc == 0) {
            for (int iproc = 0; iproc < num_procs; ++iproc) {
                const int recvSizeFromProc = recvOffsets[iproc+1]-recvOffsets[iproc];
                EXPECT_EQ(0, recvSizeFromProc);
            }
        }
        else {
            const int recvSizeFromProc0 = recvOffsets[1]-recvOffsets[0];
            EXPECT_EQ(1, recvSizeFromProc0);
            EXPECT_DOUBLE_EQ(1.0, recv_data[recvOffsets[0]]);
            for (int iproc = 1; iproc < num_procs; ++iproc) {
                const int recvSizeFromProc = recvOffsets[iproc+1]-recvOffsets[iproc];
                EXPECT_EQ(0, recvSizeFromProc);
            }
        }
    }
}

TEST(ParallelDataExchange, test_nonsym_known_sizes_from_all_other_procs_to_proc0)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int num_procs = stk::parallel_machine_size(comm);
    if (num_procs >= 2) {
        int this_proc = stk::parallel_machine_rank(comm);

        std::vector<int> sendOffsets(num_procs+1);
        std::vector<int> recvOffsets(num_procs+1);
        std::vector<double> send_data(1);
        std::vector<double> recv_data(num_procs-1);
        if (this_proc == 0) {
            recvOffsets[0] = 0;
            recvOffsets[1] = 0;  // Don't receive from self
            for (int iproc = 1; iproc < num_procs; ++iproc) {
                recvOffsets[iproc+1] = recvOffsets[iproc]+1;
            }
        }
        else {
            std::fill(sendOffsets.begin(), sendOffsets.end(), 1);
            sendOffsets[0] = 0;  // Only send to p0
            
            send_data[0] = static_cast<double>(this_proc+1);
        }

        stk::parallel_data_exchange_nonsym_known_sizes_t(sendOffsets.data(), send_data.data(),
                                                         recvOffsets.data(), recv_data.data(), comm );

        if (this_proc == 0) {
            const int recvSizeFromProc0 = recvOffsets[1]-recvOffsets[0];
            EXPECT_EQ(0, recvSizeFromProc0);
            for (int iproc = 1; iproc < num_procs; ++iproc) {
                const int recvSizeFromProc = recvOffsets[iproc+1]-recvOffsets[iproc];
                EXPECT_EQ(1, recvSizeFromProc);
                double expected_recv_value = iproc+1;
                EXPECT_EQ(expected_recv_value, recv_data[recvOffsets[iproc]]);
            }
        }
        else {
            for (int iproc = 0; iproc < num_procs; ++iproc) {
                const int recvSizeFromProc = recvOffsets[iproc+1]-recvOffsets[iproc];
                EXPECT_EQ(0, recvSizeFromProc);
            }
        }
    }
}

TEST(STKMesh_perf, parallel_sum)
{
    const bool with_ghosts = false;
    do_stk_test(with_ghosts);
}

TEST(STKMesh_perf, parallel_sum_including_ghosts)
{
    const bool with_ghosts = true;
    do_stk_test(with_ghosts);
}

} //namespace STKperf
