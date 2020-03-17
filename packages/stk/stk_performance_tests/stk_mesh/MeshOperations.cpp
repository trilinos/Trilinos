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

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <gtest/gtest.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <string>

// Globals for command-line arguments
extern int gl_argc;
extern char** gl_argv;

class MeshOperations : public stk::unit_test_util::MeshFixture {
public:
  MeshOperations()
    : p_rank(get_parallel_rank()),
      PHASE_NAMES(NUM_PHASES+1),
      timings(NUM_PHASES + 1, 0.0),
      input_base_filename("mesh_operations.g"),
      output_base_filename("mesh_operations.e")
  {
    PHASE_NAMES[INIT_META_DATA_PHASE_ID] = "Init meta data";
    PHASE_NAMES[INIT_BULK_DATA_PHASE_ID] = "Init bulk data";
    PHASE_NAMES[SKIN_MESH_PHASE_ID]      = "Skin mesh";
    PHASE_NAMES[EXODUS_CREATE_PHASE_ID]  = "Exodus file creation";
    PHASE_NAMES[PROCESS_OUTPUT_PHASE_ID] = "Process output";
    PHASE_NAMES[NUM_PHASES]              = "Total time";
  }

protected:
  void parse_filename_args() {
    const std::string input_flag  = "--heavy-test:input-file=";
    const std::string output_flag = "--heavy-test:output-file=";
    for (int argitr = 1; argitr < gl_argc; ++argitr) {
      std::string argv(gl_argv[argitr]);
      if (argv.find(input_flag) == 0) {
        input_base_filename = argv.replace(0, input_flag.size(), "");
      }
      else if (argv.find(output_flag) == 0) {
        output_base_filename = argv.replace(0, output_flag.size(), "");
      }
    }
  }

  void initialize_meta(stk::io::StkMeshIoBroker& broker) {
    double start_time = stk::wall_time();
    broker.add_mesh_database(input_base_filename, stk::io::READ_MESH);
    broker.create_input_mesh();
    timings[INIT_META_DATA_PHASE_ID] += stk::wall_dtime(start_time);
  }

  void initialize_bulk(stk::io::StkMeshIoBroker& broker) {
    double start_time = stk::wall_time();
    broker.populate_bulk_data();
    timings[INIT_BULK_DATA_PHASE_ID] += stk::wall_dtime(start_time);
  }

  void skin_bulk(stk::io::StkMeshIoBroker& broker) {
    double start_time = stk::wall_time();
    stk::mesh::skin_mesh(broker.bulk_data());
    timings[SKIN_MESH_PHASE_ID] += stk::wall_dtime(start_time);
  }

  size_t create_exodus_file(stk::io::StkMeshIoBroker& broker) {
    double start_time = stk::wall_time();
    size_t index = broker.create_output_mesh( output_base_filename, stk::io::WRITE_RESULTS);
    timings[EXODUS_CREATE_PHASE_ID] += stk::wall_dtime(start_time);
    return index;
  }

  void process_output(stk::io::StkMeshIoBroker& broker, size_t index) {
    double time_step = 0;
    double start_time = stk::wall_time();
    broker.process_output_request(index, time_step);
    timings[PROCESS_OUTPUT_PHASE_ID] += stk::wall_dtime(start_time);
  }

  void sum_timings() {
    for (unsigned i = 0; i < NUM_PHASES; ++i) {
      timings[NUM_PHASES] += timings[i];
    }
  }

  void print_mesh_stats(stk::mesh::BulkData& bulk) {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk, counts);

    if (p_rank == 0) {
      std::cout << "  mesh_stats\n";
      for (unsigned i = 0; i < counts.size(); ++i) {
        std::cout << "    entity rank=" << i << " count=" << counts[i] << "\n";
      }
      std::cout << std::endl;
    }
  }

  void print_timings() {
    stk::all_reduce(get_comm(), stk::ReduceMax<NUM_PHASES+1>(&timings[0]));

    if (p_rank == 0) {
      stk::print_timers_and_memory(&PHASE_NAMES[0], &timings[0], NUM_PHASES + 1);
    }

    stk::parallel_print_time_for_performance_compare(get_comm(), timings[NUM_PHASES]);
  }

private:
  const size_t p_rank;

  static constexpr unsigned NUM_PHASES = 5;

  const unsigned INIT_META_DATA_PHASE_ID = 0;
  const unsigned INIT_BULK_DATA_PHASE_ID = 1;
  const unsigned SKIN_MESH_PHASE_ID      = 2;
  const unsigned EXODUS_CREATE_PHASE_ID  = 3;
  const unsigned PROCESS_OUTPUT_PHASE_ID = 4;

  std::vector<std::string> PHASE_NAMES;
  std::vector<double> timings;

  std::string input_base_filename;
  std::string output_base_filename;
};

TEST_F( MeshOperations, PerformanceTimings )
{
  const unsigned NUM_RUNS = 30;

  parse_filename_args();

  for (unsigned run=0; run<NUM_RUNS; run++) {
    stk::io::StkMeshIoBroker broker(get_comm());

    initialize_meta(broker);
    initialize_bulk(broker);
    skin_bulk(broker);
    size_t index = create_exodus_file(broker);
    process_output(broker, index);

    if (run == NUM_RUNS-1) {
      sum_timings();
      print_timings();
      print_mesh_stats(broker.bulk_data());
    }
  }
}
