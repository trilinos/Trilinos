#include <stk_io/MeshReadWriteUtils.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/perf_util.hpp>

#include <Teuchos_ParameterList.hpp>

#include <string>

// Globals for command-line arguments
extern int* STKUNIT_ARGC;
extern char** STKUNIT_ARGV;

STKUNIT_UNIT_TEST( heavy, heavy )
{
  // A performance test that stresses important parts of stk_mesh
  // (such as mesh-modification in parallel, creation of relations,
  // skinning, etc) so that we can measure stk_mesh performance over
  // time to protect against degradation due to continuing code
  // development.
  //
  // This test uses IO_Fixture to conveniently get a usable mesh from
  // an exodus file.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const size_t p_rank = stk::parallel_machine_rank(pm);

  stk::io::MeshData fixture(pm);

  // Test constants:

  const unsigned NUM_PHASES = 6;

  const unsigned INIT_META_DATA_PHASE_ID = 0;
  const unsigned INIT_BULK_DATA_PHASE_ID = 1;
  const unsigned REBALANCE_PHASE_ID      = 2;
  const unsigned SKIN_MESH_PHASE_ID      = 3;
  const unsigned EXODUS_CREATE_PHASE_ID  = 4;
  const unsigned PROCESS_OUTPUT_PHASE_ID = 5;

  std::vector<std::string> PHASE_NAMES(NUM_PHASES+1);
  PHASE_NAMES[INIT_META_DATA_PHASE_ID] = "Init meta data";
  PHASE_NAMES[INIT_BULK_DATA_PHASE_ID] = "Init bulk data";
  PHASE_NAMES[REBALANCE_PHASE_ID]      = "Rebalance mesh";
  PHASE_NAMES[SKIN_MESH_PHASE_ID]      = "Skin mesh";
  PHASE_NAMES[EXODUS_CREATE_PHASE_ID]  = "Exodus file creation";
  PHASE_NAMES[PROCESS_OUTPUT_PHASE_ID] = "Process output";
  PHASE_NAMES[NUM_PHASES]              = "Total time";

  // timings[6] = sum(timings[0:5])
  std::vector<double> timings(NUM_PHASES + 1, 0.0); // leave room for sum

  // Compute input/output filename

  std::string input_base_filename = "heavy.g"; // Default
  std::string output_base_filename = "heavy.e"; // Default

  // Search cmd-line args
  const std::string input_flag  = "--heavy-test:input-file=";
  const std::string output_flag = "--heavy-test:output-file=";
  for (int argitr = 1; argitr < *STKUNIT_ARGC; ++argitr) {
    std::string argv(STKUNIT_ARGV[argitr]);
    if (argv.find(input_flag) == 0) {
      input_base_filename = argv.replace(0, input_flag.size(), "");
    }
    else if (argv.find(output_flag) == 0) {
      output_base_filename = argv.replace(0, output_flag.size(), "");
    }
  }

  // time meta_data initialize
  {
    double start_time = stk::wall_time();
    fixture.open_mesh_database(input_base_filename);
    fixture.create_input_mesh();
    timings[INIT_META_DATA_PHASE_ID] = stk::wall_dtime(start_time);
  }

  // Commit meta_data
  stk::mesh::MetaData & meta_data = fixture.meta_data();
  meta_data.commit();

  // time bulk_data initialize
  {
    double start_time = stk::wall_time();
    fixture.populate_bulk_data();
    timings[INIT_BULK_DATA_PHASE_ID] = stk::wall_dtime(start_time);
  }

  stk::mesh::BulkData & bulk_data = fixture.bulk_data();

  // time rebalance bulk_data
  {
    Teuchos::ParameterList emptyList;
    stk::rebalance::Zoltan zoltan_partition(pm, meta_data.spatial_dimension(), emptyList);
    stk::mesh::Selector selector(meta_data.locally_owned_part());

    double start_time = stk::wall_time();

    typedef stk::mesh::Field< double, stk::mesh::Cartesian> coord_field_type;
    coord_field_type *coord_field = fixture.meta_data().get_field<coord_field_type>("coordinates");
    stk::rebalance::rebalance(bulk_data,
                              selector,
                              coord_field,
                              NULL /*weight field*/,
                              zoltan_partition);
    timings[REBALANCE_PHASE_ID] = stk::wall_dtime(start_time);
  }

  // time skin bulk_data
  {
    double start_time = stk::wall_time();
    stk::mesh::skin_mesh(bulk_data);
    timings[SKIN_MESH_PHASE_ID] = stk::wall_dtime(start_time);
  }

  // time exodus file creation
  {
    double start_time = stk::wall_time();
    fixture.create_output_mesh( output_base_filename );
    timings[EXODUS_CREATE_PHASE_ID] = stk::wall_dtime(start_time);
  }

  // time process output
  {
    double time_step = 0;
    double start_time = stk::wall_time();
    fixture.process_output_request( time_step );
    timings[PROCESS_OUTPUT_PHASE_ID] = stk::wall_dtime(start_time);
  }

  // Sum times:
  for (unsigned i = 0; i < NUM_PHASES; ++i) {
    timings[NUM_PHASES] += timings[i];
  }

  ThrowRequireMsg(timings.size() == NUM_PHASES+1, "Do not push back on to timings vector");

  // Now process and print times, we print in XML to make parsing easier
  {
    stk::all_reduce(pm, stk::ReduceMax<NUM_PHASES+1>(&timings[0]));

    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( bulk_data , counts);

    if (p_rank == 0) {
      std::cout << "  mesh_stats\n";
      for (unsigned i = 0; i < counts.size(); ++i) {
        std::cout << "    entity rank=" << i << " count=" << counts[i] << "\n";
      }
      std::cout << std::endl;

      stk::print_timers_and_memory(&PHASE_NAMES[0], &timings[0], NUM_PHASES + 1);
    }
  }

  if (bulk_data.parallel_rank() == 0) {
    std::cout<<"### Total Number of Steps Taken ###: 1"<<std::endl;
    std::cout<<"### Total Wall Clock Run Time Used ###: "<< timings[NUM_PHASES] <<std::endl;
    size_t current_memory = 0, high_water_memory = 0;
    stk::get_memory_usage(current_memory, high_water_memory);
    std::cout<<"Total Memory In Use "<<high_water_memory<<std::endl;
  }
}
