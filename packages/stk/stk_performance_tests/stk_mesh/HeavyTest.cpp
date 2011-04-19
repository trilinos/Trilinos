#include <stk_io/util/IO_Fixture.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <Teuchos_ParameterList.hpp>

#include <string>

STKUNIT_UNIT_TEST( HeavyTest, heavytest )
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

  const size_t p_size = stk::parallel_machine_size(pm);
  const size_t p_rank = stk::parallel_machine_rank(pm);

  stk::io::util::IO_Fixture fixture(pm);

  const unsigned num_phases = 6;
  const char* phase_names[num_phases] = {
    "Init meta data",       // timings[0] = init meta data
    "Init bulk data",       // timings[1] = init bulk data
    "Rebalance mesh",       // timings[2] = rebalance mesh
    "Skin mesh",            // timings[3] = skin mesh
    "Exodus file creation", // timings[4] = exodus file creation
    "Process output"        // timings[5] = process output
  };
  // timings[6] = sum(timings[0:5])
  double timings[num_phases + 1]; // leave room for sum
  unsigned phase = 0;
  std::string input_base_filename = "heavy_performance_test.g";

  // time meta_data initialize
  {
    double start_time = stk::wall_time();
    fixture.initialize_meta_data( input_base_filename, "exodusii" );
    timings[phase++] = stk::wall_dtime(start_time);
  }

  // Commit meta_data
  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();
  meta_data.commit();

  // time bulk_data initialize
  {
    double start_time = stk::wall_time();
    fixture.initialize_bulk_data();
    timings[phase++] = stk::wall_dtime(start_time);
  }

  stk::mesh::BulkData & bulk_data = fixture.bulk_data();

  // time rebalance bulk_data
  {
    Teuchos::ParameterList emptyList;
    stk::rebalance::Zoltan zoltan_partition(pm, meta_data.spatial_dimension(), emptyList);
    stk::mesh::Selector selector(meta_data.locally_owned_part());

    double start_time = stk::wall_time();
    stk::rebalance::rebalance(bulk_data,
                              selector,
                              &fixture.get_coordinate_field(),
                              NULL /*weight field*/,
                              zoltan_partition);
    timings[phase++] = stk::wall_dtime(start_time);
  }

  // time skin bulk_data
  {
    double start_time = stk::wall_time();
    stk::mesh::skin_mesh( bulk_data, meta_data.spatial_dimension());
    timings[phase++] = stk::wall_dtime(start_time);
  }

  // time exodus file creation
  {
    std::string output_base_filename = "heavy_performance_test.e";
    double start_time = stk::wall_time();
    fixture.create_output_mesh( output_base_filename, "exodusii" );
    timings[phase++] = stk::wall_dtime(start_time);
  }

  // time process output
  {
    double time_step = 0;
    double start_time = stk::wall_time();
    fixture.add_timestep_to_output_mesh( time_step );
    timings[phase++] = stk::wall_dtime(start_time);
  }

  ThrowRequireMsg(phase == num_phases, phase << " != " << num_phases);

  // Sum times;
  timings[num_phases] = 0;
  for (unsigned i = 0; i < num_phases; ++i) {
    timings[num_phases] += timings[i];
  }

  // Now process and print times
  {
    stk::all_reduce(pm, stk::ReduceMax<num_phases+1>(timings));

    std::vector<size_t> counts ;
    stk::mesh::fem::comm_mesh_counts( bulk_data , counts);

    if (p_rank == 0) {
      std::cout << "\n\n";
      std::cout << "Num procs: " << p_size << "\n";
      for (unsigned i = 0; i < num_phases; ++i) {
        std::cout << "\t" << phase_names[i] << ": " << timings[i] << std::endl;
      }
      std::cout << "\tTotal time: " << timings[num_phases] << std::endl;
      std::cout << "\n\n";

      for (unsigned i = 0; i < counts.size(); ++i) {
        std::cout << "Counts for rank: " << i << ": " << counts[i] << std::endl;
      }
    }
  }
}
