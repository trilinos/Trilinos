#include <stk_io/util/IO_Fixture.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <string>

STKUNIT_UNIT_TEST( IOFixture, iofixture )
{
  // A simple test for reading and writing an exodus file using the IOFixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::io::util::IO_Fixture fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.initialize_meta_data( input_base_filename, "exodusii" );

  // Commit meta_data
  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.initialize_bulk_data();

  // exodus file creation
  std::string output_base_filename = "unit_test_output.e";
  fixture.create_output_mesh( output_base_filename, "exodusii" );

  // process output
  double time_step = 0;
  fixture.add_timestep_to_output_mesh( time_step );

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}
