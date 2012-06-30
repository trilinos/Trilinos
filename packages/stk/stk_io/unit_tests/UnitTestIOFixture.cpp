#include <stk_io/util/IO_Fixture.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <string>
#include <stdlib.h>

namespace {

void activate_entities(stk::io::util::IO_Fixture &fixture,
                       stk::mesh::Part &active_part)
{
  // Seed generator so multiple calls produce same result
  srand(999999u);
  stk::mesh::fem::FEMMetaData & meta = fixture.meta_data();
  stk::mesh::BulkData &bulk = fixture.bulk_data();

  stk::mesh::EntityRank elem_rank = meta.element_rank();

  stk::mesh::PartVector add_parts(1, &active_part);

  bulk.modification_begin();
  const stk::mesh::PartVector & all_parts = meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == elem_rank) {
      // Get all entities (elements) on this part...
      std::vector<stk::mesh::Entity*> entities;
      stk::mesh::Selector select = meta.locally_owned_part() & *part;
      stk::mesh::get_selected_entities(select, bulk.buckets(elem_rank), entities);
      for (size_t i=0; i < entities.size(); i++) {
        if (rand() > (RAND_MAX/4)*3)
          bulk.change_entity_parts(*entities[i], add_parts);
      }
    }
  }
  bulk.modification_end();
}

}

STKUNIT_UNIT_TEST( IOFixture, iofixture )
{
  // A simple test for reading and writing an exodus file using the IOFixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::io::util::IO_Fixture fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.initialize_meta_data( input_base_filename, "exodusii" );

  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();

  // Commit meta_data
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.initialize_bulk_data();

  // exodus file creation
  std::string output_base_filename = "unit_test_output.e";
  fixture.create_output_mesh( output_base_filename, "exodusii" );

  // process output
  const double time_step = 0;
  fixture.add_timestep_to_output_mesh( time_step );

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( IOFixture, active_only )
{
  // A simple test for reading and writing an exodus file using the IOFixture.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::util::IO_Fixture fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.initialize_meta_data( input_base_filename, "exodusii" );
  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", meta_data.element_rank());
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.initialize_bulk_data();

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // Set the output filter on the mesh_data...
  stk::mesh::Selector active_selector(active);
  fixture.mesh_data().m_anded_selector = &active_selector;

  // exodus file creation
  std::string output_base_filename = "unit_test_output_filtered.e";
  fixture.create_output_mesh( output_base_filename, "exodusii" );

  // process output
  const double time_step = 0;
  fixture.add_timestep_to_output_mesh( time_step );


  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( IOFixture, active_and_all )
{
  // A simple test for reading and writing two exodus files using the IOFixture.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::util::IO_Fixture fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.initialize_meta_data( input_base_filename, "exodusii" );
  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", meta_data.element_rank());
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.initialize_bulk_data();

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // Set the output filter on the mesh_data...
  stk::mesh::Selector active_selector(active);
  fixture.mesh_data().m_anded_selector = &active_selector;

  // exodus file creation
  std::string filtered_output_base_filename = "unit_test_output_first_of_two.e";
  fixture.create_output_mesh( filtered_output_base_filename, "exodusii" );

  // process output
  double time_step = 0;
  fixture.add_timestep_to_output_mesh( time_step );


  Teuchos::RCP<Ioss::Region> active_output_ioss_region = fixture.output_ioss_region();

  // Set the output filter on the mesh_data...
  stk::mesh::Selector universal_selector(meta_data.universal_part());
  fixture.mesh_data().m_anded_selector = &universal_selector;

  // exodus file creation
  std::string unfiltered_output_base_filename = "unit_test_output_second_of_two.e";
  fixture.create_output_mesh( unfiltered_output_base_filename, "exodusii" );
  fixture.add_timestep_to_output_mesh( time_step );

  Teuchos::RCP<Ioss::Region> universal_output_ioss_region = fixture.output_ioss_region();

  ++time_step;

  fixture.output_ioss_region(active_output_ioss_region);
  fixture.mesh_data().m_anded_selector = &active_selector;
  fixture.add_timestep_to_output_mesh( time_step );

  fixture.output_ioss_region(universal_output_ioss_region);
  fixture.mesh_data().m_anded_selector = &universal_selector;
  fixture.add_timestep_to_output_mesh( time_step );

  ++time_step;

  fixture.output_ioss_region(active_output_ioss_region);
  fixture.mesh_data().m_anded_selector = &active_selector;
  fixture.add_timestep_to_output_mesh( time_step );

  fixture.output_ioss_region(universal_output_ioss_region);
  fixture.mesh_data().m_anded_selector = &universal_selector;
  fixture.add_timestep_to_output_mesh( time_step );
  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

