#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <string>
#include <stdlib.h>

namespace {

void activate_entities(stk::io::MeshData &fixture,
                       stk::mesh::Part &active_part)
{
  // Seed generator so multiple calls produce same result
  srand(999999u);
  stk::mesh::MetaData & meta = fixture.meta_data();
  stk::mesh::BulkData &bulk = fixture.bulk_data();

  stk::mesh::EntityRank elem_rank = stk::mesh::MetaData::ELEMENT_RANK;

  stk::mesh::PartVector add_parts(1, &active_part);

  bulk.modification_begin();
  const stk::mesh::PartVector & all_parts = meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;
    if (stk::io::is_part_io_part(*part) && part->primary_entity_rank() == elem_rank) {
      // Get all entities (elements) on this part...
      std::vector<stk::mesh::Entity> entities;
      stk::mesh::Selector select = meta.locally_owned_part() & *part;
      stk::mesh::get_selected_entities(select, bulk.buckets(elem_rank), entities);
      for (size_t i=0; i < entities.size(); i++) {
        if (rand() > (RAND_MAX/4)*3)
          bulk.change_entity_parts(entities[i], add_parts);
      }
    }
  }
  bulk.modification_end();
}

}

STKUNIT_UNIT_TEST( MeshData, iofixture )
{
  // A simple test for reading and writing an exodus file using the MeshData

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::io::MeshData fixture;

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.create_input_mesh( "exodus", input_base_filename, pm );

  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Commit meta_data
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data();

  // exodus file creation
  std::string output_base_filename = "unit_test_output.e";
  fixture.create_output_mesh(output_base_filename);

  // process output
  const double time_step = 0;
  fixture.process_output_request( time_step );

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( MeshData, active_only )
{
  // A simple test for reading and writing an exodus file using the MeshData.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::MeshData fixture;

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.create_input_mesh("exodus", input_base_filename, pm);
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", stk::mesh::MetaData::ELEMENT_RANK);
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data();

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // Set the output filter on the mesh_data...
  stk::mesh::Selector active_selector(active);
  fixture.set_selector(&active_selector);

  // exodus file creation
  std::string output_base_filename = "unit_test_output_filtered.e";
  fixture.create_output_mesh( output_base_filename );

  // process output
  const double time_step = 0;
  fixture.process_output_request( time_step );


  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( MeshData, active_and_all )
{
  // A simple test for reading and writing two exodus files using the MeshData.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::MeshData fixture;

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  fixture.create_input_mesh("exodus", input_base_filename, pm );
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", stk::mesh::MetaData::ELEMENT_RANK);
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data();

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // Set the output filter on the mesh_data...
  stk::mesh::Selector active_selector(active);
  fixture.set_selector(&active_selector);

  // exodus file creation
  std::string filtered_output_base_filename = "unit_test_output_first_of_two.e";
  fixture.create_output_mesh( filtered_output_base_filename );

  // process output
  double time_step = 0;
  fixture.process_output_request( time_step );


  Teuchos::RCP<Ioss::Region> active_output_io_region = fixture.output_io_region();

  // Set the output filter on the mesh_data...
  stk::mesh::Selector universal_selector(meta_data.universal_part());
  fixture.set_selector(&universal_selector);

  // exodus file creation
  std::string unfiltered_output_base_filename = "unit_test_output_second_of_two.e";
  fixture.create_output_mesh( unfiltered_output_base_filename );
  fixture.process_output_request( time_step );

  Teuchos::RCP<Ioss::Region> universal_output_io_region = fixture.output_io_region();

  ++time_step;

  fixture.set_output_io_region(active_output_io_region);
  fixture.set_selector(&active_selector);
  fixture.process_output_request( time_step );

  fixture.set_output_io_region(universal_output_io_region);
  fixture.set_selector(&universal_selector);
  fixture.process_output_request( time_step );

  ++time_step;

  fixture.set_output_io_region(active_output_io_region);
  fixture.set_selector(&active_selector);
  fixture.process_output_request( time_step );

  fixture.set_output_io_region(universal_output_io_region);
  fixture.set_selector(&universal_selector);
  fixture.process_output_request( time_step );
  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( MeshData, large_mesh_test )
{
  // A simple test for reading and writing two exodus files using the MeshData.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::MeshData fixture;

  std::string input_base_filename = "1mCube_20x20x20.g";

  // Initialize meta data from exodus file
  fixture.create_input_mesh("exodus", input_base_filename, pm );
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Commit
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data();
  stk::mesh::BulkData &bulk_data = fixture.bulk_data();

  const std::vector< stk::mesh::Bucket * > & element_buckets
    = bulk_data.buckets( stk::mesh::MetaData::ELEMENT_RANK);

  // iterate elements and check num nodal relations
  for ( std::vector<stk::mesh::Bucket*>::const_iterator ib = element_buckets.begin() ;
        ib != element_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const int length   = b.size();
    for ( int k = 0 ; k < length ; ++k ) {
      // get element
      stk::mesh::Entity elem = b[k];
      stk::mesh::PairIterRelation elem_node_rels = elem.relations(stk::mesh::MetaData::NODE_RANK);
      STKUNIT_EXPECT_EQ( 8u, elem_node_rels.size());
    }
  }
}
