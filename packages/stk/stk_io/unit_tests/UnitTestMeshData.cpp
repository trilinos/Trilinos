#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <string>
#include <stdlib.h>

namespace {

void activate_entities(stk::io::StkMeshIoBroker &fixture,
                       stk::mesh::Part &active_part)
{
  // Seed generator so multiple calls produce same result
  srand(999999u);
  stk::mesh::MetaData & meta = fixture.meta_data();
  stk::mesh::BulkData &bulk = fixture.bulk_data();

  stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;

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

STKUNIT_UNIT_TEST( StkMeshIoBroker, iofixture )
{
  // A simple test for reading and writing an exodus file using the StkMeshIoBroker

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  size_t index = fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);
  fixture.create_input_mesh(index);

  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Commit meta_data
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data(index);

  // exodus file creation
  std::string output_base_filename = "unit_test_output.e";
  size_t output_index = fixture.create_output_mesh(output_base_filename, stk::io::WRITE_RESULTS);

  // process output
  const double time_step = 0;
  fixture.process_output_request(output_index, time_step);

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( StkMeshIoBroker, active_only )
{
  // A simple test for reading and writing an exodus file using the StkMeshIoBroker.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  // Initialize meta data from exodus file
  size_t input_index = fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);

  fixture.create_input_mesh(input_index);
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", stk::topology::ELEMENT_RANK);
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data(input_index);

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // exodus file creation
  std::string output_base_filename = "unit_test_output_filtered.e";
  size_t index = fixture.create_output_mesh( output_base_filename, stk::io::WRITE_RESULTS );

  // Set the output filter on the mesh_data...
  stk::mesh::Selector active_selector(active);
  fixture.set_subset_selector(index, active_selector);

  // process output
  const double time_step = 0;
  fixture.process_output_request(index, time_step);


  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

STKUNIT_UNIT_TEST( StkMeshIoBroker, active_and_all )
{
  // A simple test for reading and writing two exodus files using the StkMeshIoBroker.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  size_t input_index = fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);
  
  fixture.create_input_mesh(input_index);
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Add an "active" part...
  stk::mesh::Part &active = meta_data.declare_part("active", stk::topology::ELEMENT_RANK);
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data(input_index);

  // Put some entities into the "active" part...
  // This will be used to test the I/O filtering via a selector...
  activate_entities(fixture, active);

  // exodus file creation
  std::string filtered_output_base_filename = "unit_test_output_first_of_two.e";
  std::string unfiltered_output_base_filename = "unit_test_output_second_of_two.e";
  size_t filtered_index =  fixture.create_output_mesh( filtered_output_base_filename, stk::io::WRITE_RESULTS );
  size_t universal_index = fixture.create_output_mesh( unfiltered_output_base_filename, stk::io::WRITE_RESULTS );

  // Set the output filter on the mesh_data...
  // Only output the part declared above as "active"
  stk::mesh::Selector active_selector(active);
  fixture.set_subset_selector(filtered_index, active_selector);

  // process output
  double time_step = 0;
  fixture.process_output_request(filtered_index,  time_step);
  fixture.process_output_request(universal_index, time_step);

  ++time_step;

  fixture.process_output_request(filtered_index,  time_step);
  fixture.process_output_request(universal_index, time_step);

  ++time_step;

  fixture.process_output_request(filtered_index,  time_step);
  fixture.process_output_request(universal_index, time_step);

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
 }

STKUNIT_UNIT_TEST( StkMeshIoBroker, large_mesh_test )
{
  // A simple test for reading and writing two exodus files using the StkMeshIoBroker.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "1mCube_20x20x20.g";

  // Initialize meta data from exodus file
  size_t input_index = fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);

  fixture.create_input_mesh(input_index);
  stk::mesh::MetaData & meta_data = fixture.meta_data();

  // Commit
  meta_data.commit();

  // bulk_data initialize (from exodus file)
  fixture.populate_bulk_data(input_index);
  stk::mesh::BulkData &bulk_data = fixture.bulk_data();

  const std::vector< stk::mesh::Bucket * > & element_buckets
    = bulk_data.buckets( stk::topology::ELEMENT_RANK);

  // iterate elements and check num nodal relations
  for ( std::vector<stk::mesh::Bucket*>::const_iterator ib = element_buckets.begin() ;
        ib != element_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const int length   = b.size();
    for ( int k = 0 ; k < length ; ++k ) {
      // get element
      stk::mesh::Entity elem = b[k];
      size_t num_elem_node_rels = bulk_data.count_valid_connectivity(elem, stk::topology::NODE_RANK);
      STKUNIT_EXPECT_EQ( 8u, num_elem_node_rels);
    }
  }
}
