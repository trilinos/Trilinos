// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for rand, srand, RAND_MAX
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

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

TEST( StkMeshIoBroker, iofixture )
{
  // A simple test for reading and writing an exodus file using the StkMeshIoBroker

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  bool ok = false;
  try {
    // Initialize meta data from exodus file
    fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);
    fixture.create_input_mesh();
    ok = true;
    
    stk::mesh::MetaData & meta_data = fixture.meta_data();

    // Commit meta_data
    meta_data.commit();

    // bulk_data initialize (from exodus file)
    fixture.populate_bulk_data();

    // exodus file creation
    std::string output_base_filename = "unit_test_output.e";
    size_t output_index = fixture.create_output_mesh(output_base_filename, stk::io::WRITE_RESULTS);

    // process output
    const double time_step = 0;
    fixture.process_output_request(output_index, time_step);
  }
  catch(...) {
    ASSERT_TRUE(ok);
  }
  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

TEST( StkMeshIoBroker, testModifyTopology )
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(comm) == 1)
    {
        stk::io::StkMeshIoBroker fixture(comm);
        std::string generated_mesh_spec = "generated:1x1x2";
        fixture.add_mesh_database(generated_mesh_spec, stk::io::READ_MESH);
        fixture.create_input_mesh();

        stk::mesh::MetaData & meta_data = fixture.meta_data();
        stk::mesh::Part & side_set_part = meta_data.declare_part_with_topology("Side Set Part",stk::topology::QUAD_4);
        stk::io::put_io_part_attribute(side_set_part);

        stk::mesh::Part & inactive_part = meta_data.declare_part("Inactive Part");

        fixture.populate_bulk_data();

        {
            stk::mesh::BulkData & mesh = fixture.bulk_data();
            mesh.modification_begin();
            stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
            stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);

            stk::mesh::PartVector add_parts;
            add_parts.push_back(&inactive_part);
            mesh.change_entity_parts(element2, add_parts);

            stk::mesh::Entity side = mesh.declare_entity(stk::topology::FACE_RANK, 1, side_set_part);
            const int elem1_side_ordinal = 5;
            stk::mesh::declare_element_side(mesh, element1, side, elem1_side_ordinal);

            const int elem2_side_ordinal = 4;
            stk::mesh::declare_element_side(mesh, element2, side, elem2_side_ordinal);
            mesh.modification_end();
        }

        const std::string output_base_filename = "StkMeshIoBroker.testModifyTopology.e";
        size_t output_index = 0;
        ASSERT_NO_THROW(output_index = fixture.create_output_mesh(output_base_filename, stk::io::WRITE_RESULTS));
        stk::mesh::Selector active_selector = !inactive_part;
        fixture.set_subset_selector(output_index,active_selector);
        const double time_step = 0.1;
        ASSERT_NO_THROW(fixture.process_output_request(output_index, time_step));

        unlink(output_base_filename.c_str());
    }
}

TEST( StkMeshIoBroker, active_only )
{
  // A simple test for reading and writing an exodus file using the StkMeshIoBroker.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  bool ok = false;
  try {
    // Initialize meta data from exodus file
    fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);

    fixture.create_input_mesh();
    ok = true;
    stk::mesh::MetaData & meta_data = fixture.meta_data();

    // Add an "active" part...
    stk::mesh::Part &active = meta_data.declare_part("active", stk::topology::ELEMENT_RANK);
    meta_data.commit();

    // bulk_data initialize (from exodus file)
    fixture.populate_bulk_data();

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
  }
  catch(...) {
    ASSERT_TRUE(ok);
  }


  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

TEST( StkMeshIoBroker, active_and_all )
{
  // A simple test for reading and writing two exodus files using the StkMeshIoBroker.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 1) {
    return;
  }
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "unit_test.g";

  bool ok = false;
  try {
    fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);
    fixture.create_input_mesh();
    ok = true;
    stk::mesh::MetaData & meta_data = fixture.meta_data();

    // Add an "active" part...
    stk::mesh::Part &active = meta_data.declare_part("active", stk::topology::ELEMENT_RANK);
    meta_data.commit();

    // bulk_data initialize (from exodus file)
    fixture.populate_bulk_data();

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
  }
  catch(...) {
    ASSERT_TRUE(ok);
  }

  // Since correctness can only be established by running SEACAS tools, correctness
  // checking is left to the test XML.
}

TEST( StkMeshIoBroker, large_mesh_test )
{
  // A simple test for reading and writing two exodus files using the StkMeshIoBroker.
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(pm);
  if (numProcs != 1) {
    return;
  }
  stk::io::StkMeshIoBroker fixture(pm);

  std::string input_base_filename = "1mCube_20x20x20.g";

  bool ok = false;
  try {
    // Initialize meta data from exodus file
    fixture.add_mesh_database(input_base_filename, stk::io::READ_MESH);
    fixture.create_input_mesh();
    ok = true;
    
    stk::mesh::MetaData & meta_data = fixture.meta_data();

    // Commit
    meta_data.commit();

    // bulk_data initialize (from exodus file)
    fixture.populate_bulk_data();
    stk::mesh::BulkData &bulk_data = fixture.bulk_data();

    const stk::mesh::BucketVector & element_buckets
      = bulk_data.buckets( stk::topology::ELEMENT_RANK);

    // iterate elements and check num nodal relations
    for ( stk::mesh::BucketVector::const_iterator ib = element_buckets.begin() ;
	  ib != element_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const int length   = b.size();
      for ( int k = 0 ; k < length ; ++k ) {
	// get element
	stk::mesh::Entity elem = b[k];
	size_t num_elem_node_rels = bulk_data.count_valid_connectivity(elem, stk::topology::NODE_RANK);
	EXPECT_EQ( 8u, num_elem_node_rels);
      }
    }
  }
  catch(...) {
    ASSERT_TRUE(ok);
  }
}

