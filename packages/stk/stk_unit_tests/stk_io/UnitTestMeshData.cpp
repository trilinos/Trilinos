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
#include "stk_io/impl/StkIoImpl.hpp"


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


void call_get_or_create_face_at_element_side_and_check(stk::mesh::BulkData & mesh, stk::mesh::Entity element, unsigned side_ordinal, unsigned new_face_global_id, stk::mesh::Part & part) {
    mesh.modification_begin();
    stk::mesh::PartVector add_parts(1, &part);
    stk::mesh::Entity new_face = stk::io::impl::get_or_create_face_at_element_side(mesh, element, side_ordinal, new_face_global_id, part);
    mesh.modification_end();
    ASSERT_TRUE( mesh.is_valid(new_face) );
    EXPECT_EQ( new_face_global_id, mesh.identifier(new_face) );
    ASSERT_EQ( 1u, mesh.num_elements(new_face));
    stk::mesh::Entity attached_element = *mesh.begin_elements(new_face);
    EXPECT_EQ( element, attached_element );
    unsigned attached_side_ordinal = *mesh.begin_element_ordinals(new_face);
    EXPECT_EQ( side_ordinal, attached_side_ordinal );
    EXPECT_TRUE( mesh.bucket(new_face).member(part) );
    unsigned other_element_global_id = 2;
    stk::mesh::Entity other_element = mesh.get_entity(stk::topology::ELEMENT_RANK, other_element_global_id);
    EXPECT_EQ( 0u, mesh.num_faces(other_element) );
}

TEST( StkMeshIoBroker, test_create_face_for_sideset)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if (numProcs != 1) {
        return;
    }
    stk::io::StkMeshIoBroker fixture(pm);

    fixture.add_mesh_database("generated:1x1x2", stk::io::READ_MESH);
    fixture.create_input_mesh();
    stk::topology quad4_topology = stk::topology::QUAD_4;
    stk::mesh::Part & quad4_part = fixture.meta_data().get_topology_root_part(quad4_topology);

    stk::mesh::Part & new_topology_sub_part = fixture.meta_data().declare_part("My Fancy Part",stk::topology::FACE_RANK);
    fixture.meta_data().declare_part_subset(quad4_part, new_topology_sub_part);
    fixture.populate_bulk_data();

    stk::mesh::BulkData & mesh = fixture.bulk_data();

    unsigned elem_global_id = 1;
    stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,elem_global_id);
    unsigned side_ordinal = 5;
    unsigned new_face_global_id = 42;


    stk::mesh::Entity new_face = mesh.get_entity(stk::topology::FACE_RANK, new_face_global_id);
    EXPECT_FALSE( mesh.is_valid(new_face) );
    call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
    call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
    call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,new_topology_sub_part);
}


TEST( StkMeshIoBroker, test_connect_face_to_other_elements)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    int numProcs = stk::parallel_machine_size(pm);
    if (numProcs != 1) {
        return;
    }
    stk::io::StkMeshIoBroker fixture(pm);

    fixture.add_mesh_database("generated:1x1x2", stk::io::READ_MESH);
    fixture.create_input_mesh();
    fixture.populate_bulk_data();

    stk::mesh::BulkData & mesh = fixture.bulk_data();

    stk::topology quad4_topology = stk::topology::QUAD_4;
    stk::mesh::Part & quad4_part = fixture.meta_data().get_topology_root_part(quad4_topology);

    unsigned elem_global_id = 1;
    stk::mesh::Entity element = mesh.get_entity(stk::topology::ELEMENT_RANK,elem_global_id);
    unsigned side_ordinal = 5;
    unsigned new_face_global_id = 42;

    call_get_or_create_face_at_element_side_and_check(mesh,element,side_ordinal,new_face_global_id,quad4_part);
    stk::mesh::Entity new_face = mesh.get_entity(stk::topology::FACE_RANK, new_face_global_id);
    ASSERT_TRUE( mesh.is_valid(new_face) );
    mesh.modification_begin();
    stk::io::impl::connect_face_to_other_elements(mesh,new_face,element,side_ordinal);
    mesh.modification_end();

    ASSERT_EQ( 2u, mesh.num_elements(new_face));
    unsigned other_element_global_id = 2;
    stk::mesh::Entity other_element = mesh.get_entity(stk::topology::ELEMENT_RANK, other_element_global_id);
    ASSERT_EQ( 1u, mesh.num_faces(other_element) );
    stk::mesh::Entity attached_face = *mesh.begin_faces(other_element);
    EXPECT_EQ( new_face, attached_face );

    unsigned attached_side_ordinal = *mesh.begin_face_ordinals(other_element);
    const unsigned expected_other_element_side_ordinal = 4;
    EXPECT_EQ( expected_other_element_side_ordinal, attached_side_ordinal );
}


TEST( StkMeshIoBroker, test_create_shell_status) {
    std::vector<stk::topology> element_topology_vector;
    element_topology_vector.push_back(stk::topology(stk::topology::HEX_8));
    stk::topology original_element_topology = stk::topology::HEX_8;
    std::vector<stk::io::impl::ShellStatus> element_shell_status;
    stk::io::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
    ASSERT_EQ( 1u, element_shell_status.size() );
    EXPECT_EQ( stk::io::impl::NO_SHELLS, element_shell_status[0] );

    element_topology_vector.push_back(stk::topology(stk::topology::HEX_8));
    stk::io::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
    ASSERT_EQ( 2u, element_shell_status.size() );
    EXPECT_EQ( stk::io::impl::NO_SHELLS, element_shell_status[0] );
    EXPECT_EQ( stk::io::impl::NO_SHELLS, element_shell_status[1] );

    element_topology_vector.push_back(stk::topology(stk::topology::SHELL_QUAD_4));
    stk::io::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
    ASSERT_EQ( 3u, element_shell_status.size() );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[0] );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[1] );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[2] );

    original_element_topology = stk::topology::SHELL_QUAD_4;
    stk::io::impl::create_shell_status(element_topology_vector,original_element_topology,element_shell_status);
    ASSERT_EQ( 3u, element_shell_status.size() );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[0] );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID, element_shell_status[1] );
    EXPECT_EQ( stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS, element_shell_status[2] );

}



TEST( StkMeshIoBroker, test_connect_face_to_other_elements_2)
{
    std::vector<int> face_nodes(4);
    face_nodes[0] = 5;
    face_nodes[1] = 6;
    face_nodes[2] = 7;
    face_nodes[3] = 8;

    stk::topology element_side_topology = stk::topology::QUAD_4;

    std::vector<int> element_side_nodes(4);
    element_side_nodes[0] = 5;
    element_side_nodes[1] = 6;
    element_side_nodes[2] = 7;
    element_side_nodes[3] = 8;

    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::NO_SHELLS));
    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

    element_side_nodes[0] = 5;
    element_side_nodes[1] = 8;
    element_side_nodes[2] = 7;
    element_side_nodes[3] = 6;

    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::NO_SHELLS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

    element_side_nodes[0] = 5;
    element_side_nodes[1] = 6;
    element_side_nodes[2] = 7;
    element_side_nodes[3] = 9;

    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::NO_SHELLS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

    face_nodes[0] = 5;
    face_nodes[1] = 8;
    face_nodes[2] = 7;
    face_nodes[3] = 6;

    element_side_nodes[0] = 5;
    element_side_nodes[1] = 6;
    element_side_nodes[2] = 7;
    element_side_nodes[3] = 8;

    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::NO_SHELLS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));

    element_side_nodes[0] = 5;
    element_side_nodes[1] = 8;
    element_side_nodes[2] = 7;
    element_side_nodes[3] = 6;

    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::NO_SHELLS));
    EXPECT_TRUE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_BOTH_SHELLS_OR_BOTH_SOLIDS));
    EXPECT_FALSE(stk::io::impl::should_face_be_connected_to_element_side(face_nodes,element_side_nodes,element_side_topology,stk::io::impl::YES_SHELLS_ONE_SHELL_ONE_SOLID));


}























