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

#include <gtest/gtest.h>
#include <random> // for random_device, mt19937, etc. 
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include "stk_unit_test_utils/stk_mesh_fixtures/Gear.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/GearsFixture.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace {

/** To create a movie from this performance test, change the
  output_exodus_file flag to true and run.  Then do the following:
  > conjoin -output gears.e mesh_*
  > module load viz
  > paraview
  [ open conjoin file, filter until happy ]
  [ change viewport to largest possible (to increase resolution of image)]
  [ press play till happy ]
  [ save as animation with 60fps ]
**/

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

typedef stk::mesh::fixtures::GearsFixture::CartesianField CartesianField;
typedef stk::mesh::Field<int> IntField;

//
//-----------------------------------------------------------------------------
//

/**
 * Make a single wedge seperate itself from the cylinder and fly off on its
 * own.
 */
// if do_separate_wedge == true then wedge must be nonnull pointer
void separate_wedge(
    bool do_separate_wedge,
    stk::mesh::fixtures::GearsFixture   & fixture,
    stk::mesh::Entity wedge,
    CartesianField & velocity_field,
    stk::mesh::Part & skin_part
    )
{
  // Parallel collective call:
  fixture.bulk_data.modification_begin();

  // Request new nodes
  const size_t num_nodes_per_wedge = 6;
  const size_t spatial_dim = fixture.meta_data.spatial_dimension();
  stk::mesh::EntityVector new_nodes;

  std::vector<size_t> requests(fixture.meta_data.entity_rank_count(), 0);
  if (do_separate_wedge) {
    requests[NODE_RANK] = num_nodes_per_wedge;
  } else {
    requests[NODE_RANK] = 0;
  }

  // Parallel collective call:
  fixture.bulk_data.generate_new_entities(requests, new_nodes);

  if (do_separate_wedge) {

    // Remove wedge from cylidrical_coord_part
    stk::mesh::PartVector empty_parts, remove_parts;
    remove_parts.push_back(& fixture.cylindrical_coord_part);
    fixture.bulk_data.change_entity_parts(wedge, empty_parts, remove_parts);

    // Replace wedge's nodes with new nodes because any nodes shared with other
    // entities cannot be taken along with the separated wedge. As a
    // simplification, we simply leave all the old nodes behind.
    STK_ThrowAssert(static_cast<size_t>(fixture.bulk_data.num_nodes(wedge)) == num_nodes_per_wedge);
    stk::mesh::Entity const *rel_nodes = fixture.bulk_data.begin_nodes(wedge);

    for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
      stk::mesh::Entity old_node = rel_nodes[i];
      stk::mesh::Entity new_node = new_nodes[i];

      fixture.bulk_data.destroy_relation(wedge, old_node, i);
      fixture.bulk_data.declare_relation(wedge, new_node, i);

      fixture.bulk_data.copy_entity_fields( old_node, new_node);
    }

    // Compute the velocities of the nodes by taking the average of the
    // differences in the displacements in the last time step. Note that we need
    // all the nodes to have the same velocity; otherwise, the wedge will stretch
    std::vector<double> avg_velocity_data(spatial_dim, 0);
    for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
      stk::mesh::Entity new_node = new_nodes[i];
      const double * const new_displacement_data =
        stk::mesh::field_data( fixture.displacement_field->field_of_state(stk::mesh::StateNew), new_node);

      const double * const old_displacement_data =
        stk::mesh::field_data( fixture.displacement_field->field_of_state(stk::mesh::StateOld), new_node);

      for (size_t k=0 ; k < spatial_dim ; ++k) {
        avg_velocity_data[k] += new_displacement_data[k] - old_displacement_data[k];
      }

    }

    for (size_t k=0 ; k < spatial_dim ; ++k) {
      avg_velocity_data[k] /= 1.0*num_nodes_per_wedge;
    }

    const double detached_wedge_speedup_multiplier = 1.1;
    for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
      stk::mesh::Entity new_node = new_nodes[i];
      double * const velocity_data =
        stk::mesh::field_data( velocity_field , new_node );

      for (size_t k=0 ; k < spatial_dim ; ++k) {
        velocity_data[k] = detached_wedge_speedup_multiplier*avg_velocity_data[k];
      }

    }
  }

  // Parallel collective call:
  fixture.bulk_data.modification_end();

  // Parallel collective call:
  stk::mesh::PartVector add_parts(1,&skin_part);
  stk::mesh::skin_mesh(fixture.bulk_data, add_parts);

  // Parallel collective call:
  fixture.communicate_model_fields();
}

//
//-----------------------------------------------------------------------------
//

/**
 * Stores all the wedges still attached to the cylinder, in random order,
 * in the wedges argument.
 */
void find_and_shuffle_wedges_to_separate(
    stk::mesh::fixtures::GearsFixture & fixture,
    stk::mesh::EntityVector & wedges
    )
{
  // Get all wedges still attached to the cylinder and shuffle them.

  stk::mesh::Selector select_wedge =
    fixture.cylindrical_coord_part &
    fixture.wedge_part &
    fixture.meta_data.locally_owned_part();

  const stk::mesh::BucketVector wedge_buckets = fixture.bulk_data.buckets(stk::topology::ELEMENT_RANK);

  stk::mesh::get_selected_entities(
      select_wedge,
      wedge_buckets,
      wedges
      );

#if __cplusplus >= 201402L
  std::random_device random_number_generator;
  std::mt19937 uniform_random_number_generator(random_number_generator());
  std::shuffle(wedges.begin(), wedges.end(), uniform_random_number_generator);
#else
  std::random_shuffle(wedges.begin(),wedges.end());
#endif
}

//
//-----------------------------------------------------------------------------
//

/**
 * Make all the wedges that have already been detached from the cylinder
 * continue flying through the air.
 */
void move_detached_wedges(
    stk::mesh::fixtures::GearsFixture & fixture,
    CartesianField & velocity_field
    )
{

  // Select all detached nodes by creating a selector for things not in the
  // cylinder part.
  stk::mesh::Selector select_detached_wedges =  (! fixture.cylindrical_coord_part ) &
    (fixture.meta_data.locally_owned_part() | fixture.meta_data.globally_shared_part());

  stk::mesh::BucketVector const& node_buckets = fixture.bulk_data.get_buckets(NODE_RANK, select_detached_wedges);

  // Iterate over selected node_buckets, then iterate over each node in the
  // bucket, adjusting the node's displacement according to its velocity.
  for (stk::mesh::BucketVector::const_iterator b_itr = node_buckets.begin();
      b_itr != node_buckets.end();
      ++b_itr)
  {
    stk::mesh::Bucket & b = **b_itr;

    const CartesianField::value_type*  velocity_data = stk::mesh::field_data(velocity_field, b);
    CartesianField::value_type*  old_displacement_data = stk::mesh::field_data(fixture.displacement_field->field_of_state(stk::mesh::StateOld), b);
    CartesianField::value_type*  new_displacement_data = stk::mesh::field_data(fixture.displacement_field->field_of_state(stk::mesh::StateNew), b);
    int ndim = fixture.meta_data.spatial_dimension();

    for (size_t i = 0; i < b.size(); ++i) {
      for (size_t j = 0; j < fixture.meta_data.spatial_dimension(); ++j) {
        new_displacement_data[j+i*ndim] = old_displacement_data[j+i*ndim] + velocity_data[j+i*ndim];
      }
    }
  }
}

//
//-----------------------------------------------------------------------------
//

void populate_processor_id_field_data( stk::mesh::fixtures::GearsFixture & fixture,
    IntField & processor_field
    )
{
  const unsigned p_rank = fixture.bulk_data.parallel_rank();

  stk::mesh::Selector locally_owned_selector = fixture.meta_data.locally_owned_part();

  stk::mesh::BucketVector const& element_buckets = fixture.bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_owned_selector);

  for (stk::mesh::BucketVector::const_iterator b_itr = element_buckets.begin();
      b_itr != element_buckets.end();
      ++b_itr)
  {
    stk::mesh::Bucket & b = **b_itr;
    int* processor_data = stk::mesh::field_data(processor_field, b);
    for (size_t index = 0; index < b.size(); ++index) {
      processor_data[index] = p_rank;
    }
  }
}

} // namespace

TEST( gears_skinning, gears_skinning )
{
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  const size_t NUM_GEARS = 1;
  double start_time = stk::wall_time();

  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, NUM_GEARS,
                                            stk::mesh::fixtures::GearParams(0.025, 0.6, 1.05, -0.4, 0.4));
  const unsigned p_rank = fixture.bulk_data.parallel_rank();
  std::srand(p_rank); // Seed pseudo-random generator based on processor rank.

  stk::mesh::Part & skin_part = fixture.meta_data.declare_part("Skin_part",fixture.meta_data.side_rank());

  const unsigned ONE_STATE = 1;
  CartesianField & velocity_field = fixture.meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity",ONE_STATE);
  stk::topology::rank_t face_rank = static_cast<stk::topology::rank_t>(fixture.meta_data.side_rank());
  CartesianField & displacement = fixture.meta_data.declare_field<double>(face_rank, "face_displacement",ONE_STATE);
  IntField & processor_field = fixture.meta_data.declare_field<int>(stk::topology::ELEMENT_RANK, "processor_id",ONE_STATE);

  stk::mesh::put_field_on_mesh(
      velocity_field,
      fixture.meta_data.universal_part(),
      fixture.meta_data.spatial_dimension(),
      nullptr
      );

  stk::mesh::put_field_on_mesh(
      processor_field,
      fixture.meta_data.universal_part(),
      nullptr
      );

  // add io parts
  stk::io::put_io_part_attribute( fixture.hex_part);
  stk::io::put_io_part_attribute( fixture.wedge_part);
  stk::io::put_io_part_attribute( skin_part);

  const stk::mesh::PartVector & parts = fixture.meta_data.get_parts();
  for ( stk::mesh::PartVector::const_iterator ip = parts.begin(); ip != parts.end(); ++ip ) {
    stk::mesh::Part & topo_part = **ip;
    if(topo_part.primary_entity_rank() == stk::topology::ELEMENT_RANK-1 && std::string::npos!=topo_part.name().find("FEM_ROOT_CELL_TOPOLOGY_PART"))  {
      std::string t;
      if      (t.npos != topo_part.name().find("Triangle_3"))      t = "skin_wedge6_tri3_1";
      else if (t.npos != topo_part.name().find("Triangle_6"))      t = "skin_wedge15_tri6_2";
      else if (t.npos != topo_part.name().find("Triangle_4"))      t = "skin_wedge6_tri4_3";
      else if (t.npos != topo_part.name().find("Quadrilateral_4")) t = "skin_hex8_quad4_4";
      else if (t.npos != topo_part.name().find("Quadrilateral_8")) t = "skin_hex20_quad8_5";
      else if (t.npos != topo_part.name().find("Quadrilateral_9")) t = "skin_hex27_quad9_6";
      else {
        t = topo_part.name()+std::string("_Skin_part");
        t.erase(t.find("FEM_ROOT_CELL_TOPOLOGY_PART"), sizeof("FEM_ROOT_CELL_TOPOLOGY_PART"));
      }
      stk::mesh::Part & topo_skin_part = fixture.meta_data.declare_part(t, fixture.meta_data.side_rank());
      stk::io::put_io_part_attribute(topo_skin_part);
      fixture.meta_data.declare_part_subset(topo_part, topo_skin_part);
      fixture.meta_data.declare_part_subset(skin_part, topo_skin_part);
      if (   t == "skin_hex8_quad4_4" || t == "skin_hex20_quad8_5") {
        if  (t == "skin_hex8_quad4_4")   t = "skin_wedge6_quad4_4";
        else                             t = "skin_wedge15_quad4_8";
        stk::mesh::Part & topo_skin_part2 = fixture.meta_data.declare_part(t, fixture.meta_data.side_rank());
        stk::io::put_io_part_attribute(topo_skin_part2);
        fixture.meta_data.declare_part_subset(topo_part, topo_skin_part2);
        fixture.meta_data.declare_part_subset(skin_part, topo_skin_part2);
      }
    }
  }
  stk::mesh::Selector surface_select = fixture.meta_data.locally_owned_part();
  {
    stk::mesh::put_field_on_mesh( displacement, skin_part, 3, nullptr);
    const stk::mesh::PartVector &surf_parts = skin_part.subsets();
    for ( stk::mesh::PartVector::const_iterator ip = surf_parts.begin(); ip != surf_parts.end(); ++ip ) {
      stk::mesh::Part & surf_part = **ip;
      if (surf_part.primary_entity_rank() == stk::topology::ELEMENT_RANK-1) {
        surface_select |= surf_part;
        stk::mesh::put_field_on_mesh( displacement, surf_part, 3, nullptr);
      }
    }
  }

  fixture.meta_data.commit();

  fixture.generate_mesh();

  populate_processor_id_field_data(fixture,processor_field);

  stk::mesh::PartVector add_parts(1,&skin_part);
  stk::mesh::skin_mesh(fixture.bulk_data, add_parts);

  stk::mesh::EntityVector wedges_to_separate;

  find_and_shuffle_wedges_to_separate( fixture, wedges_to_separate );

  const size_t NUM_TIME_STEPS = 150;
  const size_t separation_interval = 30;

  const double rotation = TWO_PI*4.0/NUM_TIME_STEPS;
  const double x = 0;
  const double y = 0;
  const double z = 0;
  const stk::mesh::fixtures::GearMovement gear_movement_data(rotation,x,y,z);

  stk::mesh::fixtures::Gear & gear = fixture.get_gear(0);

  // Iterate over the time steps, updating the locations of the entities and
  // writing the current mesh state to output files.
  const bool output_exodus_file = true;

  stk::io::StkMeshIoBroker stkMeshIoBroker;
  stkMeshIoBroker.set_bulk_data(fixture.bulk_data);

  size_t vol_mesh_index = std::numeric_limits<size_t>::max();
  size_t surf_mesh_index = std::numeric_limits<size_t>::max();

  for (size_t time_step = 0; time_step < NUM_TIME_STEPS; ++time_step) {

    // Determine if it's time to separate a wedge
    const bool do_separate_wedge = !wedges_to_separate.empty() && (time_step%separation_interval == 0);

    if (time_step > 0) {

      fixture.bulk_data.update_field_data_states();

      // Move the gear; this will only rotate the gear since x,y,z are all 0
      gear.move(gear_movement_data);
      stk::mesh::Entity wedge = stk::mesh::Entity();

      // Separate the wedge if it's time
      if (do_separate_wedge) {
        wedge = wedges_to_separate.back();
        wedges_to_separate.pop_back();
      }

      separate_wedge(
          do_separate_wedge,
          fixture,
          wedge,
          velocity_field,
          skin_part);

      move_detached_wedges(
          fixture,
          velocity_field
          );

      // Parallel collective call:
      fixture.communicate_model_fields();
    }

    // update a face field
    stk::mesh::BucketVector const& face_buckets = fixture.bulk_data.get_buckets(stk::topology::FACE_RANK, surface_select);
    for (stk::mesh::BucketVector::const_iterator b_itr = face_buckets.begin(); b_itr != face_buckets.end(); ++b_itr) {
      stk::mesh::Bucket & b = **b_itr;
      for (size_t i = 0; i < b.size(); ++i) {
        stk::mesh::Entity face = b[i];
        double *elem_node_disp = stk::mesh::field_data(displacement, face);
        if (elem_node_disp) {
          const size_t num_nodes = b.num_nodes(i);

          stk::mesh::Entity const *node_rels_itr = b.begin_nodes(i);
          stk::mesh::Entity const *node_rels_end = b.end_nodes(i);
          for ( ; node_rels_itr != node_rels_end; ++node_rels_itr)
          {
            const stk::mesh::Entity node = *node_rels_itr;
            double* node_disp = stk::mesh::field_data(*fixture.displacement_field, node);
            elem_node_disp[0] = node_disp[0];
            elem_node_disp[1] = node_disp[1];
            elem_node_disp[2] = node_disp[2];
          }
          elem_node_disp[0] /= num_nodes;
          elem_node_disp[1] /= num_nodes;
          elem_node_disp[2] /= num_nodes;
        }
      }
    }

    //This section writes mesh data out to an exodus file:
    if (output_exodus_file) {
      // Write the output file at the first time step and every time the mesh is modified.
      const bool create_output_file = do_separate_wedge || !time_step;
      // Write the model to the mesh file (topology, coordinates, attributes, etc)
      if (create_output_file)  {
        std::ostringstream volume_out_filename;
        volume_out_filename << "volume_mesh_" << std::setw(7) << std::setfill('0') << time_step << ".e";

        vol_mesh_index = stkMeshIoBroker.create_output_mesh(volume_out_filename.str(), stk::io::WRITE_RESULTS);

        std::ostringstream surface_out_filename;
        surface_out_filename << "surface_mesh_" << std::setw(7) << std::setfill('0') << time_step << ".e";
        surf_mesh_index = stkMeshIoBroker.create_output_mesh(volume_out_filename.str(), stk::io::WRITE_RESULTS);

        stk::io::set_field_role(fixture.displacement_field->field_of_state(stk::mesh::StateNew), Ioss::Field::TRANSIENT);
        stk::io::set_field_role(displacement,    Ioss::Field::TRANSIENT);
        stk::io::set_field_role(processor_field, Ioss::Field::TRANSIENT);

        stkMeshIoBroker.add_field(vol_mesh_index, fixture.displacement_field->field_of_state(stk::mesh::StateNew));
        stkMeshIoBroker.add_field(vol_mesh_index, displacement);
        stkMeshIoBroker.add_field(vol_mesh_index, processor_field);
      }

      stkMeshIoBroker.process_output_request(vol_mesh_index, time_step/60.0);
      stkMeshIoBroker.process_output_request(surf_mesh_index, time_step/60.0);
    }
  }

  double total_time = stk::wall_time() - start_time;
  const char* timer_label = "Total Time";
  if (p_rank == 0) {
    stk::print_timers_and_memory(&timer_label, &total_time, 1);
  }
  stk::parallel_print_time_without_output_and_hwm(MPI_COMM_WORLD, total_time);
}
