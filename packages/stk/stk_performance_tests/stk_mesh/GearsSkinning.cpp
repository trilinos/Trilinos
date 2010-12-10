#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/fixtures/GearsFixture.hpp>
#include <stk_mesh/fixtures/Gear.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_io/util/UseCase_mesh.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

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

using stk::mesh::fem::NODE_RANK;

const double PI     = 3.14159265358979;
const double TWO_PI = 2 * PI;

typedef stk::mesh::fixtures::GearsFixture::CartesianField CartesianField;
typedef stk::mesh::Field<int> IntField;

//
//-----------------------------------------------------------------------------
//

/**
 * Make a single wedge seperate itself from the cylinder and fly off on its
 * own.
 */
void separate_wedge(
    stk::mesh::fixtures::GearsFixture & fixture,
    stk::mesh::Entity & wedge,
    CartesianField & velocity_field,
    stk::mesh::Part & skin_part
    )
{
  fixture.bulk_data.modification_begin();

  // Request new nodes
  const size_t num_nodes_per_wedge = 6;
  const size_t spatial_dim = fixture.fem.get_spatial_dimension();
  stk::mesh::EntityVector new_nodes;
  std::vector<size_t> requests(fixture.bulk_data.mesh_meta_data().entity_rank_count(), 0);
  requests[NODE_RANK] = num_nodes_per_wedge;
  fixture.bulk_data.generate_new_entities(requests, new_nodes);

  // Remove wedge from cylidrical_coord_part
  stk::mesh::PartVector empty_parts, remove_parts;
  remove_parts.push_back(& fixture.cylindrical_coord_part);
  fixture.bulk_data.change_entity_parts(wedge, empty_parts, remove_parts);

  // Replace wedge's nodes with new nodes because any nodes shared with other
  // entities cannot be taken along with the separated wedge. As a
  // simplification, we simply leave all the old nodes behind.
  stk::mesh::PairIterRelation relations = wedge.relations(NODE_RANK);
  ThrowRequire(relations.size() == num_nodes_per_wedge);
  for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
    stk::mesh::Entity & old_node = *(relations[ i ].entity());
    stk::mesh::Entity & new_node = *(new_nodes[i]);

    fixture.bulk_data.destroy_relation(wedge, old_node);
    fixture.bulk_data.declare_relation(wedge, new_node, i);

    fixture.bulk_data.copy_entity_fields( old_node, new_node);
  }

  fixture.bulk_data.modification_end();

  // Reskin the mesh
  stk::mesh::skin_mesh( fixture.bulk_data, fixture.element_rank, &skin_part);

  // Compute the velocities of the nodes by taking the average of the
  // differences in the displacements in the last time step. Note that we need
  // all the nodes to have the same velocity; otherwise, the wedge will stretch
  std::vector<double> avg_velocity_data(spatial_dim, 0);
  for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
    stk::mesh::Entity & new_node = *(new_nodes[i]);
    const double * const new_displacement_data =
      stk::mesh::field_data( fixture.displacement_field, new_node);

    const double * const old_displacement_data =
      stk::mesh::field_data( fixture.displacement_field.field_of_state(stk::mesh::StateOld), new_node);

    for (size_t k=0 ; k < spatial_dim ; ++k) {
      avg_velocity_data[k] += new_displacement_data[k] - old_displacement_data[k];
    }

  }

  for (size_t k=0 ; k < spatial_dim ; ++k) {
    avg_velocity_data[k] /= 1.0*num_nodes_per_wedge;
  }

  const double detached_wedge_speedup_multiplier = 1.1;
  for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
    stk::mesh::Entity & new_node = *(new_nodes[i]);
    double * const velocity_data =
      stk::mesh::field_data( velocity_field , new_node );

    for (size_t k=0 ; k < spatial_dim ; ++k) {
      velocity_data[k] = detached_wedge_speedup_multiplier*avg_velocity_data[k];
    }
  }
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

  stk::mesh::Selector select_wedge = fixture.cylindrical_coord_part & fixture.wedge_part;

  const stk::mesh::BucketVector wedge_buckets =
    fixture.bulk_data.buckets(fixture.element_rank);

  stk::mesh::get_selected_entities(
      select_wedge,
      wedge_buckets,
      wedges
      );

  std::random_shuffle(wedges.begin(),wedges.end());
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
  stk::mesh::Selector select_detached_wedges =  ! fixture.cylindrical_coord_part;

  const stk::mesh::BucketVector all_node_buckets =
    fixture.bulk_data.buckets(NODE_RANK);

  stk::mesh::BucketVector node_buckets;

  stk::mesh::get_buckets(
      select_detached_wedges,
      all_node_buckets,
      node_buckets
      );

  // Iterate over selected node_buckets, then iterate over each node in the
  // bucket, adjusting the node's displacement according to its velocity.
  for (stk::mesh::BucketVector::iterator b_itr = node_buckets.begin();
      b_itr != node_buckets.end();
      ++b_itr)
  {
    stk::mesh::Bucket & b = **b_itr;

    const stk::mesh::BucketArray<CartesianField>  velocity_data( velocity_field, b);
    stk::mesh::BucketArray<CartesianField>  displacement_data( fixture.displacement_field, b);

    for (size_t i = 0; i < b.size(); ++i) {
      for (size_t j = 0; j < fixture.fem.get_spatial_dimension(); ++j) {
        displacement_data(j,i) += velocity_data(j,i);
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

  stk::mesh::BucketVector all_element_buckets =
    fixture.bulk_data.buckets(fixture.element_rank);
  stk::mesh::BucketVector element_buckets;

  stk::mesh::get_buckets(
      locally_owned_selector,
      all_element_buckets,
      element_buckets
      );

  for (stk::mesh::BucketVector::iterator b_itr = element_buckets.begin();
      b_itr != element_buckets.end();
      ++b_itr)
  {
    stk::mesh::Bucket & b = **b_itr;
    stk::mesh::BucketArray<IntField> processor_data( processor_field, b);
    for (size_t index = 0; index < b.size(); ++index) {
      processor_data(index) = p_rank;
    }
  }
}

//
//-----------------------------------------------------------------------------
//

} // unnamed namespace

STKUNIT_UNIT_TEST( GearsDemo, skin_gear )
{
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  const size_t NUM_GEARS = 1;

  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, NUM_GEARS);

  stk::mesh::Part & skin_part = fixture.meta_data.declare_part("Skin_part");

  const unsigned ONE_STATE = 1;
  CartesianField & velocity_field = fixture.meta_data.declare_field<CartesianField>("velocity",ONE_STATE);
  IntField & processor_field = fixture.meta_data.declare_field<IntField>("processor_id",ONE_STATE);

  stk::mesh::put_field(
      velocity_field,
      NODE_RANK,
      fixture.meta_data.universal_part(),
      fixture.fem.get_spatial_dimension()
      );

  stk::mesh::put_field(
      processor_field,
      fixture.element_rank,
      fixture.meta_data.universal_part()
      );

  // add io parts
  stk::io::put_io_part_attribute( fixture.hex_part);
  stk::io::put_io_part_attribute( fixture.wedge_part);
  stk::io::put_io_part_attribute(skin_part);
  stk::io::set_field_role(fixture.displacement_field, Ioss::Field::TRANSIENT);
  stk::io::set_field_role(processor_field, Ioss::Field::TRANSIENT);

  fixture.meta_data.commit();

  fixture.generate_mesh();

  populate_processor_id_field_data(fixture,processor_field);

  stk::mesh::Selector select_detached_entites = fixture.meta_data.universal_part();

  for (size_t i = 0; i<NUM_GEARS; ++i) {
    select_detached_entites &= ! stk::mesh::Selector(fixture.get_gear(i).gear_part);
  }

  stk::mesh::skin_mesh( fixture.bulk_data,
                        fixture.fem.get_spatial_dimension(),
                        &skin_part );

  stk::mesh::EntityVector wedges_to_separate;

  find_and_shuffle_wedges_to_separate( fixture, wedges_to_separate );

  const size_t NUM_TIME_STEPS = 15000;
  const size_t separation_interval = 30;

  const double rotation = TWO_PI*4.0/NUM_TIME_STEPS;
  const double x = 0;
  const double y = 0;
  const double z = 0;

  Ioss::Region * out_region = NULL;

  // Iterate over the time steps, updating the locations of the entities and
  // writing the current mesh state to output files.
  const bool output_exodus_file = false;
  for (size_t time_step = 0; time_step < NUM_TIME_STEPS; ++time_step) {
    // Determine if it's time to separate a wedge
    const bool do_separate_wedge = !wedges_to_separate.empty() &&
                                   (time_step % separation_interval == 0);

    if (time_step > 0) {
      // Iterate over gears
      for ( size_t i = 0; i < NUM_GEARS; ++i) {
        stk::mesh::fixtures::Gear & gear = fixture.get_gear(i);

        fixture.bulk_data.update_field_data_states();

        // Move the gear; this will only rotate the gear since x,y,z are all 0
        const stk::mesh::fixtures::GearMovement movement(rotation,x,y,z);
        gear.move(movement);

        if (do_separate_wedge) {
          // Separate the wedge if it's time
          stk::mesh::Entity & wedge = *(wedges_to_separate.back());
          wedges_to_separate.pop_back();

          separate_wedge(
              fixture,
              wedge,
              velocity_field,
              skin_part);
        }
      }

      move_detached_wedges( fixture, velocity_field );

      fixture.communicate_model_fields();
    }

    if (output_exodus_file) {
      // Write the output file at the first time step and every time the mesh is modified.
      const bool create_output_file = do_separate_wedge || (time_step == 0);
      // Write the model to the mesh file (topology, coordinates, attributes, etc)
      if (create_output_file)  {
        delete out_region;
        out_region = NULL;

        std::ostringstream out_filename;
        out_filename << "mesh_" << std::setw(7) << std::setfill('0') << time_step << ".e";

        out_region = stk::io::util::create_output_mesh(out_filename.str(), "", "",
            MPI_COMM_WORLD,
            fixture.bulk_data,
            NULL,
            fixture.meta_data,
            true,
            false);
      }

      out_region->begin_mode(Ioss::STATE_TRANSIENT);
      int out_step = out_region->add_state(time_step / 60.0);
      stk::io::util::process_output_request(*out_region, fixture.bulk_data, out_step);
      out_region->end_mode(Ioss::STATE_TRANSIENT);
    }
  }

  delete out_region;
  out_region = NULL;
}

