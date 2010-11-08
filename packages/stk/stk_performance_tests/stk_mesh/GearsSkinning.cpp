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

#include <stk_io/IossBridge.hpp>

#include <stk_io/util/UseCase_mesh.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>


#define PI     3.14159265358979
#define TWO_PI 6.28210184121061



namespace {

typedef stk::mesh::fixtures::GearsFixture::CartesianField CartesianField;


//
//-----------------------------------------------------------------------------
//

void seperate_wedge(
    stk::mesh::fixtures::GearsFixture   & fixture,
    stk::mesh::Entity & wedge,
    CartesianField & velocity_field,
    stk::mesh::Part & skin_part
    )
{
  enum { NodeRank = 0 };

  fixture.bulk_data.modification_begin();

  std::vector<size_t> requests(fixture.bulk_data.mesh_meta_data().entity_rank_count(), 0);
  requests[NodeRank] = 6;

  stk::mesh::EntityVector new_nodes;

  fixture.bulk_data.generate_new_entities(requests, new_nodes);

  stk::mesh::PartVector empty_parts, remove_parts;
  remove_parts.push_back(& fixture.cylindrical_coord_part);

  fixture.bulk_data.change_entity_parts(wedge, empty_parts, remove_parts);

  stk::mesh::PairIterRelation relations = wedge.relations(NodeRank);

  for (size_t i = 0; i < 6; ++i) {
    stk::mesh::Entity & old_node = *(relations[ i ].entity());
    stk::mesh::Entity & new_node = *(new_nodes[i]);

    fixture.bulk_data.destroy_relation(wedge, old_node);
    fixture.bulk_data.declare_relation(wedge, new_node, i);

    fixture.bulk_data.copy_entity_fields( old_node, new_node);
  }

  fixture.bulk_data.modification_end();

  stk::mesh::skin_mesh( fixture.bulk_data, 3, &skin_part);


  double avg_velocity_data[3] = {0};
  for (size_t i = 0; i < 6; ++i) {
    stk::mesh::Entity & new_node = *(new_nodes[i]);
    const double * const new_displacement_data =
      stk::mesh::field_data( fixture.displacement_field, new_node);

    const double * const old_displacement_data =
      stk::mesh::field_data( fixture.displacement_field.field_of_state(stk::mesh::StateOld), new_node);

    avg_velocity_data[0] += new_displacement_data[0] - old_displacement_data[0];
    avg_velocity_data[1] += new_displacement_data[1] - old_displacement_data[1];
    avg_velocity_data[2] += new_displacement_data[2] - old_displacement_data[2];

  }

    avg_velocity_data[0] /= 6.0;
    avg_velocity_data[1] /= 6.0;
    avg_velocity_data[2] /= 6.0;

  for (size_t i = 0; i < 6; ++i) {
    stk::mesh::Entity & new_node = *(new_nodes[i]);
    double * const velocity_data =
      stk::mesh::field_data( velocity_field , new_node );

    velocity_data[0] = 1.1*avg_velocity_data[0];
    velocity_data[1] = 1.1*avg_velocity_data[1];
    velocity_data[2] = 1.1*avg_velocity_data[2];

  }
}

//
//-----------------------------------------------------------------------------
//

void find_wedges_to_seperate(
    stk::mesh::fixtures::GearsFixture & fixture,
    stk::mesh::EntityVector & wedges
    )
{
  stk::mesh::Selector select_wedge = fixture.meta_data.universal_part();
  select_wedge &= fixture.cylindrical_coord_part;
  select_wedge &= fixture.wedge_part;

  const stk::mesh::BucketVector wedge_buckets = fixture.bulk_data.buckets(3);

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

void move_detached_wedges(
    stk::mesh::fixtures::GearsFixture & fixture,
    CartesianField & velocity_field
    )
{

  stk::mesh::Selector select_detached_wedges = fixture.meta_data.universal_part();
  select_detached_wedges &= ! stk::mesh::Selector(fixture.cylindrical_coord_part);

  const stk::mesh::BucketVector all_node_buckets = fixture.bulk_data.buckets(0);

  stk::mesh::BucketVector node_buckets;

  stk::mesh::get_buckets(
      select_detached_wedges,
      all_node_buckets,
      node_buckets
      );

  for (stk::mesh::BucketVector::iterator b_itr = node_buckets.begin();
      b_itr != node_buckets.end();
      ++b_itr)
  {
    stk::mesh::Bucket & b = **b_itr;

    const stk::mesh::BucketArray<CartesianField>  velocity_data( velocity_field, b);
    stk::mesh::BucketArray<CartesianField>  displacement_data( fixture.displacement_field, b);

    for (size_t i = 0; i < b.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        displacement_data(j,i) += velocity_data(j,i);
      }
    }
  }

}

} // unnamed namespace

STKUNIT_UNIT_TEST( GearsDemo, skin_gear ) {

  enum {Node = 0};
  enum {SpatialDimension = 3};

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  const size_t NUM_GEARS = 1;

  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, NUM_GEARS);

  stk::mesh::Part & skin_part = fixture.meta_data.declare_part("Skin_part");

  CartesianField & velocity_field = fixture.meta_data.declare_field<CartesianField>("velocity",1);

  stk::mesh::put_field(
      velocity_field,
      Node,
      fixture.meta_data.universal_part(),
      SpatialDimension
      );


  // add io parts
  stk::io::put_io_part_attribute( fixture.hex_part);
  stk::io::put_io_part_attribute( fixture.wedge_part);
  stk::io::put_io_part_attribute(skin_part);
  stk::io::set_field_role(fixture.displacement_field, Ioss::Field::TRANSIENT);


  fixture.meta_data.commit();

  fixture.generate_mesh();


  stk::mesh::Selector select_detached_entites = fixture.meta_data.universal_part();

  for (size_t i = 0; i<NUM_GEARS; ++i) {
    select_detached_entites &= ! stk::mesh::Selector(fixture.get_gear(i).gear_part);
  }

  stk::mesh::skin_mesh( fixture.bulk_data, SpatialDimension, &skin_part);

  //count the types of node, edges, faces, and elements
  {
    std::vector<size_t> counts ;
    comm_mesh_counts( fixture.bulk_data , counts);

    if ( fixture.bulk_data.parallel_rank() == 0 ) {
      std::cout << "N_GEARS Meshing completed and verified" << std::endl ;

      std::cout << "N_GEARS Global Counts:\n{\n "
        << "\tnode = " << counts[0] << "\n"
        << "\tedge = " << counts[1] << "\n"
        << "\tface = " << counts[2] << "\n"
        << "\telem = " << counts[3] << "\n}" << std::endl;
    }
  }

  stk::mesh::EntityVector wedges_to_seperate;
  find_wedges_to_seperate(
      fixture,
      wedges_to_seperate
      );

  const size_t NUM_TIME_STEPS = 15000;

  const double rotation = TWO_PI*4.0/NUM_TIME_STEPS;
  const double x = 0;
  const double y = 0;
  const double z = 0;

  Ioss::Region * out_region = NULL;

  const bool output_exodus_file = false;
  for (size_t time_step = 0; time_step < NUM_TIME_STEPS; ++time_step) {
    //This following section writes mesh data out to an exodus file:

    std::cout << "timestep: " << time_step << std::endl;


    const bool modify_mesh = (!wedges_to_seperate.empty()) && ((time_step %30) == 0);

    if (time_step > 0) {
      for ( size_t i = 0; i < NUM_GEARS; ++i) {
        stk::mesh::fixtures::Gear & gear = fixture.get_gear(i);

        fixture.bulk_data.update_field_data_states();

        const stk::mesh::fixtures::GearData data(rotation,x,y,z);
        gear.move(data);


        if (modify_mesh) {

          stk::mesh::Entity & wedge = *(wedges_to_seperate.back());
          wedges_to_seperate.pop_back();

          seperate_wedge(
              fixture,
              wedge,
              velocity_field,
              skin_part);
        }

      }

      move_detached_wedges(
          fixture,
          velocity_field
          );

    }

    if (output_exodus_file) {
      const bool create_output_file = modify_mesh || (time_step == 0);
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
      int out_step = out_region->add_state(time_step/60.0);
      stk::io::util::process_output_request(*out_region, fixture.bulk_data, out_step);
      out_region->end_mode(Ioss::STATE_TRANSIENT);
    }

  }
  delete out_region; out_region = NULL;
}

