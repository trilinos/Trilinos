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

#include <stk_io/IossBridge.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>

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

static const size_t NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

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
    stk::mesh::Entity * wedge,
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
    fixture.bulk_data.change_entity_parts(*wedge, empty_parts, remove_parts);

    // Replace wedge's nodes with new nodes because any nodes shared with other
    // entities cannot be taken along with the separated wedge. As a
    // simplification, we simply leave all the old nodes behind.
    stk::mesh::PairIterRelation relations = wedge->relations(NODE_RANK);
    ThrowRequire(relations.size() == num_nodes_per_wedge);

    for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
      stk::mesh::Entity & old_node = *(relations[ i ].entity());
      stk::mesh::Entity & new_node = *(new_nodes[i]);

      fixture.bulk_data.destroy_relation(*wedge, old_node, i);
      fixture.bulk_data.declare_relation(*wedge, new_node, i);

      fixture.bulk_data.copy_entity_fields( old_node, new_node);
    }

    // Compute the velocities of the nodes by taking the average of the
    // differences in the displacements in the last time step. Note that we need
    // all the nodes to have the same velocity; otherwise, the wedge will stretch
    std::vector<double> avg_velocity_data(spatial_dim, 0);
    for (size_t i = 0; i < num_nodes_per_wedge; ++i) {
      stk::mesh::Entity & new_node = *(new_nodes[i]);
      const double * const new_displacement_data =
        stk::mesh::field_data( fixture.displacement_field.field_of_state(stk::mesh::StateNew), new_node);

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

  // Parallel collective call:
  fixture.bulk_data.modification_end();

  // Parallel collective call:
  stk::mesh::skin_mesh( fixture.bulk_data, fixture.element_rank, &skin_part);

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

  const stk::mesh::BucketVector wedge_buckets = fixture.bulk_data.buckets(fixture.element_rank);

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
  stk::mesh::Selector select_detached_wedges =  (! fixture.cylindrical_coord_part ) &
    (fixture.meta_data.locally_owned_part() | fixture.meta_data.globally_shared_part());

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
    stk::mesh::BucketArray<CartesianField>  old_displacement_data( fixture.displacement_field.field_of_state(stk::mesh::StateOld), b);
    stk::mesh::BucketArray<CartesianField>  new_displacement_data( fixture.displacement_field.field_of_state(stk::mesh::StateNew), b);

    for (size_t i = 0; i < b.size(); ++i) {
      for (size_t j = 0; j < fixture.meta_data.spatial_dimension(); ++j) {
        new_displacement_data(j,i) = old_displacement_data(j,i) + velocity_data(j,i);
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

namespace {

Ioss::Region *create_output_mesh(
    const std::string &mesh_filename,
    stk::mesh::BulkData &bulk_data,
    const bool skin = false)
{
  const bool add_all_fields = false;
  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(
      "exodusII",
      mesh_filename,
      Ioss::WRITE_RESULTS,
      bulk_data.parallel()
      );
  if (dbo == NULL || !dbo->ok()) {
    std::cerr << "ERROR: Could not open results database '" << mesh_filename
      << "' of type 'exodusII'\n";
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'out_region' owns 'dbo' pointer at this time...
  const std::string name = std::string("results_output_")+mesh_filename;
  Ioss::Region *out_region = new Ioss::Region(dbo, name);

  const Ioss::Region * null_in_region = NULL;
  stk::io::define_output_db(*out_region, bulk_data, null_in_region);
  stk::io::write_output_db (*out_region, bulk_data);

  out_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  // Special processing for nodeblock (all nodes in model)...
  const Ioss::Field::RoleType role_type = skin ? Ioss::Field::ATTRIBUTE : Ioss::Field::TRANSIENT;
  stk::io::ioss_add_fields(stk::mesh::MetaData::get(bulk_data).universal_part(), NODE_RANK,
                           out_region->get_node_blocks()[0], role_type, add_all_fields);

  const stk::mesh::PartVector & all_parts = stk::mesh::MetaData::get(bulk_data).get_parts();
  for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {
    stk::mesh::Part * const part = *ip; 

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = out_region->get_entity(part->name());
      if (entity != NULL) {
        if (skin) {
          const bool skin_part = part->name() == "Skin_part";
          if(skin_part) {
            Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
            if (!sset) {
              std::cerr << "ERROR: could not dynamic_cast entity.\n";
              std::exit(EXIT_FAILURE);
            }
            for (unsigned i=0; i < sset->block_count(); i++) {
              Ioss::SideBlock *fb = sset->get_block(i);
              stk::io::ioss_add_fields(*part, part->primary_entity_rank(),
                                       fb, Ioss::Field::TRANSIENT);
            }
          }
        } else if (entity->type() == Ioss::ELEMENTBLOCK) {
          stk::io::ioss_add_fields(*part, part->primary_entity_rank(), 
                                    entity, Ioss::Field::TRANSIENT, add_all_fields);
        } 
      }
    }
  }
  out_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  return out_region;
}


} // namespace

STKUNIT_UNIT_TEST( GearsDemo, skin_gear ) {


  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  const size_t NUM_GEARS = 1;

  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, NUM_GEARS);
  const unsigned p_rank = fixture.bulk_data.parallel_rank();
  std::srand(p_rank); // Seed pseudo-random generator based on processor rank.

  stk::mesh::Part & skin_part = fixture.meta_data.declare_part("Skin_part",fixture.element_rank-1);

  const unsigned ONE_STATE = 1;
  CartesianField & velocity_field = fixture.meta_data.declare_field<CartesianField>("velocity",ONE_STATE);
  CartesianField & displacement = fixture.meta_data.declare_field<CartesianField>("face_displacement",ONE_STATE);
  IntField & processor_field = fixture.meta_data.declare_field<IntField>("processor_id",ONE_STATE);

  stk::mesh::put_field(
      velocity_field,
      NODE_RANK,
      fixture.meta_data.universal_part(),
      fixture.meta_data.spatial_dimension()
      );

  stk::mesh::put_field(
      processor_field,
      fixture.element_rank,
      fixture.meta_data.universal_part()
      );

  // add io parts
  stk::io::put_io_part_attribute( fixture.hex_part);
  stk::io::put_io_part_attribute( fixture.wedge_part);
  stk::io::put_io_part_attribute( skin_part);
  stk::io::set_field_role(fixture.displacement_field.field_of_state(stk::mesh::StateNew), Ioss::Field::TRANSIENT);
  stk::io::set_field_role(displacement,    Ioss::Field::TRANSIENT);
  stk::io::set_field_role(processor_field, Ioss::Field::TRANSIENT);

  std::set<const stk::mesh::Part*> skin_io_parts; skin_io_parts.insert(&skin_part);
  const stk::mesh::PartVector & parts = fixture.meta_data.get_parts();
  for ( stk::mesh::PartVector::const_iterator ip = parts.begin(); ip != parts.end(); ++ip ) {
    stk::mesh::Part & topo_part = **ip;
    if(topo_part.primary_entity_rank() == fixture.element_rank-1 && std::string::npos!=topo_part.name().find("FEM_ROOT_CELL_TOPOLOGY_PART"))  {
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
      stk::mesh::Part & topo_skin_part = fixture.meta_data.declare_part(t, fixture.element_rank-1);
      skin_io_parts.insert(&topo_skin_part);
      stk::io::put_io_part_attribute(topo_skin_part);
      fixture.meta_data.declare_part_subset(topo_part, topo_skin_part);
      fixture.meta_data.declare_part_subset(skin_part, topo_skin_part);
      if (   t == "skin_hex8_quad4_4" || t == "skin_hex20_quad8_5") {
        if  (t == "skin_hex8_quad4_4")   t = "skin_wedge6_quad4_4";
        else                             t = "skin_wedge15_quad4_8";
        stk::mesh::Part & topo_skin_part2 = fixture.meta_data.declare_part(t, fixture.element_rank-1);
        skin_io_parts.insert(&topo_skin_part2);
        stk::io::put_io_part_attribute(topo_skin_part2);
        fixture.meta_data.declare_part_subset(topo_part, topo_skin_part2);
        fixture.meta_data.declare_part_subset(skin_part, topo_skin_part2);
      }
    } 
  }
  stk::mesh::Selector surface_select = fixture.meta_data.locally_owned_part();
  {
    stk::mesh::put_field( displacement, fixture.element_rank-1, skin_part);
    const stk::mesh::PartVector &surf_parts = skin_part.subsets();
    for ( stk::mesh::PartVector::const_iterator ip = surf_parts.begin(); ip != surf_parts.end(); ++ip ) {
      stk::mesh::Part & surf_part = **ip;
      if (surf_part.primary_entity_rank() == fixture.element_rank-1) {
        surface_select |= surf_part;
        stk::mesh::put_field( displacement, fixture.element_rank-1, surf_part);
      }
    }
  }

  fixture.meta_data.commit();

  fixture.generate_mesh();

  populate_processor_id_field_data(fixture,processor_field);

  const size_t spatial_dim = fixture.meta_data.spatial_dimension();
  stk::mesh::skin_mesh( fixture.bulk_data, spatial_dim, &skin_part);

  stk::mesh::EntityVector wedges_to_separate;

  find_and_shuffle_wedges_to_separate( fixture, wedges_to_separate );

  const size_t NUM_TIME_STEPS = 150;
  const size_t separation_interval = 30;

  const double rotation = TWO_PI*4.0/NUM_TIME_STEPS;
  const double x = 0;
  const double y = 0;
  const double z = 0;
  const stk::mesh::fixtures::GearMovement gear_movement_data(rotation,x,y,z);

  Ioss::Region * volume_out_region = NULL;
  Ioss::Region * surface_out_region = NULL;

  stk::mesh::fixtures::Gear & gear = fixture.get_gear(0);

  // Iterate over the time steps, updating the locations of the entities and
  // writing the current mesh state to output files.
  const bool output_exodus_file = true;
  for (size_t time_step = 0; time_step < NUM_TIME_STEPS; ++time_step) {

    // Determine if it's time to separate a wedge
    const bool do_separate_wedge = !wedges_to_separate.empty() && (time_step%separation_interval == 0);

    if (time_step > 0) {

      fixture.bulk_data.update_field_data_states();

      // Move the gear; this will only rotate the gear since x,y,z are all 0
      gear.move(gear_movement_data);
      stk::mesh::Entity * wedge = NULL;

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
    stk::mesh::BucketVector face_buckets;
    stk::mesh::BucketVector all_face_buckets = fixture.bulk_data.buckets(fixture.element_rank-1);
    stk::mesh::get_buckets( surface_select, all_face_buckets, face_buckets);
    for (stk::mesh::BucketVector::iterator b_itr = face_buckets.begin(); b_itr != face_buckets.end(); ++b_itr) {
      stk::mesh::Bucket & b = **b_itr;
      for (size_t i = 0; i < b.size(); ++i) {
        stk::mesh::Entity& face = b[i];
        double *elem_node_disp = field_data(displacement, face);
        if (elem_node_disp) {
          stk::mesh::PairIterRelation node_rels = face.node_relations();
          const size_t num_nodes = node_rels.size();
          for(; !node_rels.empty(); ++node_rels) {
            const stk::mesh::Entity& node = *node_rels->entity();
            double* node_disp = stk::mesh::field_data(fixture.displacement_field, node);
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
        delete volume_out_region;   volume_out_region = NULL;
        delete surface_out_region; surface_out_region = NULL;
        std::ostringstream volume_out_filename;
        volume_out_filename << "volume_mesh_" << std::setw(7) << std::setfill('0') << time_step << ".e";
        volume_out_region = create_output_mesh( volume_out_filename.str(), fixture.bulk_data,  false);
        std::ostringstream surface_out_filename;
        surface_out_filename << "surface_mesh_" << std::setw(7) << std::setfill('0') << time_step << ".e";
        surface_out_region = create_output_mesh( surface_out_filename.str(), fixture.bulk_data,  true);
      }

      stk::io::MeshData mesh;
      mesh.m_output_region=volume_out_region;
      stk::io::process_output_request(mesh, fixture.bulk_data, time_step/60.0, skin_io_parts);
      mesh.m_output_region = NULL;
      mesh.m_output_region=surface_out_region;
      stk::io::process_output_request(mesh, fixture.bulk_data, time_step/60.0);
      mesh.m_output_region = NULL;
    }
  }

  delete volume_out_region;  volume_out_region = NULL;
  delete surface_out_region; surface_out_region = NULL;
}

