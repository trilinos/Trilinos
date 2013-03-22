 /*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_io/MeshReadWriteUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_io/IossBridge.hpp>

#include <Ioss_SubSystem.h>

#include <stk_util/util/tokenize.hpp>
#include <iostream>
#include <sstream>
#include <cmath>

#include <limits>
#include <assert.h>

namespace {
void process_surface_entity(Ioss::SideSet *sset, stk::mesh::fem::FEMMetaData &fem_meta)
{
  assert(sset->type() == Ioss::SIDESET);
  const Ioss::SideBlockContainer& blocks = sset->get_side_blocks();
  stk::io::default_part_processing(blocks, fem_meta);

  stk::mesh::Part* const ss_part = fem_meta.get_part(sset->name());
  assert(ss_part != NULL);

  stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
  bool surface_df_defined = false; // Has the surface df field been defined yet?

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *sb = sset->get_block(i);
    if (stk::io::include_entity(sb)) {
      stk::mesh::Part * const sb_part = fem_meta.get_part(sb->name());
      assert(sb_part != NULL);
      fem_meta.declare_part_subset(*ss_part, *sb_part);

      if (sb->field_exists("distribution_factors")) {
        if (!surface_df_defined) {
          std::string field_name = sset->name() + "_df";
          distribution_factors_field =
            &fem_meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(field_name);
          stk::io::set_field_role(*distribution_factors_field, Ioss::Field::MESH);
          stk::io::set_distribution_factor_field(*ss_part, *distribution_factors_field);
          surface_df_defined = true;
        }
        stk::io::set_distribution_factor_field(*sb_part, *distribution_factors_field);
        int side_node_count = sb->topology()->number_nodes();
        stk::mesh::put_field(*distribution_factors_field,
                             stk::io::part_primary_entity_rank(*sb_part),
                             *sb_part, side_node_count);
      }
    }
  }
}

  size_t get_entities(stk::mesh::Part &part,
		      const stk::mesh::BulkData &bulk,
		      std::vector<stk::mesh::Entity*> &entities,
		      const stk::mesh::Selector *anded_selector)
  {
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(part);
    
    stk::mesh::Selector own = meta.locally_owned_part();
    stk::mesh::Selector selector = part & own;
    if (anded_selector) selector &= *anded_selector;
    
    get_selected_entities(selector, bulk.buckets(type), entities);
    return entities.size();
  }
}

// ========================================================================
template <typename INT>
void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk, INT /*dummy*/)
{
  assert(sset->type() == Ioss::SIDESET);

  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *block = sset->get_block(i);
    if (stk::io::include_entity(block)) {
      std::vector<INT> side_ids ;
      std::vector<INT> elem_side ;

      stk::mesh::Part * const sb_part = fem_meta.get_part(block->name());
      stk::mesh::EntityRank elem_rank = fem_meta.element_rank();

      block->get_field_data("ids", side_ids);
      block->get_field_data("element_side", elem_side);

      assert(side_ids.size() * 2 == elem_side.size());
      stk::mesh::PartVector add_parts( 1 , sb_part );

      size_t side_count = side_ids.size();
      std::vector<stk::mesh::Entity*> sides(side_count);
      for(size_t is=0; is<side_count; ++is) {
        stk::mesh::Entity* const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

        // If NULL, then the element was probably assigned to an
        // element block that appears in the database, but was
        // subsetted out of the analysis mesh. Only process if
        // non-null.
        if (elem != NULL) {
          // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
          int side_ordinal = elem_side[is*2+1] - 1;

          stk::mesh::Entity* side_ptr = NULL;
          side_ptr = &stk::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);
          stk::mesh::Entity& side = *side_ptr;

          bulk.change_entity_parts( side, add_parts );
          sides[is] = &side;
        } else {
          sides[is] = NULL;
        }
      }

      const stk::mesh::Field<double, stk::mesh::ElementNode> *df_field =
        stk::io::get_distribution_factor_field(*sb_part);
      if (df_field != NULL) {
        stk::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      block->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(field, sides, block, *I);
      }
    }
  }
}

void process_surface_entity(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(sset) == 4)
    process_surface_entity(sset, bulk, (int)0);
  else
    process_surface_entity(sset, bulk, (int64_t)0);
}

void process_nodeblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  assert(nb->field_exists("mesh_model_coordinates"));
  Ioss::Field coordinates = nb->get_field("mesh_model_coordinates");
  int spatial_dim = coordinates.transformed_storage()->component_count();

  stk::mesh::Field<double,stk::mesh::Cartesian> & coord_field =
    fem_meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

  stk::mesh::put_field( coord_field, fem_meta.node_rank(), fem_meta.universal_part(), spatial_dim);
  stk::io::define_io_fields(nb, Ioss::Field::ATTRIBUTE, fem_meta.universal_part(), 0);
}

void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  std::vector<stk::mesh::Entity*> nodes;
  stk::io::get_entity_list(nb, fem_meta.node_rank(), bulk, nodes);

  stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
    fem_meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

  stk::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");

  // Add all attributes as fields.
  // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
  // named components of the 'attribute' field, so add them instead.
  Ioss::NameList names;
  nb->field_describe(Ioss::Field::ATTRIBUTE, &names);
  for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    if(*I == "attribute" && names.size() > 1)
      continue;
    stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase> (*I);
    if (field)
      stk::io::field_data_from_ioss(field, nodes, nb, *I);
  }
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, fem_meta);
}

template <typename INT>
void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  const stk::mesh::fem::FEMMetaData& fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = fem_meta.get_part(name);
      assert(part != NULL);

      const CellTopologyData* cell_topo = stk::io::get_cell_topology(*part);
      if (cell_topo == NULL) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " returned NULL from get_cell_topology()";
        throw std::runtime_error( msg.str() );
      }

      std::vector<INT> elem_ids ;
      std::vector<INT> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = cell_topo->node_count ;

      std::vector<stk::mesh::EntityId> id_vec(nodes_per_elem);
      std::vector<stk::mesh::Entity*> elements(element_count);

      for(size_t i=0; i<element_count; ++i) {
        INT *conn = &connectivity[i*nodes_per_elem];
        std::copy(&conn[0], &conn[0+nodes_per_elem], id_vec.begin());
        elements[i] = &stk::mesh::fem::declare_element(bulk, *part, elem_ids[i], &id_vec[0]);
      }

      // Add all element attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(field, elements, entity, *I);
      }
    }
  }
}

// ========================================================================
// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, fem_meta);

  stk::mesh::Field<double> & distribution_factors_field =
    fem_meta.declare_field<stk::mesh::Field<double> >("distribution_factors");
  stk::io::set_field_role(distribution_factors_field, Ioss::Field::MESH);

  /** \todo REFACTOR How to associate distribution_factors field
   * with the nodeset part if a node is a member of multiple
   * nodesets
   */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = fem_meta.get_part(entity->name());
      assert(part != NULL);
      assert(entity->field_exists("distribution_factors"));

      stk::mesh::put_field(distribution_factors_field, fem_meta.node_rank(), *part);
    }
  }
}

// ========================================================================
// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  stk::io::default_part_processing(side_sets, fem_meta);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, fem_meta);
    }
  }
}

// ========================================================================
template <typename INT>
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk, INT /*dummy*/)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      stk::mesh::Part* const part = fem_meta.get_part(name);
      assert(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<INT> node_ids ;
      size_t node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity*> nodes(node_count);
      stk::mesh::EntityRank n_rank = fem_meta.node_rank();
      for(size_t i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity(n_rank, node_ids[i] );
        if (nodes[i] != NULL)
          bulk.declare_entity(n_rank, node_ids[i], add_parts );
      }

      stk::mesh::Field<double> *df_field =
        fem_meta.get_field<stk::mesh::Field<double> >("distribution_factors");

      if (df_field != NULL) {
        stk::io::field_data_from_ioss(df_field, nodes, entity, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(field, nodes, entity, *I);
      }
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, bulk);
    }
  }
}

// ========================================================================
void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    Ioss::Field::RoleType filter_role,
		    const stk::mesh::Selector *anded_selector=NULL)
{
  std::vector<stk::mesh::Entity*> entities;
  if (io_entity->type() == Ioss::SIDEBLOCK) {
    // Temporary Kluge to handle sideblocks which contain internally generated sides
    // where the "ids" field on the io_entity doesn't work to get the correct side...
    // NOTE: Could use this method for all entity types, but then need to correctly 
    // specify whether shared entities are included/excluded (See IossBridge version).
    size_t num_sides = get_entities(part, bulk, entities, anded_selector);
    if (num_sides != (size_t)io_entity->get_property("entity_count").get_int()) {
      std::ostringstream msg ;
      msg << " INTERNAL_ERROR: Number of sides on part " << part.name() << " (" << num_sides
	  << ") does not match number of sides in the associated Ioss SideBlock named "
	  << io_entity->name() << " (" << io_entity->get_property("entity_count").get_int()
	  << ").";
        throw std::runtime_error( msg.str() );
    }
  } else {
    stk::io::get_entity_list(io_entity, part_type, bulk, entities);
  }

  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);
  const std::vector<stk::mesh::FieldBase*> &fields = fem_meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    stk::io::field_data_to_ioss(f, entities, io_entity, f->name(), filter_role);
  }
}

void internal_process_output_request(stk::io::MeshData &mesh_data,
                                     stk::mesh::BulkData &bulk,
                                     int step,
                                     const std::set<const stk::mesh::Part*> &exclude)
{
  Ioss::Region *region = mesh_data.m_output_region;
  region->begin_state(step);
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  // Special processing for nodeblock (all nodes in model)...
  put_field_data(bulk, fem_meta.universal_part(), fem_meta.node_rank(),
                 region->get_node_blocks()[0], Ioss::Field::Field::TRANSIENT,
		 mesh_data.m_anded_selector);

  // Now handle all non-nodeblock parts...
  const stk::mesh::PartVector & all_parts = fem_meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part) && !exclude.count(part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = region->get_entity(part->name());
      if (entity != NULL && entity->type() != Ioss::SIDESET) {
        put_field_data(bulk, *part, stk::io::part_primary_entity_rank(*part),
                       entity, Ioss::Field::Field::TRANSIENT,
		       mesh_data.m_anded_selector);
      }
    }
  }
  region->end_state(step);
}

namespace stk {
namespace io {

MeshData::~MeshData()
{
  delete m_input_region;
  delete m_output_region;
}

void show_mesh_help()
{
  std::cerr << "Options are:\n"
            << "\n"
            << "filename -- specify the name of the file from which to read the\n"
            << "            mesh file. If the --directory option is specified, it will be\n"
            << "            prepended to the filename unless the filename specifies an absolute path.\n"
            << "\n"
            << "gen:NxMxL -- internally generate a hex mesh of size N by M by L\n"
            << "             intervals. See 'Generated Options' below for more options.\n"
            << "\n"
            << "Generated Options:\n"
            << "shell:xXyYzZ\n"
            << "The argument specifies whether there is a shell block\n"
            << "at the location. 'x' is minX, 'X' is maxX, etc.\n"
            << "\n"
            << "help -- no argument, shows valid options\n"
            << "\n"
            << "show -- no argument, prints out a summary of the settings used to\n"
            << "generate the mesh. The output will look similar to:\n"
            << "    \"10x12x8|shell:xX|bbox:-10,-10,-10,10,10,10|show\"\n"
            << "\n"
            << "    Mesh Parameters:\n"
            << "\tIntervals: 10 by 12 by 8\n"
            << "\tX = 2       * (0..10) + -10     Range: -10 <= X <= 10\n"
            << "\tY = 1.66667 * (0..12) + -10     Range: -10 <= Y <= 10\n"
            << "\tZ = 2.5     * (0..8)  + -10     Range: -10 <= Z <= 10\n"
            << "\tNode Count (total)    = 1287\n"
            << "\tElement Count (total) = 1152\n"
            << "\tBlock Count           = 3\n"
            << "\n"
            << "shell:xXyYzZ \n"
            << "which specifies whether there is a shell block at that\n"
            << "location. 'x' is minimum x face, 'X' is maximum x face,\n"
            << "similarly for y and z.  Note that the argument string is a\n"
            << "single multicharacter string.  You can add multiple shell blocks\n"
            << "to a face, for example, shell:xxx would add three layered shell\n"
            << "blocks on the minimum x face.  An error is output if a non\n"
            << "xXyYzZ character is found, but execution continues.\n"
            << "\n"
            << "zdecomp:n0 n1,n2,...,n#proc-1\n"
            << "which are the number of intervals in the z direction for each\n"
            << "processor in a pallel run.  If this option is specified, then\n"
            << "the total number of intervals in the z direction is the sum of\n"
            << "the n0, n1, ... An interval count must be specified for each\n"
            << "processor.  If this option is not specified, then the number of\n"
            << "intervals on each processor in the z direction is numZ/numProc\n"
            << "with the extras added to the lower numbered processors.\n"
            << "\n"
            << "scale:xs,ys,zs\n"
            << "which are the scale factors in the x, y, and z directions. All\n"
            << "three must be specified if this option is present.\n"
            << "\n"
            << "- offset -- argument = xoff, yoff, zoff which are the offsets in the\n"
            << "x, y, and z directions.  All three must be specified if this option\n"
            << "is present.\n"
            << "\n"
            << "- bbox -- argument = xmin, ymin, zmin, xmax, ymax, zmax\n"
            << "which specify the lower left and upper right corners of\n"
            << "the bounding box for the generated mesh.  This will\n"
            << "calculate the scale and offset which will fit the mesh in\n"
            << "the specified box.  All calculations are based on the currently\n"
            << "active interval settings. If scale or offset or zdecomp\n"
            << "specified later in the option list, you may not get the\n"
            << "desired bounding box.\n"
            << "\n"
            << "- rotate -- argument = axis,angle,axis,angle,...\n"
            << "where axis is 'x', 'y', or 'z' and angle is the rotation angle in\n"
            << "degrees. Multiple rotations are cumulative. The composite rotation\n"
            << "matrix is applied at the time the coordinates are retrieved after\n"
            << "scaling and offset are applied.\n"
            << "\n"
            << "The unrotated coordinate of a node at grid location i,j,k is:\n"
            << "\n"
            << "\tx = x_scale * i + x_off,\n"
            << "\ty = z_scale * j + y_off,\n"
            << "\tz = z_scale * k + z_off,\n"
            << "\n"
            << "The extent of the unrotated mesh will be:\n"
            << "\n"
            << "\tx_off <= x <= x_scale * numX + x_off\n"
            << "\ty_off <= y <= y_scale * numY + y_off\n"
            << "\tz_off <= z <= z_scale * numZ + z_off\n"
            << "\n"
            << "If an unrecognized option is specified, an error message will be\n"
            << "output and execution will continue.\n"
            << "\n"
            << "An example of valid input is:\n"
            << "\n"
            << "\t\"10x20x40|scale:1,0.5,0.25|offset:-5,-5,-5|shell:xX\"\n"
            << "\n"
            << "\n"
            << "This would create a mesh with 10 intervals in x, 20 in y, 40 in z\n"
            << "The mesh would be centered on 0,0,0 with a range of 10 in each\n"
            << "direction. There would be a shell layer on the min and max\n"
            << "x faces.\n"
            << "\n"
            << "NOTE: All options are processed in the order they appear in\n"
            << "the parameters string (except rotate which is applied at the\n"
            << "time the coordinates are generated/retrieved)\n"
            << "\n";
}

void create_input_mesh(const std::string &mesh_type,
                       const std::string &mesh_filename,
                       stk::ParallelMachine comm,
                       stk::mesh::fem::FEMMetaData &fem_meta,
                       stk::io::MeshData &mesh_data)
{
  Ioss::Region *in_region = mesh_data.m_input_region;
  if (in_region == NULL) {
	// If in_region is NULL, then open the file;
	// If in_region is non-NULL, then user has given us a valid Ioss::Region that
	// should be used.
	Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, mesh_filename,
                                                    Ioss::READ_MODEL, comm,
                                                    mesh_data.m_property_manager);
	if (dbi == NULL || !dbi->ok()) {
	  std::cerr  << "ERROR: Could not open database '" << mesh_filename
                 << "' of type '" << mesh_type << "'\n";
	  Ioss::NameList db_types;
	  Ioss::IOFactory::describe(&db_types);
	  std::cerr << "\nSupported database types:\n\t";
	  for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
	    std::cerr << *IF << "  ";
	  }
	  std::cerr << "\n\n";
	}

	// NOTE: 'in_region' owns 'dbi' pointer at this time...
	in_region = new Ioss::Region(dbi, "input_model");
	mesh_data.m_input_region = in_region;
  }

  size_t spatial_dimension = in_region->get_property("spatial_dimension").get_int();
  initialize_spatial_dimension(fem_meta, spatial_dimension, stk::mesh::fem::entity_rank_names(spatial_dimension));

  process_elementblocks(*in_region, fem_meta);
  process_nodeblocks(*in_region,    fem_meta);
  process_sidesets(*in_region,      fem_meta);
  process_nodesets(*in_region,      fem_meta);
}


void create_output_mesh(const std::string &filename,
                        stk::ParallelMachine comm,
                        stk::mesh::BulkData &bulk_data,
                        MeshData &mesh_data)
{
  Ioss::Region *out_region = NULL;

  std::string out_filename = filename;
  if (filename.empty()) {
	out_filename = "default_output_mesh";
  } else {
	// These filenames may be coming from the generated options which
	// may have forms similar to: "2x2x1|size:.05|height:-0.1,1"
	// Strip the name at the first "+:|," character:
	std::vector<std::string> tokens;
	stk::util::tokenize(out_filename, "+|:,", tokens);
	out_filename = tokens[0];
  }

  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", out_filename,
                                                  Ioss::WRITE_RESULTS,
                                                  comm, mesh_data.m_property_manager);
  if (dbo == NULL || !dbo->ok()) {
	std::cerr << "ERROR: Could not open results database '" << out_filename
              << "' of type 'exodusII'\n";
	std::exit(EXIT_FAILURE);
  }

  // NOTE: 'out_region' owns 'dbo' pointer at this time...
  out_region = new Ioss::Region(dbo, "results_output");

  stk::io::define_output_db(*out_region, bulk_data, mesh_data.m_input_region, mesh_data.m_anded_selector);
  stk::io::write_output_db(*out_region,  bulk_data, mesh_data.m_anded_selector);
  mesh_data.m_output_region = out_region;
}

// ========================================================================
int process_output_request(MeshData &mesh_data,
                           stk::mesh::BulkData &bulk,
                           double time,
                           const std::set<const stk::mesh::Part*> &exclude)
{
  Ioss::Region *region = mesh_data.m_output_region;
  region->begin_mode(Ioss::STATE_TRANSIENT);

  int out_step = region->add_state(time);
  internal_process_output_request(mesh_data, bulk, out_step,exclude);

  region->end_mode(Ioss::STATE_TRANSIENT);

  return out_step;
}

// ========================================================================
void populate_bulk_data(stk::mesh::BulkData &bulk_data,
                        MeshData &mesh_data)
{
  Ioss::Region *region = mesh_data.m_input_region;

  if (region) {

    bulk_data.modification_begin();
    process_mesh_bulk_data(region, bulk_data);
    bulk_data.modification_end();

  } else {
	std::cerr << "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_bulk_data.\n";
	std::exit(EXIT_FAILURE);
  }
}

// ========================================================================
void process_mesh_bulk_data(Ioss::Region *region, stk::mesh::BulkData &bulk_data)
{
    bool ints64bit = db_api_int_size(region) == 8;
    if (ints64bit) {
      int64_t zero = 0;
      process_elementblocks(*region, bulk_data, zero);
      process_nodeblocks(*region,    bulk_data);
      process_nodesets(*region,      bulk_data, zero);
      process_sidesets(*region,      bulk_data);
    } else {
      int zero = 0;
      process_elementblocks(*region, bulk_data, zero);
      process_nodeblocks(*region,    bulk_data);
      process_nodesets(*region,      bulk_data, zero);
      process_sidesets(*region,      bulk_data);
    }
}

namespace {
// ========================================================================
// Transfer transient field data from mesh file for io_entity to
// the corresponding stk_mesh entities If there is a stk_mesh
// field with the same name as the database field.
// Assumes that mesh is positioned at the correct state for reading.
void internal_process_input_request(Ioss::GroupingEntity *io_entity,
                                    stk::mesh::EntityRank entity_rank,
                                    stk::mesh::BulkData &bulk)
{
  assert(io_entity != NULL);
  std::vector<stk::mesh::Entity*> entity_list;
  stk::io::get_entity_list(io_entity, entity_rank, bulk, entity_list);

  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);

  Ioss::NameList names;
  io_entity->field_describe(Ioss::Field::TRANSIENT, &names);
  for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase>(*I);
    if (field) {
      stk::io::field_data_from_ioss(field, entity_list, io_entity, *I);
    }
  }
}

void input_nodeblock_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];
  internal_process_input_request(nb, fem_meta.node_rank(), bulk);
}

void input_elementblock_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(size_t i=0; i < elem_blocks.size(); i++) {
    if (stk::io::include_entity(elem_blocks[i])) {
      internal_process_input_request(elem_blocks[i], fem_meta.element_rank(), bulk);
    }
  }
}

void input_nodeset_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);
  const Ioss::NodeSetContainer& nodesets = region.get_nodesets();
  for(size_t i=0; i < nodesets.size(); i++) {
    if (stk::io::include_entity(nodesets[i])) {
      internal_process_input_request(nodesets[i], fem_meta.node_rank(), bulk);
    }
  }
}

void input_sideset_fields(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk);
  if (fem_meta.spatial_dimension() <= fem_meta.side_rank())
    return;

  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;
    if (stk::io::include_entity(entity)) {
      const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
      for(size_t i=0; i < blocks.size(); i++) {
        if (stk::io::include_entity(blocks[i])) {
          internal_process_input_request(blocks[i], fem_meta.side_rank(), bulk);
        }
      }
    }
  }
}

void define_input_nodeblock_fields(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];
  stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT,
                            fem_meta.universal_part(), fem_meta.node_rank());
}

void define_input_elementblock_fields(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(size_t i=0; i < elem_blocks.size(); i++) {
    if (stk::io::include_entity(elem_blocks[i])) {
      stk::mesh::Part* const part = fem_meta.get_part(elem_blocks[i]->name());
      assert(part != NULL);
      stk::io::define_io_fields(elem_blocks[i], Ioss::Field::TRANSIENT,
                                *part, part_primary_entity_rank(*part));
    }
  }
}

void define_input_nodeset_fields(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  const Ioss::NodeSetContainer& nodesets = region.get_nodesets();
  for(size_t i=0; i < nodesets.size(); i++) {
    if (stk::io::include_entity(nodesets[i])) {
      stk::mesh::Part* const part = fem_meta.get_part(nodesets[i]->name());
      assert(part != NULL);
      stk::io::define_io_fields(nodesets[i], Ioss::Field::TRANSIENT,
                                *part, part_primary_entity_rank(*part));
    }
  }
}

void define_input_sideset_fields(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta)
{
  if (fem_meta.spatial_dimension() <= fem_meta.side_rank())
    return;

  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;
    if (stk::io::include_entity(entity)) {
      const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
      for(size_t i=0; i < blocks.size(); i++) {
        if (stk::io::include_entity(blocks[i])) {
          stk::mesh::Part* const part = fem_meta.get_part(blocks[i]->name());
          assert(part != NULL);
          stk::io::define_io_fields(blocks[i], Ioss::Field::TRANSIENT,
                                    *part, part_primary_entity_rank(*part));
        }
      }
    }
  }
}

}

// ========================================================================
// Iterate over all Ioss entities in the input mesh database and
// define a stk_field for all transient fields found.  The stk
// field will have the same name as the field on the database.
//
// Note that all fields found on the database will have a
// corresponding stk field defined.  If you want just a selected
// subset of the defined fields, you will need to define the
// fields manually.
//
// To populate the stk field with data from the database, call
// process_input_request().
void define_input_fields(MeshData &mesh_data, stk::mesh::fem::FEMMetaData &fem_meta)
{
  Ioss::Region *region = mesh_data.m_input_region;
  if (region) {
	define_input_nodeblock_fields(*region, fem_meta);
	define_input_elementblock_fields(*region, fem_meta);
	define_input_nodeset_fields(*region, fem_meta);
	define_input_sideset_fields(*region, fem_meta);
  } else {
	std::cerr << "INTERNAL ERROR: Mesh Input Region pointer is NULL in process_input_request.\n";
	std::exit(EXIT_FAILURE);
  }
}

// ========================================================================
// Iterate over all fields defined in the stk mesh data structure.
// If the field has the io_attribute set, then define that field
// on the corresponding io entity on the output mesh database.
// The database field will have the same name as the stk field.
//
// To export the data to the database, call
// process_output_request().

void define_output_fields(const MeshData &mesh_data, const stk::mesh::fem::FEMMetaData &fem_meta,
                          bool add_all_fields)
{
  Ioss::Region *region = mesh_data.m_output_region;
  if (region) {
	region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

	// Special processing for nodeblock (all nodes in model)...
	stk::io::ioss_add_fields(fem_meta.universal_part(), fem_meta.node_rank(),
                             region->get_node_blocks()[0],
                             Ioss::Field::TRANSIENT, add_all_fields);

	const stk::mesh::PartVector & all_parts = fem_meta.get_parts();
	for ( stk::mesh::PartVector::const_iterator
            ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

	  stk::mesh::Part * const part = *ip;

	  // Check whether this part should be output to results database.
	  if (stk::io::is_part_io_part(*part)) {
	    // Get Ioss::GroupingEntity corresponding to this part...
	    Ioss::GroupingEntity *entity = region->get_entity(part->name());
	    if (entity != NULL) {
	      stk::io::ioss_add_fields(*part, part_primary_entity_rank(*part),
                                   entity, Ioss::Field::TRANSIENT, add_all_fields);
	    }
	  }
	}
	region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  } else {
	std::cerr << "INTERNAL ERROR: Mesh Input Region pointer is NULL in process_input_request.\n";
	std::exit(EXIT_FAILURE);
  }
}
// ========================================================================
void process_input_request(MeshData &mesh_data, stk::mesh::BulkData &bulk, double time)
{
  // Find the step on the database with time closest to the requested time...
  Ioss::Region *region = mesh_data.m_input_region;
  int step_count = region->get_property("state_count").get_int();
  double delta_min = 1.0e30;
  int    step_min  = 0;
  for (int istep = 0; istep < step_count; istep++) {
	double state_time = region->get_state_time(istep+1);
	double delta = state_time - time;
	if (delta < 0.0) delta = -delta;
	if (delta < delta_min) {
	  delta_min = delta;
	  step_min  = istep;
	  if (delta == 0.0) break;
	}
  }
  // Exodus steps are 1-based;
  process_input_request(mesh_data, bulk, step_min+1);
}

void process_input_request(MeshData &mesh_data,
                           stk::mesh::BulkData &bulk,
                           int step)
{
  if (step <= 0)
	return;

  Ioss::Region *region = mesh_data.m_input_region;
  if (region) {
	bulk.modification_begin();

  input_mesh_fields(region, bulk, step);

	bulk.modification_end();

  } else {
	std::cerr << "INTERNAL ERROR: Mesh Input Region pointer is NULL in process_input_request.\n";
	std::exit(EXIT_FAILURE);
  }
}

void input_mesh_fields(Ioss::Region *region, stk::mesh::BulkData &bulk,
                           double time)
{
  // Find the step on the database with time closest to the requested time...
  int step_count = region->get_property("state_count").get_int();
  double delta_min = 1.0e30;
  int    step_min  = 0;
  for (int istep = 0; istep < step_count; istep++) {
	double state_time = region->get_state_time(istep+1);
	double delta = state_time - time;
	if (delta < 0.0) delta = -delta;
	if (delta < delta_min) {
	  delta_min = delta;
	  step_min  = istep;
	  if (delta == 0.0) break;
	}
  }
  // Exodus steps are 1-based;
  input_mesh_fields(region, bulk, step_min+1);
}

void input_mesh_fields(Ioss::Region *region, stk::mesh::BulkData &bulk,
                           int step)
{
	// Pick which time index to read into solution field.
	region->begin_state(step);

	input_nodeblock_fields(*region, bulk);
	input_elementblock_fields(*region, bulk);
	input_nodeset_fields(*region, bulk);
	input_sideset_fields(*region, bulk);

	region->end_state(step);
}

// ========================================================================
template <typename INT>
void get_element_block_sizes(MeshData &mesh_data,
                             std::vector<INT>& el_blocks)
{
  Ioss::Region *io = mesh_data.m_input_region;
  const Ioss::ElementBlockContainer& elem_blocks = io->get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin(); it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;
    if (stk::io::include_entity(entity)) {
      el_blocks.push_back(entity->get_property("entity_count").get_int());
    }
  }
}
template void get_element_block_sizes(MeshData &mesh_data, std::vector<int>& el_blocks);
template void get_element_block_sizes(MeshData &mesh_data, std::vector<int64_t>& el_blocks);
} // namespace io
} // namespace stk
