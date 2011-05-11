/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/parallel/Parallel.hpp>

#include <stk_io/util/UseCase_mesh.hpp>
#include <stk_io/util/Gears.hpp>
#include <stk_io/util/Skinning.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_util/util/tokenize.hpp>
#include <iostream>
#include <sstream>
#include <cmath>

#include <limits>

namespace {
  void generate_gears(stk::ParallelMachine comm,
		      const std::string &parameters,
		      stk::mesh::fem::FEMMetaData &fem_meta,
		      std::vector<stk::io::util::Gear*> &gears);

  void generate_gears(stk::mesh::BulkData &mesh,
		      std::vector<stk::io::util::Gear*> &gears);

  void process_surface_entity(Ioss::SideSet *entity,
			      stk::mesh::MetaData &meta,
			      stk::mesh::EntityRank)
  {
    assert(entity->type() == Ioss::SIDESET);
    const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
    stk::io::default_part_processing(blocks, meta, 0);

    stk::mesh::Part* const fs_part = meta.get_part(entity->name());
    assert(fs_part != NULL);

    stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
    bool surface_df_defined = false; // Has the surface df field been defined yet?


    size_t block_count = entity->block_count();
    for (size_t i=0; i < block_count; i++) {
      Ioss::EntityBlock *fb = entity->get_block(i);
      if (stk::io::include_entity(fb)) {
	stk::mesh::Part * const fb_part = meta.get_part(fb->name());
	assert(fb_part != NULL);
	meta.declare_part_subset(*fs_part, *fb_part);

	if (fb->field_exists("distribution_factors")) {
	  if (!surface_df_defined) {
	    std::string field_name = entity->name() + "_distribution_factors";
	    distribution_factors_field =
	      &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(field_name);
	    stk::io::set_distribution_factor_field(*fs_part, *distribution_factors_field);
	    surface_df_defined = true;
	  }
	  stk::io::set_distribution_factor_field(*fb_part, *distribution_factors_field);
	  int side_node_count = fb->topology()->number_nodes();
	  stk::mesh::put_field(*distribution_factors_field,
			       stk::io::part_primary_entity_rank(*fb_part),
			       *fb_part, side_node_count);
	}

	/** \todo IMPLEMENT truly handle fields... For this case we
	 * are just defining a field for each transient field that is
	 * present in the mesh...
	 */
	stk::io::define_io_fields(fb, Ioss::Field::TRANSIENT,
				  *fb_part, stk::io::part_primary_entity_rank(*fb_part));
      }
    }
  }

  // ========================================================================
  void process_surface_entity(const Ioss::SideSet* sset ,
			      stk::mesh::BulkData & bulk)
  {
    assert(sset->type() == Ioss::SIDESET);
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

    size_t block_count = sset->block_count();
    for (size_t i=0; i < block_count; i++) {
      Ioss::EntityBlock *block = sset->get_block(i);
      if (stk::io::include_entity(block)) {
	std::vector<int> side_ids ;
	std::vector<int> elem_side ;

	stk::mesh::Part * const fb_part = meta.get_part(block->name());
	stk::mesh::EntityRank elem_rank = stk::io::element_rank(meta);

	block->get_field_data("ids", side_ids);
	block->get_field_data("element_side", elem_side);

	assert(side_ids.size() * 2 == elem_side.size());
	stk::mesh::PartVector add_parts( 1 , fb_part );

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

            stk::mesh::Entity &side = stk::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);

	    bulk.change_entity_parts( side, add_parts );
	    sides[is] = &side;
	  } else {
	    sides[is] = NULL;
	  }
	}

	const stk::mesh::Field<double, stk::mesh::ElementNode> *df_field =
	  stk::io::get_distribution_factor_field(*fb_part);
	if (df_field != NULL) {
	  stk::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
	}
      }
    }
  }


}

namespace stk {
namespace io {
namespace util {

MeshData::~MeshData()
{
  delete m_region;
  for (size_t i=0; i < m_gears.size(); i++) {
    delete m_gears[i];
  }
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
            << "gears:IxJxK -- internally generate the gears mesh. See 'Gear Options' below for more options.\n"
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
            << "\n"
            << "Gear Options:\n"
            << "gears:IxJxK -- internally generate the gears mesh with 'I' gears in\n"
            << "the X direction, 'J' gears in the Y direction and 'K' gears in the Z\n"
            << "direction. \n"
            << "\n"
            << "The parameters should be of the form:  \"IxJxK|option:param,param,param|option:a,b,c\"\n"
            << "Each \"|\" or \"+\" separated section of the parameters is a \"group\"\n"
            << "Each group is then split into options and params\n";
}

void create_input_mesh(const std::string &mesh_type,
                       const std::string &mesh_filename,
                       const std::string &working_directory,
                       stk::ParallelMachine comm,
                       stk::mesh::fem::FEMMetaData &fem_meta,
                       stk::io::util::MeshData &mesh_data,
                       bool hex_only)
{
  Ioss::Region *in_region = NULL;
  stk::mesh::MetaData &meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  if (mesh_type == "exodusii" || mesh_type == "generated" || mesh_type == "pamgen" ) {

    // Prepend the working directory onto the mesh filename iff the
    // directory was specified *and* the mesh filename does not
    // specify an absolute path *and* the type is not "generated"
    std::string filename = mesh_filename;
    if (mesh_filename[0] != '/' && !working_directory.empty() && mesh_type != "generated") {
      filename = working_directory + mesh_filename;
    }

    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, filename,
                                                    Ioss::READ_MODEL,
                                                    comm);
    if (dbi == NULL || !dbi->ok()) {
      std::cerr  << "ERROR: Could not open database '" << filename
                 << "' of type '" << mesh_type << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'in_region' owns 'dbi' pointer at this time...
    in_region = new Ioss::Region(dbi, "input_model");
    size_t spatial_dimension = in_region->get_property("spatial_dimension").get_int();
    initialize_spatial_dimension(meta_data, spatial_dimension, stk::mesh::fem::entity_rank_names(spatial_dimension));

    // Filter out all non-hex8 element blocks...
    if (hex_only) {
      const Ioss::ElementBlockContainer& elem_blocks = in_region->get_element_blocks();
      for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
          it != elem_blocks.end(); ++it) {
        Ioss::ElementBlock *entity = *it;
        std::string name = entity->topology()->name();
        if (name != "hex8") {
          entity->property_add(Ioss::Property(std::string("omitted"), 1));
        }
      }
    }
    process_elementblocks(*in_region, meta_data);
    process_nodeblocks(*in_region,    meta_data);
    process_sidesets(*in_region,      meta_data);
    process_nodesets(*in_region,      meta_data);

  }

  else if (mesh_type == "gears") {
    generate_gears(comm, mesh_filename, fem_meta, mesh_data.m_gears);
  }
  mesh_data.m_region = in_region;

  // NOTE: THIS SHOULD NOT BE USED; USE STK_MESH SKINNING INSTEAD
  // See if caller requested that the model be "skinned".  If
  // so, all exposed sides are generated and put in a part named
  // "skin"
  if (mesh_data.m_generateSkinFaces) {
    stk::mesh::Part &skin_part = fem_meta.declare_part("skin", fem_meta.side_rank());
    stk::io::put_io_part_attribute(skin_part);
    /** \todo REFACTOR Query all parts to determine topology of the skin. */
    stk::mesh::fem::set_cell_topology(skin_part, shards::getCellTopologyData<shards::Quadrilateral<4> >());
  }
}


Ioss::Region *create_output_mesh(const std::string &mesh_filename,
                                 const std::string &mesh_extension,
                                 const std::string &working_directory,
                                 MPI_Comm comm,
                                 stk::mesh::BulkData &bulk_data,
                                 const Ioss::Region *in_region,
                                 stk::mesh::fem::FEMMetaData &fem_meta,
                                 bool add_transient ,
                                 bool add_all_fields ) {
  return create_output_mesh(mesh_filename,
                            mesh_extension,
                            working_directory,
                            comm,
                            bulk_data,
                            in_region,
                            stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta),
                            add_transient,
                            add_all_fields);
}




void create_output_mesh(const std::string &mesh_filename,
                        const std::string &mesh_extension,
                        const std::string &working_directory,
                        stk::ParallelMachine comm,
                        stk::mesh::BulkData &bulk_data,
                        stk::mesh::fem::FEMMetaData &fem_meta,
                        MeshData &mesh_data,
                        bool add_transient,
                        bool add_all_fields)
{
  mesh_data.m_region = create_output_mesh(mesh_filename,
                                          mesh_extension, working_directory,
                                          comm, bulk_data,
                                          NULL, //mesh_data.m_region,
                                          fem_meta, add_transient, add_all_fields);
}
// ========================================================================

Ioss::Region *create_output_mesh(const std::string &mesh_filename,
                                 const std::string &mesh_extension,
                                 const std::string &working_directory,
                                 stk::ParallelMachine comm,
                                 stk::mesh::BulkData &bulk_data,
                                 const Ioss::Region *in_region,
                                 stk::mesh::MetaData &meta_data,
                                 bool add_transient,
                                 bool add_all_fields)
{
  std::string filename = mesh_filename;
  if (filename.empty()) {
    filename = "usecase_mesh";
  } else {
    // These filenames may be coming from the generated or gears options which
    // may have forms similar to: "2x2x1|size:.05|height:-0.1,1"
    // Strip the name at the first "+:|," character:
    std::vector<std::string> tokens;
    stk::util::tokenize(filename, "+|:,", tokens);
    filename = tokens[0];
  }

  std::string out_filename = filename;
  if (!mesh_extension.empty())
    out_filename += "." + mesh_extension;

  // Prepend the working directory onto the mesh filename iff the
  // directory was specified *and* the mesh filename does not
  // specify an absolute path.
  if (out_filename[0] != '/' && !working_directory.empty() > 0) {
    out_filename = working_directory + out_filename;
  }

  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", out_filename,
                                                  Ioss::WRITE_RESULTS,
                                                  comm);
  if (dbo == NULL || !dbo->ok()) {
    std::cerr << "ERROR: Could not open results database '" << out_filename
              << "' of type 'exodusII'\n";
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'out_region' owns 'dbo' pointer at this time...
  Ioss::Region *out_region = new Ioss::Region(dbo, "results_output");

  stk::io::define_output_db(*out_region, bulk_data, in_region);
  stk::io::write_output_db(*out_region,  bulk_data);

  if (add_transient) {
    out_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // Special processing for nodeblock (all nodes in model)...
    stk::io::ioss_add_fields(meta_data.universal_part(), node_rank(meta_data),
                             out_region->get_node_blocks()[0],
                             Ioss::Field::TRANSIENT, add_all_fields);

    const stk::mesh::PartVector & all_parts = meta_data.get_parts();
    for ( stk::mesh::PartVector::const_iterator
            ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk::mesh::Part * const part = *ip;

      // Check whether this part should be output to results database.
      if (stk::io::is_part_io_part(*part)) {
        // Get Ioss::GroupingEntity corresponding to this part...
        Ioss::GroupingEntity *entity = out_region->get_entity(part->name());
        if (entity != NULL) {
          if (entity->type() == Ioss::ELEMENTBLOCK) {
            stk::io::ioss_add_fields(*part, part_primary_entity_rank(*part),
                                     entity, Ioss::Field::TRANSIENT, add_all_fields);
          }
        }
      }
    }
    out_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }
  return out_region;
}

// ========================================================================

void process_nodeblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta) {
  process_nodeblocks(region,  stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta));
}
void process_elementblocks(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta) {
  process_elementblocks(region, stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta));
}
void process_sidesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta) {
  process_sidesets(region, stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta));
}
void process_nodesets(Ioss::Region &region, stk::mesh::fem::FEMMetaData &fem_meta) {
  process_nodesets(region, stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta));
}

void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  assert(nb->field_exists("mesh_model_coordinates"));
  Ioss::Field coordinates = nb->get_field("mesh_model_coordinates");
  int spatial_dim = coordinates.transformed_storage()->component_count();

  stk::mesh::Field<double,stk::mesh::Cartesian> & coord_field =
    meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

  stk::mesh::put_field( coord_field, node_rank(meta), meta.universal_part(), spatial_dim);
}

void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  std::vector<stk::mesh::Entity*> nodes;
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
  stk::io::get_entity_list(nb, node_rank(meta), bulk, nodes);

  /// \todo REFACTOR Application would probably store this field
  /// (and others) somewhere after the declaration instead of
  /// looking it up each time it is needed.
  stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
    meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

  stk::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");

  // Transfer any nodal "transient" fields from Ioss to stk
  // ... only if current state is set by begin_state call,
  // AND fields are in database
  int step = region.get_property("current_state").get_int();
  if (step>0) {
    Ioss::NameList names;
    nb->field_describe(Ioss::Field::TRANSIENT, &names);
    for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
      stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase>(*I);
      stk::io::field_data_from_ioss(field, nodes, nb, *I);
    }
  }
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, meta, 0);
}

void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();

  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);

      const CellTopologyData* cell_topo = stk::io::get_cell_topology(*part);
      if (cell_topo == NULL) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " returned NULL from get_cell_topology()";
        throw std::runtime_error( msg.str() );
      }

      std::vector<int> elem_ids ;
      std::vector<int> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = cell_topo->node_count ;

      std::vector<stk::mesh::Entity*> elements(element_count);
      for(size_t i=0; i<element_count; ++i) {
        /// \todo REFACTOR cast from int to unsigned is unsafe and ugly.
        /// change function to take int[] argument.
        int *conn = &connectivity[i*nodes_per_elem];
        std::vector<stk::mesh::EntityId> id_vec;

        for( int k = 0; k < nodes_per_elem; ++k )
          id_vec.push_back(conn[k]);
        elements[i] = &stk::mesh::fem::declare_element(bulk, *part, elem_ids[i], &id_vec[0]);
      }

      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase> (*I);
        if (field)
          stk::io::field_data_from_ioss(field, elements, entity, *I);
      }
    }
  }
}

// ========================================================================
// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta, 0);

  /** \todo REFACTOR should "distribution_factor" be a default field
   * that is automatically declared on all objects that it exists
   * on as is done in current framework?
   */
  stk::mesh::Field<double> & distribution_factors_field =
    meta.declare_field<stk::mesh::Field<double> >("distribution_factors");

  /** \todo REFACTOR How to associate distribution_factors field
   * with the nodeset part if a node is a member of multiple
   * nodesets
   */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());
      assert(part != NULL);
      assert(entity->field_exists("distribution_factors"));

      stk::mesh::put_field(distribution_factors_field, node_rank(meta), *part);

      /** \todo IMPLEMENT truly handle fields... For this case we
       * are just defining a field for each transient field that is
       * present in the mesh...
       */
      stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                *part, part_primary_entity_rank(*part));
    }
  }
}

// ========================================================================
// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  stk::io::default_part_processing(side_sets, meta, 0);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, meta, side_rank(meta));
    }
  }
}

// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<int> node_ids ;
      int node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity*> nodes(node_count);
      stk::mesh::EntityRank n_rank = node_rank(meta);
      for(int i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity(n_rank, node_ids[i] );
        if (nodes[i] != NULL)
          bulk.declare_entity(n_rank, node_ids[i], add_parts );
      }


      /** \todo REFACTOR Application would probably store this field
       * (and others) somewhere after the declaration instead of
       * looking it up each time it is needed.
       */
      stk::mesh::Field<double> *df_field =
        meta.get_field<stk::mesh::Field<double> >("distribution_factors");

      if (df_field != NULL) {
        stk::io::field_data_from_ioss(df_field, nodes, entity, "distribution_factors");
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
                    bool add_all)
{
  std::vector<stk::mesh::Entity*> entities;
  stk::io::get_entity_list(io_entity, part_type, bulk, entities);

  stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
  const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    if (stk::io::is_valid_part_field(f, part_type, part, meta.universal_part(),
                                     filter_role, add_all)) {
      stk::io::field_data_to_ioss(f, entities, io_entity, f->name(), filter_role);
    }
  }
}

int process_output_request(MeshData &mesh_data,
                           stk::mesh::BulkData &bulk,
                           double time, bool output_all_fields)
{
  Ioss::Region &region = *(mesh_data.m_region);
  region.begin_mode(Ioss::STATE_TRANSIENT);

  int out_step = region.add_state(time);

  process_output_request(region, bulk, out_step, output_all_fields);
  region.end_mode(Ioss::STATE_TRANSIENT);

  return out_step;
}

void process_output_request(Ioss::Region &region,
                            stk::mesh::BulkData &bulk,
                            int step, bool output_all_fields)
{
  region.begin_state(step);
  // Special processing for nodeblock (all nodes in model)...
  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);

  put_field_data(bulk, meta.universal_part(), node_rank(meta),
                 region.get_node_blocks()[0], Ioss::Field::Field::TRANSIENT,
                 output_all_fields);

  const stk::mesh::PartVector & all_parts = meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = region.get_entity(part->name());
      if (entity != NULL) {
        if (entity->type() == Ioss::ELEMENTBLOCK) {
          put_field_data(bulk, *part, part_primary_entity_rank(*part),
                         entity, Ioss::Field::Field::TRANSIENT,
                         output_all_fields);
        }
      }
    }
  }
  region.end_state(step);
}

void populate_bulk_data(stk::mesh::BulkData &bulk_data,
                        MeshData &mesh_data,
                        const std::string &mesh_type,
                        int step)
{
  if (mesh_type == "exodusii" || mesh_type == "generated" || mesh_type == "pamgen" ) {
    Ioss::Region *region = mesh_data.m_region;
    bulk_data.modification_begin();

    // Pick which time index to read into solution field.
    if (step>0) region->begin_state(step);

    stk::io::util::process_elementblocks(*region, bulk_data);
    stk::io::util::process_nodeblocks(*region, bulk_data); // solution field read here
    stk::io::util::process_nodesets(*region, bulk_data);
    stk::io::util::process_sidesets(*region, bulk_data);

    if (step>0) region->end_state(step);
    bulk_data.modification_end();
  }
  else if (mesh_type == "gears") {
    generate_gears(bulk_data, mesh_data.m_gears);
  }

  // NOTE: DO NOT USE THIS SKINNING
  if (mesh_data.m_generateSkinFaces) {
    stk::mesh::fem::FEMMetaData &fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);
    stk::mesh::Part* const skin_part = fem_meta.get_part("skin");
    stk::io::util::generate_sides(bulk_data, *skin_part, true);
  }
}
} // namespace util
} // namespace io
} // namespace stk

namespace {
void generate_gears(stk::ParallelMachine comm, const std::string &parameters, stk::mesh::fem::FEMMetaData &fem_meta,
		      std::vector<stk::io::util::Gear*> &gears)
  {
    const size_t spatial_dimension = 3;
    stk::io::initialize_spatial_dimension(stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta),
					  spatial_dimension,
					  stk::mesh::fem::entity_rank_names(spatial_dimension));
    const double TWO_PI = 2.0 * std::acos( static_cast<double>(-1.0) );

    int p_size = stk::parallel_machine_size( comm );
    int p_rank = stk::parallel_machine_rank( comm );

    // The parameters should be of the form:  "IxJxK|option:param,param,param|option:a,b,c"
    // Each "|" or "+" separated section of the parameters is a "group"
    // Each group is then split into options and params

    std::vector<std::string> groups;
    stk::util::tokenize(parameters, "|+", groups);

    // First 'group' is the interval specification -- IxJxK
    std::vector<std::string> ijktokens;
    stk::util::tokenize(groups[0], "x", ijktokens);
    assert(ijktokens.size() == 3);

    size_t max_end = std::numeric_limits<size_t>::max();
    double max_end_double = static_cast<double>(max_end);

    size_t i_end = std::strtol(ijktokens[0].c_str(), NULL, 10);
    size_t j_end = std::strtol(ijktokens[1].c_str(), NULL, 10);
    size_t k_end = std::strtol(ijktokens[2].c_str(), NULL, 10);

    if (i_end == 0 || j_end == 0 || k_end == 0) {
      std::cerr << "ERROR: Invalid parameters for gears.\n"
		<< "       All increments must be positive. Found: "
		<< "       i, j k = " << i_end << ", " << j_end << ", "<< k_end << "\n";
      std::exit(EXIT_FAILURE);
    }

    double ijk_double = static_cast<double>(i_end) * static_cast<double>(j_end) * static_cast<double>(k_end);
    if (i_end > max_end || j_end > max_end || k_end > max_end || ijk_double > max_end_double) {
      std::cerr << "ERROR: Invalid parameters for gears.\n"
		<< "       All increments must be less than " << max_end << ". Found: "
		<< "       i, j k = " << i_end << ", " << j_end << ", "<< k_end << "\n";
      std::exit(EXIT_FAILURE);
    }

    // Gear parameters.... Defaults set here...
    // Exactly touch = 1.0, force overlap by adding a little
    double rad_max = 1.0 + 0.001 ;
    double rad_min = 0.6 ;

    double z_min   = -0.4 ;
    double z_max   =  0.4 ;

    double elem_h = 0.10 ;

    for (size_t i=1; i < groups.size(); i++) {
      std::vector<std::string> option;
      stk::util::tokenize(groups[i], ":", option);
      // option[0] is the type of the option and option[1] is the argument to the option.

      if (option[0] == "radius") {
	std::vector<std::string> tokens;
        stk::util::tokenize(option[1], ",", tokens);
	assert(tokens.size() == 2);
	rad_min = std::strtod(tokens[0].c_str(), NULL);
	rad_max = std::strtod(tokens[1].c_str(), NULL);
	if (rad_max < rad_min) std::swap(rad_max, rad_min);
      }

      else if (option[0] == "height") {
	std::vector<std::string> tokens;
        stk::util::tokenize(option[1], ",", tokens);
	assert(tokens.size() == 2);
	z_min = std::strtod(tokens[0].c_str(), NULL);
	z_max = std::strtod(tokens[1].c_str(), NULL);
	if (z_max < z_min) std::swap(z_max, z_min);
      }

      else if (option[0] == "size") {
	elem_h = std::strtod(option[1].c_str(), NULL);
      }
    }

    stk::io::util::GearFields gear_fields( fem_meta );

    const size_t angle_num = static_cast<size_t>( TWO_PI / elem_h );
    const size_t rad_num   = static_cast<size_t>( 1 + ( rad_max - rad_min ) / elem_h );
    const size_t z_num     = static_cast<size_t>( 1 + ( z_max   - z_min )   / elem_h );
    const size_t elem_gear = angle_num * ( rad_num - 1 ) * ( z_num - 1 );
    const size_t num_gear  = k_end * j_end * i_end ;
    const size_t num_elem  = elem_gear * num_gear ;

    if ( p_rank == 0 ) {
      std::cout << "\n"
		<< "GEARS meshing:" << "\n"
		<< "  Number of Processors = " << p_size    << "\n"
		<< "  Number of Elements   = " << num_elem  << "\n"
		<< "  Number of Gears      = " << num_gear  << "\n"
		<< "  Number of Elem/gear  = " << elem_gear << "\n"
		<< "  Number of Angles     = " << angle_num << "\n"
		<< "  Number of Radial     = " << rad_num   << "\n"
		<< "  Number through Thick = " << z_num     << "\n"
		<< "  Gear Radii           = " << rad_max << "\t" << rad_min
		<< " (Touch = 1.0, force overlap by adding a little)\n"
		<< "  Gear Height (z-range)= " << z_min   << "\t" << z_max   << "\n"
		<< "  Element Size         = " << elem_h << "\n"
		<< "  Increments: i, j k   = " << i_end << ", " << j_end << ", "<< k_end << std::endl;
    }

    gears.resize( i_end * j_end * k_end );

    const double sqrt_3 = std::sqrt( 3.0 );
    for ( size_t k = 0 ; k < k_end ; ++k ) {
      for ( size_t j = 0 ; j < j_end ; ++j ) {
	double center[3] ;
	center[2] = k - z_min ;
	center[1] = sqrt_3 * j ;
	for ( size_t i = 0 ; i < i_end ; ++i ) {
	  int dir = i % 2 ? 1 : -1 ;

	  if ( j % 2 ) { // Odd
	    center[0] = static_cast<double>(i * 3 + i % 2) ;
	    dir = - dir ;
	  }
	  else { // Even
	    center[0] = static_cast<double>(i * 3 + ( 1 - i % 2 ));
	  }

	  std::ostringstream name ; name << "G_" << i << "_" << j << "_" << k ;

	  stk::io::util::Gear * g = new stk::io::util::Gear( fem_meta , name.str() , gear_fields ,
							     center ,
							     rad_min , rad_max , rad_num ,
							     z_min , z_max , z_num ,
							     angle_num , dir );

	  gears[ k * j_end * i_end + j * i_end + i ] = g ;

	}
      }
    }
  }

  void generate_gears(stk::mesh::BulkData &mesh,
		      std::vector<stk::io::util::Gear*> &gears)
  {
    const unsigned p_rank = mesh.parallel_rank();

    for ( std::vector<stk::io::util::Gear*>::iterator
	    i = gears.begin() ; i != gears.end() ; ++i ) {
      (*i)->mesh( mesh );
    }

    mesh.modification_end();

    // Copy coordinates to the aura nodes
    {
      std::vector< const stk::mesh::FieldBase *> fields ;
      const stk::mesh::FieldBase * ptr = NULL ;

      stk::io::util::Gear  *gear = *gears.begin();
      ptr = & gear->m_gear_coord;    fields.push_back( ptr );
      ptr = & gear->m_model_coord;   fields.push_back( ptr );

      stk::mesh::communicate_field_data(mesh.shared_aura(), fields);
    }
    //------------------------------
    {
      std::vector<size_t> counts ;
      stk::mesh::comm_mesh_counts( mesh , counts);

      if ( p_rank == 0 ) {
	std::cout << "N_GEARS Meshing completed and verified" << std::endl ;

  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(mesh);
	std::cout << "N_GEARS Global Counts { "
                  << " Node = " << counts[ stk::io::node_rank(meta)]
                  << " Edge = " << counts[ stk::io::edge_rank(meta)]
                  << " Face = " << counts[ stk::io::face_rank(meta)]
                  << " Elem = " << counts[ stk::io::element_rank(meta)] << std::endl;
      }
    }
  }
}
