/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <assert.h>

#include <stk_util/parallel/Parallel.hpp>
#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>

static const size_t spatial_dimension = 3;

namespace stk_examples {
namespace app {

  // ========================================================================
  void process_nodeblocks(Ioss::Region &region, stk_classic::mesh::fem::FEMMetaData &meta)
  {
    const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
    assert(node_blocks.size() == 1);

    Ioss::NodeBlock *nb = node_blocks[0];

    assert(nb->field_exists("mesh_model_coordinates"));
    Ioss::Field coordinates = nb->get_field("mesh_model_coordinates");
    int spatial_dim = coordinates.transformed_storage()->component_count();

    stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & coord_field =
      meta.declare_field<stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> >("coordinates");

    stk_classic::mesh::put_field( coord_field, stk_classic::mesh::fem::FEMMetaData::NODE_RANK, meta.universal_part(),
                          spatial_dim);

    /// \todo IMPLEMENT truly handle fields... For this case we
    /// are just defining a field for each transient field that is
    /// present in the mesh...
    stk_classic::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk_classic::mesh::fem::FEMMetaData::NODE_RANK);
  }

  // ========================================================================
  void process_elementblocks(Ioss::Region &region, stk_classic::mesh::fem::FEMMetaData &meta)
  {
    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
    stk_classic::io::default_part_processing(elem_blocks, meta);

    // Parts were created above, now handle element block specific
    // information (topology, attributes, ...);
    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	it != elem_blocks.end(); ++it) {
      Ioss::ElementBlock *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	stk_classic::mesh::Part* const part = meta.get_part(entity->name());
	assert(part != NULL);

        const stk_classic::mesh::EntityRank part_rank = part->primary_entity_rank();

	// Element Block attributes (if any)...
	/// \todo IMPLEMENT truly handle attribute fields... For this case we
	/// are just defining a field for each attribute field that
	/// is present in the mesh...
	stk_classic::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE,
                                  *part, part_rank );

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk_classic::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                  *part, part_rank);

	std::cout << entity->type_string() << ": " << entity->name()
		  << " , celltop = " << meta.get_cell_topology(*part).getCellTopologyData()->name
		  << std::endl ;
      }
    }
  }

  // ========================================================================
  void process_nodesets(Ioss::Region &region, stk_classic::mesh::fem::FEMMetaData &meta)
  {
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
    stk_classic::io::default_part_processing(node_sets, meta);

    /// \todo REFACTOR should "distribution_factor" be a default field
    /// that is automatically declared on all objects that it exists
    /// on as is done in current framework?
    stk_classic::mesh::Field<double> & distribution_factors_field =
      meta.declare_field<stk_classic::mesh::Field<double> >("distribution_factors");

    /// \todo REFACTOR How to associate distribution_factors field
    /// with the nodeset part if a node is a member of multiple
    /// nodesets

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
	it != node_sets.end(); ++it) {
      Ioss::NodeSet *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	stk_classic::mesh::Part* const part = meta.get_part(entity->name());
	assert(part != NULL);
	assert(entity->field_exists("distribution_factors"));

        const stk_classic::mesh::EntityRank part_rank = part->primary_entity_rank();

	stk_classic::mesh::put_field(distribution_factors_field, stk_classic::mesh::fem::FEMMetaData::NODE_RANK, *part);

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk_classic::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                  *part, part_rank);
      }
    }
  }

  // ========================================================================
  void process_surface_entity(Ioss::SideSet *sset, stk_classic::mesh::fem::FEMMetaData &meta)
  {
    assert(sset->type() == Ioss::SIDESET);
    const Ioss::SideBlockContainer& blocks = sset->get_side_blocks();
    stk_classic::io::default_part_processing(blocks, meta);

    stk_classic::mesh::Part* const sideset_part = meta.get_part(sset->name());
    assert(sideset_part != NULL);

    stk_classic::mesh::Field<double, stk_classic::mesh::ElementNode> *distribution_factors_field = NULL;
    bool surface_df_defined = false; // Has the surface df field been defined yet?


    int block_count = sset->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *sideblock = sset->get_block(i);
      if (stk_classic::io::include_entity(sideblock)) {
	std::cout << sideblock->type_string() << " " << sideblock->name() << "\n";
	stk_classic::mesh::Part * const sideblock_part = meta.get_part(sideblock->name());
	assert(sideblock_part != NULL);
	meta.declare_part_subset(*sideset_part, *sideblock_part);

        const stk_classic::mesh::EntityRank part_rank = sideblock_part->primary_entity_rank();

	/// \todo REFACTOR How to associate distribution_factors field
	/// with the sideset part if a face is a member of multiple
	/// sidesets

	if (sideblock->field_exists("distribution_factors")) {
	  if (!surface_df_defined) {
	    std::string field_name = sset->name() + "_distribution_factors";
	    distribution_factors_field =
	      &meta.declare_field<stk_classic::mesh::Field<double, stk_classic::mesh::ElementNode> >(field_name);
	    stk_classic::io::set_distribution_factor_field(*sideset_part, *distribution_factors_field);
	    surface_df_defined = true;
	  }
	  stk_classic::io::set_distribution_factor_field(*sideblock_part, *distribution_factors_field);
	  int side_node_count = sideblock->topology()->number_nodes();
	  stk_classic::mesh::put_field(*distribution_factors_field,
                               part_rank ,
                               *sideblock_part, side_node_count);
	}

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk_classic::io::define_io_fields(sideblock, Ioss::Field::TRANSIENT,
                                  *sideblock_part, part_rank );
      }
    }
  }

  // ========================================================================
  void process_sidesets(Ioss::Region &region, stk_classic::mesh::fem::FEMMetaData &meta)
  {
    const Ioss::SideSetContainer& side_sets = region.get_sidesets();
    stk_classic::io::default_part_processing(side_sets, meta);

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
	it != side_sets.end(); ++it) {
      Ioss::SideSet *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	process_surface_entity(entity, meta);
      }
    }
  }

  // ========================================================================
  // Bulk Data
  // ========================================================================
  void process_nodeblocks(Ioss::Region &region, stk_classic::mesh::BulkData &bulk)
  {
    // This must be called after the "process_element_blocks" call
    // since there may be nodes that exist in the database that are
    // not part of the analysis mesh due to subsetting of the element
    // blocks.

    const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
    assert(node_blocks.size() == 1);

    Ioss::NodeBlock *nb = node_blocks[0];

    std::vector<stk_classic::mesh::Entity*> nodes;
    stk_classic::io::get_entity_list(nb, stk_classic::mesh::fem::FEMMetaData::NODE_RANK, bulk, nodes);

    /// \todo REFACTOR Application would probably store this field
    /// (and others) somewhere after the declaration instead of
    /// looking it up each time it is needed.
    const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get(bulk);
    stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> *coord_field =
      meta.get_field<stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> >("coordinates");

    stk_classic::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");
  }

  // ========================================================================
  void process_elementblocks(Ioss::Region &region, stk_classic::mesh::BulkData &bulk)
  {
    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();

    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	it != elem_blocks.end(); ++it) {
      Ioss::ElementBlock *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	const std::string &name = entity->name();
  const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get(bulk);
	stk_classic::mesh::Part* const part = meta.get_part(name);
	assert(part != NULL);

        const CellTopologyData* cell_topo = stk_classic::mesh::fem::FEMMetaData::get(bulk).get_cell_topology(*part).getCellTopologyData();
	assert(cell_topo != NULL);

	std::vector<int> elem_ids ;
	std::vector<int> connectivity ;
	std::vector<stk_classic::mesh::EntityId> connectivity2 ;

	entity->get_field_data("ids", elem_ids);
	entity->get_field_data("connectivity", connectivity);
        connectivity2.reserve(connectivity.size());
        std::copy(connectivity.begin(), connectivity.end(), std::back_inserter(connectivity2));
          
	int element_count = elem_ids.size();
	int nodes_per_elem = cell_topo->node_count ;

	std::vector<stk_classic::mesh::Entity*> elements(element_count);
	for(int i=0; i<element_count; ++i) {
	  /// \todo REFACTOR cast from int to unsigned is unsafe and ugly.
	  /// change function to take int[] argument.
	  stk_classic::mesh::EntityId *conn = &connectivity2[i*nodes_per_elem];
	  elements[i] = &stk_classic::mesh::fem::declare_element(bulk, *part, elem_ids[i], conn);
	}

	// For this example, we are just taking all attribute fields
	// found on the io database and populating fields on the
	// corresponding mesh part.  In practice, would probably be
	// selective about which attributes to use...
	Ioss::NameList names;
	entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
	for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
	  if (*I == "attribute" && names.size() > 1)
	    continue;
	  stk_classic::mesh::FieldBase *field = meta.get_field<stk_classic::mesh::FieldBase>(*I);
	  stk_classic::io::field_data_from_ioss(field, elements, entity, *I);

	}
      }
    }
  }

  // ========================================================================
  void process_nodesets(Ioss::Region &region, stk_classic::mesh::BulkData &bulk)
  {
    // Should only process nodes that have already been defined via the element
    // blocks connectivity lists.
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
	it != node_sets.end(); ++it) {
      Ioss::NodeSet *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	const std::string & name = entity->name();
  const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get(bulk);
	stk_classic::mesh::Part* const part = meta.get_part(name);
	assert(part != NULL);
	stk_classic::mesh::PartVector add_parts( 1 , part );

	std::vector<int> node_ids ;
	int node_count = entity->get_field_data("ids", node_ids);

	std::vector<stk_classic::mesh::Entity*> nodes(node_count);
	for(int i=0; i<node_count; ++i) {
	  nodes[i] = bulk.get_entity( stk_classic::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i] );
	  if (nodes[i] != NULL)
	    bulk.declare_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i], add_parts );
	}


	/// \todo REFACTOR Application would probably store this field
	/// (and others) somewhere after the declaration instead of
	/// looking it up each time it is needed.
	stk_classic::mesh::Field<double> *df_field =
	  meta.get_field<stk_classic::mesh::Field<double> >("distribution_factors");

	if (df_field != NULL) {
	  stk_classic::io::field_data_from_ioss(df_field, nodes, entity, "distribution_factors");
	}
      }
    }
  }

  // ========================================================================
  void process_surface_entity(const Ioss::SideSet* io, stk_classic::mesh::BulkData & bulk)
  {
    assert(io->type() == Ioss::SIDESET);
    const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get(bulk);
    stk_classic::mesh::fem::FEMMetaData &fem = stk_classic::mesh::fem::FEMMetaData::get(meta);

    const stk_classic::mesh::EntityRank element_rank = fem.element_rank();

    int block_count = io->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *block = io->get_block(i);
      if (stk_classic::io::include_entity(block)) {
	std::vector<int> side_ids ;
	std::vector<int> elem_side ;

	stk_classic::mesh::Part * const fb_part = meta.get_part(block->name());

	block->get_field_data("ids", side_ids);
	block->get_field_data("element_side", elem_side);

	assert(side_ids.size() * 2 == elem_side.size());
	stk_classic::mesh::PartVector add_parts( 1 , fb_part );

	int side_count = side_ids.size();
	std::vector<stk_classic::mesh::Entity*> sides(side_count);
	for(int is=0; is<side_count; ++is) {

	  stk_classic::mesh::Entity* const elem = bulk.get_entity(element_rank, elem_side[is*2]);
	  // If NULL, then the element was probably assigned to an
	  // Ioss uses 1-based side ordinal, stk_classic::mesh uses 0-based.
	  // Hence the '-1' in the following line.
	  int side_ordinal = elem_side[is*2+1] - 1 ;

	  // element block that appears in the database, but was
	  // subsetted out of the analysis mesh. Only process if
	  // non-null.
	  if (elem != NULL) {
	    stk_classic::mesh::Entity& side =
	      stk_classic::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);
	    bulk.change_entity_parts( side, add_parts );
	    sides[is] = &side;
	  } else {
	    sides[is] = NULL;
	  }
	}

	const stk_classic::mesh::Field<double, stk_classic::mesh::ElementNode> *df_field =
	  stk_classic::io::get_distribution_factor_field(*fb_part);
	if (df_field != NULL) {
	  stk_classic::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
	}
      }
    }
  }

  // ========================================================================
  void process_sidesets(Ioss::Region &region, stk_classic::mesh::BulkData &bulk)
  {
    const Ioss::SideSetContainer& side_sets = region.get_sidesets();

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
	it != side_sets.end(); ++it) {
      Ioss::SideSet *entity = *it;

      if (stk_classic::io::include_entity(entity)) {
	process_surface_entity(entity, bulk);
      }
    }
  }

  // ========================================================================
  // ========================================================================
  void get_field_data(stk_classic::mesh::BulkData &bulk, stk_classic::mesh::Part &part,
		      stk_classic::mesh::EntityRank part_type,
		      Ioss::GroupingEntity *io_entity,
		      Ioss::Field::RoleType filter_role)
  {
    std::vector<stk_classic::mesh::Entity*> entities;
    stk_classic::io::get_entity_list(io_entity, part_type, bulk, entities);

    stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(part);
    stk_classic::mesh::Part &universal = meta.universal_part();
    const std::vector<stk_classic::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk_classic::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk_classic::mesh::FieldBase *f = *I; ++I;
      if (stk_classic::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
	stk_classic::io::field_data_from_ioss(f, entities, io_entity, f->name());
      }
    }
  }

  void process_input_request(Ioss::Region &region,
			     stk_classic::mesh::BulkData &bulk,
			     int step)
  {
    region.begin_state(step);

    // Special processing for nodeblock (all nodes in model)...
    const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(bulk);

    // ??? Get field data from nodeblock...
    app::get_field_data(bulk, meta.universal_part(), stk_classic::mesh::fem::FEMMetaData::NODE_RANK,
			region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

    const stk_classic::mesh::PartVector & all_parts = meta.get_parts();
    for ( stk_classic::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk_classic::mesh::Part * const part = *ip;

      const stk_classic::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk_classic::io::is_part_io_part(*part)) {
	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = region.get_entity(part->name());
	if (entity != NULL) {
	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    assert(sset != NULL);
	    int block_count = sset->block_count();
	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *fb = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
		app::get_field_data(bulk, *part, part_rank,
				    fb, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    app::get_field_data(bulk, *part, part_rank,
				entity, Ioss::Field::TRANSIENT);
	  }
	} else {
	  /// \todo IMPLEMENT handle error... Possibly an assert since
	  /// I think the corresponding entity should always exist...
	}
      }
    }

    region.end_state(step);
  }

  void put_field_data(stk_classic::mesh::BulkData &bulk, stk_classic::mesh::Part &part,
		      stk_classic::mesh::EntityRank part_type,
		      Ioss::GroupingEntity *io_entity,
		      Ioss::Field::RoleType filter_role)
  {
    std::vector<stk_classic::mesh::Entity*> entities;
    stk_classic::io::get_entity_list(io_entity, part_type, bulk, entities);

    stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(part);
    stk_classic::mesh::Part &universal = meta.universal_part();
    const std::vector<stk_classic::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk_classic::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk_classic::mesh::FieldBase *f = *I; ++I;
      if (stk_classic::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
	stk_classic::io::field_data_to_ioss(f, entities, io_entity, f->name(), filter_role);
      }
    }
  }

  void process_output_request(Ioss::Region &region,
			     stk_classic::mesh::BulkData &bulk,
			     int step)
  {
    region.begin_state(step);
    // Special processing for nodeblock (all nodes in model)...
    const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(bulk);

    app::put_field_data(bulk, meta.universal_part(), stk_classic::mesh::fem::FEMMetaData::NODE_RANK,
			region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

    const stk_classic::mesh::PartVector & all_parts = meta.get_parts();
    for ( stk_classic::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk_classic::mesh::Part * const part = *ip;

      const stk_classic::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk_classic::io::is_part_io_part(*part)) {

	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = region.get_entity(part->name());
	if (entity != NULL) {

	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    int block_count = sset->block_count();

	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *fb = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
	      app::put_field_data(bulk, *part, part_rank,
				  fb, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    app::put_field_data(bulk, *part, part_rank,
				entity, Ioss::Field::TRANSIENT);
	  }
	} else {
	  /// \todo IMPLEMENT handle error... Possibly an assert since
	  /// I think the corresponding entity should always exist...
	}
      }
    }
    region.end_state(step);
  }

} // namespace application

  // ========================================================================
  //------------------------------------------------------------------------
  void example_io_2( stk_classic::ParallelMachine comm,
		     const std::string& in_filename,
		     const std::string& out_filename)
  {
    // Initialize IO system.  Registers all element types and storage
    // types and the exodusII default database type.
    Ioss::Init::Initializer init_db;

    std::cout << "========================================================================\n"
	      << " Use Case 2: Subsetting with df and attribute field input/output        \n"
	      << "             (Also contains functionality of use case 1.5)              \n"
	      << "========================================================================\n";

    std::string dbtype("exodusII");
    MPI_Comm mpicomm = comm;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(dbtype, in_filename, Ioss::READ_MODEL,
						    mpicomm);
    if (dbi == NULL || !dbi->ok()) {
      std::cerr  << "ERROR: Could not open database '" << in_filename
		 << "' of type '" << dbtype << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'in_region' owns 'dbi' pointer at this time...
    Ioss::Region in_region(dbi, "input_model");

    // SUBSETTING PARSING/PREPROCESSING...
    // Just an example of how application could control whether an
    // entity is subsetted or not...


    // Example command line in current code corresponding to behavior below:
    std::cout << "\nWhen processing file multi-block.g for use case 2, the blocks below will be omitted:\n";
    std::cout << "\tOMIT BLOCK Cblock Eblock I1 I2\n\n";
    Ioss::ElementBlock *eb = in_region.get_element_block("cblock");
    if (eb != NULL)
      eb->property_add(Ioss::Property(std::string("omitted"), 1));

    eb = in_region.get_element_block("eblock");
    if (eb != NULL)
      eb->property_add(Ioss::Property(std::string("omitted"), 1));

    eb = in_region.get_element_block("i1");
    if (eb != NULL)
      eb->property_add(Ioss::Property(std::string("omitted"), 1));

    eb = in_region.get_element_block("i2");
    if (eb != NULL)
      eb->property_add(Ioss::Property(std::string("omitted"), 1));

#if 0
    // Example for subsetting -- omit "odd" blocks
    if (entity->type() == Ioss::ELEMENTBLOCK) {
      int id = entity->get_property("id").get_int();
      if (id % 2) {
	entity->property_add(Ioss::Property(std::string("omitted"), 1));
	std::cout << "Skipping " << entity->type_string() << ": "  << entity->name() << "\n";
      }
    }
#endif

    //----------------------------------
    // Process Entity Types. Subsetting is possible.
    stk_classic::mesh::fem::FEMMetaData meta_data( spatial_dimension );
    app::process_elementblocks(in_region, meta_data);
    app::process_nodeblocks(in_region,    meta_data);
    app::process_sidesets(in_region,      meta_data);
    app::process_nodesets(in_region,      meta_data);

    //----------------------------------
    // Done populating meta data, commit and create bulk data
    meta_data.commit();

    //----------------------------------
    // Process Bulkdata for all Entity Types. Subsetting is possible.
    stk_classic::mesh::BulkData bulk_data(meta_data.get_meta_data(meta_data), comm);
    bulk_data.modification_begin();
    app::process_elementblocks(in_region, bulk_data);
    app::process_nodeblocks(in_region,    bulk_data);
    app::process_sidesets(in_region,      bulk_data);
    app::process_nodesets(in_region,      bulk_data);
    bulk_data.modification_end();

    //----------------------------------
    // OUTPUT...Create the output "mesh" portion

    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, out_filename,
						    Ioss::WRITE_RESULTS,
						    mpicomm);
    if (dbo == NULL || !dbo->ok()) {
      std::cerr << "ERROR: Could not open results database '" << out_filename
		<< "' of type '" << dbtype << "'\n";
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'out_region' owns 'dbo' pointer at this time...
    Ioss::Region out_region(dbo, "results_output");

    /// \todo REFACTOR Need better interface for passing application
    /// fields to stk_classic::io.  This maybe works with only a single field,
    /// but adding attributes, distribution factors and perhaps other
    /// fields will make this unwieldy.
    stk_classic::io::define_output_db(out_region, bulk_data, &in_region);
    stk_classic::io::write_output_db(out_region,  bulk_data);

    // ------------------------------------------------------------------------
    /// \todo REFACTOR A real app would register a subset of the
    /// fields on the mesh database as fields that the app would want
    // read at one or all or specified steps.  In this example, all
    // fields existing on the input mesh database are defined on the
    // parts in the stk_classic::mesh.
    //
    // The real app would also only register a subset of the stk_classic::mesh
    // fields as output fields and would probably have a mapping from
    // the internally used name to some name picked by the user. In
    // this example, all TRANSIENT fields defined on the stk_classic::mesh are
    // output to the results database and the internal stk_classic::mesh field
    // name is used as the name on the database....

    out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // Special processing for nodeblock (all nodes in model)...
    stk_classic::io::ioss_add_fields(meta_data.universal_part(), stk_classic::mesh::fem::FEMMetaData::NODE_RANK,
			     out_region.get_node_blocks()[0],
			     Ioss::Field::TRANSIENT);

    const stk_classic::mesh::PartVector & all_parts = meta_data.get_parts();
    for ( stk_classic::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk_classic::mesh::Part * const part = *ip;

      const stk_classic::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk_classic::io::is_part_io_part(*part)) {
	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = out_region.get_entity(part->name());
	if (entity != NULL) {
	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    int block_count = sset->block_count();
	    for (int i=0; i < block_count; i++) {
	      Ioss::EntityBlock *fb = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
		stk_classic::io::ioss_add_fields(*part, part_rank,
					 fb, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    stk_classic::io::ioss_add_fields(*part, part_rank,
				     entity, Ioss::Field::TRANSIENT);
	  }
	} else {
	  /// \todo IMPLEMENT handle error... Possibly an assert since
	  /// I think the corresponding entity should always exist...
	}
      }
    }
    out_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    // ------------------------------------------------------------------------

    // Read and Write transient fields...
    out_region.begin_mode(Ioss::STATE_TRANSIENT);
    int timestep_count = in_region.get_property("state_count").get_int();
    for (int step = 1; step <= timestep_count; step++) {
      double time = in_region.get_state_time(step);

      // Read data from the io input mesh database into stk_classic::mesh fields...
      app::process_input_request(in_region, bulk_data, step);

      // app::execute()

      // Write data from the stk_classic::mesh fields out to the output database.a
      int out_step = out_region.add_state(time);
      app::process_output_request(out_region, bulk_data, out_step);
    }
    out_region.end_mode(Ioss::STATE_TRANSIENT);

  }
} // namespace stk_examples

