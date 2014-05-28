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
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>

/** \addtogroup stk_io_module
 * \{
 */

/**
 * Example code showing a basic, but complete, mesh to results output
 * coding including subsetting and periodic field input and output.
 * Includes handling of nodeblocks, element blocks, nodesets, and
 * sidesets.  Attribute fields and distribution factor fields are also
 * supported.
 *
 * This example can serve as the basis for adding binary IO support to
 * an application.  The code here uses the Ioss to/from stk::mesh
 * bridge functions in the stk::io namespace defined in IossBridge.hpp
 * include file.
 */
namespace stk_example_io {

  /// Declare "coordinates" field and put it on the universal part. This
  /// example also defines all Ioss::Field::TRANSIENT fields that exist on the
  /// Ioss::Nodeblock as fields on the universal part.
  void process_nodeblocks    (Ioss::Region &region, stk::mesh::MetaData &meta);

  /// Declare a part for each element block on the Ioss::Region
  /// 'region' unless the element block has the "omitted" property set
  /// to the value 1. The example then iterates each element block and
  /// defines any Ioss::Field::ATTRIBUTE and Ioss::Field::TRANSIENT fields that exist on the
  /// Ioss::ElementBlock as fields on the corresponding part.
  void process_elementblocks (Ioss::Region &region, stk::mesh::MetaData &meta);

  /// Declare a part for each Ioss::NodeSet on the Ioss::Region
  /// 'region' unless the nodeset has the "omitted" property set
  /// to the value 1. The example then iterates each nodeset and
  /// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
  /// exist on the Ioss::NodeSet as fields on the corresponding
  /// part.
  void process_nodesets      (Ioss::Region &region, stk::mesh::MetaData &meta);

  /// Declare a part for each Ioss::SideSet on the Ioss::Region
  /// 'region' unless the sideset has the "omitted" property set
  /// to the value 1. The example then iterates each sideset and
  /// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
  /// exist on the Ioss::SideSet as fields on the corresponding
  /// part.
  ///
  /// Each sideblock in the active sidesets is then processed by
  /// defining a part for each Ioss::SideBlock on the Ioss::SideSet
  /// unless the sideblock has the "omitted" property set to the value
  /// 1. The example then iterates each sideblock and defines any
  /// "distribution factor" and Ioss::Field::TRANSIENT fields that exist on the
  /// Ioss::SideBlock as fields on the corresponding part.
  void process_sidesets      (Ioss::Region &region, stk::mesh::MetaData &meta);

  /// NOTE: This must be called after the process_elementblocks() call
  /// since there may be nodes that exist in the database that are
  /// not part of the analysis mesh due to subsetting of the element
  /// blocks.
  ///
  /// Populates  the "coordinates" field for all active nodes in the model.
  void process_nodeblocks    (Ioss::Region &region, stk::mesh::BulkData &bulk);

  /// NOTE: This should be the first function called of any of the
  /// "process_X" type functions that take an stk::mesh::BulkData
  /// argument, especially if the input Ioss::Region mesh is going to
  /// be subsetted (have element blocks omitted).
  ///
  /// This function iterates all non-omitted element blocks and
  /// declares each element (and the corresponding nodes) in the
  /// element block. If there are any Ioss::Field::ATTRIBUTE fields on the element
  /// block (for example, shell thickness or particle radius), then
  /// that field data is alse read and the corresponding
  /// stk::mesh::Field populated.
  void process_elementblocks (Ioss::Region &region, stk::mesh::BulkData &bulk);

  /// Iterates each non-omitted Ioss::NodeSet and then iterates each
  /// node in the Ioss::NodeSet.  If the node exists (that is, it is
  /// connected to a non-omitted Ioss::ElementBlock), then that node
  /// is associated with the part corresponding to this
  /// Ioss::NodeSet. If the "distribution_factor" field exists, then
  /// that data is also associated with the field.
  void process_nodesets      (Ioss::Region &region, stk::mesh::BulkData &bulk);

  /// Process each non-omitted Ioss::SideSet and the contained
  /// non-omitted Ioss::SideBlock and associate each element-side pair with
  /// the corresponding part if the underlying element is active.  If
  /// the "distribution_factor" field exists, then that data is also
  /// associated with the corresponding field.
  void process_sidesets      (Ioss::Region &region, stk::mesh::BulkData &bulk);

  /// A minimal example function showing how field data on the
  /// Ioss::Region entities can be periodically transferred to the
  /// corresponding field(s) on the stk::mesh entities. This would be
  /// used to bring in initial condition data or interpolation data or
  /// any other scenario in which data on the mesh file needs to be
  /// transferred to the stk::mesh fields.
  void process_input_request (Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

  /// A minimal example function showing how stk::mesh field data can
  /// periodically be output to a results, history, heartbeat, or
  /// restart database.  The scheduling would be done either in this
  /// function or at a higher level and is not shown here. The
  /// function iterates all parts and if there is a corresponding Ioss
  /// part on the Ioss::Region, all fields defined to be output are
  /// iterated and their data output to the corresponding
  /// Ioss::Field. The function calls the
  /// stk::io::is_valid_part_field() function to determine whether the
  /// field should be output and then calls the
  /// stk::io::field_data_to_ioss() function to do the actual output
  /// of the field.
  void process_output_request(Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

  /// This function shows the basic calls needed to perform definition
  /// and input of the mesh model and definition and periodic output
  /// of a results database. The function is given the mesh filename
  /// and the output filename and goes through all steps of
  /// associating the filename with an Ioss::DatabaseIO object of the
  /// correct type ("exodusII" in this example); creating an
  /// Ioss::Region and then defining an stk::mesh corresponding to
  /// this mesh.  The function also provides an example of how
  /// specific element blocks existing in the mesh database could be
  /// omitted from the analysis model.
  ///
  /// The example then shows how to define a results database
  /// corresponding to the analysis model and periodically output the
  /// results in an execute loop.
  ///
  /// A true application would have to provide additional
  /// functionality and robustness, but the example shows how the
  /// basic functionality can be provided by an application.
  ///
  /// Note that the paradigm illustrated here is different than the
  /// mesh input and output paradigm provided in the current
  /// framework.  In this case, the application is responsible for the
  /// majority of the IO behavior and the toolkit only provides some
  /// helper functions to bridge between the Ioss and the stk::mesh.
  /// It is hoped that this paradigm will result in more functionality
  /// for the application with less complication and overhead.
  void io_example( stk::ParallelMachine comm,
		   const std::string& in_filename,
		   const std::string& out_filename,
		   const std::string& decomp_method)
  {
    // Initialize IO system.  Registers all element types and storage
    // types and the exodusII default database type.
    Ioss::Init::Initializer init_db;

    std::cout << "========================================================================\n"
	      << " Copy input mesh to output mesh.                                        \n"
	      << "========================================================================\n";

    std::string dbtype("exodusII");
    Ioss::PropertyManager properties;
    if (!decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", Ioss::Utils::uppercase(decomp_method)));
    }
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(dbtype, in_filename, Ioss::READ_MODEL,
						    comm, properties);
    if (dbi == NULL || !dbi->ok()) {
      std::cerr  << "ERROR: Could not open database '" << in_filename
		 << "' of type '" << dbtype << "'\n";
      std::exit(EXIT_FAILURE);
    }

    std::cout << "Reading input file:   " << in_filename << "\n";
    // NOTE: 'in_region' owns 'dbi' pointer at this time...
    Ioss::Region in_region(dbi, "input_model");

    // SUBSETTING PARSING/PREPROCESSING...
    // Just an example of how application could control whether an
    // entity is subsetted or not...


#if 0
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
#endif

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

    static size_t spatial_dimension = in_region.get_property("spatial_dimension").get_int();

    stk::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension );
    stk::mesh::MetaData &meta_data = fem_meta_data.get_meta_data(fem_meta_data);
    process_elementblocks(in_region, meta_data);
    process_nodeblocks(in_region,    meta_data);
    process_sidesets(in_region,      meta_data);
    process_nodesets(in_region,      meta_data);

    //----------------------------------
    // Done populating meta data, commit and create bulk data
    meta_data.commit();

    //----------------------------------
    // Process Bulkdata for all Entity Types. Subsetting is possible.
    stk::mesh::BulkData bulk_data(meta_data, comm);

    bulk_data.modification_begin();
    process_elementblocks(in_region, bulk_data);
    process_nodeblocks(in_region,    bulk_data);
    process_sidesets(in_region,      bulk_data);
    process_nodesets(in_region,      bulk_data);
    bulk_data.modification_end();

    //----------------------------------
    // OUTPUT...Create the output "mesh" portion

    std::cout << "Creating output file: " << out_filename << "\n";
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, out_filename,
						    Ioss::WRITE_RESULTS,
						    comm);
    if (dbo == NULL || !dbo->ok()) {
      std::cerr << "ERROR: Could not open results database '" << out_filename
		<< "' of type '" << dbtype << "'\n";
      std::exit(EXIT_FAILURE);
    }

#if 0
    {
      // Code to test the remove_io_part_attribute functionality.
      // Hook this up to a command line option at some point to test nightly...
      const stk::mesh::PartVector & all_parts = meta_data.get_parts();
      for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {
	stk::mesh::Part * const part = *ip;
	const stk::mesh::EntityRank part_rank = part->primary_entity_rank();
	
	if (stk::io::is_part_io_part(*part) && part_rank == 2) {
	  std::cout << "Removing part attribute from " << part->name() << "\n";
	  stk::io::remove_io_part_attribute(*part);
	}
      }
    }
#endif

    // NOTE: 'out_region' owns 'dbo' pointer at this time...
    Ioss::Region out_region(dbo, "results_output");

    stk::io::define_output_db(out_region, bulk_data, &in_region);
    stk::io::write_output_db(out_region,  bulk_data);

    // ------------------------------------------------------------------------
    /** \todo REFACTOR A real app would register a subset of the
     * fields on the mesh database as fields that the app would want
     * read at one or all or specified steps.  In this example, all
     * fields existing on the input mesh database are defined on the
     * parts in the stk::mesh.
     *
     * The real app would also only register a subset of the stk::mesh
     * fields as output fields and would probably have a mapping from
     * the internally used name to some name picked by the user. In
     * this example, all Ioss::Field::TRANSIENT fields defined on the stk::mesh are
     * output to the results database and the internal stk::mesh field
     * name is used as the name on the database....
     */

    out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // Special processing for nodeblock (all nodes in model)...
    stk::io::ioss_add_fields(meta_data.universal_part(), stk::mesh::fem::FEMMetaData::NODE_RANK,
			     out_region.get_node_blocks()[0],
			     Ioss::Field::TRANSIENT);

    const stk::mesh::PartVector & all_parts = meta_data.get_parts();
    for ( stk::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk::mesh::Part * const part = *ip;

      const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk::io::is_part_io_part(*part)) {
	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = out_region.get_entity(part->name());
	if (entity != NULL) {
	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    assert(sset != NULL);
	    int block_count = sset->block_count();
	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *fb = sset->get_block(i);
	      stk::io::ioss_add_fields(*part, part_rank,
				       fb, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    stk::io::ioss_add_fields(*part, part_rank,
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

      // Read data from the io input mesh database into stk::mesh fields...
      process_input_request(in_region, bulk_data, step);

      // execute()

      // Write data from the stk::mesh fields out to the output database.a
      int out_step = out_region.add_state(time);
      process_output_request(out_region, bulk_data, out_step);
    }
    out_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  // ========================================================================
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

    stk::mesh::put_field( coord_field, stk::mesh::fem::FEMMetaData::NODE_RANK, meta.universal_part(),
                          spatial_dim);

    /** \todo IMPLEMENT truly handle fields... For this case we are
     * just defining a field for each transient field that is present
     * in the mesh...
     */
    stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk::mesh::fem::FEMMetaData::NODE_RANK);
  }

  // ========================================================================
  void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    const stk::mesh::fem::FEMMetaData &fem_meta_data = stk::mesh::fem::FEMMetaData::get(meta);
    const stk::mesh::EntityRank element_rank = fem_meta_data.element_rank();

    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
    stk::io::default_part_processing(elem_blocks, meta, element_rank);

    // Parts were created above, now handle element block specific
    // information (topology, attributes, ...);
    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	it != elem_blocks.end(); ++it) {
      Ioss::ElementBlock *entity = *it;

      if (stk::io::include_entity(entity)) {
	stk::mesh::Part* const part = meta.get_part(entity->name());
	assert(part != NULL);

      const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

	// Element Block attributes (if any)...
	/** \todo IMPLEMENT truly handle attribute fields... For this
	 * case we are just defining a field for each attribute field
	 * that is present in the mesh...
	 */
	stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE,
				  *part,
                                  part_rank);

	/** \todo IMPLEMENT truly handle fields... For this case we
	 * are just defining a field for each transient field that is
	 * present in the mesh...
	 */
	stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
				  *part,
                                  part_rank);

	const CellTopologyData* cell_topo = fem_meta_data.get_cell_topology(*part).getCellTopologyData();
	std::string cell_topo_name = "UNKNOWN";
	if (cell_topo != NULL)
	  cell_topo_name = cell_topo->name;
      }
    }
  }

  // ========================================================================
  void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
    stk::io::default_part_processing(node_sets, meta, stk::mesh::fem::FEMMetaData::NODE_RANK);

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

	stk::mesh::put_field(distribution_factors_field, stk::mesh::fem::FEMMetaData::NODE_RANK, *part);

	/** \todo IMPLEMENT truly handle fields... For this case we
	 * are just defining a field for each transient field that is
	 * present in the mesh...
	 */
	stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
				  *part,
                                  part->primary_entity_rank() ) ;
      }
    }
  }

  // ========================================================================
  void process_surface_entity(Ioss::SideSet *sset, stk::mesh::MetaData &meta,
			      stk::mesh::EntityRank sset_rank)
  {
    assert(sset->type() == Ioss::SIDESET);
    Ioss::SideSet *fs = dynamic_cast<Ioss::SideSet *>(sset);
    assert(fs != NULL);
    const Ioss::SideBlockContainer& blocks = fs->get_side_blocks();
    stk::io::default_part_processing(blocks, meta, sset_rank);

    stk::mesh::Part* const fs_part = meta.get_part(sset->name());
    assert(fs_part != NULL);

    stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
    bool surface_df_defined = false; // Has the surface df field been defined yet?


    int block_count = sset->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *side_block = sset->get_block(i);
      if (stk::io::include_entity(side_block)) {
	stk::mesh::Part * const side_block_part = meta.get_part(side_block->name());
	assert(side_block_part != NULL);
	meta.declare_part_subset(*fs_part, *side_block_part);

        const stk::mesh::EntityRank part_rank = side_block_part->primary_entity_rank();

	if (side_block->field_exists("distribution_factors")) {
	  if (!surface_df_defined) {
	    std::string field_name = sset->name() + "_distribution_factors";
	    distribution_factors_field =
	      &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(field_name);
	    stk::io::set_distribution_factor_field(*fs_part, *distribution_factors_field);
	    surface_df_defined = true;
	  }
	  stk::io::set_distribution_factor_field(*side_block_part, *distribution_factors_field);
	  int side_node_count = side_block->topology()->number_nodes();
	  stk::mesh::put_field(*distribution_factors_field,
                               part_rank,
                               *side_block_part, side_node_count);
	}

	/** \todo IMPLEMENT truly handle fields... For this case we
	 * are just defining a field for each transient field that is
	 * present in the mesh...
	 */
	stk::io::define_io_fields(side_block, Ioss::Field::TRANSIENT,
				  *side_block_part,
                                  part_rank);
      }
    }
  }

  // ========================================================================
  void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    stk::mesh::fem::FEMMetaData &fem = stk::mesh::fem::FEMMetaData::get(meta);
    const stk::mesh::EntityRank side_rank = fem.side_rank();

    const Ioss::SideSetContainer& side_sets = region.get_sidesets();
    stk::io::default_part_processing(side_sets, meta, side_rank);

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
	it != side_sets.end(); ++it) {
      Ioss::SideSet *entity = *it;

      if (stk::io::include_entity(entity)) {
	process_surface_entity(entity, meta, side_rank);
      }
    }
  }

  // ========================================================================
  // Bulk Data
  // ========================================================================
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
    stk::io::get_entity_list(nb, stk::mesh::fem::FEMMetaData::NODE_RANK, bulk, nodes);

    /** \todo REFACTOR Application would probably store this field
     * (and others) somewhere after the declaration instead of
     * looking it up each time it is needed.
     */
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
    stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
      meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

    stk::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");
  }

  // ========================================================================
  void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
  {
    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();

    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	it != elem_blocks.end(); ++it) {
      Ioss::ElementBlock *entity = *it;

      if (stk::io::include_entity(entity)) {
	const std::string &name = entity->name();
        const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
        const stk::mesh::fem::FEMMetaData &fem_meta_data = stk::mesh::fem::FEMMetaData::get(meta);
	stk::mesh::Part* const part = meta.get_part(name);
	assert(part != NULL);

	const CellTopologyData* cell_topo = fem_meta_data.get_cell_topology(*part).getCellTopologyData();
	if (cell_topo == NULL) {
          std::ostringstream msg ;
	  msg << " INTERNAL_ERROR: Part " << part->name() << " returned NULL from get_cell_topology()";
	  throw std::runtime_error( msg.str() );
	}

	std::vector<int> elem_ids ;
	std::vector<int> connectivity ;
	std::vector<stk::mesh::EntityId> connectivity2 ;

	entity->get_field_data("ids", elem_ids);
	entity->get_field_data("connectivity", connectivity);
        connectivity2.reserve(connectivity.size());
        std::copy(connectivity.begin(), connectivity.end(), std::back_inserter(connectivity2));

	size_t element_count = elem_ids.size();
	int nodes_per_elem = cell_topo->node_count ;

	std::vector<stk::mesh::Entity*> elements(element_count);
	for(size_t i=0; i<element_count; ++i) {
	  stk::mesh::EntityId *conn = &connectivity2[i*nodes_per_elem];
	  elements[i] = &stk::mesh::fem::declare_element(bulk, *part, elem_ids[i], conn);
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
	  stk::mesh::FieldBase *field = meta.get_field<stk::mesh::FieldBase>(*I);
	  stk::io::field_data_from_ioss(field, elements, entity, *I);

	}
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
	for(int i=0; i<node_count; ++i) {
	  nodes[i] = bulk.get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i] );
	  if (nodes[i] != NULL)
	    bulk.declare_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node_ids[i], add_parts );
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
  void process_surface_entity(const Ioss::SideSet* sset ,
			      stk::mesh::BulkData & bulk)
  {
    assert(sset->type() == Ioss::SIDESET);

    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
    const stk::mesh::fem::FEMMetaData &fem_meta_data = stk::mesh::fem::FEMMetaData::get(meta);
    const stk::mesh::EntityRank element_rank = fem_meta_data.element_rank();

    int block_count = sset->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *block = sset->get_block(i);
      if (stk::io::include_entity(block)) {
	std::vector<int> side_ids ;
	std::vector<int> elem_side ;

	stk::mesh::Part * const side_block_part = meta.get_part(block->name());
	stk::mesh::EntityRank side_rank = side_block_part->primary_entity_rank();

	block->get_field_data("ids", side_ids);
	block->get_field_data("element_side", elem_side);

	assert(side_ids.size() * 2 == elem_side.size());
	stk::mesh::PartVector add_parts( 1 , side_block_part );

	size_t side_count = side_ids.size();
	std::vector<stk::mesh::Entity*> sides(side_count);
	for(size_t is=0; is<side_count; ++is) {

	  stk::mesh::Entity* const elem = bulk.get_entity(element_rank, elem_side[is*2]);

	  // If NULL, then the element was probably assigned to an
	  // element block that appears in the database, but was
	  // subsetted out of the analysis mesh. Only process if
	  // non-null.
	  if (elem != NULL) {
	    // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
	    // Hence the '-1' in the following line.
	    int side_ordinal = elem_side[is*2+1] - 1 ;

	    stk::mesh::Entity *side = NULL;
	    if (side_rank == 2) {
	      side = &stk::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);
	    } else {
	      side = &stk::mesh::fem::declare_element_edge(bulk, side_ids[is], *elem, side_ordinal);
	    }
	    bulk.change_entity_parts( *side, add_parts );
	    sides[is] = side;
	  } else {
	    sides[is] = NULL;
	  }
	}

	const stk::mesh::Field<double, stk::mesh::ElementNode> *df_field =
	  stk::io::get_distribution_factor_field(*side_block_part);
	if (df_field != NULL) {
	  stk::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
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
  // ========================================================================
  void get_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
		      stk::mesh::EntityRank part_type,
		      Ioss::GroupingEntity *io_entity,
		      Ioss::Field::RoleType filter_role)
  {
    std::vector<stk::mesh::Entity*> entities;
    stk::io::get_entity_list(io_entity, part_type, bulk, entities);

    stk::mesh::MetaData& meta = stk::mesh::MetaData::get(part);
    stk::mesh::Part &universal = meta.universal_part();
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (stk::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
	stk::io::field_data_from_ioss(f, entities, io_entity, f->name());
      }
    }
  }

  void process_input_request(Ioss::Region &region,
			     stk::mesh::BulkData &bulk,
			     int step)
  {
    region.begin_state(step);

    // Special processing for nodeblock (all nodes in model)...
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

    // ??? Get field data from nodeblock...
    get_field_data(bulk, meta.universal_part(), stk::mesh::fem::FEMMetaData::NODE_RANK,
		   region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

    const stk::mesh::PartVector & all_parts = meta.get_parts();
    for ( stk::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk::mesh::Part * const part = *ip;

      const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk::io::is_part_io_part(*part)) {
	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = region.get_entity(part->name());
	if (entity != NULL) {
	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    assert(sset != NULL);
	    int block_count = sset->block_count();
	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *side_block = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
		get_field_data(bulk, *part,
                               part_rank,
			       side_block, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    get_field_data(bulk, *part,
                           part_rank,
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

  void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
		      stk::mesh::EntityRank part_type,
		      Ioss::GroupingEntity *io_entity,
		      Ioss::Field::RoleType filter_role)
  {
    std::vector<stk::mesh::Entity*> entities;
    stk::io::get_entity_list(io_entity, part_type, bulk, entities);

    stk::mesh::MetaData& meta = stk::mesh::MetaData::get(part);
    stk::mesh::Part &universal = meta.universal_part();
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (stk::io::is_valid_part_field(f, part_type, part, universal, filter_role)) {
	stk::io::field_data_to_ioss(f, entities, io_entity, f->name(), filter_role);
      }
    }
  }

  void process_output_request(Ioss::Region &region,
			      stk::mesh::BulkData &bulk,
			      int step)
  {
    region.begin_state(step);
    // Special processing for nodeblock (all nodes in model)...
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

    put_field_data(bulk, meta.universal_part(), stk::mesh::fem::FEMMetaData::NODE_RANK,
		   region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

    const stk::mesh::PartVector & all_parts = meta.get_parts();
    for ( stk::mesh::PartVector::const_iterator
	    ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

      stk::mesh::Part * const part = *ip;

      const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Check whether this part should be output to results database.
      if (stk::io::is_part_io_part(*part)) {

	// Get Ioss::GroupingEntity corresponding to this part...
	Ioss::GroupingEntity *entity = region.get_entity(part->name());
	if (entity != NULL) {

	  if (entity->type() == Ioss::SIDESET) {
	    Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
	    assert(sset != NULL);
	    int block_count = sset->block_count();

	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *side_block = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
		put_field_data(bulk, *part, part_rank,
			       side_block, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    put_field_data(bulk, *part, part_rank,
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
}

// ========================================================================
#include <boost/program_options.hpp>

#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/ProgramOptions.hpp>

  namespace bopt = boost::program_options;
  int main(int argc, char** argv)
  {
    //----------------------------------
    // Broadcast argc and argv to all processors.

    stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

    stk::BroadcastArg b_arg(comm, argc, argv);

    //----------------------------------
    // Process the broadcast command line arguments

    bopt::options_description desc("options");

    desc.add_options()
      ("help,h",        "produce help message")
      ("mesh",         bopt::value<std::string>(), "mesh file" )
      ("decomposition,D", bopt::value<std::string>(), "decomposition method" )
      ("directory,d",  bopt::value<std::string>(), "working directory" )
      ("output-log,o", bopt::value<std::string>(), "output log path" )
      ("runtest,r",    bopt::value<std::string>(), "runtest pid file" );

    stk::get_options_description().add(desc);

    bopt::variables_map &vm = stk::get_variables_map();
    try {
      bopt::store(bopt::parse_command_line(b_arg.m_argc, b_arg.m_argv, desc), vm);
      bopt::notify(vm);
    }
    catch (std::exception & /* x */) {
      std::exit(1);
    }

    if (vm.count("help")) {
      std::cout << desc << "\n";
      std::exit(EXIT_SUCCESS);
    }

    //----------------------------------

    if ( vm.count("mesh") ) {
      std::string in_filename = boost::any_cast<std::string>(vm["mesh"].value());
      std::string out_filename = in_filename + ".out";
      std::string decomp_method;
      if (vm.count("decomposition")) {
	decomp_method = boost::any_cast<std::string>(vm["decomposition"].value());
      }
      stk_example_io::io_example(comm, in_filename, out_filename, decomp_method );
    } else {
      std::cout << "OPTION ERROR: The '--mesh <filename>' option is required!\n";
      std::exit(EXIT_FAILURE);
    }
    stk::parallel_machine_finalize();

    return 0;
  }

/**
 * \}
 */
