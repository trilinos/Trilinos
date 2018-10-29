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

#include <iostream>
#include <assert.h>

#include <stk_util/parallel/Parallel.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>

static const size_t spatial_dimension = 3;

namespace stk_examples {
namespace app {

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
      meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::put_field_on_mesh( coord_field, meta.universal_part(),
                          spatial_dim, nullptr);

    /// \todo IMPLEMENT truly handle fields... For this case we
    /// are just defining a field for each transient field that is
    /// present in the mesh...
    stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk::topology::NODE_RANK);
  }

  // ========================================================================
  void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
    stk::io::default_part_processing(elem_blocks, meta);

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
	/// \todo IMPLEMENT truly handle attribute fields... For this case we
	/// are just defining a field for each attribute field that
	/// is present in the mesh...
	stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE,
                                  *part, part_rank );

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                  *part, part_rank);

	std::cout << entity->type_string() << ": " << entity->name()
		  << " , celltop = " << meta.get_cell_topology(*part).getCellTopologyData()->name
		  << std::endl ;
      }
    }
  }

  // ========================================================================
  void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
    stk::io::default_part_processing(node_sets, meta);

    /// \todo REFACTOR should "distribution_factor" be a default field
    /// that is automatically declared on all objects that it exists
    /// on as is done in current framework?
    stk::mesh::Field<double> & distribution_factors_field =
      meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors");

    /// \todo REFACTOR How to associate distribution_factors field
    /// with the nodeset part if a node is a member of multiple
    /// nodesets

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
	it != node_sets.end(); ++it) {
      Ioss::NodeSet *entity = *it;

      if (stk::io::include_entity(entity)) {
	stk::mesh::Part* const part = meta.get_part(entity->name());
	assert(part != NULL);
	assert(entity->field_exists("distribution_factors"));

        const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

	stk::mesh::put_field_on_mesh(distribution_factors_field, *part, nullptr);

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                  *part, part_rank);
      }
    }
  }

  // ========================================================================
  void process_surface_entity(Ioss::SideSet *sset, stk::mesh::MetaData &meta)
  {
    assert(sset->type() == Ioss::SIDESET);
    const Ioss::SideBlockContainer& blocks = sset->get_side_blocks();
    stk::io::default_part_processing(blocks, meta);

    stk::mesh::Part* const sideset_part = meta.get_part(sset->name());
    assert(sideset_part != NULL);

    stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
    bool surface_df_defined = false; // Has the surface df field been defined yet?


    int block_count = sset->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *sideblock = sset->get_block(i);
      if (stk::io::include_entity(sideblock)) {
	std::cout << sideblock->type_string() << " " << sideblock->name() << "\n";
	stk::mesh::Part * const sideblock_part = meta.get_part(sideblock->name());
	assert(sideblock_part != NULL);
	meta.declare_part_subset(*sideset_part, *sideblock_part);

        const stk::mesh::EntityRank part_rank = sideblock_part->primary_entity_rank();

	/// \todo REFACTOR How to associate distribution_factors field
	/// with the sideset part if a face is a member of multiple
	/// sidesets

	if (sideblock->field_exists("distribution_factors")) {
	  if (!surface_df_defined) {
	    std::string field_name = sset->name() + "_distribution_factors";
	    distribution_factors_field =
	      &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(static_cast<stk::topology::rank_t>(part_rank), field_name);
	    stk::io::set_distribution_factor_field(*sideset_part, *distribution_factors_field);
	    surface_df_defined = true;
	  }
	  stk::io::set_distribution_factor_field(*sideblock_part, *distribution_factors_field);
	  int side_node_count = sideblock->topology()->number_nodes();
	  stk::mesh::put_field_on_mesh(*distribution_factors_field,
                               *sideblock_part, side_node_count, nullptr);
	}

	/// \todo IMPLEMENT truly handle fields... For this case we
	/// are just defining a field for each transient field that is
	/// present in the mesh...
	stk::io::define_io_fields(sideblock, Ioss::Field::TRANSIENT,
                                  *sideblock_part, part_rank );
      }
    }
  }

  // ========================================================================
  void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
  {
    const Ioss::SideSetContainer& side_sets = region.get_sidesets();
    stk::io::default_part_processing(side_sets, meta);

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
	it != side_sets.end(); ++it) {
      Ioss::SideSet *entity = *it;

      if (stk::io::include_entity(entity)) {
	process_surface_entity(entity, meta);
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

    std::vector<stk::mesh::Entity> nodes;
    stk::io::get_input_entity_list(nb, stk::topology::NODE_RANK, bulk, nodes);

    /// \todo REFACTOR Application would probably store this field
    /// (and others) somewhere after the declaration instead of
    /// looking it up each time it is needed.
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
    stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
      meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

    stk::io::field_data_from_ioss(bulk, coord_field, nodes, nb, "mesh_model_coordinates");
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
	stk::mesh::Part* const part = meta.get_part(name);
	assert(part != NULL);

        const CellTopologyData* cell_topo = stk::mesh::MetaData::get(bulk).get_cell_topology(*part).getCellTopologyData();
	assert(cell_topo != NULL);

	std::vector<int> elem_ids ;
	std::vector<int> connectivity ;

	entity->get_field_data("ids", elem_ids);
	entity->get_field_data("connectivity", connectivity);

	int element_count = elem_ids.size();
	int nodes_per_elem = cell_topo->node_count ;

	std::vector<stk::mesh::EntityId> connectivity2(nodes_per_elem);

        std::vector<int>::const_iterator connBegin = connectivity.begin();
	std::vector<stk::mesh::Entity> elements(element_count);
	for(int i=0; i<element_count; ++i, connBegin += nodes_per_elem) {
          std::copy(connBegin, connBegin + nodes_per_elem, connectivity2.begin());
	  /// \todo REFACTOR cast from int to unsigned is unsafe and ugly.
	  /// change function to take int[] argument.
	  elements[i] = stk::mesh::declare_element(bulk, *part, elem_ids[i], connectivity2);
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
	  stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEMENT_RANK, *I);
	  stk::io::field_data_from_ioss(bulk, field, elements, entity, *I);

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

	std::vector<stk::mesh::Entity> nodes(node_count);
	for(int i=0; i<node_count; ++i) {
	  nodes[i] = bulk.get_entity( stk::topology::NODE_RANK, node_ids[i] );
	  if (bulk.is_valid(nodes[i]))
	    bulk.declare_node(node_ids[i], add_parts);
	}


	/// \todo REFACTOR Application would probably store this field
	/// (and others) somewhere after the declaration instead of
	/// looking it up each time it is needed.
	stk::mesh::Field<double> *df_field =
	  meta.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors");

	if (df_field != NULL) {
	  stk::io::field_data_from_ioss(bulk, df_field, nodes, entity, "distribution_factors");
	}
      }
    }
  }

  // ========================================================================
  void process_surface_entity(const Ioss::SideSet* io, stk::mesh::BulkData & bulk)
  {
    assert(io->type() == Ioss::SIDESET);
    const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

    const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

    int block_count = io->block_count();
    for (int i=0; i < block_count; i++) {
      Ioss::SideBlock *block = io->get_block(i);
      if (stk::io::include_entity(block)) {

	stk::mesh::Part * const fb_part = meta.get_part(block->name());

	std::vector<int> elem_side ;
	block->get_field_data("element_side", elem_side);

	stk::mesh::PartVector add_parts( 1 , fb_part );

        int side_count = elem_side.size() / 2;
	std::vector<stk::mesh::Entity> sides(side_count);
	for(int is=0; is<side_count; ++is) {

	  stk::mesh::Entity const elem = bulk.get_entity(element_rank, elem_side[is*2]);
	  // If NULL, then the element was probably assigned to an
	  // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
	  // Hence the '-1' in the following line.
	  int side_ordinal = elem_side[is*2+1] - 1 ;

	  // element block that appears in the database, but was
	  // subsetted out of the analysis mesh. Only process if
	  // non-null.
	  if (bulk.is_valid(elem)) {
	    stk::mesh::Entity side = bulk.declare_element_side(elem, side_ordinal, add_parts);
	    bulk.change_entity_parts( side, add_parts );
	    sides[is] = side;
	  } else {
	    sides[is] = stk::mesh::Entity();
	  }
	}

	const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*fb_part);
	if (df_field != NULL) {
	  stk::io::field_data_from_ioss(bulk, df_field, sides, block, "distribution_factors");
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
    std::vector<stk::mesh::Entity> entities;
    stk::io::get_input_entity_list(io_entity, part_type, bulk, entities);

    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
	stk::io::field_data_from_ioss(bulk, f, entities, io_entity, f->name());
      }
    }
  }

  void process_input_request(Ioss::Region &region,
			     stk::mesh::BulkData &bulk,
			     int step)
  {
    region.begin_state(step);

    // Special processing for nodeblock (all nodes in model)...
    const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);

    // ??? Get field data from nodeblock...
    app::get_field_data(bulk, meta.universal_part(), stk::topology::NODE_RANK,
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

  void put_field_data(stk::io::OutputParams& params, stk::mesh::Part &part,
		      stk::mesh::EntityRank part_type,
		      Ioss::GroupingEntity *io_entity,
		      Ioss::Field::RoleType filter_role)
  {
    std::vector<stk::mesh::Entity> entities;
    stk::io::get_output_entity_list(io_entity, part_type, params, entities);

    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I; ++I;
      if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
	stk::io::field_data_to_ioss(params.bulk_data(), f, entities, io_entity, f->name(), filter_role);
      }
    }
  }

  void process_output_request(Ioss::Region &region,
			     stk::mesh::BulkData &bulk,
			     int step)
  {
    region.begin_state(step);
    // Special processing for nodeblock (all nodes in model)...
    const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);

    stk::io::OutputParams params(region, bulk);
    app::put_field_data(params, meta.universal_part(), stk::topology::NODE_RANK,
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
	    int block_count = sset->block_count();

	    for (int i=0; i < block_count; i++) {
	      Ioss::SideBlock *fb = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
	      app::put_field_data(params, *part, part_rank,
				  fb, Ioss::Field::TRANSIENT);
	    }
	  } else {
	    app::put_field_data(params, *part, part_rank,
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
  void example_io_2( stk::ParallelMachine comm,
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

    //----------------------------------
    // Process Entity Types. Subsetting is possible.
    stk::mesh::MetaData meta_data( spatial_dimension );
    app::process_elementblocks(in_region, meta_data);
    app::process_nodeblocks(in_region,    meta_data);
    app::process_sidesets(in_region,      meta_data);
    app::process_nodesets(in_region,      meta_data);

    //----------------------------------
    // Done populating meta data, commit and create bulk data
    meta_data.commit();

    //----------------------------------
    // Process Bulkdata for all Entity Types. Subsetting is possible.
    stk::mesh::BulkData bulk_data(meta_data, comm);
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
    /// fields to stk::io.  This maybe works with only a single field,
    /// but adding attributes, distribution factors and perhaps other
    /// fields will make this unwieldy.
    stk::io::OutputParams params(out_region, bulk_data);
    stk::io::define_output_db(params, {}, &in_region);
    stk::io::write_output_db(params);

    // ------------------------------------------------------------------------
    /// \todo REFACTOR A real app would register a subset of the
    /// fields on the mesh database as fields that the app would want
    // read at one or all or specified steps.  In this example, all
    // fields existing on the input mesh database are defined on the
    // parts in the stk::mesh.
    //
    // The real app would also only register a subset of the stk::mesh
    // fields as output fields and would probably have a mapping from
    // the internally used name to some name picked by the user. In
    // this example, all TRANSIENT fields defined on the stk::mesh are
    // output to the results database and the internal stk::mesh field
    // name is used as the name on the database....

    out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // Special processing for nodeblock (all nodes in model)...
    stk::io::ioss_add_fields(meta_data.universal_part(), stk::topology::NODE_RANK,
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
	    int block_count = sset->block_count();
	    for (int i=0; i < block_count; i++) {
	      Ioss::EntityBlock *fb = sset->get_block(i);
	      /// \todo REFACTOR Need filtering mechanism.
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
      app::process_input_request(in_region, bulk_data, step);

      // app::execute()

      // Write data from the stk::mesh fields out to the output database.a
      int out_step = out_region.add_state(time);
      app::process_output_request(out_region, bulk_data, out_step);
    }
    out_region.end_mode(Ioss::STATE_TRANSIENT);
    stk::io::delete_selector_property(out_region);
  }
} // namespace stk_examples

