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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <Ionit_Initializer.h>                          // for Initializer
#include <cstdlib>                                      // for exit, size_t
#include <algorithm>                                    // for copy, max
#include <cassert>                                      // for assert
#include <cstdint>                                      // for int64_t
#include <iostream>                                     // for operator<<
#include <stdexcept>                                    // for runtime_error
#include <stk_io/IossBridge.hpp>                        // for include_entity
#include <stk_mesh/base/BulkData.hpp>                   // for BulkData
#include <stk_mesh/base/MeshBuilder.hpp>                // for MeshBuilder
#include <stk_mesh/base/FEMHelpers.hpp>                 // for declare_element
#include <stk_mesh/base/Field.hpp>                      // for Field
#include <stk_mesh/base/MetaData.hpp>                   // for MetaData, put...
#include <stk_util/command_line/CommandLineParser.hpp>  // for CommandLinePa...
#include <stk_util/parallel/Parallel.hpp>               // for parallel_mach...
#include <string>                                       // for string, opera...
#include <vector>                                       // for vector, vecto...
#include "Ioss_DBUsage.h"                               // for READ_MODEL
#include "Ioss_DatabaseIO.h"                            // for DatabaseIO
#include "Ioss_ElementBlock.h"                          // for ElementBlock
#include "Ioss_ElementTopology.h"                       // for ElementTopology
#include "Ioss_EntityType.h"                            // for SIDESET
#include "Ioss_Field.h"                                 // for Field, Field:...
#include "Ioss_GroupingEntity.h"                        // for GroupingEntity
#include "Ioss_IOFactory.h"                             // for NameList, IOF...
#include "Ioss_NodeBlock.h"                             // for NodeBlock
#include "Ioss_NodeSet.h"                               // for NodeSet
#include "Ioss_Property.h"                              // for Property
#include "Ioss_PropertyManager.h"                       // for PropertyManager
#include "Ioss_Region.h"                                // for Region, NodeB...
#include "Ioss_SideBlock.h"                             // for SideBlock
#include "Ioss_SideSet.h"                               // for SideSet, Side...
#include "Ioss_State.h"                                 // for STATE_DEFINE_...
#include "Ioss_Utils.h"                                 // for Utils
#include "Ioss_VariableType.h"                          // for VariableType
#include "mpi.h"                                        // for MPI_COMM_WORLD
#include "stk_io/OutputParams.hpp"                      // for OutputParams
#include "stk_mesh/base/Entity.hpp"                     // for Entity
#include "stk_mesh/base/FieldBase.hpp"                  // for FieldBase
#include "stk_mesh/base/Part.hpp"                       // for Part
#include "stk_mesh/base/Types.hpp"                      // for PartVector
#include "stk_topology/topology.hpp"                    // for topology, top...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

/** \addtogroup stk_io_module
 * \{
 */

/**
 * NOTE NOTE NOTE
 *
 * This code demonstrates the low-level IossBridge capabilities.
 * For a higher-level interface, see io_mesh_read_write_example.cpp which uses
 * the StkMeshIoBroker class defined in StkMeshIoBroker.{hpp|cpp}.
 *
 * NOTE NOTE NOTE
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
 *
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
void io_example( const std::string& in_filename,
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
                                                  MPI_COMM_WORLD, properties);
  if (dbi == nullptr || !dbi->ok()) {
    std::cerr  << "ERROR: Could not open database '" << in_filename << "' of type '" << dbtype << "'\n";
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Reading input file:   " << in_filename << "\n";
  // NOTE: 'in_region' owns 'dbi' pointer at this time...
  Ioss::Region in_region(dbi, "input_model");

  // SUBSETTING PARSING/PREPROCESSING...
  // Just an example of how application could control whether an
  // entity is subsetted or not...

  //----------------------------------
  // Process Entity Types. Subsetting is possible.

  static size_t spatial_dimension = in_region.get_property("spatial_dimension").get_int();

  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatial_dimension);
  std::shared_ptr<stk::mesh::BulkData> bulk_data = builder.create();

  stk::mesh::MetaData& fem_meta_data = bulk_data->mesh_meta_data();
  process_elementblocks(in_region, fem_meta_data);
  process_nodeblocks(in_region,    fem_meta_data);
  process_sidesets(in_region,      fem_meta_data);
  process_nodesets(in_region,      fem_meta_data);

  //----------------------------------
  // Done populating meta data, commit and create bulk data
  fem_meta_data.commit();

  //----------------------------------
  // Process Bulkdata for all Entity Types. Subsetting is possible.
  bulk_data->modification_begin();
  process_elementblocks(in_region, *bulk_data);
  process_nodeblocks(in_region,    *bulk_data);
  process_sidesets(in_region,      *bulk_data);
  process_nodesets(in_region,      *bulk_data);
  bulk_data->modification_end();

  //----------------------------------
  // OUTPUT...Create the output "mesh" portion

  std::cout << "Creating output file: " << out_filename << "\n";
  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, out_filename,
                                                  Ioss::WRITE_RESULTS,
                                                  MPI_COMM_WORLD);
  if (dbo == nullptr || !dbo->ok()) {
    std::cerr << "ERROR: Could not open results database '" << out_filename << "' of type '" << dbtype << "'\n";
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'out_region' owns 'dbo' pointer at this time...
  Ioss::Region out_region(dbo, "results_output");

  stk::io::OutputParams params(out_region, *bulk_data);
  stk::io::define_output_db(params, {}, &in_region);
  stk::io::write_output_db(params);

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
  stk::io::ioss_add_fields(fem_meta_data.universal_part(), stk::topology::NODE_RANK,
                           out_region.get_node_blocks()[0],
      Ioss::Field::TRANSIENT);

  const stk::mesh::PartVector & all_parts = fem_meta_data.get_parts();
  for ( stk::mesh::PartVector::const_iterator
        ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = out_region.get_entity(part->name());
      if (entity != nullptr) {
        if (entity->type() == Ioss::SIDESET) {
          Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          assert(sset != nullptr);
          int block_count = sset->block_count();
          for (int i=0; i < block_count; i++) {
            Ioss::SideBlock *fb = sset->get_block(i);
            stk::io::ioss_add_fields(*part, part_rank, fb, Ioss::Field::TRANSIENT);
          }
        } else {
          stk::io::ioss_add_fields(*part, part_rank, entity, Ioss::Field::TRANSIENT);
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
    process_input_request(in_region, *bulk_data, step);

    // execute()

    // Write data from the stk::mesh fields out to the output database.a
    int out_step = out_region.add_state(time);
    process_output_request(out_region, *bulk_data, out_step);
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

  stk::mesh::Field<double> & coord_field = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates");

  stk::mesh::put_field_on_mesh(coord_field, meta.universal_part(), spatial_dim, nullptr);

  /** \todo IMPLEMENT truly handle fields... For this case we are
     * just defining a field for each transient field that is present
     * in the mesh...
     */
  stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(), stk::topology::NODE_RANK);
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  stk::io::default_part_processing(elem_blocks, meta);

  // Parts were created above, now handle element block specific
  // information (topology, attributes, ...);
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin(); it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());
      STKIORequire(part != nullptr);

      const stk::mesh::EntityRank part_rank = part->primary_entity_rank();

      // Element Block attributes (if any)...
      /** \todo IMPLEMENT truly handle attribute fields... For this
         * case we are just defining a field for each attribute field
         * that is present in the mesh...
         */
      stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE, *part, part_rank);

      /** \todo IMPLEMENT truly handle fields... For this case we
         * are just defining a field for each transient field that is
         * present in the mesh...
         */
      stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT, *part, part_rank);
    }
  }
}

// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta);

  /** \todo REFACTOR should "distribution_factor" be a default field
     * that is automatically declared on all objects that it exists
     * on as is done in current framework?
     */
  stk::mesh::Field<double> & distribution_factors_field = meta.declare_field<double>(stk::topology::NODE_RANK,
                                                                                     "distribution_factors");

  /** \todo REFACTOR How to associate distribution_factors field
     * with the nodeset part if a node is a member of multiple
     * nodesets
     */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin(); it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());
      STKIORequire(part != nullptr);
      STKIORequire(entity->field_exists("distribution_factors"));

      stk::mesh::put_field_on_mesh(distribution_factors_field, *part, nullptr);

      /** \todo IMPLEMENT truly handle fields... For this case we
         * are just defining a field for each transient field that is
         * present in the mesh...
         */
      stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT, *part, part->primary_entity_rank());
    }
  }
}

// ========================================================================
void process_surface_entity(Ioss::SideSet *sset, stk::mesh::MetaData &meta,
                            stk::mesh::EntityRank sset_rank)
{
  assert(sset->type() == Ioss::SIDESET);
  Ioss::SideSet *fs = dynamic_cast<Ioss::SideSet *>(sset);
  assert(fs != nullptr);
  const Ioss::SideBlockContainer& blocks = fs->get_side_blocks();
  stk::io::default_part_processing(blocks, meta);

  stk::mesh::Part* const fs_part = meta.get_part(sset->name());
  STKIORequire(fs_part != nullptr);

  stk::mesh::Field<double> *distribution_factors_field = nullptr;
  bool surface_df_defined = false; // Has the surface df field been defined yet?


  int block_count = sset->block_count();
  for (int i=0; i < block_count; i++) {
    Ioss::SideBlock *side_block = sset->get_block(i);
    if (stk::io::include_entity(side_block)) {
      stk::mesh::Part * const side_block_part = meta.get_part(side_block->name());
      STKIORequire(side_block_part != nullptr);
      meta.declare_part_subset(*fs_part, *side_block_part);

      const stk::mesh::EntityRank part_rank = side_block_part->primary_entity_rank();

      if (side_block->field_exists("distribution_factors")) {
        if (!surface_df_defined) {
          std::string field_name = sset->name() + "_distribution_factors";
          distribution_factors_field =
              &meta.declare_field<double>(static_cast<stk::topology::rank_t>(part_rank), field_name);
          stk::io::set_distribution_factor_field(*fs_part, *distribution_factors_field);
          surface_df_defined = true;
        }
        stk::io::set_distribution_factor_field(*side_block_part, *distribution_factors_field);
        int side_node_count = side_block->topology()->number_nodes();
        stk::mesh::put_field_on_mesh(*distribution_factors_field, *side_block_part, side_node_count, nullptr);
      }

      /** \todo IMPLEMENT truly handle fields... For this case we
         * are just defining a field for each transient field that is
         * present in the mesh...
         */
      stk::io::define_io_fields(side_block, Ioss::Field::TRANSIENT, *side_block_part, part_rank);
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const stk::mesh::EntityRank side_rank = meta.side_rank();

  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  stk::io::default_part_processing(side_sets, meta);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin(); it != side_sets.end(); ++it) {
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

  std::vector<stk::mesh::Entity> nodes = stk::io::get_input_entity_list(nb, stk::topology::NODE_RANK, bulk);

  /** \todo REFACTOR Application would probably store this field
     * (and others) somewhere after the declaration instead of
     * looking it up each time it is needed.
     */
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Field<double> *coord_field = meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");

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
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != nullptr);

      const stk::topology topo = part->topology();
      if (topo == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
        throw std::runtime_error( msg.str() );
      }

      std::vector<int> elem_ids ;
      std::vector<int> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = topo.num_nodes();

      stk::mesh::EntityIdVector connectivity2(nodes_per_elem);

      std::vector<int>::const_iterator connBegin = connectivity.begin();
      std::vector<stk::mesh::Entity> elements(element_count);
      for(size_t i=0; i<element_count; ++i, connBegin += nodes_per_elem) {
        std::copy(connBegin, connBegin + nodes_per_elem, connectivity2.begin());
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
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != nullptr);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<int> node_ids ;
      int node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity> nodes(node_count);
      for(int i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity( stk::topology::NODE_RANK, node_ids[i] );
        if (bulk.is_valid(nodes[i])) {
          bulk.declare_node(node_ids[i], add_parts );
        }
      }

      /** \todo REFACTOR Application would probably store this field
         * (and others) somewhere after the declaration instead of
         * looking it up each time it is needed.
         */
      stk::mesh::Field<double> *df_field = meta.get_field<double>(stk::topology::NODE_RANK, "distribution_factors");

      if (df_field != nullptr) {
        stk::io::field_data_from_ioss(bulk, df_field, nodes, entity, "distribution_factors");
      }
    }
  }
}

// ========================================================================
void process_surface_entity(const Ioss::SideSet* sset ,
                            stk::mesh::BulkData & bulk)
{
  assert(sset->type() == Ioss::SIDESET);

  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;
  const int64_t ten = 10;

  int block_count = sset->block_count();
  for (int i=0; i < block_count; i++) {
    Ioss::SideBlock *block = sset->get_block(i);
    if (stk::io::include_entity(block)) {

      stk::mesh::Part * const side_block_part = meta.get_part(block->name());
      stk::mesh::EntityRank side_rank = side_block_part->primary_entity_rank();

      std::vector<int> elem_side ;
      block->get_field_data("element_side", elem_side);

      stk::mesh::PartVector add_parts( 1 , side_block_part );

      size_t side_count = elem_side.size() / 2;
      std::vector<stk::mesh::Entity> sides(side_count);
      for(size_t is=0; is<side_count; ++is) {

        stk::mesh::Entity const elem = bulk.get_entity(element_rank, elem_side[is*2]);

        // If NULL, then the element was probably assigned to an
        // element block that appears in the database, but was
        // subsetted out of the analysis mesh. Only process if
        // non-null.
        if (bulk.is_valid(elem)) {
          // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
          // Hence the '-1' in the following line.
          int side_ordinal = elem_side[is*2+1] - 1 ;

          stk::mesh::Entity side = stk::mesh::Entity();
          if (side_rank == 2) {
            side = bulk.declare_element_side(elem, side_ordinal, add_parts);
          } else {
            int64_t side_id = ten * elem_side[is*2+0] + elem_side[is*2+1];
            side = stk::mesh::declare_element_edge(bulk, side_id, elem, side_ordinal);
            bulk.change_entity_parts( side, add_parts );
          }
          sides[is] = side;
        } else {
          sides[is] = stk::mesh::Entity();
        }
      }

      const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*side_block_part);

      if (df_field != nullptr) {
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
  std::vector<stk::mesh::Entity> entities = stk::io::get_input_entity_list(io_entity, part_type, bulk);

  stk::mesh::MetaData& meta = stk::mesh::MetaData::get(part);
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
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  // ??? Get field data from nodeblock...
  get_field_data(bulk, meta.universal_part(), stk::topology::NODE_RANK,
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
      if (entity != nullptr) {
        if (entity->type() == Ioss::SIDESET) {
          Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          assert(sset != nullptr);
          int block_count = sset->block_count();
          for (int i=0; i < block_count; i++) {
            Ioss::SideBlock *side_block = sset->get_block(i);
            /// \todo REFACTOR Need filtering mechanism.
            get_field_data(bulk, *part, part_rank, side_block, Ioss::Field::TRANSIENT);
          }
        } else {
          get_field_data(bulk, *part, part_rank, entity, Ioss::Field::TRANSIENT);
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
  std::vector<stk::mesh::Entity> entities;
  stk::io::OutputParams params(bulk);
  stk::io::get_output_entity_list(io_entity, part_type, params, entities);

  stk::mesh::MetaData& meta = stk::mesh::MetaData::get(part);
  const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
      stk::io::field_data_to_ioss(bulk, f, entities, io_entity, f->name(), filter_role);
    }
  }
}

void process_output_request(Ioss::Region &region,
                            stk::mesh::BulkData &bulk,
                            int step)
{
  region.begin_state(step);
  // Special processing for nodeblock (all nodes in model)...
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  put_field_data(bulk, meta.universal_part(), stk::topology::NODE_RANK,
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
      if (entity != nullptr) {

        if (entity->type() == Ioss::SIDESET) {
          Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          assert(sset != nullptr);
          int block_count = sset->block_count();

          for (int i=0; i < block_count; i++) {
            Ioss::SideBlock *side_block = sset->get_block(i);
            /// \todo REFACTOR Need filtering mechanism.
            put_field_data(bulk, *part, part_rank, side_block, Ioss::Field::TRANSIENT);
          }
        } else {
          put_field_data(bulk, *part, part_rank, entity, Ioss::Field::TRANSIENT);
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

int main(int argc, const char** argv)
{
  stk::parallel_machine_init(&argc, const_cast<char***>(&argv));

  stk::CommandLineParser cmdLine;
  cmdLine.add_required<std::string>({"mesh", "m", "mesh file"});
  cmdLine.add_optional<std::string>({"decomposition", "D", "decomposition method"}, "not provided");
  cmdLine.add_optional<std::string>({"directory", "d", "working directory"}, "not provided");
  cmdLine.add_optional<std::string>({"output-log", "o", "output log path"}, "not provided");
  cmdLine.add_optional<std::string>({"runtest", "r", "runtest pid file"}, "not provided");

  stk::CommandLineParser::ParseState parseState = cmdLine.parse(argc, argv);

  if (stk::CommandLineParser::ParseComplete == parseState) {
    std::string in_filename = cmdLine.get_option_value<std::string>("mesh");
    std::string out_filename = in_filename + ".out";
    std::string decomp_method;
    if (cmdLine.is_option_provided("decomposition")) {
      decomp_method = cmdLine.get_option_value<std::string>("decomposition");
    }

    stk_example_io::io_example(in_filename, out_filename, decomp_method );
  }
  else if (stk::CommandLineParser::ParseHelpOnly == parseState ||
           stk::CommandLineParser::ParseError == parseState) {
    std::cout << cmdLine.get_usage() << "\n";
  }
  else if (stk::CommandLineParser::ParseVersionOnly == parseState) {
    std::cout << "This program doesn't have a version number." << "\n";
  }

  stk::parallel_machine_finalize();

  return 0;
}

/**
 * \}
 */
