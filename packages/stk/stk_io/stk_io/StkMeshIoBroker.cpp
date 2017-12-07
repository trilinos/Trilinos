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

#include <stk_io/StkMeshIoBroker.hpp>
#include <assert.h>                     // for assert
#include <Ionit_Initializer.h>     // for Initializer
#include <stdint.h>                     // for int64_t
#include <stdlib.h>                     // for exit, EXIT_FAILURE
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <iterator>                     // for back_insert_iterator, etc
#include <limits>
#include <stdexcept>                    // for runtime_error
#include <stk_io/InputFile.hpp>         // for InputFile
#include <stk_io/IossBridge.hpp>        // for FieldAndName, STKIORequire, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element_edge, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/Comm.hpp> 
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowErrorMsgIf
#include <utility>                      // for make_pair, pair
#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"               // for DatabaseUsage, etc
#include "Ioss_DatabaseIO.h"            // for DatabaseIO
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_ElementTopology.h"       // for ElementTopology
#include "Ioss_EntityType.h"            // for EntityType::SIDEBLOCK, etc
#include "Ioss_Field.h"                 // for Field, Field::BasicType, etc
#include "Ioss_GroupingEntity.h"        // for GroupingEntity
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_NodeSet.h"               // for NodeSet
#include "Ioss_ParallelUtils.h"         // for ParallelUtils
#include "Ioss_Property.h"              // for Property
#include "Ioss_PropertyManager.h"       // for PropertyManager
#include "Ioss_Region.h"                // for Region, NodeSetContainer, etc
#include "Ioss_SideBlock.h"             // for SideBlock
#include "Ioss_SideSet.h"               // for SideSet, SideBlockContainer
#include "Ioss_State.h"
#include "Ioss_VariableType.h"          // for NameList, VariableType

#include "ProcessSetsOrBlocks.hpp"
#include "SidesetTranslator.hpp"
#include "StkIoUtils.hpp"

#include "Teuchos_RCP.hpp"              // for RCP::operator->, etc
#include "boost/any.hpp"                // for any_cast, any
#include "boost/cstdint.hpp"            // for int64_t
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose, etc
#include "stk_io/MeshField.hpp"         // for MeshField, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/FieldState.hpp"  // for FieldState
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/TopologyDimensions.hpp"  // for ElementNode
#include "stk_mesh/base/Types.hpp"      // for EntityId, PartVector, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.hpp"    // for topology::num_nodes
#include "stk_util/util/ParameterList.hpp"  // for Type, Type::DOUBLE, etc
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_util/parallel/CommSparse.hpp"


namespace stk {
namespace io {

  template <typename T>
  bool is_index_valid(const std::vector<T> &file_vector, size_t input_file_index)
  {
    bool invalid = file_vector.empty() ||
      input_file_index >= file_vector.size() ||
      Teuchos::is_null(file_vector[input_file_index]);
    return !invalid;
  }

  bool fieldOrdinalSort(const stk::io::FieldAndName& f1, const stk::io::FieldAndName &f2) {
    return f1.field()->mesh_meta_data_ordinal() < f2.field()->mesh_meta_data_ordinal();
  }

  std::pair<size_t, Ioss::Field::BasicType> get_io_parameter_size_and_type(const stk::util::ParameterType::Type type,
                                                                           const boost::any &value)
  {
    try {
      switch(type)  {
      case stk::util::ParameterType::INTEGER: {
        return std::make_pair(1, Ioss::Field::INTEGER);
      }

      case stk::util::ParameterType::INT64: {
        return std::make_pair(1, Ioss::Field::INT64);
      }

      case stk::util::ParameterType::DOUBLE: {
        return std::make_pair(1, Ioss::Field::REAL);
      }

      case stk::util::ParameterType::DOUBLEVECTOR: {
        std::vector<double> vec = boost::any_cast<std::vector<double> >(value);
        return std::make_pair(vec.size(), Ioss::Field::REAL);
      }

      case stk::util::ParameterType::INTEGERVECTOR: {
        std::vector<int> vec = boost::any_cast<std::vector<int> >(value);
        return std::make_pair(vec.size(), Ioss::Field::INTEGER);
      }

      case stk::util::ParameterType::INT64VECTOR: {
        std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(value);
        return std::make_pair(vec.size(), Ioss::Field::INT64);
      }

      default: {
        return std::make_pair(0, Ioss::Field::INVALID);
      }
      }
    }
    catch(...) {
      std::cerr << "ERROR: Actual type of parameter does not match the declared type. Something went wrong."
                << " Maybe you need to call add_global_ref() instead of add_global() or vice-versa.";
      throw;
    }
  }

  template <typename DataType>
  void internal_write_global(Teuchos::RCP<Ioss::Region> output_region, const std::string &globalVarName,
                             DataType globalVarData)
  {
    ThrowErrorMsgIf (Teuchos::is_null(output_region),
                     "There is no Output mesh region associated with this Mesh Data.");
    ThrowErrorMsgIf (output_region->get_state() != Ioss::STATE_TRANSIENT,
                     "The output region " << output_region->name() <<
                     " is not in the correct state for outputting data at this time.");
    ThrowErrorMsgIf (!output_region->field_exists(globalVarName),
                     "The field named '" << globalVarName << "' does not exist.");
    output_region->put_field_data(globalVarName, &globalVarData, sizeof(DataType));
  }

  template <typename DataType>
  void internal_write_global(Teuchos::RCP<Ioss::Region> output_region, const std::string &globalVarName,
                             std::vector<DataType> &globalVarData)
  {
    ThrowErrorMsgIf (Teuchos::is_null(output_region),
                     "There is no Output mesh region associated with this Mesh Data.");
    ThrowErrorMsgIf (output_region->get_state() != Ioss::STATE_TRANSIENT,
                     "The output region " << output_region->name() <<
                     " is not in the correct state for outputting data at this time.");
    ThrowErrorMsgIf (!output_region->field_exists(globalVarName),
                     "The field named '" << globalVarName << "' does not exist "
                     "on output region "  << output_region->name());
    size_t comp_count = output_region->get_fieldref(globalVarName).raw_storage()->component_count();
    ThrowErrorMsgIf (comp_count != globalVarData.size(),
                     "On output region "  << output_region->name() <<
                     ", the field named '" << globalVarName << "' was registered with size "
                     << comp_count
                     << " but the output size is " << globalVarData.size());

    output_region->put_field_data(globalVarName, globalVarData);
  }

  template <typename DataType>
  bool internal_read_global(Teuchos::RCP<Ioss::Region> input_region, const std::string &globalVarName,
                            DataType &globalVarData, Ioss::Field::BasicType iossType,
                            bool abort_if_not_found)
  {
    ThrowErrorMsgIf (Teuchos::is_null(input_region),
                     "There is no Input mesh region associated with this Mesh Data.");

    if (input_region->field_exists(globalVarName)) {
      input_region->get_fieldref(globalVarName).check_type(iossType);
      input_region->get_field_data(globalVarName, &globalVarData, sizeof(DataType));
      return true;
    }
    else {
      if (abort_if_not_found) {
        std::ostringstream msg;
        msg << "ERROR: The field named '" << globalVarName << "' does not exist "
            << "on input region "  << input_region->name();
        throw std::runtime_error( msg.str() );
      }
      return false;
    }
  }

  template <typename DataType>
  bool internal_read_global(Teuchos::RCP<Ioss::Region> input_region, const std::string &globalVarName,
                            std::vector<DataType> &globalVarData, Ioss::Field::BasicType iossType,
                            bool abort_if_not_found)
  {
    ThrowErrorMsgIf (Teuchos::is_null(input_region),
                     "There is no Input mesh region associated with this Mesh Data.");

    if (input_region->field_exists(globalVarName)) {
      input_region->get_fieldref(globalVarName).check_type(iossType);
      input_region->get_field_data(globalVarName, globalVarData);
      return true;
    }
    else {
      if (abort_if_not_found) {
        std::ostringstream msg;
        msg << "ERROR: The field named '" << globalVarName << "' does not exist "
            << "on input region "  << input_region->name();
        throw std::runtime_error( msg.str() );
      }
      return false;
    }
  }

  void internal_write_parameter(Teuchos::RCP<Ioss::Region> output_region,
                                const std::string &name, const boost::any &any_value,
                                stk::util::ParameterType::Type type)
  {
    try {
      switch(type)
        {
        case stk::util::ParameterType::INTEGER: {
          int value = boost::any_cast<int>(any_value);
          internal_write_global(output_region, name, value);
          break;
        }

        case stk::util::ParameterType::INT64: {
          int64_t value = boost::any_cast<int64_t>(any_value);
          internal_write_global(output_region, name, value);
          break;
        }

        case stk::util::ParameterType::DOUBLE: {
          double value = boost::any_cast<double>(any_value);
          internal_write_global(output_region, name, value);
          break;
        }

        case stk::util::ParameterType::DOUBLEVECTOR: {
          std::vector<double> vec = boost::any_cast<std::vector<double> >(any_value);
          internal_write_global(output_region, name, vec);
          break;
        }

        case stk::util::ParameterType::INTEGERVECTOR: {
          std::vector<int> vec = boost::any_cast<std::vector<int> >(any_value);
          internal_write_global(output_region, name, vec);
          break;
        }

        case stk::util::ParameterType::INT64VECTOR: {
          std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(any_value);
          internal_write_global(output_region, name, vec);
          break;
        }

        default: {
          std::cerr << "WARNING: '" << name << "' is not a supported type. It's value cannot be output." << std::endl;
          break;
        }
        }
    }
    catch(...) {
      std::cerr << "ERROR: Actual type of parameter named '" << name
                << "' does not match the declared type. Something went wrong."
                << " Maybe you need to call add_global_ref() instead of add_global() or vice-versa.";
      throw;
    }
  }

  void write_defined_global_any_fields(Teuchos::RCP<Ioss::Region> region,
                                       std::vector<stk::io::GlobalAnyVariable> &global_any_fields)
  {
    for (size_t i=0; i < global_any_fields.size(); i++) {
      const std::string &name = global_any_fields[i].m_name;
      const boost::any *value = global_any_fields[i].m_value;
      stk::util::ParameterType::Type type = global_any_fields[i].m_type;
      internal_write_parameter(region, name, *value, type);
    }
  }

  bool internal_read_parameter(Teuchos::RCP<Ioss::Region> input_region,
                               const std::string &globalVarName,
                               boost::any &any_value, stk::util::ParameterType::Type type,
                               bool abort_if_not_found)
  {
    bool success = false;
    switch(type)
      {
      case stk::util::ParameterType::INTEGER: {
        int value = 0;
        success = internal_read_global(input_region, globalVarName, value, Ioss::Field::INTEGER,
                                       abort_if_not_found);
        any_value = value;
        break;
      }

      case stk::util::ParameterType::INT64: {
        int64_t value = 0;
        success = internal_read_global(input_region, globalVarName, value, Ioss::Field::INT64,
                                       abort_if_not_found);
        any_value = value;
        break;
      }

      case stk::util::ParameterType::DOUBLE: {
        double value = 0;
        success = internal_read_global(input_region, globalVarName, value, Ioss::Field::REAL,
                                       abort_if_not_found);
        any_value = value;
        break;
      }

      case stk::util::ParameterType::DOUBLEVECTOR: {
        std::vector<double> vec;
        success = internal_read_global(input_region, globalVarName, vec, Ioss::Field::REAL,
                                       abort_if_not_found);
        any_value = vec;
        break;
      }

      case stk::util::ParameterType::INTEGERVECTOR: {
        std::vector<int> vec;
        success = internal_read_global(input_region, globalVarName, vec, Ioss::Field::INTEGER,
                                       abort_if_not_found);
        any_value = vec;
        break;
      }

      case stk::util::ParameterType::INT64VECTOR: {
        std::vector<int64_t> vec;
        success = internal_read_global(input_region, globalVarName, vec, Ioss::Field::INT64,
                                       abort_if_not_found);
        any_value = vec;
        break;
      }

      default: {
        std::cerr << "WARNING: '" << globalVarName << "' is not a supported type. It's value cannot be input." << std::endl;
        break;
      }
      }
    return success;
  }

  void internal_add_global(Teuchos::RCP<Ioss::Region> region,
                           const std::string &globalVarName, const std::string &storage,
                           Ioss::Field::BasicType dataType)
  {
    ThrowErrorMsgIf (region->field_exists(globalVarName),
                     "On region named " << region->name() <<
                     " Attempt to add global variable '" << globalVarName << "' twice.");

    region->field_add(Ioss::Field(globalVarName, dataType, storage, Ioss::Field::TRANSIENT, 1));
  }

  void internal_add_global(Teuchos::RCP<Ioss::Region> region,
                           const std::string &globalVarName, int component_count, Ioss::Field::BasicType dataType)
  {
    if (component_count == 1) {
      internal_add_global(region, globalVarName, "scalar", dataType);
    } else {
      std::ostringstream type;
      type << "Real[" << component_count << "]";
      internal_add_global(region, globalVarName, type.str(), dataType);
    }
  }

  size_t get_entities(stk::mesh::Part &part,
                      const stk::mesh::BulkData &bulk,
                      std::vector<stk::mesh::Entity> &entities,
                      const stk::mesh::Selector *subset_selector)
  {
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(part);

    stk::mesh::Selector own = meta.locally_owned_part();
    stk::mesh::Selector selector = part & own;
    if (subset_selector) selector &= *subset_selector;

    get_selected_entities(selector, bulk.buckets(type), entities);
    return entities.size();
  }

  bool is_skipped_attribute_field(const std::string &name, size_t numAttrFields)
  {
      return name == "attribute" && numAttrFields > 1;
  }

// ========================================================================
template <typename INT>
void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
  assert(sset->type() == Ioss::SIDESET);

  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  size_t block_count = sset->block_count();
  for (size_t i=0; i < block_count; i++) {
    Ioss::SideBlock *block = sset->get_block(i);
    if (stk::io::include_entity(block)) {

      stk::mesh::Part * const sb_part = meta.get_part(block->name());
      // Get topology of the sides being defined to see if they
      // are 'faces' or 'edges'.  This is needed since for shell-type
      // elements, (and actually all elements) a sideset can specify either a face or an edge...
      // For a quad shell, sides 1,2 are faces and 3,4,5,6 are edges.

      // NOTE: This assumes that the sides within a side block are homogenous.  If the side_set
      //       is not split into homogenous side_blocks, then the topology will not necessarily
      //       be the same and this could fail (a sideset of mixed edges and faces)
      int par_dimen = block->topology()->parametric_dimension();
      STKIORequire(par_dimen == 1 || par_dimen == 2);

      stk::mesh::EntityRank side_rank = par_dimen == 1 ? stk::topology::EDGE_RANK : stk::topology::FACE_RANK;

      // Would be nice to do:
      //    std::vector<stk::mesh::Entity> sides ;
      //    get_entities(*sb_part, bulk, sides, NULL);
      // But, we need the entities in the exact same order and count  as they appear on the
      // mesh file.  The get_entities can give them in different order and will ignore duplicated
      // "sides" that occur in some exodus files.

      std::vector<INT> elemSidePairs;
      block->get_field_data("element_side", elemSidePairs);
      size_t side_count = elemSidePairs.size()/2;

      std::vector<stk::mesh::Entity> sides;
      sides.reserve(side_count);

      for(size_t is = 0; is < side_count; ++is)
          sides.push_back(stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk, elemSidePairs[is*2], elemSidePairs[is*2+1]-1, side_rank));

      const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*sb_part);
      if (df_field != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field, sides, block, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      block->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(is_skipped_attribute_field(*I, names.size()))
          continue;
        stk::mesh::FieldBase *field = meta.get_field(side_rank, *I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, sides, block, *I);
      }
    }
  }
}



void process_surface_entity_df(const Ioss::SideSet* sset, stk::mesh::BulkData & bulk)
{
  if (stk::io::db_api_int_size(sset) == 4) {
    process_surface_entity_df<int>(sset, bulk);
  }
  else {
    process_surface_entity_df<int64_t>(sset, bulk);
  }
}



template <typename INT>
void process_node_coords_and_attributes(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  // Currently, all nodes found in the finite element mesh are defined
  // as nodes in the stk_mesh database. If some of the element blocks
  // are omitted, then there will be disconnected nodes defined.
  // However, if we only define nodes that are connected to elements,
  // then we risk missing "free" nodes that the user may want to have
  // existing in the model.
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  size_t node_count = nb->get_property("entity_count").get_int();

  std::vector<INT> ids;
  nb->get_field_data("ids", ids);

  std::vector<stk::mesh::Entity> nodes;
  nodes.reserve(node_count);
  for (size_t i=0; i < ids.size(); i++) {
    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, ids[i]);
    nodes.push_back(node);
  }

  // Temporary (2013/04/02) kluge for Salinas porting to stk-based mesh.
  // Salinas uses the "implicit id" which is the ordering of the nodes
  // in the "undecomposed" or "serial" mesh as user-visible ids
  // instead of the "global" ids. If there exists a stk-field with the
  // name "implicit_node_ids", then populate the field with the correct
  // data.
  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);
  stk::mesh::FieldBase *implicit_node_id_field = meta.get_field(stk::topology::NODE_RANK, "implicit_node_ids");
  if (implicit_node_id_field) {
    stk::io::field_data_from_ioss(bulk, implicit_node_id_field, nodes, nb, "implicit_ids");
  }


  stk::mesh::FieldBase const* coord_field = meta.coordinate_field();
  stk::io::field_data_from_ioss(bulk, coord_field, nodes, nb, "mesh_model_coordinates");

  // Add all attributes as fields.
  // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
  // named components of the 'attribute' field, so add them instead.
  Ioss::NameList names;
  nb->field_describe(Ioss::Field::ATTRIBUTE, &names);
  for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    if(is_skipped_attribute_field(*I, names.size()))
      continue;
    stk::mesh::FieldBase *field = meta.get_field(stk::topology::NODE_RANK, *I);
    if (field)
      stk::io::field_data_from_ioss(bulk, field, nodes, nb, *I);
  }
}

// ========================================================================

template <typename INT>
void process_elem_attributes_and_implicit_ids(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != NULL);

      stk::topology topo = part->topology();
      if (topo == stk::topology::INVALID_TOPOLOGY) {
        std::ostringstream msg ;
        msg << " INTERNAL_ERROR: Part " << part->name() << " has invalid topology";
        throw std::runtime_error( msg.str() );
      }

      // See if we need to get the list of elements...
      bool elements_needed = false;
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);

      stk::mesh::FieldBase *implicit_elem_id_field = meta.get_field(stk::topology::ELEMENT_RANK, "implicit_element_ids");
      if (implicit_elem_id_field) {
        elements_needed = true;
      } else {
        for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
          if(is_skipped_attribute_field(*I, names.size()))
            continue;
          stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEMENT_RANK, *I);
          if (field) {
            elements_needed = true;
            break;
          }
        }
      }

      if (!elements_needed)
        continue;

      std::vector<INT> elem_ids ;
      entity->get_field_data("ids", elem_ids);

      size_t element_count = elem_ids.size();
      std::vector<stk::mesh::Entity> elements;
      elements.reserve(element_count);

      for(size_t i=0; i<element_count; ++i) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elem_ids[i]);
        elements.push_back(elem);
      }

      // Temporary (2013/04/17) kluge for Salinas porting to stk-based mesh.
      // Salinas uses the "implicit id" which is the ordering of the nodes
      // in the "undecomposed" or "serial" mesh as user-visible ids
      // instead of the "global" ids. If there exists a stk-field with the
      // name "implicit_element_ids", then populate the field with the correct
      // data.
      if (implicit_elem_id_field) {
        stk::io::field_data_from_ioss(bulk, implicit_elem_id_field, elements, entity, "implicit_ids");
      }

      // Add all element attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(is_skipped_attribute_field(*I, names.size()))
          continue;
        stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEMENT_RANK, *I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, elements, entity, *I);
      }
    }
  }
}

// ========================================================================
// ========================================================================


// ========================================================================
template <typename INT>
void process_nodesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(bulk);

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      stk::mesh::Part* const part = meta.get_part(name);
      STKIORequire(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<INT> node_ids ;
      size_t node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity> nodes(node_count);
      for(size_t i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity(stk::topology::NODE_RANK, node_ids[i] );
        if (!bulk.is_valid(nodes[i])) {
          bulk.declare_entity(stk::topology::NODE_RANK, node_ids[i], add_parts );
        }
      }

      std::string distributionFactorsPerNodesetFieldName = "distribution_factors_" + part->name();

      stk::mesh::Field<double> *df_field_per_nodeset =
                meta.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, distributionFactorsPerNodesetFieldName);

      if (df_field_per_nodeset != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field_per_nodeset, nodes, entity, "distribution_factors");
      }

      // Add all attributes as fields.
      // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
      // named components of the 'attribute' field, so add them instead.
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if(is_skipped_attribute_field(*I, names.size()))
          continue;
        stk::mesh::FieldBase *field = meta.get_field(stk::topology::NODE_RANK, *I);
        if (field)
          stk::io::field_data_from_ioss(bulk, field, nodes, entity, *I);
      }
    }
  }
}

// ========================================================================

// ========================================================================
void process_sidesets_df(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity_df(entity, bulk);
    }
  }
}

void put_field_data(const stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    const std::vector<stk::io::FieldAndName> &namedFields,
                    Ioss::Field::RoleType filter_role,
                    const stk::mesh::Selector *subset_selector)
{
  std::vector<stk::mesh::Entity> entities;
  if(io_entity->type() == Ioss::SIDEBLOCK)
  {
    // Temporary Kluge to handle sideblocks which contain internally generated sides
    // where the "ids" field on the io_entity doesn't work to get the correct side...
    // NOTE: Could use this method for all entity types, but then need to correctly
    // specify whether shared entities are included/excluded (See IossBridge version).

    Ioss::SideBlock *block = dynamic_cast<Ioss::SideBlock *>(io_entity);
    if(block==nullptr)
    {
        std::ostringstream msg;
        msg << " INTERNAL ERROR: grouping entity failed to a SideBlock. Contact sierra-help@sandia.gov for assistance.\n";
        throw std::runtime_error(msg.str());
    }

    std::vector<int> elem_side_ids;
    stk::topology stk_elem_topology = map_ioss_topology_to_stk(block->parent_element_topology(), bulk.mesh_meta_data().spatial_dimension());

    stk::io::fill_element_and_side_ids(*io_entity, &part, bulk,
                                   stk_elem_topology,
                                   subset_selector,
                                   entities,
                                   elem_side_ids);

    static int counter = 0;
    counter++;
    size_t num_sides = entities.size();
    if(num_sides != static_cast<size_t>(io_entity->get_property("entity_count").get_int()))
    {
        std::ostringstream msg;
        msg << " INTERNAL_ERROR: Number of sides on part " << part.name() << " (" << num_sides
                << ") does not match number of sides in the associated Ioss SideBlock named "
                << io_entity->name() << " (" << io_entity->get_property("entity_count").get_int()
                << ").";
        throw std::runtime_error(msg.str());
    }
  }
  else
  {
    stk::io::get_entity_list(io_entity, part_type, bulk, entities);
  }

  for (size_t i=0;i<namedFields.size();i++)
  {
    const stk::mesh::FieldBase *f = namedFields[i].field();
    std::string field_name = namedFields[i].db_name();
    // There is ugliness here to deal with the output of fields on the nodes of a part,
    // either on the nodeblock or on a nodeset part created from the part nodes.
    if ((namedFields[i].m_forceNodeblockOutput && io_entity->type() == Ioss::NODEBLOCK) ||
        stk::io::is_field_on_part(f, part_type, part) ||
        (io_entity->type() != Ioss::NODEBLOCK && io_entity->field_exists(field_name))) {
      // NOTE: Multi-state issues are currently handled at field_add time
      stk::io::field_data_to_ioss(bulk, f, entities, io_entity, field_name, filter_role);
    }
  }
}



// ========================================================================
void put_field_data(const stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    const std::vector<stk::io::FieldAndName> &namedFields,
                    const stk::mesh::Selector *subset_selector=NULL)
{
  put_field_data(bulk, part, part_type, io_entity, namedFields, Ioss::Field::Field::TRANSIENT, subset_selector);
}

    StkMeshIoBroker::StkMeshIoBroker()
      : m_communicator(MPI_COMM_NULL), m_connectivity_map(NULL), m_active_mesh_index(0), m_sideset_face_creation_behavior(STK_IO_SIDE_CREATION_USING_GRAPH_TEST)
    {
      Ioss::Init::Initializer::initialize_ioss();
    }

    StkMeshIoBroker::StkMeshIoBroker(stk::ParallelMachine comm, const stk::mesh::ConnectivityMap * connectivity_map)
      : m_communicator(comm), m_connectivity_map(connectivity_map), m_active_mesh_index(0), m_sideset_face_creation_behavior(STK_IO_SIDE_CREATION_USING_GRAPH_TEST)
    {
      Ioss::Init::Initializer::initialize_ioss();
    }

    StkMeshIoBroker::~StkMeshIoBroker()
    { }

    void StkMeshIoBroker::property_add(const Ioss::Property &property)
    {
      m_property_manager.add(property);
      //In case there are already input/output files, put the property on them too.
      if (get_input_io_region().get() != NULL)
      {
          get_input_io_region()->property_add(property);
      }
      for(size_t i=0; i<m_output_files.size(); ++i)
      {
          get_output_io_region(i)->property_add(property);
      }
    }

    void StkMeshIoBroker::remove_property_if_exists(const std::string &property_name)
    {
      m_property_manager.erase(property_name);
    }

    stk::mesh::FieldBase const& StkMeshIoBroker::get_coordinate_field()
    {
      stk::mesh::FieldBase const* coord_field = meta_data().coordinate_field();
      STKIORequire( coord_field != NULL);
      return * coord_field;
    }

    size_t StkMeshIoBroker::add_mesh_database(Teuchos::RCP<Ioss::Region> ioss_input_region)
    {
      auto input_file = Teuchos::rcp(new InputFile(ioss_input_region));
      m_input_files.push_back(input_file);

      size_t index_of_input_file = m_input_files.size()-1;
      return index_of_input_file;
    }

    void StkMeshIoBroker::set_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
    {
      ThrowErrorMsgIf( !Teuchos::is_null(m_bulk_data),
                       "Bulk data already initialized" );
      m_bulk_data = arg_bulk_data;

      if (Teuchos::is_null(m_meta_data)) {
        m_meta_data = Teuchos::rcpFromRef(bulk_data().mesh_meta_data());
      }

#ifdef STK_BUILT_IN_SIERRA
      m_communicator = m_bulk_data->parallel();
#endif
    }

    void StkMeshIoBroker::replace_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
    {
      ThrowErrorMsgIf( Teuchos::is_null(m_bulk_data),
                       "There is  no bulk data to replace." );
      ThrowErrorMsgIf( Teuchos::is_null(m_meta_data),
                       "Meta data must be non-null when calling StkMeshIoBroker::replace_bulk_data." );

      stk::mesh::MetaData &new_meta_data = arg_bulk_data->mesh_meta_data();
      ThrowErrorMsgIf( &(*m_meta_data) != &new_meta_data, 
                       "Meta data for both new and old bulk data must be the same." );

      m_bulk_data = arg_bulk_data;
    }

    size_t StkMeshIoBroker::add_mesh_database(std::string filename, DatabasePurpose purpose)
    {
      std::string type = "exodus";

      // See if filename contains a ":" at the beginning of the filename
      // and if the text preceding that filename specifies a valid
      // database type.  If so, set that as the file type and extract
      // the portion following the colon as the filename.
      // If no colon in name, use default type.

      size_t colon = filename.find(':');
      if (colon != std::string::npos && colon > 0) {
        type = filename.substr(0, colon);
        filename = filename.substr(colon+1);
      }
      return add_mesh_database(filename, type, purpose);
    }


    size_t StkMeshIoBroker::add_mesh_database(const std::string &filename,
                                              const std::string &type,
                                              DatabasePurpose purpose)
    {
      auto input_file = Teuchos::rcp(new InputFile(filename, m_communicator, type, purpose, m_property_manager));
      m_input_files.push_back(input_file);

      size_t index_of_input_file = m_input_files.size()-1;
      return index_of_input_file;
    }

    Teuchos::RCP<Ioss::Region> StkMeshIoBroker::get_input_io_region()
    {
      if (is_index_valid(m_input_files, m_active_mesh_index)) {
        return m_input_files[m_active_mesh_index]->get_input_io_region();
      }
      else {
        return Teuchos::RCP<Ioss::Region>();
      }
    }

    InputFile &StkMeshIoBroker::get_mesh_database(size_t input_file_index)
    {
      validate_input_file_index(input_file_index);
      return *m_input_files[input_file_index];
    }

    void StkMeshIoBroker::remove_mesh_database(size_t input_file_index)
    {
      validate_input_file_index(input_file_index);
      // It would be nice to be able to just delete the database, but
      // its io_region may be being used by one or more output files
      // in the Ioss::Region::synchronize... function.  Therefore, we
      // need to keep it around (Could also use some shared pointers,
      // but not allowed to use C++11 yet).  But, we want this database to be
      // inaccessible and close the file associated with the database.  Therefore,
      // we add an empty InputFile to the end of 'm_input_files' and then swap it with
      // this one--That way the 'input_file_index' points to an invalid InputFile class as
      // it should.
      m_input_files[input_file_index]->get_input_io_region()->get_database()->closeDatabase();
      m_input_files.push_back(m_input_files[input_file_index]);
      m_input_files[input_file_index] = Teuchos::RCP<InputFile>();
      assert(Teuchos::is_null(m_input_files[input_file_index]));
  }

    size_t StkMeshIoBroker::set_active_mesh(size_t input_file_index)
    {
      validate_input_file_index(input_file_index);
      size_t old = m_active_mesh_index;
      m_active_mesh_index = input_file_index;
      m_input_files[m_active_mesh_index]->create_ioss_region();
      return old;
    }

      void StkMeshIoBroker::create_ioss_region()
      {
        validate_input_file_index(m_active_mesh_index);
        m_input_files[m_active_mesh_index]->create_ioss_region();
      }

      void StkMeshIoBroker::set_rank_name_vector(const std::vector<std::string> &rank_names)
      {
        ThrowErrorMsgIf(!Teuchos::is_null(m_meta_data),
                        "There meta data associated with this StkMeshIoBroker has already been created. "
                        "It is not permissible to set the rank_name_vector() at this time.");

        m_rank_names.clear();
        std::copy(rank_names.begin(), rank_names.end(), std::back_inserter(m_rank_names));
      }

      std::vector<std::string> get_ordered_attribute_field_names(const Ioss::GroupingEntity &iossElementBlock)
      {
          std::vector<std::string> names;
          iossElementBlock.field_describe(Ioss::Field::ATTRIBUTE, &names);
          std::vector<std::pair<int, std::string>> indexAndName;
          for(size_t i = 0; i < names.size(); i++)
          {
              int attributeIndex = iossElementBlock.get_field(names[i]).get_index();
              if(!is_skipped_attribute_field(names[i], names.size()))
                  indexAndName.push_back(std::make_pair(attributeIndex, names[i]));
          }
          std::sort(indexAndName.begin(), indexAndName.end());

          names.resize(indexAndName.size());
          for(size_t i=0;i<names.size();++i)
              names[i] = indexAndName[i].second;
          return names;
      }

      stk::mesh::FieldVector get_fields_by_name(const stk::mesh::MetaData &meta, const std::vector<std::string> &names)
      {
          stk::mesh::FieldVector attrFields(names.size());
          for(size_t i=0;i<attrFields.size();++i)
          {
              attrFields[i] = meta.get_field(stk::topology::ELEM_RANK, names[i]);
              ThrowRequireMsg(attrFields[i] != nullptr, "Can't find field named " << names[i]);
          }
          return attrFields;
      }

      void StkMeshIoBroker::store_attribute_field_ordering()
      {
          const stk::mesh::PartVector& parts = meta_data().get_parts();
          attributeFieldOrderingByPartOrdinal.clear();
          attributeFieldOrderingByPartOrdinal.resize(parts.size());
          for(const stk::mesh::Part* part : parts)
          {
              const Ioss::GroupingEntity* iossGroupingEntity = part->attribute<Ioss::GroupingEntity>();
              if(iossGroupingEntity != nullptr)
              {
                  std::vector<std::string> names = get_ordered_attribute_field_names(*iossGroupingEntity);
                  const stk::mesh::FieldVector attrFields = get_fields_by_name(*m_meta_data, names);
                  int partOrd = part->mesh_meta_data_ordinal();
                  attributeFieldOrderingByPartOrdinal[partOrd].resize(attrFields.size());
                  for(size_t i = 0; i < attrFields.size(); i++)
                  {
                      attributeFieldOrderingByPartOrdinal[partOrd][i] = attrFields[i]->mesh_meta_data_ordinal();
                  }
              }
          }
      }

      void StkMeshIoBroker::create_input_mesh()
      {
        validate_input_file_index(m_active_mesh_index);
        if (Teuchos::is_null(m_input_files[m_active_mesh_index]->get_input_io_region())) {
          m_input_files[m_active_mesh_index]->create_ioss_region();
        }

        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        ThrowErrorMsgIf (region==NULL,
                         "INTERNAL ERROR: Mesh Input Region pointer is NULL in create_input_mesh.");

        // See if meta data is null, if so, create a new one...
        if (Teuchos::is_null(m_meta_data)) {
          m_meta_data = Teuchos::rcp(new stk::mesh::MetaData());
        }

        size_t spatial_dimension = region->get_property("spatial_dimension").get_int();
        if (m_rank_names.empty()) {
          initialize_spatial_dimension(meta_data(), spatial_dimension, stk::mesh::entity_rank_names());
        } else {
          initialize_spatial_dimension(meta_data(), spatial_dimension, m_rank_names);
        }

        process_nodeblocks(*region,    meta_data());
        process_elementblocks(*region, meta_data());
        process_sidesets(*region,      meta_data());
        process_nodesets(*region,      meta_data());

        IossBlockMembership blockMemberships = get_block_memberships(*this);
        for(IossBlockMembership::iterator iter = blockMemberships.begin(); iter != blockMemberships.end(); iter++)
        {
          stk::mesh::Part* sidesetPart = meta_data().get_part(iter->first);
          if(sidesetPart != nullptr && sidesetPart->primary_entity_rank() == meta_data().side_rank())
          {
            std::vector<const stk::mesh::Part*> blocks;
            fill_block_parts_given_names(iter->second, meta_data(), blocks);
            meta_data().set_surface_to_block_mapping(sidesetPart, blocks);
          }
        }

        store_attribute_field_ordering();
      }


      size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type, 
                                                 char const* type)
      {
        return create_output_mesh(filename, db_type, m_property_manager, type);
      }


      size_t StkMeshIoBroker::create_output_mesh(const std::string &filename, DatabasePurpose db_type,
                                                 Ioss::PropertyManager &properties,
                                                 char const* type)
      {
        if(db_type == stk::io::APPEND_RESULTS)
        {
            properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND));
            db_type = stk::io::WRITE_RESULTS;
        }

        std::string out_filename = filename;
        stk::util::filename_substitution(out_filename);
        Ioss::Region *input_region = NULL;
        if (is_index_valid(m_input_files, m_active_mesh_index)) {
          input_region = get_input_io_region().get();
        }

        // Determine whether 64-bit integers are required for the output mesh...
        if (!properties.exists("INTEGER_SIZE_DB")){
          bool requires_64bit = check_integer_size_requirements() == 8;
          if (requires_64bit) {
            properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
          }
        }

        auto output_file = Teuchos::rcp(new impl::OutputFile(out_filename, m_communicator, db_type,
                                                                           properties, input_region, type));
        m_output_files.push_back(output_file);

        size_t index_of_output_file = m_output_files.size()-1;
        return index_of_output_file;
      }


      void StkMeshIoBroker::write_output_mesh(size_t output_file_index)
      {
        validate_output_file_index(output_file_index);
        stk::io::create_bulkdata_sidesets(bulk_data());
        m_output_files[output_file_index]->write_output_mesh(*m_bulk_data, attributeFieldOrderingByPartOrdinal);
        stk::io::clear_bulkdata_sidesets(bulk_data());
      }

      void StkMeshIoBroker::flush_output() const
      {
	for (const auto& out_file : m_output_files) {
          out_file->flush_output();
	}
	for (const auto& hb : m_heartbeat) {
          hb->flush_output();
	}
      }

      int StkMeshIoBroker::write_defined_output_fields(size_t output_file_index)
      {
        validate_output_file_index(output_file_index);
        int current_output_step = m_output_files[output_file_index]->write_defined_output_fields(*m_bulk_data);
        return current_output_step;
      }

      int StkMeshIoBroker::process_output_request(size_t output_file_index, double time)
      {
        validate_output_file_index(output_file_index);
        stk::io::create_bulkdata_sidesets(bulk_data());
        int current_output_step = m_output_files[output_file_index]->process_output_request(time, *m_bulk_data, attributeFieldOrderingByPartOrdinal);
        stk::io::clear_bulkdata_sidesets(bulk_data());
        return current_output_step;
      }

      void StkMeshIoBroker::begin_output_step(size_t output_file_index, double time)
      {
        validate_output_file_index(output_file_index);
        stk::io::create_bulkdata_sidesets(bulk_data());
        m_output_files[output_file_index]->begin_output_step(time, *m_bulk_data, attributeFieldOrderingByPartOrdinal);
      }

      void StkMeshIoBroker::end_output_step(size_t output_file_index)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->end_output_step();
        stk::io::clear_bulkdata_sidesets(bulk_data());
      }

      void StkMeshIoBroker::populate_mesh(bool delay_field_data_allocation)
      {
        validate_input_file_index(m_active_mesh_index);

        create_bulk_data();

        if (delay_field_data_allocation) {
          bulk_data().deactivate_field_updating();
        }

        bool i_started_modification_cycle = bulk_data().modification_begin("Mesh Read");

        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        bool ints64bit = db_api_int_size(region) == 8;
        if (ints64bit) {
#ifdef STK_BUILT_IN_SIERRA
          process_nodeblocks<int64_t>(*region,    bulk_data());
#else
          process_nodeblocks<int64_t>(*region,    bulk_data(), m_communicator);
#endif
          process_elementblocks<int64_t>(*region, bulk_data());
          process_nodesets<int64_t>(*region,      bulk_data());
        } else {
#ifdef STK_BUILT_IN_SIERRA
          process_nodeblocks<int>(*region,    bulk_data());
#else
          process_nodeblocks<int>(*region,    bulk_data(), m_communicator);
#endif
          process_elementblocks<int>(*region, bulk_data());
          process_nodesets<int>(*region,      bulk_data());
        }

        stk_mesh_resolve_node_sharing();

        stk::mesh::EntityIdProcMap elemIdMovedToProc;

        if(m_sideset_face_creation_behavior!=STK_IO_SIDE_CREATION_USING_GRAPH_TEST)
        {
            process_sidesets(*region, bulk_data(), elemIdMovedToProc, m_sideset_face_creation_behavior);
            bool saveOption = bulk_data().use_entity_ids_for_resolving_sharing();
            bulk_data().set_use_entity_ids_for_resolving_sharing(true);
            stk_mesh_modification_end_after_node_sharing_resolution();
            bulk_data().set_use_entity_ids_for_resolving_sharing(saveOption);
        }
        else
        {
            bulk_data().initialize_face_adjacent_element_graph();
            process_sidesets(*region, bulk_data(), elemIdMovedToProc, m_sideset_face_creation_behavior);
            stk_mesh_modification_end_after_node_sharing_resolution();
        }

        // Not sure if this is needed anymore. Don't think it'll be called with a nested modification cycle
        if(!i_started_modification_cycle)
            bulk_data().modification_begin();
      }

      void StkMeshIoBroker::populate_field_data()
      {
        validate_input_file_index(m_active_mesh_index);

        //if field-data has already been allocated, then the allocate_field_data() method
        //is a harmless no-op.
        bulk_data().allocate_field_data();

        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        ThrowErrorMsgIf (region==NULL,
                         "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_field_data.");

        bool ints64bit = db_api_int_size(region) == 8;
        if (ints64bit) {
          process_node_coords_and_attributes<int64_t>(*region, bulk_data());
          process_elem_attributes_and_implicit_ids<int64_t>(*region, bulk_data());
          process_nodesets_df<int64_t>(*region,      bulk_data());
          process_sidesets_df(*region,      bulk_data());
        }
        else {
          process_node_coords_and_attributes<int>(*region, bulk_data());
          process_elem_attributes_and_implicit_ids<int>(*region, bulk_data());
          process_nodesets_df<int>(*region,      bulk_data());
          process_sidesets_df(*region,      bulk_data());
        }

        if (region->get_property("state_count").get_int() == 0) {
          region->get_database()->release_memory();
        }
      }

      void StkMeshIoBroker::create_bulk_data()
      {
        if (!meta_data().is_commit())
          meta_data().commit();

        validate_input_file_index(m_active_mesh_index);
        ThrowErrorMsgIf (Teuchos::is_null(m_input_files[m_active_mesh_index]->get_input_io_region()),
                         "There is no Input mesh region associated with this Mesh Data.");

        Ioss::Region *region = m_input_files[m_active_mesh_index]->get_input_io_region().get();
        ThrowErrorMsgIf (region==NULL,
                         "INTERNAL ERROR: Mesh Input Region pointer is NULL in populate_mesh.");

        // Check if bulk_data is null; if so, create a new one...
        if (Teuchos::is_null(m_bulk_data)) {
          set_bulk_data(Teuchos::rcp( new stk::mesh::BulkData(   meta_data()
                                                                 , region->get_database()->util().communicator()
                                                                 , stk::mesh::BulkData::AUTO_AURA
#ifdef SIERRA_MIGRATION
                                                                 , false
#endif
                                                                 , m_connectivity_map
                                                                 )));
        }
      }

      // ========================================================================
      void StkMeshIoBroker::populate_bulk_data()
      {
        validate_input_file_index(m_active_mesh_index);

        create_bulk_data();

        // to preserve behavior for callers of this method, don't do the
        // delay-field-data-allocation optimization.
        // If want the optimization, call the population_mesh/populate_field_data methods separately.

        bool delay_field_data_allocation = false;
        populate_mesh(delay_field_data_allocation);

        populate_field_data();

        if(m_bulk_data->is_automatic_aura_on())
        {
            std::vector< const stk::mesh::FieldBase *> fields(m_meta_data->get_fields().begin(), m_meta_data->get_fields().end());
            stk::mesh::communicate_field_data(m_bulk_data->aura_ghosting(), fields);
        }
      }

      void StkMeshIoBroker::add_input_field(const stk::io::MeshField &mesh_field)
      {
        add_input_field(m_active_mesh_index, mesh_field);
      }

      void StkMeshIoBroker::add_input_field(size_t mesh_index, const stk::io::MeshField &mesh_field)
      {
        validate_input_file_index(mesh_index);
        m_input_files[mesh_index]->add_input_field(mesh_field);
      }

      void StkMeshIoBroker::validate_output_file_index(size_t output_file_index) const
      {
        ThrowErrorMsgIf(!is_index_valid(m_output_files, output_file_index),
                        "StkMeshIoBroker::validate_output_file_index: invalid output file index of "
                        << output_file_index << ".");

        ThrowErrorMsgIf (Teuchos::is_null(m_output_files[output_file_index]->get_output_io_region()),
                         "StkMeshIoBroker::validate_output_file_index: There is no Output mesh region associated with this output file index: " << output_file_index << ".");
      }

      void StkMeshIoBroker::validate_input_file_index(size_t input_file_index) const
      {
        ThrowErrorMsgIf(!is_index_valid(m_input_files, input_file_index),
                        "StkMeshIoBroker::validate_input_file_index: invalid input file index of "
                        << input_file_index << ".");
      }

      void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field)
      {
        add_field(output_file_index, field, field.name());
      }

      void StkMeshIoBroker::add_field(size_t output_file_index, stk::mesh::FieldBase &field, const std::string &alternate_name)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_field(field, alternate_name);
      }

      void StkMeshIoBroker::get_global_variable_names(std::vector<std::string> &names)
      {
        validate_input_file_index(m_active_mesh_index);
        m_input_files[m_active_mesh_index]->get_global_variable_names(names);
      }

      bool StkMeshIoBroker::get_global(const std::string &globalVarName,
                                       boost::any &value, stk::util::ParameterType::Type type,
                                       bool abort_if_not_found)
      {
        validate_input_file_index(m_active_mesh_index);
        auto region = m_input_files[m_active_mesh_index]->get_input_io_region();
        return internal_read_parameter(region, globalVarName, value, type, abort_if_not_found);
      }

      size_t StkMeshIoBroker::get_global_variable_length(const std::string& globalVarName)
      {
          validate_input_file_index(m_active_mesh_index);
          auto region = m_input_files[m_active_mesh_index]->get_input_io_region();

          size_t length = 0;
          if (region->field_exists(globalVarName)) {
              Ioss::Field field = region->get_field(globalVarName);
              length = field.raw_count() * field.raw_storage()->component_count();
          }
          return length;
      }

      bool StkMeshIoBroker::get_global(const std::string &globalVarName, std::vector<double> &globalVar,
                                       bool abort_if_not_found)
      {
        validate_input_file_index(m_active_mesh_index);
        auto region = m_input_files[m_active_mesh_index]->get_input_io_region();
        return internal_read_global(region, globalVarName, globalVar, Ioss::Field::REAL,
                                    abort_if_not_found);
      }

      bool StkMeshIoBroker::get_global(const std::string &globalVarName, std::vector<int> &globalVar,
                                       bool abort_if_not_found)
      {
        validate_input_file_index(m_active_mesh_index);
        auto region = m_input_files[m_active_mesh_index]->get_input_io_region();
        return internal_read_global(region, globalVarName, globalVar, Ioss::Field::INTEGER,
                                    abort_if_not_found);
      }

      bool StkMeshIoBroker::get_global(const std::string &globalVarName, int &globalVar,
                                       bool abort_if_not_found)
      {
        validate_input_file_index(m_active_mesh_index);
        auto region = m_input_files[m_active_mesh_index]->get_input_io_region();
        return internal_read_global(region, globalVarName, globalVar, Ioss::Field::INTEGER,
                                    abort_if_not_found);
      }

      bool StkMeshIoBroker::get_global(const std::string &globalVarName, double &globalVar,
                                       bool abort_if_not_found)
      {
        validate_input_file_index(m_active_mesh_index);
        auto region = m_input_files[m_active_mesh_index]->get_input_io_region();
        return internal_read_global(region, globalVarName, globalVar, Ioss::Field::REAL,
                                    abort_if_not_found);
      }

      bool StkMeshIoBroker::has_global(size_t output_file_index, const std::string &globalVarName) const
      {
        validate_output_file_index(output_file_index);
        return m_output_files[output_file_index]->has_global(globalVarName);
      }

      void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &name,
                                       const boost::any &value, stk::util::ParameterType::Type type)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_global(name, value, type);
      }

      void StkMeshIoBroker::add_global_ref(size_t output_file_index, const std::string &name,
                                           const boost::any *value, stk::util::ParameterType::Type type)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_global_ref(name, value, type);
      }

      void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, Ioss::Field::BasicType dataType)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_global(globalVarName, dataType);
      }

      void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, int component_count, Ioss::Field::BasicType dataType)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_global(globalVarName, component_count, dataType);
      }

      void StkMeshIoBroker::add_global(size_t output_file_index, const std::string &globalVarName, const std::string &storage, Ioss::Field::BasicType dataType)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->add_global(globalVarName, storage, dataType);
      }

      void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName,
                                         const boost::any &value, stk::util::ParameterType::Type type)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->write_global(globalVarName, value, type);
      }

      void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, double globalVarData)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->write_global(globalVarName, globalVarData);
      }

      void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, int globalVarData)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->write_global(globalVarName, globalVarData);
      }

      void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, std::vector<double>& globalVarData)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->write_global(globalVarName, globalVarData);
      }

      void StkMeshIoBroker::write_global(size_t output_file_index, const std::string &globalVarName, std::vector<int>& globalVarData)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->write_global(globalVarName, globalVarData);
      }

      FieldNameToPartVector StkMeshIoBroker::get_nodal_var_names()
      {
          validate_input_file_index(m_active_mesh_index);
          return m_input_files[m_active_mesh_index]->get_var_names(Ioss::NODEBLOCK, meta_data());
      }

      FieldNameToPartVector StkMeshIoBroker::get_elem_var_names()
      {
          validate_input_file_index(m_active_mesh_index);
          return m_input_files[m_active_mesh_index]->get_var_names(Ioss::ELEMENTBLOCK, meta_data());
      }

      FieldNameToPartVector StkMeshIoBroker::get_nodeset_var_names()
      {
          validate_input_file_index(m_active_mesh_index);
          return m_input_files[m_active_mesh_index]->get_var_names(Ioss::NODESET, meta_data());
      }

      FieldNameToPartVector StkMeshIoBroker::get_sideset_var_names()
      {
          validate_input_file_index(m_active_mesh_index);
          return m_input_files[m_active_mesh_index]->get_var_names(Ioss::SIDESET, meta_data());
      }

      void StkMeshIoBroker::add_all_mesh_fields_as_input_fields(MeshField::TimeMatchOption tmo)
      {
        validate_input_file_index(m_active_mesh_index);
        m_input_files[m_active_mesh_index]->add_all_mesh_fields_as_input_fields(meta_data(), tmo);
      }

      bool StkMeshIoBroker::read_input_field(stk::io::MeshField &mf)
      {
        validate_input_file_index(m_active_mesh_index);
        return m_input_files[m_active_mesh_index]->read_input_field(mf, bulk_data());
      }

      double StkMeshIoBroker::read_defined_input_fields(double time,
                                                        std::vector<stk::io::MeshField> *missingFields)
      {
        validate_input_file_index(m_active_mesh_index);
        return m_input_files[m_active_mesh_index]->read_defined_input_fields(time, missingFields, bulk_data());
      }

      double StkMeshIoBroker::read_defined_input_fields(int step,
                                                        std::vector<stk::io::MeshField> *missing)
      {
        if (step <= 0)
          return 0.0;

        validate_input_file_index(m_active_mesh_index);
        return m_input_files[m_active_mesh_index]->read_defined_input_fields(step, missing, bulk_data());
      }

      double StkMeshIoBroker::read_defined_input_fields_at_step(int step,
                                                        std::vector<stk::io::MeshField> *missing)
      {
        if (step <= 0)
          return 0.0;

        validate_input_file_index(m_active_mesh_index);
        return m_input_files[m_active_mesh_index]->read_defined_input_fields_at_step(step, missing, bulk_data());
      }

      bool StkMeshIoBroker::use_nodeset_for_part_nodes_fields(size_t output_file_index) const
      {
        validate_output_file_index(output_file_index);
        return m_output_files[output_file_index]->use_nodeset_for_part_nodes_fields();
      }

      void StkMeshIoBroker::use_nodeset_for_part_nodes_fields(size_t output_file_index, bool true_false)
      {
        validate_output_file_index(output_file_index);
        m_output_files[output_file_index]->use_nodeset_for_part_nodes_fields(true_false);
      }

      void StkMeshIoBroker::set_option_to_not_collapse_sequenced_fields()
      {
          char tmp[2] = {'\1', '\0'};
          property_add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", tmp));
      }

      int StkMeshIoBroker::get_num_time_steps()
      {
          int numTimeSteps = 0;
          Ioss::Region *ioRegion = get_input_io_region().get();
          if(ioRegion != nullptr)
          {
              Ioss::Property stateCount = ioRegion->get_implicit_property("state_count");
              numTimeSteps = stateCount.get_int();
          }
          return numTimeSteps;
      }

      std::vector<double> StkMeshIoBroker::get_time_steps()
      {
          int numTimeSteps = get_num_time_steps();
          std::vector<double> timeSteps;

          Ioss::Region *ioRegion = get_input_io_region().get();
          if(ioRegion != nullptr)
          {
              for(int istep = 0; istep < numTimeSteps; istep++)
              {
                  double state_time = ioRegion->get_state_time(istep + 1);
                  timeSteps.push_back(state_time);
              }
          }
          return timeSteps;
      }

      double StkMeshIoBroker::get_max_time()
      {
          return get_input_io_region()->get_max_time().second;
      }

      void StkMeshIoBroker::set_max_num_steps_before_overwrite(size_t outputFileIndex, int maxNumStepsInFile)
      {
          get_output_io_region(outputFileIndex)->get_database()->set_cycle_count(maxNumStepsInFile);
      }

      size_t StkMeshIoBroker::add_heartbeat_output(const std::string &filename, HeartbeatType hb_type,
                                                   const Ioss::PropertyManager &properties)
      {
        std::string out_filename = filename;
        stk::util::filename_substitution(out_filename);
        auto heartbeat = Teuchos::rcp(new impl::Heartbeat(out_filename, hb_type,
                                                                       properties, m_communicator));
        m_heartbeat.push_back(heartbeat);
        return m_heartbeat.size()-1;
      }

      int StkMeshIoBroker::check_integer_size_requirements()
      {
	// 1. If the INGEGER_SIZE_DB or _API property exists, then use its value no matter what...
	if (m_property_manager.exists("INTEGER_SIZE_DB")) {
	  return m_property_manager.get("INTEGER_SIZE_DB").get_int();
	}

	if (m_property_manager.exists("INTEGER_SIZE_API")) {
	  return m_property_manager.get("INTEGER_SIZE_API").get_int();
	}
	
	// 2. If input_region exists, then if it is using 64-bit integers, the output should
	//    use those also.
	Ioss::Region *input_region = NULL;
        if (is_index_valid(m_input_files, m_active_mesh_index)) {
          input_region = get_input_io_region().get();
        }
	if (input_region != NULL) {
	  // Get the integer size setting for the database associated with the region.
	  int int_size = db_api_int_size(input_region);
	  if (int_size == 8) {
	    return int_size;
	  }
	}

	// 3. If any entity count exceeds INT_MAX, then use 64-bit integers.
        if ( !Teuchos::is_null(m_bulk_data) ) {
          std::vector<size_t> entityCounts;
          stk::mesh::comm_mesh_counts(*m_bulk_data, entityCounts);
          for (size_t i=0; i < entityCounts.size(); i++) {
            if (entityCounts[i] > (size_t)std::numeric_limits<int>::max()) {
              return 8;
            }
          }
        }
	
        // 4. Should also check if the maximum node or element id exceeds INT_MAX.

	// 5. Default to 4-byte integers...
	return 4;
      }

      void StkMeshIoBroker::set_name_and_version_for_qa_record(size_t outputFileIndex, const std::string &codeName, const std::string &codeVersion)
      {
          Ioss::Region *region = get_output_io_region(outputFileIndex).get();
          region->property_add(Ioss::Property(std::string("code_name"), codeName));
          region->property_add(Ioss::Property(std::string("code_version"), codeVersion));
      }
      void StkMeshIoBroker::add_qa_records(size_t outputFileIndex, const std::vector<QaRecord> &qaRecords)
      {
          Ioss::Region *region = get_output_io_region(outputFileIndex).get();
          for(const QaRecord &qaRec : qaRecords)
            region->add_qa_record(qaRec.name, qaRec.version, qaRec.date, qaRec.time);
      }
      void StkMeshIoBroker::add_info_records(size_t outputFileIndex, const std::vector<std::string> &infoRecords)
      {
          Ioss::Region *region = get_output_io_region(outputFileIndex).get();
          region->add_information_records(infoRecords);
      }
      std::vector<QaRecord> StkMeshIoBroker::get_qa_records()
      {
          std::vector<QaRecord> qaRecords;
          Ioss::Region *region = get_input_io_region().get();
          const std::vector<std::string> &qa = region->get_qa_records();
          for (size_t i = 0; i < qa.size(); i += 4)
              qaRecords.push_back({qa[i + 0], qa[i + 1], qa[i + 2], qa[i + 3]});
          return qaRecords;
      }
      std::vector<std::string> StkMeshIoBroker::get_info_records()
      {
          Ioss::Region *region = get_input_io_region().get();
          return region->get_information_records();
      }

      stk::mesh::FieldVector StkMeshIoBroker::get_ordered_attribute_fields(const stk::mesh::Part *blockPart) const
      {
          stk::mesh::FieldVector attrFields;
          if(blockPart->mesh_meta_data_ordinal() < attributeFieldOrderingByPartOrdinal.size())
          {
              const std::vector<int> &fieldOrds = attributeFieldOrderingByPartOrdinal[blockPart->mesh_meta_data_ordinal()];
              attrFields.resize(fieldOrds.size());
              const stk::mesh::FieldVector &allFields = m_meta_data->get_fields();
              for(size_t i=0; i<fieldOrds.size(); i++)
              {
                  attrFields[i] = allFields[fieldOrds[i]];
              }
          }
          return attrFields;
      }

      const std::vector<std::vector<int>> & StkMeshIoBroker::get_attribute_field_ordering_stored_by_part_ordinal() const
      {
          return attributeFieldOrderingByPartOrdinal;
      }

      void StkMeshIoBroker::set_attribute_field_ordering_stored_by_part_ordinal(const std::vector<std::vector<int>> &ordering)
      {
          attributeFieldOrderingByPartOrdinal = ordering;
      }

      impl::Heartbeat::Heartbeat(const std::string &filename, HeartbeatType hb_type,
                           Ioss::PropertyManager properties, stk::ParallelMachine comm)
        : m_current_step(0), m_processor(0)
        {
          if (comm != MPI_COMM_NULL) {
            m_processor = stk::parallel_machine_rank(comm);
          }

          if (m_processor == 0) {
            std::string db_io_type = "exodus";
            Ioss::DatabaseUsage db_usage = Ioss::WRITE_HISTORY;

            if (hb_type != BINARY) {
              db_io_type = "heartbeat";
              db_usage = Ioss::WRITE_HEARTBEAT;

              // Always add the "time" field to all heartbeat outputs...
              if (!properties.exists("SHOW_TIME_FIELD")) {
                properties.add(Ioss::Property("SHOW_TIME_FIELD", true));
              }

              if (hb_type == SPYHIS) {
                if (!properties.exists("FILE_FORMAT")) {
                  properties.add(Ioss::Property("FILE_FORMAT", "spyhis"));
                }
              }
              else if (hb_type == CSV) {
                if (!properties.exists("SHOW_TIME_STAMP")) {
                  properties.add(Ioss::Property("SHOW_TIME_STAMP", false));
                }
                if (!properties.exists("FIELD_SEPARATOR")) {
                  properties.add(Ioss::Property("FIELD_SEPARATOR", ", "));
                }
              }
              else if (hb_type == TS_CSV) {
                if (!properties.exists("SHOW_TIME_STAMP")) {
                  properties.add(Ioss::Property("SHOW_TIME_STAMP", true));
                }
                if (!properties.exists("FIELD_SEPARATOR")) {
                  properties.add(Ioss::Property("FIELD_SEPARATOR", ", "));
                }
              }
              else if (hb_type == TEXT) {
                if (!properties.exists("SHOW_TIME_STAMP")) {
                  properties.add(Ioss::Property("SHOW_TIME_STAMP", false));
                }
                if (!properties.exists("FIELD_SEPARATOR")) {
                  properties.add(Ioss::Property("FIELD_SEPARATOR", "\t"));
                }
              }
              else if (hb_type == TS_TEXT) {
                if (!properties.exists("SHOW_TIME_STAMP")) {
                  properties.add(Ioss::Property("SHOW_TIME_STAMP", true));
                }
                if (!properties.exists("FIELD_SEPARATOR")) {
                  properties.add(Ioss::Property("FIELD_SEPARATOR", "\t"));
                }
              }
            }

            Ioss::DatabaseIO *db = Ioss::IOFactory::create(db_io_type, filename,
                                                           db_usage, comm, properties);
            if (db == NULL || !db->ok()) {
              std::cerr << "ERROR: Could not open history/heartbeat database '" << filename << "'\n";
              return;
            }

            // NOTE: 'region' owns 'db' pointer at this time...
            m_region = Teuchos::rcp(new Ioss::Region(db, filename));
          }
        }

        void impl::Heartbeat::add_global_ref(const std::string &name, const boost::any *value,
                                       stk::util::ParameterType::Type type)
        {
          if (m_processor == 0) {
            ThrowErrorMsgIf (m_current_step != 0,
                             "At least one output step has been written to the history/heartbeat file. "
                             "Variables cannot be added anymore.");

            Ioss::State currentState = m_region->get_state();
            if(currentState != Ioss::STATE_DEFINE_TRANSIENT) {
              m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
            }

            // Determine name and type of parameter...
            std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, *value);
            internal_add_global(m_region, name, parameter_type.first, parameter_type.second);
            m_fields.push_back(GlobalAnyVariable(name, value, type));
          }
        }

        void impl::Heartbeat::process_output(int step, double time)
        {
          if (m_processor == 0) {
            Ioss::State currentState = m_region->get_state();
            if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
              m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
            }

            m_region->begin_mode(Ioss::STATE_TRANSIENT);
            m_current_step = m_region->add_state(time);
            m_region->begin_state(m_current_step);

            write_defined_global_any_fields(m_region, m_fields);

            m_region->end_state(m_current_step);
            m_region->end_mode(Ioss::STATE_TRANSIENT);
          }
        }

        void impl::Heartbeat::flush_output() const
        {
          if (m_processor == 0) {
	    m_region->get_database()->flush_database();
	  }
        }

        void impl::OutputFile::write_output_mesh(const stk::mesh::BulkData& bulk_data,
                                                 const std::vector<std::vector<int>> &attributeOrdering)
        {
          if ( m_mesh_defined == false )
            {
              m_mesh_defined = true;

              // If using hdf5 as the underlying file type for exodus/netcdf,
              // it is more picky about overwriting an existing file -- if the
              // file is open, then it will abort; it will only overwrite an existing
              // file if it is not open.  Since overwriting restart files (input/output)
              // is a common usecase, we need to check at this point whether there are
              // any existing input files with the same name as the file we are attempting
              // to create here. However, due to symbolic links and other junk, it is often
              // difficult to determine that the files are the same, so..., If m_input_region
              // refers to a file, just close it since we should be done with it at this time...
              if (m_input_region) {
                m_input_region->get_database()->closeDatabase();
              }

              // used in stk_adapt/stk_percept
              bool sort_stk_parts_by_name = m_region->property_exists("sort_stk_parts");

              stk::io::define_output_db(*m_region, bulk_data, attributeOrdering, m_input_region, m_subset_selector.get(),
                                        sort_stk_parts_by_name, m_use_nodeset_for_part_nodes_fields);

              if(!m_appending_to_mesh)
                  stk::io::write_output_db(*m_region, bulk_data, m_subset_selector.get());

              //Attempt to avoid putting state change into the interface.  We'll see . . .
              m_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
            }
        }

        void impl::OutputFile::flush_output() const
        {
	    m_region->get_database()->flush_database();
        }

        void impl::OutputFile::add_field(stk::mesh::FieldBase &field, const std::string &alternate_name)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "Attempting to add fields after fields have already been written to the database.");
          ThrowErrorMsgIf (alternate_name.empty(),
                           "Attempting to output results field " << field.name() << " with no name.");

          bool fieldAlreadyExists=false;
          for (size_t i=0;i<m_named_fields.size();i++) {
            if ( &field == m_named_fields[i].field() ) {
              m_named_fields[i].set_db_name(alternate_name);
              fieldAlreadyExists = true;
              break;
            }
          }

          if (!fieldAlreadyExists) {
            if (m_db_purpose == stk::io::WRITE_RESTART) {
              int state_count = field.number_of_states();
              STKIORequire(state_count < 7);
              int num_states_to_write = std::max(state_count-1, 1);
              for(int state=0; state < num_states_to_write; state++) {
                stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
                stk::mesh::FieldBase *statedField = field.field_state(state_identifier);
                std::string field_name_with_suffix = stk::io::get_stated_field_name(alternate_name, state_identifier);
                stk::io::FieldAndName namedField(statedField, field_name_with_suffix);
                m_named_fields.push_back(namedField);
                stk::io::set_field_role(*statedField, Ioss::Field::TRANSIENT);
              }
            } else {
              stk::io::FieldAndName namedField(&field, alternate_name);
              m_named_fields.push_back(namedField);
              stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
            }
          }
        }

        void impl::OutputFile::add_global_ref(const std::string &name, const boost::any *value, stk::util::ParameterType::Type type)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "On region named " << m_region->name() <<
                           " Attempting to add global variable after data has already been written to the database.");
          std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, *value);
          internal_add_global(m_region, name, parameter_type.first, parameter_type.second);
          m_global_any_fields.push_back(GlobalAnyVariable(name, value, type));
        }

        bool impl::OutputFile::has_global(const std::string &globalVarName) const
        {
            return m_region->field_exists(globalVarName);
        }

        void impl::OutputFile::add_global(const std::string &name, const boost::any &value, stk::util::ParameterType::Type type)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "On region named " << m_region->name() <<
                           " Attempting to add global variable after data has already been written to the database.");
          std::pair<size_t, Ioss::Field::BasicType> parameter_type = get_io_parameter_size_and_type(type, value);
          m_non_any_global_variables_defined = true;  // This output file has at least 1 global variable.
          internal_add_global(m_region, name, parameter_type.first, parameter_type.second);
        }

        void impl::OutputFile::add_global(const std::string &globalVarName, Ioss::Field::BasicType dataType)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "On region named " << m_region->name() <<
                           " Attempting to add global variable after data has already been written to the database.");
          m_non_any_global_variables_defined = true;  // This output file has at least 1 global variable.
          internal_add_global(m_region, globalVarName, "scalar", dataType);
        }

        void impl::OutputFile::add_global(const std::string &globalVarName, int component_count, Ioss::Field::BasicType dataType)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "On region named " << m_region->name() <<
                           " Attempting to add global variable after data has already been written to the database.");
          m_non_any_global_variables_defined = true;  // This output file has at least 1 global variable.
          internal_add_global(m_region, globalVarName, component_count, dataType);
        }

        void impl::OutputFile::add_global(const std::string &globalVarName, const std::string &storage, Ioss::Field::BasicType dataType)
        {
          ThrowErrorMsgIf (m_fields_defined,
                           "On region named " << m_region->name() <<
                           " Attempting to add global variable after data has already been written to the database.");
          m_non_any_global_variables_defined = true;  // This output file has at least 1 global variable.
          internal_add_global(m_region, globalVarName, storage, dataType);
        }

        void impl::OutputFile::write_global(const std::string &globalVarName,
                                      const boost::any &value, stk::util::ParameterType::Type type)
        {
          internal_write_parameter(m_region, globalVarName, value, type);
        }

        void impl::OutputFile::write_global(const std::string &globalVarName, std::vector<double>& globalVarData)
        {
          internal_write_global(m_region, globalVarName, globalVarData);
        }

        void impl::OutputFile::write_global(const std::string &globalVarName, std::vector<int>& globalVarData)
        {
          internal_write_global(m_region, globalVarName, globalVarData);
        }

        void impl::OutputFile::write_global(const std::string &globalVarName, int globalVarData)
        {
          internal_write_global(m_region, globalVarName, globalVarData);
        }

        void impl::OutputFile::write_global(const std::string &globalVarName, double globalVarData)
        {
          internal_write_global(m_region, globalVarName, globalVarData);
        }

        void impl::OutputFile::setup_output_file(const std::string &filename, stk::ParallelMachine communicator,
                                           Ioss::PropertyManager &property_manager, char const* type)
        {
          ThrowErrorMsgIf (filename.empty(),
                           "No filename was specified for the output file creation.");
          Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(type, filename,
                                                          Ioss::WRITE_RESULTS,
                                                          communicator,
                                                          property_manager);

          if (dbo == NULL || !dbo->ok()) {
            std::cerr << "ERROR: Could not open output database '" << filename << "' of type 'exodus'\n";
            std::exit(EXIT_FAILURE);
          }

          // There is no "label" for the output region; just use the filename for now so
          // Ioss messages will specify a unique or identifiable instance.
          m_region = Teuchos::rcp(new Ioss::Region(dbo, filename));

          if(property_manager.exists("APPEND_OUTPUT") && property_manager.get("APPEND_OUTPUT").get_int() == Ioss::DB_APPEND)
              m_appending_to_mesh = true;
        }

        void impl::OutputFile::begin_output_step(double time, const stk::mesh::BulkData& bulk_data,
                                                 const std::vector<std::vector<int>> &attributeOrdering)
        {
          if (!m_fields_defined) {
            define_output_fields(bulk_data, attributeOrdering);
          }

          //Attempt to avoid putting state change into the interface.  We'll see . . .
          Ioss::State currentState = m_region->get_state();
          if(currentState == Ioss::STATE_DEFINE_TRANSIENT) {
            m_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
          }

          m_region->begin_mode(Ioss::STATE_TRANSIENT);
          m_current_output_step = m_region->add_state(time);
          m_region->begin_state(m_current_output_step);
        }

        // ========================================================================
        // Iterate over all fields defined in the stk mesh data structure.
        // If the field has the io_attribute set, then define that field
        // on the corresponding io entity on the output mesh database.
        // The database field will have the same name as the stk field.
        //
        // To export the data to the database, call
        // process_output_request().
        void impl::OutputFile::define_output_fields(const stk::mesh::BulkData& bulk_data,
                                                    const std::vector<std::vector<int>> &attributeOrdering)
        {
          if(m_fields_defined) {
            return;
          }

          write_output_mesh(bulk_data, attributeOrdering);

          Ioss::Region *region = m_region.get();

          // Sort the fields in m_named_fields based on the field ordinal.
          // This is needed so all processors will process the fields in the same order
          // Not guaranteed that the application has added the fields in order, so make
          // the guarantee here...
          std::sort(m_named_fields.begin(), m_named_fields.end(), fieldOrdinalSort);

          const stk::mesh::MetaData &meta_data = bulk_data.mesh_meta_data();
          // Special processing for nodeblock (all nodes in model)...
          stk::io::ioss_add_fields(meta_data.universal_part(), stk::topology::NODE_RANK,
                                   region->get_node_blocks()[0], m_named_fields);

          const stk::mesh::PartVector &all_parts = meta_data.get_parts();
          for(stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip) {
            stk::mesh::Part * const part = *ip;

            // Check whether this part should be output to database.
            if(stk::io::is_part_io_part(*part)) {
              stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
              // Get Ioss::GroupingEntity corresponding to this part...
              Ioss::GroupingEntity *entity = region->get_entity(part->name());

              if(entity!= nullptr)
              {
                if (entity->type() == Ioss::SIDESET)
                {
                  Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
                  int block_count = sset->block_count();

                  for(int i = 0; i < block_count; i++)
                  {
                      Ioss::SideBlock *fb = sset->get_block(i);
                      stk::io::ioss_add_fields(*part, rank, fb, m_named_fields);
                  }
                }
                else
                {
                  stk::io::ioss_add_fields(*part, rank, entity, m_named_fields);
                }
              }
              // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
              // (should probably do edges and faces also...)
              // Get Ioss::GroupingEntity corresponding to the nodes on this part...
              if(rank != stk::topology::NODE_RANK) {
                Ioss::GroupingEntity *node_entity = NULL;
                if (m_use_nodeset_for_part_nodes_fields) {
                  std::string nodes_name = part->name() + "_nodes";
                  node_entity = region->get_entity(nodes_name);
                } else {
                  node_entity = region->get_entity("nodeblock_1");
                }
                if(node_entity != NULL) {
                  stk::io::ioss_add_fields(*part, stk::topology::NODE_RANK, node_entity, m_named_fields);
                }
              }
            }
          }
          m_fields_defined = true;
        }

        int impl::OutputFile::process_output_request(double time, const stk::mesh::BulkData& bulk_data,
                                                     const std::vector<std::vector<int>> &attributeOrdering)
        {
          ThrowErrorMsgIf(m_non_any_global_variables_defined,
                          "The output database " << m_region->name() << " has defined global variables, "
                          "but is calling the process_output_request() function which does not output global "
                          "variables.  Call begin_output_step() instead.");

          begin_output_step(time, bulk_data, attributeOrdering);
          write_defined_output_fields(bulk_data);
          write_defined_global_any_fields(m_region, m_global_any_fields);
          end_output_step();

          return m_current_output_step;
        }

        int impl::OutputFile::write_defined_output_fields(const stk::mesh::BulkData& bulk_data)
        {
          Ioss::Region *region = m_region.get();
          ThrowErrorMsgIf (region==NULL, "INTERNAL ERROR: Mesh Output Region pointer is NULL in write_defined_output_fields.");

          const stk::mesh::MetaData& meta_data = bulk_data.mesh_meta_data();
          // Special processing for nodeblock (all nodes in model)...
          put_field_data(bulk_data, meta_data.universal_part(), stk::topology::NODE_RANK,
                         region->get_node_blocks()[0], m_named_fields, m_subset_selector.get());

          // Now handle all non-nodeblock parts...
          const stk::mesh::PartVector &all_parts = meta_data.get_parts();
          for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {
            stk::mesh::Part * const part = *ip;

            // Check whether this part should be output to database.
            if (stk::io::is_part_io_part(*part)) {
              stk::mesh::EntityRank rank = part_primary_entity_rank(*part);
              // Get Ioss::GroupingEntity corresponding to this part...
              Ioss::GroupingEntity *entity = region->get_entity(part->name());

              // If the sideset has only a single sideblock and it
              // shares the same name as the parent sideset, then the
              // entity that we want to output on at this point is the
              // sideblock and not the sideset.  If there are multiple
              // sideblocks in the sideset, then they will be output
              // separately...
              if (entity != NULL && entity->type() == Ioss::SIDESET) {
                Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
                size_t block_count = sset->block_count();
                if (block_count == 1) {
                  Ioss::GroupingEntity *ssblock = sset->get_side_block(part->name());
                  if (ssblock) {
                    // NOTE: 'entity' is reset at this point.
                    entity = ssblock;
                  }
                }
              }

              if (entity != NULL && entity->type() != Ioss::SIDESET) {
                put_field_data(bulk_data, *part, rank, entity, m_named_fields, m_subset_selector.get());
              }

              // If rank is != NODE_RANK, then see if any fields are defined on the nodes of this part
              // (should probably do edges and faces also...)
              // Get Ioss::GroupingEntity corresponding to the nodes on this part...
              if (rank != stk::topology::NODE_RANK && m_use_nodeset_for_part_nodes_fields) {
                std::string nodes_name = part->name() + "_nodes";
                Ioss::GroupingEntity *node_entity = region->get_entity(nodes_name);
                if (node_entity != NULL) {
                  put_field_data(bulk_data, *part, stk::topology::NODE_RANK, node_entity, m_named_fields,
                                 m_subset_selector.get());
                }
              }
            }
          }
          return m_current_output_step;
        }

        void impl::OutputFile::end_output_step()
        {
          m_region->end_state(m_current_output_step);
          m_region->end_mode(Ioss::STATE_TRANSIENT);
        }

        void impl::OutputFile::set_subset_selector(Teuchos::RCP<stk::mesh::Selector> my_selector)
        {
          ThrowErrorMsgIf(m_mesh_defined,
                          "ERROR: On region named " << m_region->name() <<
                          " the subset_selector cannot be changed after the mesh has already been written.");
          m_subset_selector = my_selector;
        }

        bool impl::OutputFile::use_nodeset_for_part_nodes_fields() const
        {
          return m_use_nodeset_for_part_nodes_fields;
        }

        void impl::OutputFile::use_nodeset_for_part_nodes_fields(bool true_false)
        {
          ThrowErrorMsgIf(m_mesh_defined,
                          "ERROR: The use_nodeset_for_part_nodes_fields setting cannot be changed after "
                          "the mesh has already been written.");
          m_use_nodeset_for_part_nodes_fields = true_false;
        }

    } // namespace io
  } // namespace stk
