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
#include <stk_io/IossBridge.hpp>
#include <Ioss_IOFactory.h>                         // for NameList, IOFactory
#include <cassert>                                  // for assert
#include <algorithm>                                // for min, sort, max
#include <cstdint>                                  // for int64_t, uint64_t
#include <iostream>                                 // for operator<<, basic...
#include <memory>                                   // for allocator_traits<...
#include <stdexcept>                                // for runtime_error
#include <unordered_map>
#include <tuple>
#include <stk_mesh/base/BulkData.hpp>               // for BulkData
#include <stk_mesh/base/Comm.hpp>                   // for comm_mesh_counts
#include <stk_mesh/base/FEMHelpers.hpp>             // for get_side_entity_f...
#include <stk_mesh/base/Field.hpp>                  // for Field
#include <stk_mesh/base/FindRestriction.hpp>        // for find_restriction
#include <stk_mesh/base/GetEntities.hpp>            // for count_selected_en...
#include <stk_mesh/base/MetaData.hpp>               // for MetaData, put_fie...
#include <stk_mesh/base/Types.hpp>                  // for PartVector, Entit...
#include <stk_util/diag/StringUtil.hpp>             // for make_lower, to_st...
#include <stk_util/environment/RuntimeWarning.hpp>  // for RuntimeWarning
#include <stk_util/parallel/ParallelReduce.hpp>     // for all_reduce_sum
#include <stk_util/util/SortAndUnique.hpp>          // for sort_and_unique
#include <stk_util/util/string_case_compare.hpp>
#include <stk_util/util/tokenize.hpp>               // for tokenize
#include <typeinfo>                                 // for type_info
#include "Ioss_Assembly.h"                          // for Assembly
#include "Ioss_CommSet.h"                           // for CommSet
#include "Ioss_CompositeVariableType.h"             // for CompositeVariable...
#include "Ioss_ConcreteVariableType.h"
#include "Ioss_DBUsage.h"                           // for WRITE_RESULTS
#include "Ioss_DataSize.h"                          // for USE_INT64_API
#include "Ioss_DatabaseIO.h"                        // for DatabaseIO
#include "Ioss_EdgeBlock.h"                         // for EdgeBlock
#include "Ioss_ElementBlock.h"                      // for ElementBlock
#include "Ioss_ElementTopology.h"                   // for ElementTopology
#include "Ioss_EntityBlock.h"                       // for EntityBlock
#include "Ioss_EntityType.h"                        // for EntityType, NODEB...
#include "Ioss_FaceBlock.h"                         // for FaceBlock
#include "Ioss_Field.h"                             // for Field, Field::Rol...
#include "Ioss_GroupingEntity.h"                    // for GroupingEntity
#include "Ioss_NodeBlock.h"                         // for NodeBlock
#include "Ioss_NodeSet.h"                           // for NodeSet
#include "Ioss_Property.h"                          // for Property
#include "Ioss_Region.h"                            // for Region, SideSetCo...
#include "Ioss_SideBlock.h"                         // for SideBlock
#include "Ioss_SideSet.h"                           // for SideSet, SideBloc...
#include "Ioss_State.h"                             // for STATE_DEFINE_MODEL
#include "Ioss_SurfaceSplit.h"                      // for SPLIT_BY_ELEMENT_...
#include "Ioss_VariableType.h"                      // for VariableType
#include "Ioss_Utils.h"
#include "SidesetTranslator.hpp"                    // for get_number_sides_...
#include "StkIoUtils.hpp"                           // for part_primary_enti...
#include "mpi.h"                                    // for MPI_COMM_SELF
#include "stk_io/FieldAndName.hpp"                  // for FieldAndName
#include "stk_mesh/base/Bucket.hpp"                 // for Bucket
#include "stk_mesh/base/Entity.hpp"                 // for Entity
#include "stk_mesh/base/EntityKey.hpp"              // for operator<<
#include "stk_mesh/base/FieldBase.hpp"              // for FieldBase, FieldB...
#include "stk_mesh/base/FieldRestriction.hpp"       // for FieldRestriction
#include "stk_mesh/base/Part.hpp"                   // for Part, Part::INVAL...
#include "stk_mesh/base/Selector.hpp"               // for Selector, operator&
#include "stk_mesh/baseImpl/PartAttribute.hpp"
#include "stk_topology/topology.hpp"                // for topology, topolog...
#include "stk_util/util/ReportHandler.hpp"          // for ThrowRequireMsg
#include "stk_util/util/string_utils.hpp"           // for string_starts_with

namespace stk { namespace mesh { class Bucket; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk {
  namespace io {
    bool is_field_on_part(const stk::mesh::FieldBase *field,
                          const stk::mesh::EntityRank partType,
                          const stk::mesh::Part &part);

    stk::mesh::EntityRank get_entity_rank(const Ioss::GroupingEntity *entity,
                                          const stk::mesh::MetaData &meta)
    {
      switch (entity->type()) {
      case Ioss::NODEBLOCK:
        return stk::topology::NODE_RANK;

      case Ioss::NODESET:
        return stk::topology::NODE_RANK;

      case Ioss::ELEMENTBLOCK:
        return stk::topology::ELEMENT_RANK;

      case Ioss::SUPERELEMENT:
        return stk::topology::ELEMENT_RANK;

      case Ioss::SIDESET:
        {
          const Ioss::SideSet *sset = dynamic_cast<const Ioss::SideSet*>(entity);
          assert(sset != nullptr);
          int my_rank = sset->max_parametric_dimension();
          if (my_rank == 2)
            return stk::topology::FACE_RANK;
          if (my_rank == 1)
            return stk::topology::EDGE_RANK;
          if (my_rank == 0)
            return stk::topology::NODE_RANK;
          else
            return stk::mesh::InvalidEntityRank;
        }

      case Ioss::SIDEBLOCK:
        {
          const Ioss::SideBlock *sblk = dynamic_cast<const Ioss::SideBlock*>(entity);
          assert(sblk != nullptr);
          int rank = sblk->topology()->parametric_dimension();
          if (sblk->topology()->shape() == Ioss::ElementShape::UNKNOWN) {
            rank = sblk->owner()->max_parametric_dimension();
          }
          if (rank == 2)
            return stk::topology::FACE_RANK;
          if (rank == 1)
            return stk::topology::EDGE_RANK;
          if (rank == 0)
            return stk::topology::NODE_RANK;
          else
            return stk::mesh::InvalidEntityRank;
        }

      case Ioss::FACEBLOCK:
        return stk::topology::FACE_RANK;

      case Ioss::EDGEBLOCK:
        return stk::topology::EDGE_RANK;

      default:
        return stk::mesh::InvalidEntityRank;
      }
    }

  }
}

void STKIORequire(bool cond)
{
  if (!cond) throw std::runtime_error("");
}

namespace {

  static const std::string invalid("invalid");
  static const std::string scalar("scalar");
  static const std::string vector_2d("vector_2d");
  static const std::string vector_3d("vector_3d");
  static const std::string full_tensor_36("full_tensor_36");
  static const std::string full_tensor_32("full_tensor_32");
  static const std::string full_tensor_22("full_tensor_22");
  static const std::string full_tensor_16("full_tensor_16");
  static const std::string full_tensor_12("full_tensor_12");
  static const std::string sym_tensor_33("sym_tensor_33");
  static const std::string sym_tensor_31("sym_tensor_31");
  static const std::string sym_tensor_21("sym_tensor_21");
  static const std::string matrix_22("matrix_22");
  static const std::string matrix_33("matrix_33");

  const std::string base_stk_part_name = "_base_stk_part_name";

  const stk::mesh::FieldBase *declare_stk_field(stk::mesh::MetaData &meta,
                                                stk::mesh::EntityRank type,
                                                stk::mesh::Part &part,
                                                const Ioss::Field &ioField)
  {
    Ioss::Field::BasicType ioFieldType = ioField.get_type();
    const bool ioFieldTypeIsRecognized = (ioFieldType == Ioss::Field::INTEGER) || (ioFieldType == Ioss::Field::INT64)
                                      || (ioFieldType == Ioss::Field::REAL)    || (ioFieldType == Ioss::Field::COMPLEX);
    STK_ThrowRequireMsg(ioFieldTypeIsRecognized, "Unrecognized field type for IO field '"<<ioField.get_name()<<"'");

    return stk::io::impl::declare_stk_field_internal(meta, type, part, ioField);
  }

  template <typename Tfield, typename Tio>
  void internal_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                     const Ioss::Field &ioField,
                                     const stk::mesh::FieldBase *field,
                                     const std::vector<stk::mesh::Entity> &entities,
                                     Ioss::GroupingEntity *ioEntity)
  {
    size_t iossNumFieldComponents = ioField.transformed_storage()->component_count();

    std::vector<Tio> ioFieldData;
    size_t ioEntityCount = ioEntity->get_field_data(ioField.get_name(), ioFieldData);
    assert(ioFieldData.size() == entities.size() * iossNumFieldComponents);

    size_t entityCount = entities.size();

    if (ioEntityCount != entityCount) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << ioField.get_name()
             << "' on " << ioEntity->type_string() << " " << ioEntity->name()
             << ". The IO system has " << ioEntityCount
             << " entries, but the stk:mesh system has " << entityCount
             << " entries. The two counts must match.";
      throw std::runtime_error(errmsg.str());
    }

    field->sync_to_host();
    field->modify_on_host();
    for (size_t i=0; i < entityCount; ++i) {
      if (mesh.is_valid(entities[i])) {
        Tfield *fldData = static_cast<Tfield*>(stk::mesh::field_data(*field, entities[i]));
        if (fldData != nullptr) {
          const size_t stkNumFieldComponents = stk::mesh::field_scalars_per_entity(*field, entities[i]);
          const size_t len = std::min(stkNumFieldComponents, iossNumFieldComponents);
          for(size_t j=0; j<len; ++j) {
            fldData[j] = ioFieldData[i*iossNumFieldComponents+j];
          }
        }
      }
    }
  }

  template <typename Tfield, typename Tio>
  void internal_subsetted_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                               const Ioss::Field &ioField,
                                               const stk::mesh::FieldBase *field,
                                               const std::vector<stk::mesh::Entity> &entities,
                                               Ioss::GroupingEntity *ioEntity,
                                               const stk::mesh::Part *stkPart)
  {
    size_t field_componentCount = ioField.transformed_storage()->component_count();
    std::vector<Tio> ioFieldData;
    size_t ioEntityCount = ioEntity->get_field_data(ioField.get_name(), ioFieldData);
    assert(ioFieldData.size() == entities.size() * field_componentCount);
    size_t entityCount = entities.size();
    if (ioEntityCount != entityCount) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << ioField.get_name()
             << "' on " << ioEntity->type_string() << " " << ioEntity->name()
             << ". The IO system has " << ioEntityCount
             << " entries, but the stk:mesh system has " << entityCount
             << " entries. The two counts must match.";
      throw std::runtime_error(errmsg.str());
    }

    stk::mesh::MetaData &meta = stk::mesh::MetaData::get(*stkPart);
    stk::mesh::Selector selector = (meta.globally_shared_part() | meta.locally_owned_part()) & *stkPart;

    field->sync_to_host();
    field->modify_on_host();
    for (size_t i=0; i < entityCount; ++i) {
      if (mesh.is_valid(entities[i])) {
        const stk::mesh::Bucket &bucket = mesh.bucket(entities[i]);
        if (selector(bucket)) {
          Tfield *fldData = static_cast<Tfield*>(stk::mesh::field_data(*field, entities[i]));
          if (fldData !=nullptr) {
            for(size_t j=0; j<field_componentCount; ++j) {
              fldData[j] = ioFieldData[i*field_componentCount+j];
            }
          }
        }
      }
    }
  }

  template <typename T>
  void internal_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                                   const Ioss::Field &ioField,
                                   const stk::mesh::FieldBase *field,
                                   std::vector<stk::mesh::Entity> &entities,
                                   Ioss::GroupingEntity *ioEntity)
  {
    auto io_db = ioEntity->get_database();
    const auto supports = io_db->entity_field_support();
    
    if (!(ioEntity->type() & supports)) {
      return;
    }
    int iossFieldLength = ioField.transformed_storage()->component_count();
    size_t entityCount = entities.size();

    std::vector<T> ioFieldData(entityCount*iossFieldLength);

    field->sync_to_host();
    const stk::mesh::Bucket* prevBkt = nullptr;
    int stkFieldLength = 0;
    int length = 0;
    for (size_t i=0; i < entityCount; ++i) {
      if (mesh.is_valid(entities[i]) && mesh.entity_rank(entities[i]) == field->entity_rank()) {
        const T *fldData = static_cast<T*>(stk::mesh::field_data(*field, entities[i]));
        if (fldData != nullptr) {
          const stk::mesh::Bucket* curBkt = mesh.bucket_ptr(entities[i]);
          if (curBkt != prevBkt) {
            prevBkt = curBkt;
            stkFieldLength = stk::mesh::field_scalars_per_entity(*field, *curBkt);
            STK_ThrowRequireMsg((iossFieldLength >= stkFieldLength), "Field "<<field->name()<<" scalars-per-entity="<<stkFieldLength<<" doesn't match Ioss iossFieldLength(="<<iossFieldLength<<") for io_entity "<<ioEntity->name());
            length = std::min(iossFieldLength, stkFieldLength);
          }
          T* ioFieldDataPtr = ioFieldData.data()+i*iossFieldLength;
          for(int j=0; j<length; ++j) {
            ioFieldDataPtr[j] = fldData[j];
          }
        }
      }
    }

    size_t ioEntityCount = ioEntity->put_field_data(ioField.get_name(), ioFieldData);
    assert(ioFieldData.size() == entities.size() * iossFieldLength);

    if (ioEntityCount != entityCount) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << ioField.get_name()
             << "' on " << ioEntity->type_string() << " " << ioEntity->name()
             << ". The IO system has " << ioEntityCount
             << " entries, but the stk:mesh system has " << entityCount
             << " entries. The two counts must match.";
      throw std::runtime_error(errmsg.str());
    }
  }

  bool will_output_lower_rank_fields(const stk::mesh::Part &part, stk::mesh::EntityRank rank)
  {
    // Determine whether the part 'part' needs to output any fields of rank 'rank'
    // to the output database.

    const std::vector<stk::mesh::FieldBase *> &fields = stk::mesh::MetaData::get(part).get_fields();
    std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
    while (I != fields.end()) {
      const stk::mesh::FieldBase *f = *I ; ++I ;

      bool validPartField = stk::io::is_valid_part_field(f, rank, part, Ioss::Field::TRANSIENT);

      if (validPartField) {
        return true;
      }
    }
    return false;
  }

  void add_canonical_name_property(Ioss::GroupingEntity* ge, stk::mesh::Part& part)
  {
    if(stk::io::has_alternate_part_name(part)) {
      std::string canonName = stk::io::get_alternate_part_name(part);
      if(canonName != ge->name()) {
        ge->property_add(Ioss::Property("db_name", canonName));
      }
    }
  }

  void add_original_topology_property(Ioss::GroupingEntity* ge, stk::mesh::Part& part)
  {
    std::string topoString("original_topology_type");

    if(stk::io::has_original_topology_type(part)) {
      std::string origTopology = stk::io::get_original_topology_type(part);
      if(!ge->property_exists(topoString) || (origTopology != ge->get_property(topoString).get_string())) {
        ge->property_add(Ioss::Property(topoString, origTopology));
      }
    }
  }

  void process_element_attributes_for_define(stk::io::OutputParams &params, stk::mesh::Part &part)
  {
      stk::mesh::EntityRank rank = stk::mesh::EntityRank::ELEM_RANK;

      STK_ThrowRequireMsg(part.primary_entity_rank() == rank, "Input part is not ELEM_RANK");

      const std::vector<stk::io::FieldAndName>& additionalFields = params.get_additional_attribute_fields();

      Ioss::Region & ioRegion = params.io_region();
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      Ioss::ElementBlock* ioBlock = ioRegion.get_element_block(stk::io::getPartName(part));

      for(const stk::io::FieldAndName& attribute : additionalFields) {
          if(attribute.apply_to_entity(ioBlock)) {
              const stk::mesh::FieldBase *stkField = attribute.field();

              STK_ThrowRequireMsg(stkField->entity_rank() == rank, "Input attribute field: " + stkField->name() + " is not ELEM_RANK");
              stk::mesh::PartVector relevantParts;
              relevantParts.push_back(&part);
              stk::io::superset_mesh_parts(part, relevantParts);
              relevantParts.push_back(&meta.universal_part());

              if(stkField->defined_on_any(relevantParts)) {
                  const std::string dbName = attribute.db_name();
                  if(!ioBlock->field_exists(dbName)) {
                      int ebSize = ioBlock->get_property("entity_count").get_int();

                      const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*stkField, rank, relevantParts);
                      STK_ThrowRequireMsg(res.num_scalars_per_entity() != 0,
                                      "Could not find a restriction for field: " + stkField->name() + " on part: " + part.name());
                      stk::io::FieldType fieldType;
                      stk::io::get_io_field_type(stkField, res, &fieldType);
                      ioBlock->field_add(Ioss::Field(dbName, fieldType.type, fieldType.name, Ioss::Field::ATTRIBUTE, ebSize));
                  }
              }
          }
      }
  }

  void process_element_attributes_for_output(stk::io::OutputParams &params, stk::mesh::Part &part)
  {
      stk::mesh::EntityRank rank = stk::mesh::EntityRank::ELEM_RANK;

      STK_ThrowRequireMsg(part.primary_entity_rank() == rank, "Input part is not ELEM_RANK");

      const std::vector<stk::io::FieldAndName>& additionalFields = params.get_additional_attribute_fields();

      Ioss::Region & ioRegion = params.io_region();
      Ioss::ElementBlock* ioBlock = ioRegion.get_element_block(stk::io::getPartName(part));

      for(const stk::io::FieldAndName& attribute : additionalFields) {
          if(attribute.apply_to_entity(ioBlock)) {
              const stk::mesh::FieldBase *stkField = attribute.field();

              STK_ThrowRequireMsg(stkField->entity_rank() == rank, "Input attribute field: " + stkField->name() + " is not ELEM_RANK");

              const std::string dbName = attribute.db_name();

              std::vector<stk::mesh::Entity> entities;
              stk::io::get_output_entity_list(ioBlock, rank, params, entities);
              stk::io::field_data_to_ioss(params.bulk_data(), stkField, entities, ioBlock, dbName, Ioss::Field::ATTRIBUTE);
          }
      }
  }

  bool contain(const stk::mesh::BulkData& stkmesh, stk::mesh::Entity elem, const stk::mesh::Part* parentBlock)
  {
    const stk::mesh::PartVector& parts = stkmesh.bucket(elem).supersets();

    unsigned int partId = parentBlock->mesh_meta_data_ordinal();
    auto i = parts.begin();
    for(; i != parts.end() && (*i)->mesh_meta_data_ordinal() != partId; ++i)
      ;

    return (i != parts.end());
  }

}//namespace <empty>


namespace stk {
namespace io {

namespace impl {

const Ioss::VariableType * get_variable_type_from_factory(const std::string & typeName)
{
  static Ioss::StorageInitializer initializeStorage;

  const Ioss::VariableType * variableType = nullptr;
  try {
    variableType = Ioss::VariableType::factory(typeName);
  }
  catch (...) {
  }

  return variableType;
}

void set_field_output_type(stk::mesh::FieldBase & field, const Ioss::VariableType * type)
{
  mesh::MetaData & meta = mesh::MetaData::get(field);
  const Ioss::VariableType * oldVariableType = field.attribute<Ioss::VariableType>();
  if (not oldVariableType) {
    meta.declare_attribute_no_delete(field, type);
  }
  else {
    if (oldVariableType->name() != type->name()) {
      const bool success = meta.remove_attribute(field, oldVariableType);
      STK_ThrowRequireMsg(success, "stk::io::impl::set_field_output_type(): Failed to remove old attribute " +
                      oldVariableType->name() + " from field " + field.name());
      meta.declare_attribute_no_delete(field, type);
    }
  }
}

void set_field_output_type(stk::mesh::FieldBase & field, const std::string & typeName)
{
  const Ioss::VariableType * variableType = get_variable_type_from_factory(typeName);

  if (not variableType) {
    STK_ThrowErrorMsg("Unrecognized Field output type '" + typeName + "'.  Valid choices with output subscripts are:\n"
                  "  - scalar\n"
                  "  - vector_2d      [x, y]\n"
                  "  - vector_3d      [x, y, z]\n"
                  "  - full_tensor_36 [xx, yy, zz, xy, yz, zx, yx, zy, xz]\n"
                  "  - full_tensor_32 [xx, yy, zz, xy, yx]\n"
                  "  - full_tensor_22 [xx, yy, xy, yx]\n"
                  "  - full_tensor_16 [xx, xy, yz, zx, yx, zy, xz]\n"
                  "  - full_tensor_12 [xx, xy, yx]\n"
                  "  - sym_tensor_33  [xx, yy, zz, xy, yz, zx]\n"
                  "  - sym_tensor_31  [xx, yy, zz, xy]\n"
                  "  - sym_tensor_21  [xx, yy, xy]\n"
                  "  - sym_tensor_13  [xx, xy, yz, zx]\n"
                  "  - sym_tensor_11  [xx, xy]\n"
                  "  - sym_tensor_10  [xx]\n"
                  "  - asym_tensor_03 [xy, yz, zx]\n"
                  "  - asym_tensor_02 [xy, yz]\n"
                  "  - asym_tensor_01 [xy]\n"
                  "  - matrix_22      [xx, xy, yx, yy]\n"
                  "  - matrix_33      [xx, xy, xz, yx, yy, yz, zx, zy, zz]\n"
                  "  - quaternion_2d  [s, q]\n"
                  "  - quaternion_3d  [x, y, z, q]\n"
                  "  - Custom named-suffix output type [user-defined]\n"
                  "Default if unspecified: Scalar or generic array [1, 2, 3, ...]");
  }

  impl::set_field_output_type(field, variableType);
}

const Ioss::VariableType * get_field_output_variable_type(const stk::mesh::FieldBase & field)
{
  return field.attribute<Ioss::VariableType>();
}

const stk::mesh::FieldBase *declare_stk_field_internal(stk::mesh::MetaData &meta,
                                                       stk::mesh::EntityRank type,
                                                       stk::mesh::Part &part,
                                                       const Ioss::Field &io_field)
{
  std::string name = io_field.get_name();
  stk::mesh::FieldBase *field = meta.get_field(type, name);
  // If the field has already been declared, don't redeclare it.
  if (field != nullptr && stk::io::is_field_on_part(field, type, part)) {
    return field;
  }

  stk::topology::rank_t entityRank = static_cast<stk::topology::rank_t>(type);

  const Ioss::VariableType* varType = io_field.transformed_storage();
  size_t numComponents = varType->component_count();
  size_t numCopies = 1;

  const Ioss::CompositeVariableType* compositeVarType = dynamic_cast<const Ioss::CompositeVariableType*>(varType);
  if (compositeVarType != nullptr) {
    const Ioss::VariableType * baseVarType = compositeVarType->get_base_type();
    numComponents = baseVarType->component_count();
    numCopies = compositeVarType->get_num_copies();
    varType = baseVarType;
  }
  std::string field_type = varType->name();

  field = &meta.declare_field<double>(entityRank, name);
  stk::mesh::put_field_on_mesh(*field, part, numComponents, numCopies, nullptr);

  const int oldVarTypeSize = has_field_output_type(*field) ? get_field_output_variable_type(*field)->component_count() : 0;
  const int newVarTypeSize = varType->component_count();

  if (newVarTypeSize > oldVarTypeSize) {
    set_field_output_type(*field, varType);
  }

  stk::io::set_field_role(*field, io_field.get_role());

  return field;
}

} //namespace impl


    struct IossDerivedNodesetAttribute
    {
      bool value;
    };

    void set_derived_nodeset_attribute(stk::mesh::Part& part, const bool hasAttribute)
    {
      stk::mesh::impl::set_part_attribute<IossDerivedNodesetAttribute>(part, hasAttribute);
    }

    bool has_derived_nodeset_attribute(stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossDerivedNodesetAttribute>(part);
    }

    bool get_derived_nodeset_attribute(stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossDerivedNodesetAttribute>(part);
    }


    struct IossOriginalPartId
    {
      int64_t value;
    };

    void set_original_part_id(stk::mesh::Part& part, const int64_t originalId)
    {
      stk::mesh::impl::set_part_attribute<IossOriginalPartId>(part, originalId);
    }

    bool has_original_part_id(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossOriginalPartId>(part);
    }

    int64_t get_original_part_id(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossOriginalPartId>(part);
    }


    struct IossOriginalBlockOrder
    {
      int64_t value;
    };

    void set_original_block_order(stk::mesh::Part& part, const int64_t originalBlockOrder)
    {
      stk::mesh::impl::set_part_attribute<IossOriginalBlockOrder>(part, originalBlockOrder);
    }

    bool has_original_block_order(stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossOriginalBlockOrder>(part);
    }

    int64_t get_original_block_order(stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossOriginalBlockOrder>(part);
    }


    struct IossOriginalTopologyType
    {
      std::string value;
    };

    void set_original_topology_type(stk::mesh::Part& part)
    {
      stk::mesh::impl::set_part_attribute<IossOriginalTopologyType>(part, part.topology().name());
    }

    void set_original_topology_type(stk::mesh::Part& part, const std::string& origTopo)
    {
      stk::mesh::impl::set_part_attribute<IossOriginalTopologyType>(part, origTopo);
    }

    std::string get_original_topology_type(stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossOriginalTopologyType>(part);
    }

    bool has_original_topology_type(stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossOriginalTopologyType>(part);
    }

    struct IossTopologyType
    {
      std::string value;
    };

    void set_topology_type(stk::mesh::Part& part)
    {
      stk::mesh::impl::set_part_attribute<IossTopologyType>(part, part.topology().name());
    }

    void set_topology_type(stk::mesh::Part& part, const std::string& topo)
    {
      stk::mesh::impl::set_part_attribute<IossTopologyType>(part, topo);
    }

    std::string get_topology_type(stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossTopologyType>(part);
    }

    bool has_topology_type(stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossTopologyType>(part);
    }

    struct IossAlternatePartName
    {
      std::string value;
    };

    void set_alternate_part_name(stk::mesh::Part& part, const std::string& altPartName)
    {
      stk::mesh::impl::set_unique_part_attribute<IossAlternatePartName>(part, altPartName);

      mesh::MetaData & meta = mesh::MetaData::get(part);
      meta.add_part_alias(part, altPartName);
    }

    bool has_alternate_part_name(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossAlternatePartName>(part);
    }

    std::string get_alternate_part_name(const stk::mesh::Part& part)
    {
        std::string name = "";

        if(has_alternate_part_name(part)) {
          name = stk::mesh::impl::get_part_attribute<IossAlternatePartName>(part);
        }
        return name;
    }

    struct IossFaceBlockPartAttribute
    {
      bool value;
    };

    void set_face_block_part_attribute(stk::mesh::Part& part, const bool isFaceBlockPart)
    {
      stk::mesh::impl::set_part_attribute<IossFaceBlockPartAttribute>(part, isFaceBlockPart);
    }

    bool has_face_block_part_attribute(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossFaceBlockPartAttribute>(part);
    }

    bool get_face_block_part_attribute(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossFaceBlockPartAttribute>(part);
    }

    struct IossEdgeBlockPartAttribute
    {
      bool value;
    };

    void set_edge_block_part_attribute(stk::mesh::Part& part, const bool isEdgeBlockPart)
    {
      stk::mesh::impl::set_part_attribute<IossEdgeBlockPartAttribute>(part, isEdgeBlockPart);
    }

    bool has_edge_block_part_attribute(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::has_part_attribute<IossEdgeBlockPartAttribute>(part);
    }

    bool get_edge_block_part_attribute(const stk::mesh::Part& part)
    {
      return stk::mesh::impl::get_part_attribute<IossEdgeBlockPartAttribute>(part);
    }

    std::string getPartName(const stk::mesh::Part& part)
    {
      std::string apn = get_alternate_part_name(part);
      if (apn.length())
        {
          return apn;
        }
      else
        return part.name();
    }

    Ioss::GroupingEntity* get_grouping_entity(const Ioss::Region& region, const stk::mesh::Part& part)
    {
      if(!stk::io::is_part_io_part(part)) { return nullptr; }

      Ioss::GroupingEntity* entity = nullptr;
      std::string partName = stk::io::getPartName(part);
      std::vector<Ioss::EntityType> types = get_ioss_entity_types(part);

      for(Ioss::EntityType type : types) {
        entity = region.get_entity(partName, type);
        if(entity != nullptr) {
          return entity;
        }
      }
      return nullptr;
    }

    std::vector<Ioss::EntityType> get_ioss_entity_types(const stk::mesh::MetaData& meta, stk::mesh::EntityRank rank)
    {
      std::vector<Ioss::EntityType> types;
      switch(rank) {
        case stk::topology::NODE_RANK:
          types.push_back(Ioss::NODEBLOCK);
          types.push_back(Ioss::NODESET);
          break;
        case stk::topology::ELEMENT_RANK:
          types.push_back(Ioss::ELEMENTBLOCK);
          types.push_back(Ioss::SUPERELEMENT);
          break;
        case stk::topology::EDGE_RANK:
          if(meta.spatial_dimension() == 2) {
            types.push_back(Ioss::SIDESET);
            types.push_back(Ioss::SIDEBLOCK);
          }
          break;
        case stk::topology::FACE_RANK:
          if(meta.spatial_dimension() == 3) {
            types.push_back(Ioss::SIDESET);
            types.push_back(Ioss::SIDEBLOCK);
          }
          break;
        default:
          break;
      }
      return types;
    }

    std::vector<Ioss::EntityType> get_ioss_entity_types(const stk::mesh::Part& part)
    {
      return get_ioss_entity_types(part.mesh_meta_data(), part.primary_entity_rank());
    }

    size_t db_api_int_size(const Ioss::GroupingEntity *entity)
    {
      return entity->get_database()->int_byte_size_api();
    }

    void initialize_spatial_dimension(stk::mesh::MetaData & meta, size_t spatialDimension,
                                      const std::vector<std::string> &entityRankNames)
    {
      if (!meta.is_initialized() ) {
        meta.initialize(spatialDimension, entityRankNames);
      }
    }

    bool is_field_on_part(const stk::mesh::FieldBase *field,
                          const stk::mesh::EntityRank partType,
                          const stk::mesh::Part &part)
    {
      const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(part);
      const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*field, partType, part);
      if (res.num_scalars_per_entity() > 0) {
        // The field exists on the current 'part'.  Now check (for
        // node types only) whether the 'part' is *either* the
        // universal_part() *or* the field *doesn't* exist on the
        // universal part...
        // Note that for "node" type parts, the IO database has a part
        // (nodeblock) that corresponds to the universal_part(), so
        // fields defined on all nodes are output on the nodeblock and
        // fields defined on only a subset of the parts should be
        // output on that subset which maps to a nodeset.  For other
        // part types, there is no IO entity that corresponds to a
        // universal part, so fields must be output on the part they
        // exist on.  There may be a problem if we start using element
        // sets ..., but wait until we get to that point; current code
        // works with current entity set.
        if (partType != stk::topology::NODE_RANK || part == meta.universal_part()) {
          return true;
        }

        const stk::mesh::FieldBase::Restriction &universalRes = stk::mesh::find_restriction(*field, partType, meta.universal_part());
        if (universalRes.num_scalars_per_entity() <= 0) {
          // Field exists on current part, but not on the universal
          // set (and this part is not the universal part)
          return true;
        }
      }
      return false;
    }

    bool is_valid_part_field(const stk::mesh::FieldBase *field,
                             const stk::mesh::EntityRank partType,
                             const stk::mesh::Part &part,
                             const Ioss::Field::RoleType filterRole)
    {
      const Ioss::Field::RoleType *role = stk::io::get_field_role(*field);

      if (role == nullptr) {
        return false;
      }

      if (role != nullptr && *role != filterRole)
        return false;

      return is_field_on_part(field, partType, part);
    }

    void assign_generic_field_type(const stk::mesh::FieldRestriction &res, FieldType *result)
    {
      int scalarsPerEntity = res.num_scalars_per_entity();
      int firstDimension = res.dimension();

      if (firstDimension == 1) {
        result->name = scalar;
      }
      else {
        result->name = "Real[" + std::to_string(firstDimension) + "]";
      }

      result->copies = scalarsPerEntity / firstDimension;
    }

    void get_io_field_type(const stk::mesh::FieldBase *field,
                           const stk::mesh::FieldRestriction &res,
                           FieldType *result)
    {
      result->type = Ioss::Field::INVALID;

      if ( field->type_is<double>() ) {
        result->type = Ioss::Field::REAL;
      }
      else if ( field->type_is<int>() || field->type_is<uint32_t>()) {
        result->type = Ioss::Field::INTEGER;
      }
      else if ( field->type_is<int64_t>() || field->type_is<uint64_t>()) {
        result->type = Ioss::Field::INT64;
      }

      const int scalarsPerEntity = res.num_scalars_per_entity();
      const int firstDimension = res.dimension();

      result->copies = 1;

      if (has_field_output_type(*field)) {
        const Ioss::VariableType * variableType = get_field_output_variable_type(*field);
        const std::string variableTypeName = Ioss::Utils::lowercase(variableType->name());

        if (variableTypeName == vector_3d || variableTypeName == vector_2d) {
          if (firstDimension == 3) {
            result->name = vector_3d;
            result->copies = scalarsPerEntity / firstDimension;
          }
          else if (firstDimension == 2) {
            result->name = vector_2d;
            result->copies = scalarsPerEntity / firstDimension;
          }
          else {
            assign_generic_field_type(res, result);
          }
        }
        else {
          if (firstDimension == variableType->component_count()) {
            result->name = variableType->name();
            result->copies = scalarsPerEntity / firstDimension;
          }
          else {
            assign_generic_field_type(res, result);
          }
        }
      }
      else {
        assign_generic_field_type(res, result);
      }
    }

    void create_named_suffix_field_output_type(const std::string & typeName, const std::vector<std::string> & suffices)
    {
      Ioss::VariableType::create_named_suffix_type(typeName, suffices);
    }

    void set_named_suffix_field_output_type(stk::mesh::FieldBase & field, const std::string & typeName)
    {
      const Ioss::VariableType * variableType = impl::get_variable_type_from_factory(typeName);

      if (not variableType) {
        STK_ThrowErrorMsg("Unrecognized custom named suffix Field output type '" + typeName + "'.\n"
                      "Please be sure to pre-register your custom output type with a call to:\n"
                      "  stk::io::create_named_suffix_field_output_type()");
      }

      for (unsigned i = 0; i < field.number_of_states(); ++i) {
        stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(i);
        stk::mesh::FieldBase * fieldOfState = field.field_state(state);

        impl::set_field_output_type(*fieldOfState, variableType);
      }
    }

    void set_field_output_type(stk::mesh::FieldBase & field, FieldOutputType fieldOutputType)
    {
      for (unsigned i = 0; i < field.number_of_states(); ++i) {
        stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(i);
        stk::mesh::FieldBase * fieldOfState = field.field_state(state);

        switch (fieldOutputType) {
          case (FieldOutputType::SCALAR)         : impl::set_field_output_type(*fieldOfState, "scalar"); break;
          case (FieldOutputType::VECTOR_2D)      : impl::set_field_output_type(*fieldOfState, "vector_2d"); break;
          case (FieldOutputType::VECTOR_3D)      : impl::set_field_output_type(*fieldOfState, "vector_3d"); break;
          case (FieldOutputType::FULL_TENSOR_36) : impl::set_field_output_type(*fieldOfState, "full_tensor_36"); break;
          case (FieldOutputType::FULL_TENSOR_32) : impl::set_field_output_type(*fieldOfState, "full_tensor_32"); break;
          case (FieldOutputType::FULL_TENSOR_22) : impl::set_field_output_type(*fieldOfState, "full_tensor_22"); break;
          case (FieldOutputType::FULL_TENSOR_16) : impl::set_field_output_type(*fieldOfState, "full_tensor_16"); break;
          case (FieldOutputType::FULL_TENSOR_12) : impl::set_field_output_type(*fieldOfState, "full_tensor_12"); break;
          case (FieldOutputType::SYM_TENSOR_33)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_33"); break;
          case (FieldOutputType::SYM_TENSOR_31)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_31"); break;
          case (FieldOutputType::SYM_TENSOR_21)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_21"); break;
          case (FieldOutputType::SYM_TENSOR_13)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_13"); break;
          case (FieldOutputType::SYM_TENSOR_11)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_11"); break;
          case (FieldOutputType::SYM_TENSOR_10)  : impl::set_field_output_type(*fieldOfState, "sym_tensor_10"); break;
          case (FieldOutputType::ASYM_TENSOR_03) : impl::set_field_output_type(*fieldOfState, "asym_tensor_03"); break;
          case (FieldOutputType::ASYM_TENSOR_02) : impl::set_field_output_type(*fieldOfState, "asym_tensor_02"); break;
          case (FieldOutputType::ASYM_TENSOR_01) : impl::set_field_output_type(*fieldOfState, "asym_tensor_01"); break;
          case (FieldOutputType::MATRIX_22)      : impl::set_field_output_type(*fieldOfState, "matrix_22"); break;
          case (FieldOutputType::MATRIX_33)      : impl::set_field_output_type(*fieldOfState, "matrix_33"); break;
          case (FieldOutputType::QUATERNION_2D)  : impl::set_field_output_type(*fieldOfState, "quaternion_2d"); break;
          case (FieldOutputType::QUATERNION_3D)  : impl::set_field_output_type(*fieldOfState, "quaternion_3d"); break;
          case (FieldOutputType::CUSTOM) : {
            STK_ThrowErrorMsg("To set a custom FieldOutputType, please call stk::io::create_named_suffix_field_output_type()"
                          " and stk::io::set_named_suffix_field_output_type() instead");
            break;
          }
          default:
            STK_ThrowErrorMsg("Unsupported FieldOutputType for Field " << field.name() << ": "
                          << static_cast<int>(fieldOutputType));
        }
      }
    }

    FieldOutputType get_field_output_type(const stk::mesh::FieldBase & field)
    {
      const Ioss::VariableType * variableType = field.attribute<Ioss::VariableType>();

      if (variableType == nullptr) {
        return FieldOutputType::SCALAR;
      }

      const std::string variableTypeName = variableType->name();

      if (variableTypeName == "scalar")              return FieldOutputType::SCALAR;
      else if (variableTypeName == "vector_2d")      return FieldOutputType::VECTOR_2D;
      else if (variableTypeName == "vector_3d")      return FieldOutputType::VECTOR_3D;
      else if (variableTypeName == "full_tensor_36") return FieldOutputType::FULL_TENSOR_36;
      else if (variableTypeName == "full_tensor_32") return FieldOutputType::FULL_TENSOR_32;
      else if (variableTypeName == "full_tensor_22") return FieldOutputType::FULL_TENSOR_22;
      else if (variableTypeName == "full_tensor_16") return FieldOutputType::FULL_TENSOR_16;
      else if (variableTypeName == "full_tensor_12") return FieldOutputType::FULL_TENSOR_12;
      else if (variableTypeName == "sym_tensor_33")  return FieldOutputType::SYM_TENSOR_33;
      else if (variableTypeName == "sym_tensor_31")  return FieldOutputType::SYM_TENSOR_31;
      else if (variableTypeName == "sym_tensor_21")  return FieldOutputType::SYM_TENSOR_21;
      else if (variableTypeName == "sym_tensor_13")  return FieldOutputType::SYM_TENSOR_13;
      else if (variableTypeName == "sym_tensor_11")  return FieldOutputType::SYM_TENSOR_11;
      else if (variableTypeName == "sym_tensor_10")  return FieldOutputType::SYM_TENSOR_10;
      else if (variableTypeName == "asym_tensor_03") return FieldOutputType::ASYM_TENSOR_03;
      else if (variableTypeName == "asym_tensor_02") return FieldOutputType::ASYM_TENSOR_02;
      else if (variableTypeName == "asym_tensor_01") return FieldOutputType::ASYM_TENSOR_01;
      else if (variableTypeName == "matrix_22")      return FieldOutputType::MATRIX_22;
      else if (variableTypeName == "matrix_33")      return FieldOutputType::MATRIX_33;
      else if (variableTypeName == "quaternion_2d")  return FieldOutputType::QUATERNION_2D;
      else if (variableTypeName == "quaternion_3d")  return FieldOutputType::QUATERNION_3D;
      else {
        return FieldOutputType::CUSTOM;
      }
    }

    bool has_field_output_type(const stk::mesh::FieldBase & field)
    {
      return (field.attribute<Ioss::VariableType>() != nullptr);
    }

    void set_field_output_variable_type(stk::mesh::FieldBase & field, const Ioss::VariableType * type)
    {
      if (type != nullptr) {
        for (unsigned i = 0; i < field.number_of_states(); ++i) {
          stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(i);
          stk::mesh::FieldBase * fieldOfState = field.field_state(state);

          impl::set_field_output_type(*fieldOfState, type);
        }
      }
    }

    const Ioss::VariableType * get_field_output_variable_type(const stk::mesh::FieldBase & field) {
      return field.attribute<Ioss::VariableType>();
    }

    bool has_io_part_attribute(mesh::Part & part)
    {
        return is_part_io_part(part);
    }

    struct IossPartAttribute
    {
      bool value;
    };

    void put_io_part_attribute(mesh::Part & part)
    {
      if (!is_part_io_part(part)) {
        stk::mesh::impl::set_part_attribute<IossPartAttribute>(part, true);
      }
    }

    bool is_part_face_block_io_part(const mesh::Part &part)
    {
      const bool isIoPart = is_part_io_part(part);
      const bool isFaceBlockPart = isIoPart ?
           (has_face_block_part_attribute(part) && get_face_block_part_attribute(part)) : false;
      return isFaceBlockPart;
    }

    void put_face_block_io_part_attribute(mesh::Part & part)
    {
      STK_ThrowRequireMsg(part.mesh_meta_data().spatial_dimension() == 3, "face-block IO attribute can not be used for 2D.");
      STK_ThrowRequireMsg(part.primary_entity_rank() == stk::topology::FACE_RANK, "face-block IO attribute must only be used for face-rank parts.");

      if (!is_part_face_block_io_part(part)) {
        put_io_part_attribute(part);
        set_face_block_part_attribute(part, true);
      }
    }

    bool is_part_edge_block_io_part(const mesh::Part &part)
    {
      const bool isIoPart = is_part_io_part(part);
      const bool isEdgeBlockPart = isIoPart ?
           (has_edge_block_part_attribute(part) && get_edge_block_part_attribute(part)) : false;
      return isEdgeBlockPart;
    }

    void put_edge_block_io_part_attribute(mesh::Part & part)
    {
      STK_ThrowRequireMsg(part.primary_entity_rank() == stk::topology::EDGE_RANK, "edge-block IO attribute must only be used for edge-rank parts.");

      if (!is_part_edge_block_io_part(part)) {
        put_io_part_attribute(part);
        set_edge_block_part_attribute(part, true);
      }
    }

    struct IossAssemblyPartAttribute
    {
      bool value;
    };

    bool is_part_assembly_io_part(const stk::mesh::Part &part)
    {
      return stk::mesh::impl::has_part_attribute<IossAssemblyPartAttribute>(part);
    }

    void put_assembly_io_part_attribute(mesh::Part & part)
    {
      if (!is_part_assembly_io_part(part)) {
        put_io_part_attribute(part);
        stk::mesh::impl::set_part_attribute<IossAssemblyPartAttribute>(part, true);
      }
    }

    std::vector<std::string> get_assembly_names(const stk::mesh::MetaData& meta)
    {
      std::vector<std::string> assemblyNames;

      for(const stk::mesh::Part* part : meta.get_parts()) {
        if (is_part_assembly_io_part(*part)) {
          assemblyNames.push_back(part->name());
        }
      }

      return assemblyNames;
    }

    bool is_in_subsets_of_parts(const stk::mesh::Part& part,
                                const stk::mesh::PartVector& parts)
    {
      for(const stk::mesh::Part* thisPart : parts) {
        if (stk::mesh::contains(thisPart->subsets(), part)) {
          return true;
        }
      }
      return false;
    }

    std::vector<std::string> get_sub_assembly_names(const stk::mesh::MetaData& meta,
                                                    const std::string& assemblyName)
    {
      const stk::mesh::Part* part = meta.get_part(assemblyName);
      std::vector<std::string> assemblyNames;
      const stk::mesh::PartVector& subsets = part->subsets();

      for(const stk::mesh::Part* subset : subsets) {
        if (is_part_assembly_io_part(*subset) &&
            !is_in_subsets_of_parts(*subset, subsets))
        {
          assemblyNames.push_back(subset->name());
        }
      }

      return assemblyNames;
    }

    bool has_sub_assemblies(const stk::mesh::MetaData& meta,
                            const std::string& assemblyName)
    {
      const stk::mesh::Part* part = meta.get_part(assemblyName);
      for(const stk::mesh::Part* subset : part->subsets()) {
        if (is_part_assembly_io_part(*subset)) {
          return true;
        }
      }
      return false;
    }

    stk::mesh::PartVector get_unique_leaf_parts(const stk::mesh::Part &assemblyPart)
    {
      stk::mesh::PartVector leafParts;

      if (is_part_assembly_io_part(assemblyPart)) {
        for (stk::mesh::Part *subset : assemblyPart.subsets()) {
          if (!is_part_assembly_io_part(*subset)) {
            leafParts.push_back(subset);
          }
        }
      }
      return leafParts;
    }

    stk::mesh::PartVector get_unique_leaf_parts(const stk::mesh::MetaData &meta, const std::string &assemblyName)
    {
      const stk::mesh::Part *part = meta.get_part(assemblyName);
      return get_unique_leaf_parts(*part);
    }

    void remove_io_part_attribute(mesh::Part & part)
    {
      const IossPartAttribute* ioPartAttr = part.attribute<IossPartAttribute>();
      if (ioPartAttr != nullptr) {
        mesh::MetaData & meta = mesh::MetaData::get(part);
        bool success = meta.remove_attribute(part, ioPartAttr);
        STK_ThrowRequireMsg(success, "stk::io::remove_io_part_attribute(" << part.name() << ") FAILED.");
        delete ioPartAttr;
      }
    }

    bool is_part_element_block_io_part(const stk::mesh::Part &part)
    {
      return (part.primary_entity_rank() == stk::topology::ELEM_RANK) && is_part_io_part(part) &&
             !is_part_assembly_io_part(part) && (part.subsets().size() == 0);
    }

    bool is_part_surface_io_part(const stk::mesh::Part &part)
    {
      return (part.primary_entity_rank() == part.mesh_meta_data().side_rank()) && is_part_io_part(part) &&
             !is_part_assembly_io_part(part);
    }

    stk::topology get_start_topology(const Ioss::ElementTopology* topology, unsigned meshSpatialDimension)
    {
        if (topology->is_element() && topology->spatial_dimension() == (int)meshSpatialDimension)
        {
            return stk::topology::BEGIN_ELEMENT_RANK;
        }
        return stk::topology::BEGIN_TOPOLOGY;
    }

    stk::topology map_ioss_topology_to_stk(const Ioss::ElementTopology *topology,
                                           unsigned meshSpatialDimension)
    {
      stk::topology beginTopo = get_start_topology(topology, meshSpatialDimension);
      for (stk::topology topo=beginTopo; topo < stk::topology::END_TOPOLOGY; ++topo) {
        if (topology->is_alias(topo.name()))
        {
           bool bothAreElements = topology->is_element() && topo.rank()==stk::topology::ELEM_RANK;
           bool dimensionsMatch = topology->spatial_dimension()==(int)topo.dimension();
           bool iossNotElemButParametricDimMatchesDim = !topology->is_element() && topology->parametric_dimension() == (int)topo.dimension();
           if (bothAreElements || dimensionsMatch || iossNotElemButParametricDimMatchesDim) {
               return topo;
           }
        }
      }
      std::string tmpCopy = sierra::make_lower(topology->name().substr(0,5));
      if (tmpCopy == "super")
      {
        return stk::create_superelement_topology(topology->number_nodes());
      }

      return stk::topology::INVALID_TOPOLOGY;
    }

    std::string map_stk_topology_to_ioss(stk::topology topo)
    {
      std::string name = topo.name();

      // FIXME SHELL SIDE TOPOLOGY
      if (topo == stk::topology::SHELL_SIDE_BEAM_2) { name = "edge3d2"; }
      if (topo == stk::topology::SHELL_SIDE_BEAM_3) { name = "edge3d3"; }

      Ioss::ElementTopology *iossTopo = Ioss::ElementTopology::factory(name, true);
      return iossTopo != nullptr ? iossTopo->name() : "invalid";
    }

    template<typename ENTITY>
    void set_io_part_attribute(const ENTITY* entity, stk::mesh::Part& part)
    {
      if (entity->type() == Ioss::FACEBLOCK) {
        stk::io::put_face_block_io_part_attribute(part);
      }
      else if (entity->type() == Ioss::EDGEBLOCK) {
        stk::io::put_edge_block_io_part_attribute(part);
      }
      else if (entity->type() == Ioss::ASSEMBLY) {
        stk::io::put_assembly_io_part_attribute(part);
      }
      else {
        stk::io::put_io_part_attribute(part);
      }
    }

    template<typename ENTITY>
    void set_original_topology_type_from_ioss(const ENTITY* entity, stk::mesh::Part& part)
    {
      const std::string origTopoStr("original_topology_type");
      if (entity->property_exists(origTopoStr)) {
        set_original_topology_type(part, entity->get_property(origTopoStr).get_string());
      }
    }

    void declare_stk_aliases(stk::mesh::Part& part, Ioss::GroupingEntity *ge, stk::mesh::MetaData &meta)
    {
      meta.add_part_alias(part, part.name());

      if(nullptr != ge && ge->get_database() != nullptr) {
        Ioss::Region* region = ge->get_database()->get_region();

        if(ge->property_exists("db_name")) {
          std::string canonName = ge->get_property("db_name").get_string();
          meta.add_part_alias(part, canonName);
        }

        const Ioss::AliasMap& ioss_alias_map = region->get_alias_map(ge->type());
        for(auto&& alias : ioss_alias_map) {
          if(stk::equal_case(alias.second, part.name())) {
            meta.add_part_alias(part, alias.first);
          }
        }
      }
    }

    stk::mesh::Part& declare_stk_part(Ioss::GroupingEntity* entity, stk::mesh::MetaData& meta)
    {
      if (entity->type() == Ioss::ASSEMBLY) {
        return meta.declare_part(entity->name());
      }
      else {
        return meta.declare_part(entity->name(), get_entity_rank(entity, meta));
      }
    }

    void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta, TopologyErrorHandler handler)
    {
      if (include_entity(entity)) {
        stk::mesh::Part & part = declare_stk_part(entity, meta);
        declare_stk_aliases(part, entity, meta);
        if (entity->property_exists("id")) {
          meta.set_part_id(part, entity->get_property("id").get_int());
        }
        set_io_part_attribute(entity, part);
      }
    }

    void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta, TopologyErrorHandler handler)
    {
      if (include_entity(entity)) {
        mesh::EntityRank type = get_entity_rank(entity, meta);
        stk::mesh::Part * part = nullptr;
        part = &meta.declare_part(entity->name(), type);
        declare_stk_aliases(*part, entity, meta);
        if (entity->property_exists("id")) {
            meta.set_part_id(*part, entity->get_property("id").get_int());
        }
        set_io_part_attribute(entity, *part);

        const Ioss::ElementTopology *topology = entity->topology();
        // Check spatial dimension of the element topology here so we can
        // issue a more meaningful error message.  If the dimension is bad
        // and we continue to the following calls, there is an exception
        // and we get unintelligible (to the user) error messages.  Could
        // also do a catch...

        if (entity->type() == Ioss::ELEMENTBLOCK) {
          assert(topology != nullptr);
          if (topology->spatial_dimension() < static_cast<int>(meta.spatial_dimension())) {
            // NOTE: The comparison is '<' and not '!=' since a 2D mesh
            // can contain a "3d" element -- a Beam is both a 2D and
            // 3D element...

            std::ostringstream msg ;
            msg << "\n\nERROR: Element Block " << entity->name()
                << " contains " << topology->name() << " elements with spatial dimension "
                << topology->spatial_dimension()
                << "\n       which does not match the spatial dimension of the model which is "
                << meta.spatial_dimension() << "\n\n";
            throw std::runtime_error( msg.str() );
          }

          set_original_topology_type_from_ioss(entity, *part);
        }

        stk::topology stkTopology = map_ioss_topology_to_stk(topology, meta.spatial_dimension());
        if (stkTopology != stk::topology::INVALID_TOPOLOGY) {
          if (stkTopology.rank() != part->primary_entity_rank() && entity->entity_count() == 0) {
            std::ostringstream os;
            os<<"stk_io WARNING: failed to obtain sensible topology for Ioss::GroupingEntity: " << entity->name()<<", iossTopology: "<<topology->name()<<", stk-part: "<<part->name()<<", rank: "<<part->primary_entity_rank()<<", stk-topology: "<<stkTopology<<". Probably because this GroupingEntity is empty on this MPI rank. Unable to set correct stk topology hierarchy. Proceeding, but beware of unexpected behavior."<<std::endl;
            std::cerr<<os.str();
          }
          else {
            stk::mesh::set_topology(*part, stkTopology);
          }
        } else {
          handler(*part);
        }
        stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE, *part, type);
      }
    }

    //----------------------------------------------------------------------
    /** Add all stk::Fields on the entities of the specified part_type
     *  on the specified part of the specified role * to the specified
     *  Ioss::GroupingEntity
     */

    void ioss_add_field_to_entity(const stk::mesh::FieldBase *f,
                                  const stk::mesh::FieldBase::Restriction &res,
                                  Ioss::GroupingEntity *entity,
                                  FieldAndName &namedField,
                                  const Ioss::Field::RoleType filterRole)
    {
        FieldType fieldType;
        get_io_field_type(f, res, &fieldType);
        if ((fieldType.type != Ioss::Field::INVALID) && namedField.apply_to_entity(entity)) {
          size_t entitySize = entity->get_property("entity_count").get_int();
          std::string name = namedField.db_name();
          std::string storage = fieldType.name;

          if (namedField.get_use_alias()) {
              Ioss::VariableType::get_field_type_mapping(f->name(), &storage);
          }

          entity->field_add(Ioss::Field(name, fieldType.type, storage,
                                        fieldType.copies, filterRole, entitySize));
          if (entity->type() == Ioss::NODEBLOCK) {
            namedField.m_forceNodeblockOutput = true;
          }
        }
    }

    bool is_valid_nodeset_field(const stk::mesh::Part &part,
                                       const stk::mesh::EntityRank partType,
                                       Ioss::GroupingEntity *entity,
                                       FieldAndName &namedField,
                                       const Ioss::Field::RoleType filterRole)
    {
        bool isValid = false;

        const stk::mesh::FieldBase *f = namedField.field();
        const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);

        bool isNodeset = (entity != nullptr) && (entity->type() == Ioss::NODESET);
        bool hasMatchingFieldRole = (role != nullptr) ? (*role == filterRole) : false;
        bool hasMatchingEntityRank = f->entity_rank() == partType;
        bool isNodesetField = namedField.is_nodeset_variable();

        if(isNodeset && hasMatchingFieldRole && hasMatchingEntityRank && isNodesetField) {

              if(namedField.apply_to_entity(entity)) {
                  const stk::mesh::EntityRank nodeRank = stk::topology::NODE_RANK;

                  const std::vector<stk::mesh::FieldBase::Restriction> & restrictions = f->restrictions();
                  if (restrictions.size() > 0 && f->entity_rank() == nodeRank) {
                      isValid = true;
                  }
              }
        }

        return isValid;
    }

    void ioss_add_field_to_derived_nodeset(const stk::mesh::Part &part,
                                           const stk::mesh::EntityRank partType,
                                           Ioss::GroupingEntity *entity,
                                           FieldAndName &namedField,
                                           const Ioss::Field::RoleType filterRole)
    {
        const bool isValid = is_valid_nodeset_field(part, partType, entity, namedField, filterRole);

        if(isValid) {
            const stk::mesh::FieldBase *f = namedField.field();
            const stk::mesh::FieldBase::Restriction *res = nullptr;

            const std::vector<stk::mesh::FieldBase::Restriction> & restrictions = f->restrictions();
            if (restrictions.size() > 0 && f->entity_rank() == stk::topology::NODE_RANK) {
              res = restrictions.data();
            }

            if(res != nullptr) {
                ioss_add_field_to_entity(f, *res, entity, namedField, filterRole);
            }
        }
    }

    bool field_should_be_added(const std::string& fieldDbName,
                               unsigned numScalarsPerEntity,
                               Ioss::GroupingEntity* iossEntity)
    {
      if (!iossEntity->field_exists(fieldDbName)) {
        return true;
      }
      Ioss::Field iossField = iossEntity->get_field(fieldDbName);
      size_t fieldComponentCount = iossField.transformed_storage()->component_count();
      if (fieldComponentCount < numScalarsPerEntity) {
        return true;
      }
      return false;
    }

    void ioss_add_fields_for_subpart(const stk::mesh::Part &part,
                                     const stk::mesh::EntityRank partType,
                                     Ioss::GroupingEntity *entity,
                                     FieldAndName &namedField,
                                     const Ioss::Field::RoleType filterRole)
    {
        stk::mesh::EntityRank partRank = part_primary_entity_rank(part);
        stk::mesh::PartVector blocks = part.subsets();
        const stk::mesh::FieldBase *f = namedField.field();

        for (size_t j = 0; j < blocks.size(); j++) {
            mesh::Part & sideBlockPart = *blocks[j];
            bool validSubsetPartField = stk::io::is_valid_part_field(f, partType, sideBlockPart, filterRole);
            Ioss::GroupingEntity* subEntity = entity;

            if (validSubsetPartField) {
                const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, partType, sideBlockPart);
                if (partRank < stk::topology::ELEM_RANK) {
                    Ioss::Region* region = entity->get_database()->get_region();
                    if (nullptr != region) {
                        Ioss::GroupingEntity* tempEntity = region->get_entity(sideBlockPart.name());
                        if (nullptr != tempEntity) {
                          const bool isEntityNodeRankOrSideSetBlock =
                            (tempEntity->type() == Ioss::NODESET ||
                             tempEntity->type() == Ioss::NODEBLOCK ||
                             (tempEntity->type() == Ioss::SIDEBLOCK && entity->type() == Ioss::SIDESET));
                          if (isEntityNodeRankOrSideSetBlock) {
                            subEntity = tempEntity;
                          }

                          if(tempEntity->type() == Ioss::SIDEBLOCK && entity->type() == Ioss::SIDEBLOCK && tempEntity != entity) {
                            subEntity = nullptr;
                          }
                        }
                    }
                }

                bool validIossField = (subEntity == nullptr) ? false :
                                      (namedField.is_nodeset_variable() ? (subEntity->type() == Ioss::NODESET) : true);

                if(validIossField) {
                  if(subEntity->type() == Ioss::SIDEBLOCK && subEntity != entity) {
                      const bool shouldAddFieldToParent =
                          field_should_be_added(namedField.db_name(), res.num_scalars_per_entity(), entity);
                      if (shouldAddFieldToParent) {
                        ioss_add_field_to_entity(f, res, entity, namedField, filterRole);
                      }
                    }

                    ioss_add_field_to_entity(f, res, subEntity, namedField, filterRole);
                }
            }
        }
    }

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank partType,
                         Ioss::GroupingEntity *entity,
                         std::vector<FieldAndName> &namedFields,
                         const Ioss::Field::RoleType filterRole)
    {
        stk::mesh::EntityRank partRank = part_primary_entity_rank(part);
        const stk::mesh::PartVector &blocks = part.subsets();
        bool checkSubparts =  (partRank == stk::topology::NODE_RANK ||
                               partRank == stk::topology::EDGE_RANK ||
                               partRank == stk::topology::FACE_RANK) &&
                              (blocks.size() > 0);
        for (size_t i=0; i<namedFields.size(); i++) {
            const stk::mesh::FieldBase *f = namedFields[i].field();

            if (stk::io::is_valid_part_field(f, partType, part, filterRole)) {
                const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, partType, part);
                ioss_add_field_to_entity(f, res, entity, namedFields[i], filterRole);
            } else if (partRank == namedFields[i].type()) {
                ioss_add_field_to_derived_nodeset(part, partType, entity, namedFields[i], filterRole);
            }

            // If this is a sideset, check the subset parts for the field also...
            if (checkSubparts) {
                ioss_add_fields_for_subpart(part, partType, entity, namedFields[i], filterRole);
            }
        }
    }

    void getNamedFields(const stk::mesh::MetaData &meta,
                        Ioss::GroupingEntity *ioEntity,
                        const Ioss::Field::RoleType filterRole,
                        std::vector<FieldAndName> &namedFields)
    {
      const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();
      namedFields.reserve(fields.size());
      std::vector<stk::mesh::FieldBase *>::const_iterator fieldIterator = fields.begin();
      for(;fieldIterator != fields.end();++fieldIterator) {
          const Ioss::Field::RoleType *role = stk::io::get_field_role(**fieldIterator);
          if (role != nullptr && *role == filterRole) {
              namedFields.emplace_back(*fieldIterator, (*fieldIterator)->name());
          }
      }
    }

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank partType,
                         Ioss::GroupingEntity *entity,
                         const Ioss::Field::RoleType filterRole)
    {
      std::vector<FieldAndName> namedFields;
      stk::io::getNamedFields(mesh::MetaData::get(part), entity, filterRole, namedFields);

      ioss_add_fields(part, partType, entity, namedFields, filterRole);
    }

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank partType,
                         Ioss::GroupingEntity *entity,
                         std::vector<FieldAndName> &namedFields)
    {
      ioss_add_fields(part, partType, entity, namedFields, Ioss::Field::Field::TRANSIENT);
    }


    /**
     * For the given Ioss::GroupingEntity "entity", find all fields that
     * exist on the input database of type "role" and declare them on
     * the give stk::mesh::Part "part". The "part_type" argument
     * specifies the entity type (node, element, ...) that the field
     * should be declared on for this "part"
     *
     * The "role" will typically be either "ATTRIBUTE" or "TRANSIENT"
     */
    void define_io_fields(Ioss::GroupingEntity *entity,
                          Ioss::Field::RoleType role,
                          stk::mesh::Part &part,
                          stk::mesh::EntityRank partType)
    {
      stk::mesh::MetaData &meta = mesh::MetaData::get(part);

      Ioss::NameList names;
      entity->field_describe(role, &names);

      for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        // \todo IMPLEMENT Need a field selection mechanism and a field naming
        // (ioss_name -> stk::name)  For now, select all and give the
        // stk field the same name as the ioss field.

        // Skip the attribute field that is named "attribute"
        if (*I == "attribute" && names.size() > 1)
          continue;

        // \todo IMPLEMENT Need to determine whether these are
        // multi-state fields or constant, or interpolated, or ...
        Ioss::Field ioField = entity->get_field(*I);
        declare_stk_field(meta, partType, part, ioField);
      }
    }

    template <typename T>
    void delete_selector_property(const std::vector<T> &entities)
    {
      for(size_t i=0; i < entities.size(); i++) {
        delete_selector_property(entities[i]);
      }
    }

    void delete_selector_property(Ioss::Region &region)
    {
      // Iterate all Ioss::GroupingEntity types on the io_region and
      // if the have a property named 'selector' of type 'pointer',
      // delete the pointer and remove the property.
      delete_selector_property(region.get_node_blocks());
      delete_selector_property(region.get_element_blocks());
      delete_selector_property(region.get_edge_blocks());
      delete_selector_property(region.get_face_blocks());
      delete_selector_property(region.get_nodesets());
      delete_selector_property(region.get_commsets());

      const Ioss::SideSetContainer& sideSets = region.get_sidesets();
      for(Ioss::SideSetContainer::const_iterator it = sideSets.begin();
          it != sideSets.end(); ++it) {
        Ioss::SideSet *sset = *it;
        delete_selector_property(*it);
        delete_selector_property(sset->get_side_blocks());
      }
    }

    void delete_selector_property(Ioss::GroupingEntity *ioEntity)
    {
      // If the Ioss::GroupingEntity has a property named 'selector' of
      // type 'pointer', delete the pointer and remove the property.
      if (ioEntity->property_exists(s_internalSelectorName)) {
        mesh::Selector *select = reinterpret_cast<mesh::Selector*>(ioEntity->get_property(s_internalSelectorName).get_pointer());
        delete select;
        ioEntity->property_erase(s_internalSelectorName);
      }
    }

    template <typename INT>
    std::vector<stk::mesh::Entity> get_entity_list(Ioss::GroupingEntity *ioEntity,
                                                   stk::mesh::EntityRank partType,
                                                   const stk::mesh::BulkData &bulk)
    {
      std::vector<stk::mesh::Entity> entities;

      if (ioEntity->type() == Ioss::SIDEBLOCK) {
        std::vector<INT> elemSide ;
        ioEntity->get_field_data("element_side", elemSide);
        size_t sideCount = elemSide.size() / 2;
        for(size_t is=0; is<sideCount; ++is)
          entities.push_back(stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk, elemSide[is*2], elemSide[is*2+1]-1, partType));
      }
      else {
        std::vector<INT> ids ;
        ioEntity->get_field_data("ids", ids);

        size_t count = ids.size();
        entities.reserve(count);

        for(size_t i=0; i<count; ++i) {
          entities.push_back(bulk.get_entity( partType, ids[i] ));
        }
      }

      return entities;
    }

    std::vector<stk::mesh::Entity> get_input_entity_list(Ioss::GroupingEntity *ioEntity,
                                                         stk::mesh::EntityRank partType,
                                                         const stk::mesh::BulkData &bulk)
    {
      STK_ThrowRequireMsg(ioEntity->get_database()->is_input(), "Database is output type");
      if (db_api_int_size(ioEntity) == 4) {
          return get_entity_list<int>(ioEntity, partType, bulk);
      } else {
          return get_entity_list<int64_t>(ioEntity, partType, bulk);
      }
    }

    void get_output_entity_list(Ioss::GroupingEntity *ioEntity,
                                stk::mesh::EntityRank partType,
                                OutputParams &params,
                                std::vector<stk::mesh::Entity> &entities)
    {
      const stk::mesh::BulkData &bulk = params.bulk_data();
      STK_ThrowRequireMsg(!ioEntity->get_database()->is_input(), "Database is input type");
      assert(ioEntity->property_exists(s_internalSelectorName));

      mesh::Selector *select = reinterpret_cast<mesh::Selector*>(ioEntity->get_property(s_internalSelectorName).get_pointer());

      if(ioEntity->type() == Ioss::NODEBLOCK) {
          get_selected_nodes(params, *select, entities);
      } else {
          const bool sortById = true;
          stk::mesh::get_entities(bulk, partType, *select, entities, sortById);
      }
    }

    const std::string get_suffix_for_field_at_state(enum stk::mesh::FieldState fieldState, std::vector<std::string>* multiStateSuffixes)
    {
      if(nullptr != multiStateSuffixes) {
          STK_ThrowRequireMsg((multiStateSuffixes->size() >= fieldState),
                          "Invalid field state index '" << fieldState << "'");
          return (*multiStateSuffixes)[fieldState];
      }

      std::string suffix = "";
      switch(fieldState)
        {
        case stk::mesh::StateN:
          suffix = ".N";
          break;
        case stk::mesh::StateNM1:
          suffix = ".NM1";
          break;
        case stk::mesh::StateNM2:
          suffix = ".NM2";
          break;
        case stk::mesh::StateNM3:
          suffix = ".NM3";
          break;
        case stk::mesh::StateNM4:
          suffix = ".NM4";
          break;
        case stk::mesh::StateNP1:
          break;
        default:
          STK_ThrowRequireMsg(false, "Internal Error: Unsupported stk::mesh::FieldState: " << fieldState << ".\n");
        }
      return suffix;
    }

    std::string get_stated_field_name(const std::string &fieldBaseName, stk::mesh::FieldState stateIdentifier,
                                      std::vector<std::string>* multiStateSuffixes)
    {
      std::string field_name_with_suffix = fieldBaseName + get_suffix_for_field_at_state(stateIdentifier, multiStateSuffixes);
      return field_name_with_suffix;
    }

    bool field_state_exists_on_io_entity(const std::string& dbName, const stk::mesh::FieldBase* field, stk::mesh::FieldState stateIdentifier,
                                         Ioss::GroupingEntity *ioEntity, std::vector<std::string>* multiStateSuffixes)
    {
        std::string fieldNameWithSuffix = get_stated_field_name(dbName, stateIdentifier, multiStateSuffixes);
        return ioEntity->field_exists(fieldNameWithSuffix);
    }

    bool all_field_states_exist_on_io_entity(const std::string& dbName, const stk::mesh::FieldBase* field, Ioss::GroupingEntity *ioEntity,
                                             std::vector<stk::mesh::FieldState> &missingStates, std::vector<std::string>* inputMultiStateSuffixes)
    {
        bool allStatesExist = true;
        size_t stateCount = field->number_of_states();

        std::vector<std::string>* multiStateSuffixes = stateCount > 2 ? inputMultiStateSuffixes : nullptr;

        if(nullptr != multiStateSuffixes) {
            STK_ThrowRequire(multiStateSuffixes->size() >= stateCount);
        }

        for(size_t state = 0; state < stateCount - 1; state++) {
            stk::mesh::FieldState stateIdentifier = static_cast<stk::mesh::FieldState>(state);
            if (!field_state_exists_on_io_entity(dbName, field, stateIdentifier, ioEntity, multiStateSuffixes)) {
                allStatesExist = false;
                missingStates.push_back(stateIdentifier);
            }
        }

        return allStatesExist;
    }

    void multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                         const stk::mesh::FieldBase *field,
                                         const std::vector<stk::mesh::Entity> &entityList,
                                         Ioss::GroupingEntity *ioEntity,
                                         const std::string &name,
                                         const size_t stateCount,
                                         bool ignoreMissingFields,
                                         std::vector<std::string>* inputMultiStateSuffixes)
    {
        std::vector<std::string>* multiStateSuffixes = stateCount > 2 ? inputMultiStateSuffixes : nullptr;

        if(nullptr != multiStateSuffixes) {
            STK_ThrowRequire(multiStateSuffixes->size() >= stateCount);
        }

        for(size_t state = 0; state < stateCount - 1; state++)
        {
            stk::mesh::FieldState stateIdentifier = static_cast<stk::mesh::FieldState>(state);
            bool fieldExists = field_state_exists_on_io_entity(name, field, stateIdentifier, ioEntity, multiStateSuffixes);
            if (!ignoreMissingFields) {
              const sierra::String s = multiStateSuffixes != nullptr ? (*multiStateSuffixes)[state] : std::to_string(state);
              STK_ThrowRequireMsg(fieldExists, "Field " << field->name() << s << " does not exist in input database");
            }
            if (fieldExists) {
                stk::mesh::FieldBase *statedField = field->field_state(stateIdentifier);
                std::string fieldNameWithSuffix = get_stated_field_name(name, stateIdentifier, multiStateSuffixes);
                stk::io::field_data_from_ioss(mesh, statedField, entityList, ioEntity, fieldNameWithSuffix);
            }
        }
    }

    void subsetted_multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                                   const stk::mesh::FieldBase *field,
                                                   const std::vector<stk::mesh::Entity> &entityList,
                                                   Ioss::GroupingEntity *ioEntity,
                                                   const stk::mesh::Part *stkPart,
                                                   const std::string &name,
                                                   const size_t stateCount,
                                                   bool ignoreMissingFields,
                                                   std::vector<std::string>* inputMultiStateSuffixes)
    {
        std::vector<std::string>* multiStateSuffixes = stateCount > 2 ? inputMultiStateSuffixes : nullptr;

        if(nullptr != multiStateSuffixes) {
            STK_ThrowRequire(multiStateSuffixes->size() >= stateCount);
        }

        for(size_t state = 0; state < stateCount - 1; state++)
        {
            stk::mesh::FieldState stateIdentifier = static_cast<stk::mesh::FieldState>(state);
            bool fieldExists = field_state_exists_on_io_entity(name, field, stateIdentifier, ioEntity, multiStateSuffixes);
            if (!fieldExists && !ignoreMissingFields) {
                STKIORequire(fieldExists);
            }
            if (fieldExists) {
                stk::mesh::FieldBase *statedField = field->field_state(stateIdentifier);
                std::string fieldNameWithSuffix = get_stated_field_name(name, stateIdentifier, multiStateSuffixes);
                stk::io::subsetted_field_data_from_ioss(mesh, statedField, entityList,
                                                      ioEntity, stkPart, fieldNameWithSuffix);
            }
        }
    }

    void field_data_from_ioss(const stk::mesh::BulkData& mesh,
                              const stk::mesh::FieldBase *field,
                              const std::vector<stk::mesh::Entity> &entities,
                              Ioss::GroupingEntity *ioEntity,
                              const std::string &ioFieldName)
    {
        /// \todo REFACTOR Need some additional compatibility checks between
        /// Ioss field and stk::mesh::Field; better error messages...

        if (field != nullptr && ioEntity->field_exists(ioFieldName)) {
            // Make sure the IO field type matches the STK field type.
            // By default, all IO fields are created of type 'double'
            const Ioss::Field &ioField = ioEntity->get_fieldref(ioFieldName);
            if (field->type_is<double>()) {
                internal_field_data_from_ioss<double, double>(mesh, ioField, field, entities, ioEntity);
            }
            else if (field->type_is<int>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_field_data_from_ioss<int, int>(mesh, ioField, field, entities, ioEntity);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_field_data_from_ioss<int, int64_t>(mesh, ioField, field, entities, ioEntity);
                }
            }
            else if (field->type_is<int64_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_field_data_from_ioss<int64_t, int>(mesh, ioField, field, entities, ioEntity);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_field_data_from_ioss<int64_t, int64_t>(mesh, ioField, field, entities, ioEntity);
                }
            }
            else if (field->type_is<uint32_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_field_data_from_ioss<uint32_t, int>(mesh, ioField, field, entities, ioEntity);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_field_data_from_ioss<uint32_t, int64_t>(mesh, ioField, field, entities, ioEntity);
                }
            }
            else if (field->type_is<uint64_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_field_data_from_ioss<uint64_t, int>(mesh, ioField, field, entities, ioEntity);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_field_data_from_ioss<uint64_t, int64_t>(mesh, ioField, field, entities, ioEntity);
                }
            }

        }
    }

    void subsetted_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                        const stk::mesh::FieldBase *field,
                                        const std::vector<stk::mesh::Entity> &entities,
                                        Ioss::GroupingEntity *ioEntity,
                                        const stk::mesh::Part *stkPart,
                                        const std::string &ioFieldName)
    {
        /// \todo REFACTOR Need some additional compatibility checks between
        /// Ioss field and stk::mesh::Field; better error messages...
        if (field != nullptr && ioEntity->field_exists(ioFieldName)) {
            // Make sure the IO field type matches the STK field type.
            // By default, all IO fields are created of type 'double'
            const Ioss::Field &ioField = ioEntity->get_fieldref(ioFieldName);
            if (field->type_is<double>()) {
                internal_subsetted_field_data_from_ioss<double, double>(mesh, ioField, field, entities, ioEntity, stkPart);
            }
            else if (field->type_is<int>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_subsetted_field_data_from_ioss<int, int>(mesh, ioField, field, entities, ioEntity, stkPart);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_subsetted_field_data_from_ioss<int, int64_t>(mesh, ioField, field, entities, ioEntity, stkPart);
                }
            }
            else if (field->type_is<int64_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_subsetted_field_data_from_ioss<int64_t, int>(mesh, ioField, field, entities, ioEntity, stkPart);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_subsetted_field_data_from_ioss<int64_t, int64_t>(mesh, ioField, field, entities, ioEntity, stkPart);
                }
            }
            else if (field->type_is<uint32_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_subsetted_field_data_from_ioss<uint32_t, int>(mesh, ioField, field, entities, ioEntity, stkPart);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_subsetted_field_data_from_ioss<uint32_t, int64_t>(mesh, ioField, field, entities, ioEntity, stkPart);
                }
            }
            else if (field->type_is<uint64_t>()) {
                if (ioField.get_type() == Ioss::Field::INTEGER) {
                    ioField.check_type(Ioss::Field::INTEGER);
                    internal_subsetted_field_data_from_ioss<uint64_t, int>(mesh, ioField, field, entities, ioEntity, stkPart);
                } else {
                    ioField.check_type(Ioss::Field::INT64);
                    internal_subsetted_field_data_from_ioss<uint64_t, int64_t>(mesh, ioField, field, entities, ioEntity, stkPart);
                }
            }
        }
    }

    void multistate_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                                       const stk::mesh::FieldBase *field,
                                       std::vector<stk::mesh::Entity> &entities,
                                       Ioss::GroupingEntity *ioEntity,
                                       const std::string &ioFieldName,
                                       Ioss::Field::RoleType filterRole,
                                       const size_t stateCount)
    {
      for(size_t state = 0; state < stateCount - 1; state++)
      {
        stk::mesh::FieldState stateIdentifier = static_cast<stk::mesh::FieldState>(state);
        std::string fieldNameWithSuffix = get_stated_field_name(ioFieldName, stateIdentifier);
        stk::mesh::FieldBase *statedField = field->field_state(stateIdentifier);
        stk::io::field_data_to_ioss(mesh, statedField, entities, ioEntity, fieldNameWithSuffix, filterRole);
      }
    }

    void field_data_to_ioss(const stk::mesh::BulkData& mesh,
                            const stk::mesh::FieldBase *field,
                            std::vector<stk::mesh::Entity> &entities,
                            Ioss::GroupingEntity *ioEntity,
                            const std::string &ioFieldName,
                            Ioss::Field::RoleType filterRole)
    {
      /// \todo REFACTOR Need some additional compatibility checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != nullptr && ioEntity->field_exists(ioFieldName)) {
        const Ioss::Field &ioField = ioEntity->get_fieldref(ioFieldName);
        if (ioField.get_role() == filterRole) {
          if (field->type_is<double>()) {
            internal_field_data_to_ioss<double>(mesh, ioField, field, entities, ioEntity);
          } else if (field->type_is<int>()) {
            ioField.check_type(Ioss::Field::INTEGER);
            internal_field_data_to_ioss<int>(mesh, ioField, field, entities, ioEntity);
          } else if (field->type_is<int64_t>()) {
            ioField.check_type(Ioss::Field::INT64);
            internal_field_data_to_ioss<int64_t>(mesh, ioField, field, entities, ioEntity);
          } else if (field->type_is<uint32_t>()) {
            ioField.check_type(Ioss::Field::INT32);
            internal_field_data_to_ioss<uint32_t>(mesh, ioField, field, entities, ioEntity);
          } else if (field->type_is<uint64_t>()) {
            ioField.check_type(Ioss::Field::INT64);
            internal_field_data_to_ioss<uint64_t>(mesh, ioField, field, entities, ioEntity);
          }
        }
      }
    }
    //----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // Returns true if 'entity' should be a 'part' in the analysis mesh.
    // Returns false if the application is only using a subset of the
    // database entities and this entity is not to be used. The
    // "omitted" property is set by the application during parsing or
    // pre-mesh reading time.
    bool include_entity(const Ioss::GroupingEntity *entity)
    {
      assert(entity);

      // Check whether entity has "omitted" property...
      bool omitted = (entity->property_exists("omitted")) &&
        (entity->get_property("omitted").get_int() == 1);

      return !omitted;
    }

      int64_t get_side_offset(const Ioss::ElementTopology* sideTopo,
                              const Ioss::ElementTopology* parentTopo)
      {
        int64_t sideOffset = 0;
        if ((  sideTopo != nullptr) && (  sideTopo->name() != "unknown") &&
	    (parentTopo != nullptr) && (parentTopo->name() != "unknown")) {
          int sideTopoDim = sideTopo->parametric_dimension();
          int elemTopoDim = parentTopo->parametric_dimension();
          int elemSpatDim = parentTopo->spatial_dimension();

          if (sideTopoDim + 1 < elemSpatDim && sideTopoDim < elemTopoDim) {
            sideOffset = parentTopo->number_faces();
          }
        }
        return sideOffset;
      }

      int64_t get_side_offset(const Ioss::SideBlock* sb)
      {
	if(nullptr != sb) {
	  const Ioss::ElementTopology *sideTopo   = sb->topology();
	  const Ioss::ElementTopology *parentTopo = sb->parent_element_topology();
	  return get_side_offset(sideTopo, parentTopo);
	}

	return 0;	
      }
      

    namespace {

    stk::mesh::EntityRank get_output_rank(stk::io::OutputParams& params)
    {
      return params.has_skin_mesh_selector() ? params.bulk_data().mesh_meta_data().side_rank() : stk::topology::ELEMENT_RANK;
    }

    //----------------------------------------------------------------------
    void define_node_block(stk::io::OutputParams &params, stk::mesh::Part &part)
    {
      //--------------------------------
      // Set the spatial dimension:
      mesh::MetaData & meta = mesh::MetaData::get(part);

      //We now get spatial-dim from meta.spatial_dimension() rather than getting
      //it from the coordinate-field's restriction onto the universal part.
      //This is because some codes (sierra framework) don't put the coordinate
      //field on the universal part. (framework puts it on active and inactive parts)
      const int spatialDim = meta.spatial_dimension();
      stk::mesh::EntityRank rank = get_output_rank(params);
      //--------------------------------
      // Create the special universal node block:
      mesh::Selector sharedSelector = params.has_shared_selector() ? *(params.get_shared_selector())
                                                                   : meta.globally_shared_part();

      mesh::Selector allSelector = meta.globally_shared_part() | meta.locally_owned_part();
      if (params.get_subset_selector(    )) allSelector &= *params.get_subset_selector();
      if (params.get_output_selector(rank)) allSelector &= *params.get_output_selector(rank);

      mesh::Selector ownSelector = meta.locally_owned_part();
      if (params.get_subset_selector(    )) ownSelector &= *params.get_subset_selector();
      if (params.get_output_selector(rank)) ownSelector &= *params.get_output_selector(rank);

      int64_t allNodes = count_selected_nodes(params, allSelector);
      int64_t ownNodes = count_selected_nodes(params, ownSelector);

      const std::string name("nodeblock_1");

      Ioss::NodeBlock * nb = params.io_region().get_node_block(name);
      if(nb == nullptr)
      {
          nb = new Ioss::NodeBlock(params.io_region().get_database(),
                                   name, allNodes, spatialDim);
          params.io_region().add( nb );
      }

      delete_selector_property(nb);
      mesh::Selector *nodeSelect = new mesh::Selector(allSelector);
      nb->property_add(Ioss::Property(s_internalSelectorName, nodeSelect));
      nb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

      // Add locally-owned property...
      nb->property_add(Ioss::Property("locally_owned_count", ownNodes));
      // Add the attribute fields.
      ioss_add_fields(part, part_primary_entity_rank(part), nb, Ioss::Field::ATTRIBUTE);
    }

    void define_node_set(stk::io::OutputParams &params,
                         stk::mesh::Part &part,
                         const std::string &name,
                         const bool isDerivedNodeset = false)
    {
      mesh::EntityRank rank = get_output_rank(params);
      mesh::MetaData & meta = mesh::MetaData::get(part);
      Ioss::Region & ioRegion = params.io_region();

      mesh::Selector sharedSelector = params.has_shared_selector() ? *(params.get_shared_selector())
                                                                   : meta.globally_shared_part();

      mesh::Selector allSelector = (meta.globally_shared_part() | meta.locally_owned_part()) & part;
      if (params.get_subset_selector(    )) allSelector &= *params.get_subset_selector();
      if (params.get_output_selector(rank)) allSelector &= *params.get_output_selector(rank);

      mesh::Selector ownSelector = meta.locally_owned_part() & part;
      if (params.get_subset_selector(    )) ownSelector &= *params.get_subset_selector();
      if (params.get_output_selector(rank)) ownSelector &= *params.get_output_selector(rank);

      int64_t allNodes = count_selected_nodes(params, allSelector);
      int64_t ownNodes = count_selected_nodes(params, ownSelector);

      Ioss::NodeSet *ns = ioRegion.get_nodeset(name);
      if(ns == nullptr)
      {
          ns = new Ioss::NodeSet( ioRegion.get_database(), name, allNodes);
          ioRegion.add(ns);

          bool use_generic_canonical_name = ioRegion.get_database()->get_use_generic_canonical_name();
          if(use_generic_canonical_name) {
            add_canonical_name_property(ns, part);
          }
      }

      ns->property_add(Ioss::Property("locally_owned_count", ownNodes));

      delete_selector_property(ns);
      mesh::Selector *select = new mesh::Selector(allSelector);
      ns->property_add(Ioss::Property(s_internalSelectorName, select));
      ns->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

      if(!isDerivedNodeset) {
          if(has_original_part_id(part)) {
              int64_t id = get_original_part_id(part);
              ns->property_add(Ioss::Property("id", id));
          } else if (params.get_use_part_id_for_output() && (part.id() != stk::mesh::Part::INVALID_ID)) {
              ns->property_add(Ioss::Property("id", part.id()));
          }
      }

      // Add the attribute fields.
      ioss_add_fields(part, stk::topology::NODE_RANK, ns, Ioss::Field::ATTRIBUTE);
    }

      std::tuple<std::string, const Ioss::ElementTopology *, stk::topology>
      get_touching_element_block_topology_from_side_block_by_tokenization(stk::io::OutputParams &params, const stk::mesh::Part& part)
      {
        const stk::mesh::BulkData &bulk = params.bulk_data();

        std::string elementTopoName = "unknown";
        const Ioss::ElementTopology *elementTopo = nullptr;
        stk::topology invalidTopology = stk::topology::INVALID_TOPOLOGY;
        stk::topology stkElementTopology = invalidTopology;

        // Try to decode from part name...not always consistent since it is client dependent
        std::vector<std::string> tokens;
        std::string sideBlockName = getPartName(part);
        stk::util::tokenize(sideBlockName, "_", tokens);

        if (tokens.size() >= 3) {
          // If the sideset has a "canonical" name as in "surface_{id}",
          // Then the sideblock name will be of the form:
          //  * "surface_eltopo_sidetopo_id" or
          //  * "surface_block_id_sidetopo_id"
          // If the sideset does *not* have a canonical name, then
          // the sideblock name will be of the form:
          //  * "{sideset_name}_eltopo_sidetopo" or
          //  * "{sideset_name}_block_id_sidetopo"

          // Check the last token and see if it is an integer...
          const bool allDigits = tokens.back().find_first_not_of("0123456789") == std::string::npos;
          if (allDigits) {
            elementTopo = Ioss::ElementTopology::factory(tokens[1], true);
          } else {
            const std::string& elemTopoToken = tokens[tokens.size()-2];
            const bool subAllDigits = elemTopoToken.find_first_not_of("0123456789") == std::string::npos;
            if (!subAllDigits) {
              elementTopo = Ioss::ElementTopology::factory(elemTopoToken, true);
            }
          }

	  const stk::mesh::MetaData& meta = params.bulk_data().mesh_meta_data();
	  std::vector<const stk::mesh::Part*> touchingParts = meta.get_blocks_touching_surface(&part);
	
	  if(touchingParts.size() == 1u) {
	    stk::topology stkTouchingTopology = touchingParts[0]->topology();
	    std::string touchingTopoName = map_stk_topology_to_ioss(stkTouchingTopology);
	    const Ioss::ElementTopology *touchingTopo = Ioss::ElementTopology::factory(touchingTopoName, true);
	    if(touchingTopo != elementTopo) {
	      elementTopo = nullptr;	      
	    }
	  }
	  
          if (elementTopo != nullptr) {
            elementTopoName = elementTopo->name();
            stkElementTopology = map_ioss_topology_to_stk(elementTopo, bulk.mesh_meta_data().spatial_dimension());
          }
        }

        return std::make_tuple(elementTopoName, elementTopo, stkElementTopology);
      }

      std::tuple<std::string, const Ioss::ElementTopology *, stk::topology>
      get_touching_element_block_topology_from_side_block(stk::io::OutputParams &params, const stk::mesh::Part& part)
      {
        return get_touching_element_block_topology_from_side_block_by_tokenization(params, part);
      }

      void define_side_block(stk::io::OutputParams &params,
                             stk::mesh::Selector selector,
                             stk::mesh::Part &part,
                             Ioss::SideSet *sset,
                             int spatialDimension,
                             bool createNodeset)
      {
        stk::mesh::EntityRank type = part.primary_entity_rank();
        const stk::mesh::EntityRank siderank = stk::mesh::MetaData::get(part).side_rank();
        const stk::mesh::EntityRank edgerank = stk::topology::EDGE_RANK;
        STKIORequire(type == siderank || type == edgerank);

        stk::topology sideTopology = part.topology();
        std::string ioTopo = map_stk_topology_to_ioss(sideTopology);
        std::string elementTopoName = "unknown";
        const Ioss::ElementTopology *elementTopo = nullptr;
        stk::topology stkElementTopology = stk::topology::INVALID_TOPOLOGY;

        std::tie(elementTopoName, elementTopo, stkElementTopology) =
            get_touching_element_block_topology_from_side_block(params, part);

        const stk::mesh::BulkData &bulk = params.bulk_data();

        const stk::mesh::Part *parentElementBlock = get_parent_element_block(bulk, params.io_region(), part.name());

        if (elementTopo != nullptr && !elementTopo->is_element()) {
          if(parentElementBlock != nullptr) {
            stkElementTopology = parentElementBlock->topology();
            elementTopoName = map_stk_topology_to_ioss(stkElementTopology);
          } else {
            elementTopoName = "unknown";
            elementTopo = nullptr;
            stkElementTopology = stk::topology::INVALID_TOPOLOGY;
          }
        }

        const Ioss::ElementTopology* iossSideTopo = Ioss::ElementTopology::factory(ioTopo, true);
        int64_t sideOffset = get_side_offset(iossSideTopo, elementTopo);
        size_t sideCount = get_number_sides_in_sideset(params, part, stkElementTopology, parentElementBlock, sideOffset);

        std::string name = getPartName(part);
        Ioss::SideBlock *sideBlock = sset->get_side_block(name);
        if(sideBlock == nullptr)
        {
	  if("unknown" == ioTopo) {
	    // Special case of heterogenous sideset encapsulated in one side block
	    // There is no side topology so we cannot specify an element topology
	    // due to internal IOSS code that gives wrong SideBlock offset
	    elementTopoName = "unknown";
	  }
	  
            sideBlock = new Ioss::SideBlock(sset->get_database(), name, ioTopo, elementTopoName, sideCount);
            sset->add(sideBlock);
        }

        const mesh::FieldBase *df = get_distribution_factor_field(part);
        if (df != nullptr) {
          int nodesPerSide = sideTopology.num_nodes();
          std::string storageType = "Real[";
          storageType += sierra::to_string(nodesPerSide);
          storageType += "]";
          sideBlock->field_add(Ioss::Field(s_distributionFactors, Ioss::Field::REAL, storageType,
                                            Ioss::Field::MESH, sideCount));
        }

        selector &= bulk.mesh_meta_data().locally_owned_part();
        delete_selector_property(sideBlock);
        mesh::Selector *select = new mesh::Selector(selector);
        sideBlock->property_add(Ioss::Property(s_internalSelectorName, select));
        sideBlock->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), sideBlock, Ioss::Field::ATTRIBUTE);

        if(createNodeset) {
            std::string nodes_name = getPartName(part) + s_entityNodesSuffix;
            bool isDerivedNodeset = true;
            define_node_set(params, part, nodes_name, isDerivedNodeset);
        }
      }

      bool should_create_nodeset_from_sideset(stk::mesh::Part &part,
                                              bool useNodesetForNodalFields,
                                              bool checkFieldExistence)
      {
          STKIORequire(part.primary_entity_rank() == stk::topology::FACE_RANK || stk::topology::EDGE_RANK);

          bool createNodesets = false;

          if (useNodesetForNodalFields) {
              if(checkFieldExistence) {
                  bool lowerRankFields = will_output_lower_rank_fields(part, stk::topology::NODE_RANK);

                  if (!lowerRankFields) {
                      // See if lower rank fields are defined on sideblock parts of this sideset...
                      const stk::mesh::PartVector &blocks = part.subsets();
                      for (size_t j = 0; j < blocks.size() && !lowerRankFields; j++) {
                          mesh::Part & side_block_part = *blocks[j];
                          lowerRankFields |= will_output_lower_rank_fields(side_block_part, stk::topology::NODE_RANK);
                      }
                  }
                  if (lowerRankFields) {
                      createNodesets = true;
                  }
              } else {
                  createNodesets = true;
              }
          }

          if(has_derived_nodeset_attribute(part)) {
              createNodesets = get_derived_nodeset_attribute(part);
          }

          return createNodesets;
      }

      void define_side_blocks(stk::io::OutputParams &params,
                              stk::mesh::Part &part,
                              Ioss::SideSet *sset,
                              stk::mesh::EntityRank type,
                              int spatialDimension)
      {
        STKIORequire(type == stk::topology::FACE_RANK || stk::topology::EDGE_RANK);

        bool createNodesets = should_create_nodeset_from_sideset(part,
                                                                 params.get_use_nodeset_for_sideset_node_fields(),
                                                                 params.check_field_existence_when_creating_nodesets());

        stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
        const stk::mesh::PartVector &blocks = part.subsets();
        if (blocks.size() > 0) {
          for (size_t j = 0; j < blocks.size(); j++) {
            mesh::Part & sideBlockPart = *blocks[j];
            mesh::Selector selector = sideBlockPart;
            if (params.get_subset_selector(    )) selector &= *params.get_subset_selector();
            if (params.get_output_selector(rank)) selector &= *params.get_output_selector(rank);
            define_side_block(params, selector,
                              sideBlockPart, sset, spatialDimension,
                              createNodesets);
          }
        } else {
          mesh::Selector selector = part;
          if (params.get_subset_selector(    )) selector &= *params.get_subset_selector();
          if (params.get_output_selector(rank)) selector &= *params.get_output_selector(rank);
          define_side_block(params, selector,
                            part, sset, spatialDimension,
                            createNodesets);
        }
      }

      void set_attribute_field_order(const stk::mesh::BulkData& bulk,
                                     const std::vector<std::vector<int> >& attributeOrdering,
                                     stk::mesh::Part& part,
                                     Ioss::ElementBlock* eb)
      {
          int attributeIndex = 1;
          const stk::mesh::FieldVector& allFields = bulk.mesh_meta_data().get_fields();
          if(part.mesh_meta_data_ordinal() < attributeOrdering.size())
          {
              for(int fieldOrd : attributeOrdering[part.mesh_meta_data_ordinal()])
              {
                  stk::mesh::FieldBase* stkField = allFields[fieldOrd];
                  const Ioss::Field& iossField = eb->get_fieldref(stkField->name());
                  iossField.set_index(attributeIndex);
                  attributeIndex += iossField.raw_storage()->component_count();
              }
          }
      }

      void set_id_property(stk::io::OutputParams &params,
                           const stk::mesh::Part &part,
                           Ioss::GroupingEntity* groupingEntity)
      {
        if (has_original_part_id(part)) {
          int64_t id = get_original_part_id(part);
          groupingEntity->property_update("id", id);
        }
        else if (params.get_use_part_id_for_output() && (part.id() != stk::mesh::Part::INVALID_ID)) {
          groupingEntity->property_update("id", part.id());
        }
      }

      bool assembly_has_valid_io_leaf_part(stk::io::OutputParams &params,
                                           const stk::mesh::Part& assemblyPart)
      {
        const stk::mesh::MetaData & meta = mesh::MetaData::get(assemblyPart);
        stk::mesh::PartVector leafParts = get_unique_leaf_parts(meta, assemblyPart.name());
        for (stk::mesh::Part* leafPart : leafParts) {
          if (is_in_subsets_of_parts(*leafPart, leafParts)) {continue;}
          if (is_valid_for_output(params, *leafPart)) {
            return true;
          }
        }
        return false;
      }

      void define_assembly(stk::io::OutputParams &params,
                           const stk::mesh::Part &part)
      {
        if (!assembly_has_valid_io_leaf_part(params, part)) {
          return;
        }
        Ioss::Region &ioRegion = params.io_region();

        std::string name = getPartName(part);
        Ioss::Assembly *assembly = ioRegion.get_assembly(name);
        if (assembly == nullptr) {
          assembly = new Ioss::Assembly(ioRegion.get_database(), name);
          set_id_property(params, part, assembly);
          ioRegion.add(assembly);
        }
      }

      bool is_valid_assembly_member_type(const Ioss::Assembly *assem, const Ioss::GroupingEntity* member)
      {
        if(nullptr == member) return false;

        if((member->type() != Ioss::ELEMENTBLOCK) && (member->type() != Ioss::SIDESET) &&
           (member->type() != Ioss::NODESET)      && (member->type() != Ioss::ASSEMBLY))   {
          std::string filename = assem->get_database()->get_filename();
          stk::RuntimeWarningP0() << "The entity type of '" << member->name() << "' (" << member->type_string() <<
                                     ") is not a valid assembly member type for "
                                     "assembly '" << assem->name() << "' (" << assem->contains_string() <<
                                     ").\n\t In the database file '" << filename << "'.\n";
          return false;
        }

        return true;
      }

      bool is_empty_element_block(stk::io::OutputParams &params, const stk::mesh::Part* leafPart)
      {
        bool isEmptyElementBlock = false;
        const std::unordered_map<unsigned, size_t>& blockSizes = params.get_block_sizes();

        if(leafPart != nullptr && is_part_element_block_io_part(*leafPart)) {
          auto iter = blockSizes.find(leafPart->mesh_meta_data_ordinal());
          STK_ThrowRequireMsg(iter != blockSizes.end(), "Could not find element block in block size list: " << leafPart->name());
          isEmptyElementBlock = (iter->second == 0);
        }

        return isEmptyElementBlock;
      }

      bool can_add_to_assembly(stk::io::OutputParams &params, const Ioss::Assembly *assembly,
                               const Ioss::GroupingEntity* leafEntity, const stk::mesh::Part* leafPart)
      {
        bool isNotCurrentMember = (leafEntity != nullptr) && (assembly->get_member(leafEntity->name()) == nullptr);
        bool isValidMemberType = is_valid_assembly_member_type(assembly, leafEntity);
        bool isEmptyElementBlock = false;

        bool filterEmptyBlocks = params.get_filter_empty_entity_blocks() ||
                                 params.get_filter_empty_assembly_entity_blocks();

        if(filterEmptyBlocks) {
          isEmptyElementBlock = is_empty_element_block(params, leafPart);
        }

        bool isValid = isNotCurrentMember && isValidMemberType && !isEmptyElementBlock;

        return isValid;
      }

      void define_assembly_hierarchy(stk::io::OutputParams &params,
                                     const stk::mesh::Part &part)
      {
        if (!assembly_has_valid_io_leaf_part(params, part)) {
          return;
        }
        const stk::mesh::MetaData & meta = mesh::MetaData::get(part);
        Ioss::Region &ioRegion = params.io_region();

        std::string name = getPartName(part);

        Ioss::Assembly *assembly = ioRegion.get_assembly(name);
        STK_ThrowRequireMsg(assembly != nullptr, "Failed to find assembly "<<name);

        if (has_sub_assemblies(meta, part.name())) {
          std::vector<std::string> subAssemblyNames = get_sub_assembly_names(meta, part.name());
          for(const std::string& subAssemblyName : subAssemblyNames) {
            const stk::mesh::Part* subAssemblyPart = meta.get_part(subAssemblyName);

            std::string iossSubAssemblyName = getPartName(*subAssemblyPart);
            const Ioss::Assembly* subAssembly = ioRegion.get_assembly(iossSubAssemblyName);
            if(subAssembly == nullptr) {
              stk::RuntimeWarning() << "Failed to find subAssembly "<<iossSubAssemblyName;
              continue;
            }
            if(assembly->get_member(subAssembly->name())==nullptr) {
                assembly->add(subAssembly);
            }
          }
        }
        else {
          stk::mesh::PartVector leafParts = get_unique_leaf_parts(meta, part.name());
          for(stk::mesh::Part* leafPart : leafParts) {
            if(is_in_subsets_of_parts(*leafPart, leafParts)) {continue;}
            std::string iossLeafPartName = getPartName(*leafPart);
            const Ioss::GroupingEntity* leafEntity = ioRegion.get_entity(iossLeafPartName);
            if (leafEntity == nullptr) {
              stk::RuntimeWarning() << "Failed to find ioss entity: '" << iossLeafPartName << "' in assembly: '" << name
                                    << "'";
            }
            if (can_add_to_assembly(params, assembly, leafEntity, leafPart)) {
              assembly->add(leafEntity);
            }
          }
        }

        if (assembly->member_count() == 0) {
          ioRegion.remove(assembly);
          delete assembly;
        }
      }

      void define_face_block(stk::io::OutputParams &params,
                             stk::mesh::Part &part)
      {
        mesh::MetaData & meta = mesh::MetaData::get(part);
        const stk::mesh::BulkData &bulk = params.bulk_data();
        Ioss::Region &ioRegion = params.io_region();

        stk::topology topo = part.topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR when defining output for region '"<<ioRegion.name()<<"': Part " << part.name()
              << " returned INVALID from get_topology(). Please contact sierra-help@sandia.gov";
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank rank = stk::topology::FACE_RANK;

        mesh::Selector selector = meta.locally_owned_part() & part;
        if (params.get_subset_selector()) selector &= *params.get_subset_selector();
        if (params.get_output_selector(rank)) selector &= *params.get_output_selector(rank);

        std::string topologyName = map_stk_topology_to_ioss(topo);
        const size_t numFaces = stk::mesh::count_entities(bulk, rank, selector);

        // Defer the counting of attributes until after we define the
        // element block so we can count them as we add them as fields to
        // the element block
        std::string name = getPartName(part);
        Ioss::FaceBlock *fb = ioRegion.get_face_block(name);
        if(fb == nullptr)
        {
            fb = new Ioss::FaceBlock(ioRegion.get_database() ,
                                     name,
                                     topologyName,
                                     numFaces);
            ioRegion.add(fb);

            bool useGenericCanonicalName = ioRegion.get_database()->get_use_generic_canonical_name();
            if(useGenericCanonicalName) {
              add_canonical_name_property(fb, part);
            }

            bool useOriginalTopology = has_original_topology_type(part);
            if(useOriginalTopology) {
                add_original_topology_property(fb, part);
            }
        }

        set_id_property(params, part, fb);

        delete_selector_property(fb);
        mesh::Selector *select = new mesh::Selector(selector);
        fb->property_add(Ioss::Property(s_internalSelectorName, select));
        fb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), fb, Ioss::Field::ATTRIBUTE);
      }

      void define_edge_block(stk::io::OutputParams &params,
                             stk::mesh::Part &part)
      {
        mesh::MetaData & meta = mesh::MetaData::get(part);
        const stk::mesh::BulkData &bulk = params.bulk_data();
        Ioss::Region &ioRegion = params.io_region();

        stk::topology topo = part.topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR when defining output for region '"<<ioRegion.name()<<"': Part " << part.name()
              << " returned INVALID from get_topology(). Please contact sierra-help@sandia.gov";
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank rank = stk::topology::EDGE_RANK;

        mesh::Selector selector = meta.locally_owned_part() & part;
        if (params.get_subset_selector()) selector &= *params.get_subset_selector();
        if (params.get_output_selector(rank)) selector &= *params.get_output_selector(rank);

        std::string topologyName = map_stk_topology_to_ioss(topo);
        const size_t numEdges = stk::mesh::count_entities(bulk, rank, selector);

        // Defer the counting of attributes until after we define the
        // element block so we can count them as we add them as fields to
        // the element block
        std::string name = getPartName(part);
        Ioss::EdgeBlock *eb = ioRegion.get_edge_block(name);
        if(eb == nullptr)
        {
            eb = new Ioss::EdgeBlock(ioRegion.get_database() ,
                                     name,
                                     topologyName,
                                     numEdges);
            ioRegion.add(eb);

            bool use_generic_canonical_name = ioRegion.get_database()->get_use_generic_canonical_name();
            if(use_generic_canonical_name) {
              add_canonical_name_property(eb, part);
            }

            bool use_original_topology = has_original_topology_type(part);
            if(use_original_topology) {
                add_original_topology_property(eb, part);
            }
        }

        set_id_property(params, part, eb);

        delete_selector_property(eb);
        mesh::Selector *select = new mesh::Selector(selector);
        eb->property_add(Ioss::Property(s_internalSelectorName, select));
        eb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), eb, Ioss::Field::ATTRIBUTE);
      }

      void define_element_block(stk::io::OutputParams &params,
                                stk::mesh::Part &part,
                                const std::vector<std::vector<int>> &attributeOrdering,
                                bool orderBlocksByCreationOrder)
      {
        mesh::MetaData & meta = mesh::MetaData::get(part);
        const stk::mesh::BulkData &bulk = params.bulk_data();
        Ioss::Region &ioRegion = params.io_region();

        stk::mesh::EntityRank rank = get_output_rank(params);

        mesh::Selector selector = impl::internal_build_selector(params.get_subset_selector(),
                                                                params.get_output_selector(rank),
                                                                nullptr, part, false);

        const size_t numElems = stk::mesh::count_entities(bulk, rank, selector);

        stk::topology topo = part.topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR when defining output for region '"<<ioRegion.name()<<"': Part " << part.name()
              << " returned INVALID from get_topology(). Please contact sierra-help@sandia.gov";
          std::cerr << msg.str()<<std::endl;
        }

        std::string topologyName = map_stk_topology_to_ioss(topo);

        if( params.has_skin_mesh_selector() ) {
          // FIX THIS...Assumes homogenous faces... <- From Frio (tookusa)
          topologyName = map_stk_topology_to_ioss(topo.face_topology(0));

          if (meta.get_topology(part) == stk::topology::PARTICLE){
            stk::RuntimeWarning() << "When using the skin output option, element topologies of type PARTICLE are not supported."
                                  << "Blocks that use this topology will not be output.";
            return;
          }
        }

        // Defer the counting of attributes until after we define the
        // element block so we can count them as we add them as fields to
        // the element block
        std::string name = getPartName(part);
        Ioss::ElementBlock *eb = ioRegion.get_element_block(name);
        if(eb == nullptr)
        {
            eb = new Ioss::ElementBlock(ioRegion.get_database() ,
                                        name,
                                        topologyName,
                                        numElems);
            ioRegion.add(eb);

            bool useGenericCanonicalName = ioRegion.get_database()->get_use_generic_canonical_name();
            if(useGenericCanonicalName) {
              add_canonical_name_property(eb, part);
            }

            bool useOriginalTopology = has_original_topology_type(part);
            if(useOriginalTopology && !params.has_skin_mesh_selector()) {
                add_original_topology_property(eb, part);
            }
        }

        if (orderBlocksByCreationOrder)
        {
            int ordinal = part.mesh_meta_data_ordinal();
            eb->property_update("original_block_order", ordinal);
        }

        set_id_property(params, part, eb);

        delete_selector_property(eb);
        mesh::Selector *select = new mesh::Selector(selector);
        eb->property_add(Ioss::Property(s_internalSelectorName, select));
        eb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), eb, Ioss::Field::ATTRIBUTE);
        process_element_attributes_for_define(params, part);

        set_attribute_field_order(bulk, attributeOrdering, part, eb);

        // Check whether there are any transient fields defined on the nodes of this elementblock
        // that are to be output.  If so, create a nodeset named "part.name()"+block_nodes_suffix
        // and output the fields on that nodeset...
        if (params.get_use_nodeset_for_block_node_fields() &&
            will_output_lower_rank_fields(part, stk::topology::NODE_RANK)) {
          std::string nodesName = getPartName(part) + s_entityNodesSuffix;
          bool isDerivedNodeset = true;
          define_node_set(params, part, nodesName, isDerivedNodeset);
        }
      }

      void define_communication_maps(stk::io::OutputParams &params)
      {
        const mesh::BulkData & bulk = params.bulk_data();
        Ioss::Region & ioRegion = params.io_region();
        mesh::EntityRank rank = get_output_rank(params);
        const stk::mesh::Selector *subsetSelector = params.get_subset_selector();
        const stk::mesh::Selector *outputSelector = params.get_output_selector(rank);

        if (bulk.parallel_size() > 1) {
          const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
          const std::string csName("node_symm_comm_spec");

          mesh::Selector selector = meta.globally_shared_part();
          if (subsetSelector) selector &= *subsetSelector;
          if (outputSelector) selector &= *outputSelector;

          std::vector<mesh::Entity> entities;
          get_selected_nodes(params, selector, entities);

          std::vector<int> sharingProcs;
          size_t size = 0;
          for (size_t i=0; i < entities.size(); i++) {
            bulk.comm_shared_procs(bulk.entity_key(entities[i]), sharingProcs);
            size+=sharingProcs.size();
          }

          Ioss::DatabaseIO *dbo = ioRegion.get_database();
          Ioss::CommSet *ioCs = new Ioss::CommSet(dbo, csName, "node", size);
          ioRegion.add(ioCs);

          delete_selector_property(ioCs);
          mesh::Selector *select = new mesh::Selector(selector);
          ioCs->property_add(Ioss::Property(s_internalSelectorName, select));

          // Update global node and element count...
          if (!ioRegion.property_exists("global_node_count") || !ioRegion.property_exists("global_element_count")) {
            std::vector<size_t> entityCounts;
            stk::mesh::comm_mesh_counts(bulk, entityCounts);

            ioRegion.property_add(Ioss::Property("global_node_count",    static_cast<int64_t>(entityCounts[stk::topology::NODE_RANK])));
            ioRegion.property_add(Ioss::Property("global_element_count", static_cast<int64_t>(entityCounts[stk::topology::ELEMENT_RANK])));
          }
        }
      }

      void define_side_set(stk::io::OutputParams &params, stk::mesh::Part &part)
      {
        const stk::mesh::EntityRank sideRank = mesh::MetaData::get(part).side_rank();

        bool createSideset = ! params.has_skin_mesh_selector();
        if (part.subsets().empty()) {
          // Only define a sideset for this part if its superset part is
          // not a side-containing part..  (i.e., this part is not a subset part
          // in a surface...)
          const stk::mesh::PartVector &supersets = part.supersets();
          for (size_t i=0; i < supersets.size(); i++) {
            if (is_part_surface_io_part(*supersets[i])) {
              createSideset = false;
              break;
            }
          }
        }

        if (createSideset) {
          std::string name = getPartName(part);
          Ioss::Region & ioRegion = params.io_region();
          Ioss::SideSet *ss = ioRegion.get_sideset(name);
          if(ss == nullptr)
          {
              ss = new Ioss::SideSet(ioRegion.get_database(), name);
              ioRegion.add(ss);

              bool useGenericCanonicalName = ioRegion.get_database()->get_use_generic_canonical_name();
              if(useGenericCanonicalName) {
                add_canonical_name_property(ss, part);
              }
          }

          if(has_original_part_id(part)) {
              int64_t id = get_original_part_id(part);
              ss->property_add(Ioss::Property("id", id));
          } else if (params.get_use_part_id_for_output() && (part.id() != stk::mesh::Part::INVALID_ID)) {
              ss->property_add(Ioss::Property("id", part.id()));
          }

          int spatialDim = ioRegion.get_property("spatial_dimension").get_int();
          define_side_blocks(params, part, ss, sideRank, spatialDim);
        }
      }

    } // namespace <blank>

    void set_element_block_order(const mesh::PartVector *parts, Ioss::Region & ioRegion)
    {
        int64_t offset=0;

        for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {
            mesh::Part * const part = *i ;

            if (is_part_io_part(*part) && (part->primary_entity_rank() == stk::topology::ELEMENT_RANK)) {
                if(has_original_block_order(*part)) {
                    int64_t order = get_original_block_order(*part);
                    Ioss::GroupingEntity *elementBlock = ioRegion.get_entity(getPartName(*part));
                    if (elementBlock) {
                        elementBlock->property_update("original_block_order", order);
                        offset = std::max(offset, order);
                    }
                }
            }
        }

        offset += 1;

        for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {
            mesh::Part * const part = *i ;

            if (is_part_io_part(*part) && (part->primary_entity_rank() == stk::topology::ELEMENT_RANK)) {
                Ioss::GroupingEntity *elementBlock = ioRegion.get_entity(getPartName(*part));
                if (elementBlock) {
                    if (!elementBlock->property_exists("original_block_order")) {
                        elementBlock->property_add(Ioss::Property("original_block_order", offset));
                        ++offset;
                    }
                }
            }
        }
    }

    struct part_compare_by_name {
      bool operator() (const stk::mesh::Part * const i, const stk::mesh::Part * const j) {
          if(i == j) return false;
          if(nullptr == i) return true;
          if(nullptr == j) return false;
          return (i->name() < j->name());
      }
    };

    bool has_io_subset_but_no_non_assembly_io_superset(const stk::mesh::Part& part)
    {
      for(const stk::mesh::Part* superset : part.supersets()) {
        if (stk::io::is_part_io_part(*superset) &&
           !stk::io::is_part_assembly_io_part(*superset)) {
          return false;
        }
      }
      for(const stk::mesh::Part* subset : part.subsets()) {
        if (stk::io::is_part_io_part(*subset)) {
          return true;
        }
      }
      return false;
    }

    bool is_edge_rank_sideset_part(const stk::mesh::Part& part)
    {
      if (!stk::io::is_part_edge_block_io_part(part)) {
        const unsigned spatialDim = part.mesh_meta_data().spatial_dimension();
        if (part.primary_entity_rank() == stk::topology::EDGE_RANK) {
          if (spatialDim == 2) {
            return true;
          }
          if (spatialDim == 3 && has_io_subset_but_no_non_assembly_io_superset(part)) {
            return true;
          }
        }
      }
      return false;
    }

    bool is_face_rank_sideset_part(const stk::mesh::Part& part)
    {
      return part.primary_entity_rank() == stk::topology::FACE_RANK
          && part.mesh_meta_data().spatial_dimension() == 3
          && !is_part_face_block_io_part(part);
    }

    bool is_sideset_part(const stk::mesh::Part& part)
    {
      return is_face_rank_sideset_part(part) || is_edge_rank_sideset_part(part);
    }

    void define_output_db_within_state_define(stk::io::OutputParams &params,
                                              const std::vector<std::vector<int>> &attributeOrdering,
                                              const Ioss::Region *inputRegion = nullptr)
    {
       Ioss::Region & ioRegion = params.io_region();
       const mesh::BulkData &bulkData = params.bulk_data();
       const bool sortStkPartsByName = params.get_sort_stk_parts_by_name();

       const mesh::MetaData & metaData = bulkData.mesh_meta_data();
       define_node_block(params, metaData.universal_part());

       // All parts of the meta data:
       const mesh::PartVector *parts = nullptr;
       mesh::PartVector allPartsSorted;

       const mesh::PartVector & allParts = metaData.get_parts();
       // sort parts so they go out the same on all processors (srk: this was induced by streaming refine)
       if (sortStkPartsByName) {
         allPartsSorted = allParts;
         std::sort(allPartsSorted.begin(), allPartsSorted.end(), part_compare_by_name());
         parts = &allPartsSorted;
       } else {
         parts = &allParts;
       }

       const bool orderBlocksByCreationOrder = (inputRegion == nullptr) && !sortStkPartsByName;
       const int spatialDim = metaData.spatial_dimension();

       for (stk::mesh::Part* const part : *parts) {
         const stk::mesh::EntityRank rank = part->primary_entity_rank();

         if (is_part_io_part(*part)) {
           bool isValidForOutput = is_valid_for_output(params, *part);

           if (is_part_assembly_io_part(*part)) {
             define_assembly(params, *part);
           }
           else if (rank == mesh::InvalidEntityRank) {
             continue;
           }
           else if ((rank == stk::topology::NODE_RANK) && isValidForOutput) {
             define_node_set(params, *part, getPartName(*part));
           }
           else if ((rank == stk::topology::ELEMENT_RANK) && isValidForOutput) {
             define_element_block(params, *part, attributeOrdering, orderBlocksByCreationOrder);
           }
           else if (is_part_face_block_io_part(*part)) {
             define_face_block(params, *part);
           }
           else if (is_part_edge_block_io_part(*part)) {
             define_edge_block(params, *part);
           } else if (is_sideset_part(*part) && isValidForOutput) {
             define_side_set(params, *part);
           } else if ((rank == stk::topology::EDGE_RANK) && spatialDim == 3 && params.get_enable_edge_io()) {
             define_edge_block(params, *part);
           }
         }
       }

       for (const stk::mesh::Part* part : *parts) {
         if (is_part_assembly_io_part(*part)) {
           define_assembly_hierarchy(params, *part);
         }
       }

       define_communication_maps(params);

       if (inputRegion != nullptr)
         ioRegion.synchronize_id_and_name(inputRegion, true);

       set_element_block_order(parts, ioRegion);
    }

    void define_output_db(stk::io::OutputParams &params,
                          const std::vector<std::vector<int>> &attributeOrdering,
                          const Ioss::Region *inputRegion)
    {
      params.io_region().begin_mode( Ioss::STATE_DEFINE_MODEL );
      define_output_db_within_state_define(params, attributeOrdering, inputRegion);
      params.io_region().end_mode( Ioss::STATE_DEFINE_MODEL );
    }

    //----------------------------------------------------------------------

    namespace {
      template <typename INT>
      void write_side_data_to_ioss( stk::io::OutputParams &params,
                                    Ioss::GroupingEntity & io ,
                                    mesh::Part * const part ,
                                    const Ioss::ElementTopology *elementTopology)
      {
        std::vector<INT> elemSideIds;
        stk::mesh::EntityVector sides;

        fill_data_for_side_block(params, io, part, elementTopology, elemSideIds, sides);
        size_t numSides = sides.size();

        const size_t numSideWritten = io.put_field_data("element_side",elemSideIds);

        if ( numSides != numSideWritten ) {
          std::ostringstream msg ;

          msg << "stk::io::write_side_data_to_ioss FAILED for " ;
          msg << io.name();
          msg << " in Ioss::GroupingEntity::put_field_data:" ;
          msg << " numSides = " << numSides ;
          msg << " , num_side_written = " << numSideWritten ;
          throw std::runtime_error( msg.str() );
        }

        const mesh::FieldBase *df = get_distribution_factor_field(*part);
        if (df != nullptr) {
          field_data_to_ioss(params.bulk_data(), df, sides, &io, s_distributionFactors, Ioss::Field::MESH);
        }

        const mesh::MetaData & metaData = mesh::MetaData::get(*part);

        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != nullptr && *role == Ioss::Field::ATTRIBUTE) {
            stk::io::field_data_to_ioss(params.bulk_data(), f, sides, &io, f->name(), Ioss::Field::ATTRIBUTE);
          }
        }
      }

      //----------------------------------------------------------------------
      template <typename INT>
      void output_node_block(stk::io::OutputParams &params,
                             Ioss::NodeBlock &nb,
                             stk::mesh::Part &part)
      {
        //----------------------------------
        // Exactly one node block to obtain the nodal coordinates and ids:
        // Note that the "ids" field of the nodes needs to be written
        // before any other bulk data that uses node ids since it sets up
        // the global->local mapping of nodes for the output database.
        // Similarly for the element "ids" field related to bulk data
        // using element ids.
        const stk::mesh::BulkData &bulk = params.bulk_data();
        stk::mesh::EntityRank rank = get_output_rank(params);

        std::vector<mesh::Entity> nodes;
        size_t numNodes = get_entities_for_nodeblock(params, part, rank,
                                                      nodes, true);

        std::vector<INT> nodeIds;
        nodeIds.reserve(numNodes);
        for(size_t i=0; i<numNodes; ++i) {
          const mesh::Entity node = nodes[i] ;
          nodeIds.push_back(bulk.identifier(node));
        }

        size_t numIdsWritten = nb.put_field_data("ids", nodeIds);
        if ( numNodes != numIdsWritten) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::NodeBlock::put_field_data:" ;
          msg << " numNodes = " << numNodes ;
          msg << " , num_ids_written = " << numIdsWritten ;
          throw std::runtime_error( msg.str() );
        }

        if (nb.get_database()->needs_shared_node_information()) {
          std::vector<int> owningProcessor;
          owningProcessor.reserve(numNodes);
          for(size_t i=0; i<numNodes; ++i) {
            owningProcessor.push_back(bulk.parallel_owner_rank(nodes[i]));
          }
          nb.put_field_data("owning_processor", owningProcessor);
        }

        const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
        const mesh::FieldBase *coordField = metaData.coordinate_field();
        assert(coordField != nullptr);
        field_data_to_ioss(bulk, coordField, nodes, &nb, "mesh_model_coordinates", Ioss::Field::MESH);

        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          if (stk::io::is_valid_part_field(f, part_primary_entity_rank(part), part,
                                           Ioss::Field::ATTRIBUTE)) {
            stk::io::field_data_to_ioss(bulk, f, nodes, &nb, f->name(), Ioss::Field::ATTRIBUTE);
          }
        }
      }

      std::pair<stk::mesh::Entity, unsigned>
      get_parent_element(stk::io::OutputParams &params, stk::mesh::Entity obj, const stk::mesh::Part* parentBlock = nullptr)
      {
        std::pair<stk::mesh::Entity, unsigned> parent(stk::mesh::Entity(), 0U);

        const stk::mesh::BulkData& stkmesh = params.bulk_data();
        const stk::topology objTopology = stkmesh.bucket(obj).topology();
        const stk::mesh::Entity* elems = stkmesh.begin_elements(obj);
        const stk::mesh::ConnectivityOrdinal* elemOrdinals = stkmesh.begin_element_ordinals(obj);
        const stk::mesh::Permutation* elemPermutations = stkmesh.begin_element_permutations(obj);

        const stk::mesh::Selector* subsetSelector = params.get_subset_selector();
        bool activeOnly = subsetSelector != nullptr;

        for(unsigned ielem = 0, e = stkmesh.num_elements(obj); ielem < e; ++ielem) {
          stk::mesh::Entity elem = elems[ielem];
          unsigned elemSideOrdinal = elemOrdinals[ielem];

          stk::mesh::Bucket &elemBucket = stkmesh.bucket(elem);

          if(stkmesh.bucket(elem).owned() && (!activeOnly || (activeOnly && (*subsetSelector)(elemBucket)))) {
            if((parentBlock == nullptr && objTopology.is_positive_polarity(elemPermutations[ielem])) ||
               (parentBlock != nullptr && contain(stkmesh, elem, parentBlock))) {
              if(params.has_output_selector(stk::topology::ELEMENT_RANK) && !params.get_is_restart()) {
                // See if elem is a member of any of the includedMeshBlocks.
                const stk::mesh::Selector* outputSelector = params.get_output_selector(stk::topology::ELEMENT_RANK);
                if((*outputSelector)(elemBucket)) {
                  parent.first = elem;
                  parent.second = elemSideOrdinal;
                  return parent;
                }
                return parent;
              }
              else {
                parent.first = elem;
                parent.second = elemSideOrdinal;
              }
              return parent;
            }
          }
        }
        return parent;
      }

      template <typename INT>
      void output_element_block_skin_map(stk::io::OutputParams &params, Ioss::ElementBlock *block,
                                         const std::vector<mesh::Entity>& meshObjects)
      {
        const stk::mesh::BulkData& stkmesh = params.bulk_data();
        bool skinMesh = params.has_skin_mesh_selector();
        if(!skinMesh) return; // This map only supported for skinning the mesh.

        size_t entitySize = block->get_property("entity_count").get_int();
        if(!block->field_exists("skin")) {
          block->field_add(Ioss::Field("skin", Ioss::Field::INTEGER, "Real[2]", Ioss::Field::MESH, entitySize));
        }

        size_t count = block->get_field("skin").raw_count();
        int mapSize = block->get_field("skin").get_size();
        std::vector<INT> elemFace(mapSize);

        if(count > 0) {
          // global element id + local face of that element.

          size_t i = 0;
          size_t faceCount = meshObjects.size();
          assert(faceCount == count);
          for(size_t j = 0; j < faceCount; j++) {
            stk::mesh::Entity face = meshObjects[j];
            std::pair<stk::mesh::Entity, unsigned> elemFacePair = get_parent_element(params, face);
            if(stkmesh.is_valid(elemFacePair.first)) {
              elemFace[i++] = stkmesh.identifier(elemFacePair.first);
              elemFace[i++] = elemFacePair.second + 1;
            }
          }

          assert(i == 2 * count);
        }
        block->put_field_data("skin", elemFace.data(), mapSize);
      }

      template <typename INT>
      void output_element_block(stk::io::OutputParams &params, Ioss::ElementBlock *block)
      {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
        const std::string& name = block->name();
        mesh::Part* part = metaData.get_part(name);
        assert(part != nullptr);

        stk::topology topo = part->topology();
        if (params.has_skin_mesh_selector() && topo == stk::topology::PARTICLE) {
          return;
        }

        std::vector<mesh::Entity> elements;
        stk::mesh::EntityRank type = part_primary_entity_rank(*part);
        if (params.has_skin_mesh_selector()) {
          type = metaData.side_rank();
        }
        size_t numElems = get_entities(params, *part, type, elements, false);

        if (numElems >  0 && topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
          throw std::runtime_error( msg.str() );
        }

        size_t nodesPerElem = block->get_property("topology_node_count").get_int();

        std::vector<INT> elemIds;
        elemIds.reserve(numElems == 0 ? 1 : numElems);
        std::vector<INT> connectivity;
        connectivity.reserve( (numElems*nodesPerElem) == 0 ? 1 : (numElems*nodesPerElem));

        for (size_t i = 0; i < numElems; ++i) {
          elemIds.push_back(bulk.identifier(elements[i]));
          stk::mesh::Entity const * elemNodes = bulk.begin_nodes(elements[i]);

          for (size_t j = 0; j < nodesPerElem; ++j) {
            connectivity.push_back(bulk.identifier(elemNodes[j]));
          }
        }

        const size_t numIdsWritten = block->put_field_data("ids", elemIds);
        const size_t numConWritten = block->put_field_data("connectivity", connectivity);

        if ( numElems != numIdsWritten || numElems != numConWritten ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::ElementBlock::put_field_data:" << std::endl ;
          msg << "  numElems = " << numElems << std::endl ;
          msg << "  numIdsWritten = " << numIdsWritten << std::endl ;
          msg << "  num_connectivity_written = " << numConWritten << std::endl ;
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank elemRank = stk::topology::ELEMENT_RANK;
        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != nullptr && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, elemRank, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, elements, block, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }

        process_element_attributes_for_output(params, *part);

        output_element_block_skin_map<INT>(params, block, elements);
      }

      template<typename INT>
      void output_nodeset_distribution_factor(const stk::mesh::BulkData& bulk,
                                              Ioss::NodeSet* ns,
                                              stk::mesh::Part* part,
                                              std::vector<stk::mesh::Entity>& nodes)
      {
          const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
          const std::string& name = ns->name();
          const std::string dfName = s_distributionFactors + "_" + name;
          stk::mesh::Field<double>* dfField = metaData.get_field<double>(stk::topology::NODE_RANK, dfName);

          if(dfField != nullptr) {
              const stk::mesh::FieldBase::Restriction& res = stk::mesh::find_restriction(*dfField, stk::topology::NODE_RANK, *part);
              if(res.num_scalars_per_entity() > 0) {
                  stk::io::field_data_to_ioss(bulk, dfField, nodes, ns, s_distributionFactors, Ioss::Field::MESH);
              }
          } else {
              assert(ns->field_exists(s_distributionFactors));
              size_t dfSize = ns->get_field(s_distributionFactors).raw_count();
              std::vector<double> df;
              df.reserve(dfSize);
              const auto* const nodeFactorVar = get_distribution_factor_field(*part);
              if((nodeFactorVar != nullptr) && (nodeFactorVar->entity_rank() == stk::topology::NODE_RANK)) {
                  nodeFactorVar->sync_to_host();
                  for(auto& node : nodes) {
                      df.push_back(*(double*) (stk::mesh::field_data(*nodeFactorVar, node)));
                  }
              } else {
                  size_t count = nodes.size();
                  for(size_t i = 0; i < count; ++i) {
                      df.push_back(1.0);
                  }
              }
              ns->put_field_data(s_distributionFactors, df);
          }
      }

      template <typename INT>
      void output_node_set(stk::io::OutputParams &params, Ioss::NodeSet *ns)
      {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
        const std::string& name = ns->name();
        mesh::Part* part = metaData.get_part(name);

        // If part is null, then it is possible that this nodeset is a "viz nodeset" which
        // means that it is a nodeset containing the nodes of an element block.
        // See if there is a property base_stk_part_name and if so, get the part with
        // that name.
        if (part == nullptr) {
          if (ns->property_exists(base_stk_part_name)) {
            std::string baseName = ns->get_property(base_stk_part_name).get_string();
            part = metaData.get_part(baseName);
          }
          if (part == nullptr) {
            std::ostringstream msg ;
            msg << " FAILED in Ioss::NodeSet::output_node_set:"
                << " Could not find stk part corresponding to nodeset named '"
                << name << "'";
            throw std::runtime_error( msg.str() );
          }
        }

        std::vector<stk::mesh::Entity> nodes;
        mesh::EntityRank rank = get_output_rank(params);
        size_t numNodes = get_entities_for_nodeblock(params, *part, rank, nodes, true);

        std::vector<INT> node_ids;
        node_ids.reserve(numNodes);
        for(size_t i=0; i<numNodes; ++i) {
          const stk::mesh::Entity node = nodes[i] ;
          node_ids.push_back(bulk.identifier(node));
        }

        size_t numIdsWritten = ns->put_field_data("ids", node_ids);
        if ( numNodes != numIdsWritten ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::NodeSet::output_node_set:"
              << " numNodes = " << numNodes
              << ", numIdsWritten = " << numIdsWritten;
          throw std::runtime_error( msg.str() );
        }

        output_nodeset_distribution_factor<INT>(bulk, ns, part, nodes);

        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != nullptr && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, stk::topology::NODE_RANK, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, nodes, ns, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }
      }

      template <typename INT>
      void output_communication_maps(stk::io::OutputParams &params)
      {
        Ioss::Region &ioRegion = params.io_region();
        const stk::mesh::BulkData &bulk = params.bulk_data();
        mesh::EntityRank rank = get_output_rank(params);
        const stk::mesh::Selector *subsetSelector = params.get_subset_selector();
        const stk::mesh::Selector *outputSelector = params.get_output_selector(rank);

        if (bulk.parallel_size() > 1) {
          const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
          mesh::Selector selector = meta.globally_shared_part();
          if (subsetSelector) selector &= *subsetSelector;
          if (outputSelector) selector &= *outputSelector;

          std::vector<mesh::Entity> entities;
          get_selected_nodes(params, selector, entities);

          const std::string csName("node_symm_comm_spec");
          Ioss::CommSet * ioCs = ioRegion.get_commset(csName);
          STKIORequire(ioCs != nullptr);

          // Allocate data space to store <id, processor> pair
          assert(ioCs->field_exists("entity_processor"));
          size_t size = ioCs->get_field("entity_processor").raw_count();

          std::vector<INT> ep;
          ep.reserve(size*2);

          std::vector<int> sharingProcs;
          for (size_t i=0; i < entities.size(); i++) {
            bulk.comm_shared_procs(bulk.entity_key(entities[i]), sharingProcs);
            for ( size_t j=0; j<sharingProcs.size(); j++ ) {
              ep.push_back(bulk.identifier(entities[i]));
              ep.push_back(sharingProcs[j]);
            }
          }
          assert(size*2 == ep.size());
          ioCs->put_field_data("entity_processor", ep);
        }
      }

      template <typename INT>
      void output_side_set(stk::io::OutputParams &params, Ioss::SideSet *ss)
      {
        const stk::mesh::MetaData & meta = params.bulk_data().mesh_meta_data();

        size_t blockCount = ss->block_count();
        for (size_t i=0; i < blockCount; i++) {
          Ioss::SideBlock *block = ss->get_block(i);
          if (stk::io::include_entity(block)) {
            stk::mesh::Part * part = meta.get_part(block->name());
            const Ioss::ElementTopology *parent_topology = block->parent_element_topology();
            stk::io::write_side_data_to_ioss<INT>(params, *block, part, parent_topology);
          }
        }
      }

      template <typename INT>
      void output_face_block(stk::io::OutputParams &params, Ioss::FaceBlock *fb)
      {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
        const std::string& name = fb->name();
        mesh::Part* part = metaData.get_part(name);
        assert(part != nullptr);

        stk::topology topo = part->topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
          throw std::runtime_error( msg.str() );
        }

        std::vector<mesh::Entity> faces;
        stk::mesh::EntityRank type = part_primary_entity_rank(*part);
        size_t numFaces = get_entities(params, *part, type, faces, false);

        size_t nodesPerFace = fb->get_property("topology_node_count").get_int();

        std::vector<INT> faceIds;
        faceIds.reserve(numFaces == 0 ? 1 : numFaces);
        std::vector<INT> connectivity; connectivity.reserve( (numFaces*nodesPerFace) == 0 ? 1 : (numFaces*nodesPerFace));

        for (size_t i = 0; i < numFaces; ++i) {
          faceIds.push_back(bulk.identifier(faces[i]));
          stk::mesh::Entity const * faceNodes = bulk.begin_nodes(faces[i]);

          for (size_t j = 0; j < nodesPerFace; ++j) {
            connectivity.push_back(bulk.identifier(faceNodes[j]));
          }
        }

        const size_t numIdsWritten = fb->put_field_data("ids", faceIds);
        const size_t numConWritten = fb->put_field_data("connectivity", connectivity);

        if ( numFaces != numIdsWritten || numFaces != numConWritten ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::FaceBlock::put_field_data:" << std::endl ;
          msg << "  numFaces = " << numFaces << std::endl ;
          msg << "  numIdsWritten = " << numIdsWritten << std::endl ;
          msg << "  num_connectivity_written = " << numConWritten << std::endl ;
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank faceRank = stk::topology::FACE_RANK;
        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        for(const mesh::FieldBase* f : fields) {
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != nullptr && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, faceRank, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, faces, fb, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }
      }

      template <typename INT>
      void output_edge_block(stk::io::OutputParams &params, Ioss::EdgeBlock *eb)
      {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        const stk::mesh::MetaData & metaData = bulk.mesh_meta_data();
        const std::string& name = eb->name();
        mesh::Part* part = metaData.get_part(name);
        assert(part != nullptr);

        stk::topology topo = part->topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
          throw std::runtime_error( msg.str() );
        }

        std::vector<mesh::Entity> edges;
        stk::mesh::EntityRank type = part_primary_entity_rank(*part);
        size_t numEdges = get_entities(params, *part, type, edges, false);

        size_t nodesPerEdge = eb->get_property("topology_node_count").get_int();

        std::vector<INT> edgeIds;
        edgeIds.reserve(numEdges == 0 ? 1 : numEdges);
        std::vector<INT> connectivity;
        connectivity.reserve( (numEdges*nodesPerEdge) == 0 ? 1 : (numEdges*nodesPerEdge));

        for (size_t i = 0; i < numEdges; ++i) {
          edgeIds.push_back(bulk.identifier(edges[i]));
          stk::mesh::Entity const * edgeNodes = bulk.begin_nodes(edges[i]);

          for (size_t j = 0; j < nodesPerEdge; ++j) {
            connectivity.push_back(bulk.identifier(edgeNodes[j]));
          }
        }

        const size_t numIdsWritten = eb->put_field_data("ids", edgeIds);
        const size_t numConWritten = eb->put_field_data("connectivity", connectivity);

        if ( numEdges != numIdsWritten || numEdges != numConWritten ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::EdgeBlock::put_field_data:" << std::endl ;
          msg << "  numEdges = " << numEdges << std::endl ;
          msg << "  numIdsWritten = " << numIdsWritten << std::endl ;
          msg << "  num_connectivity_written = " << numConWritten << std::endl ;
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
        const std::vector<mesh::FieldBase *> &fields = metaData.get_fields();
        for(const mesh::FieldBase* f : fields) {
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != nullptr && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, edge_rank, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, edges, eb, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }
      }
    }

    void write_output_db_node_block(stk::io::OutputParams &params)
    {
        const stk::mesh::MetaData & meta = params.bulk_data().mesh_meta_data();
        Ioss::Region &ioRegion = params.io_region();

        bool ints64bit = db_api_int_size(&ioRegion) == 8;

        Ioss::NodeBlock & nb = *ioRegion.get_node_blocks()[0];

        if (ints64bit)
          output_node_block<int64_t>(params, nb, meta.universal_part());
        else
          output_node_block<int>(params, nb, meta.universal_part());
    }

    void write_output_db_element_blocks(stk::io::OutputParams &params)
    {
      Ioss::Region &ioRegion = params.io_region();
      bool ints64bit = db_api_int_size(&ioRegion) == 8;

      //----------------------------------
      const Ioss::ElementBlockContainer& elemBlocks = ioRegion.get_element_blocks();
      for(Ioss::ElementBlockContainer::const_iterator it = elemBlocks.begin();
          it != elemBlocks.end(); ++it) {
        if (ints64bit)
          output_element_block<int64_t>(params, *it);
        else
          output_element_block<int>(params, *it);
      }
    }

    template <typename T>
    void write_output_db_for_entitysets_and_comm_map(stk::io::OutputParams &params)
    {
        Ioss::Region &ioRegion = params.io_region();
        auto *dbo = ioRegion.get_database();
        const auto supports = dbo->entity_field_support();

        if (supports & Ioss::NODESET) {
          for(Ioss::NodeSet *ns : ioRegion.get_nodesets()) {
            output_node_set<T>(params, ns);
          }
        }

        if (supports & Ioss::SIDESET) {
          for(Ioss::SideSet *ss : ioRegion.get_sidesets()) {
            output_side_set<T>(params, ss);
          }
        }

        if (supports & Ioss::EDGEBLOCK) {
          for(Ioss::EdgeBlock *eb: ioRegion.get_edge_blocks()) {
            output_edge_block<T>(params, eb);
          }
        }

        if (supports & Ioss::FACEBLOCK) {
          for(Ioss::FaceBlock *fb: ioRegion.get_face_blocks()) {
            output_face_block<T>(params, fb);
          }
        }

        output_communication_maps<T>(params);
    }

    void write_output_db_rest_of_mesh(stk::io::OutputParams &params)
    {
      Ioss::Region &ioRegion = params.io_region();

      write_output_db_element_blocks(params);

      bool ints64bit = db_api_int_size(&ioRegion) == 8;

        if (ints64bit) {
            write_output_db_for_entitysets_and_comm_map<int64_t>(params);
        } else {
            write_output_db_for_entitysets_and_comm_map<int>(params);
        }
    }

    void write_output_db(stk::io::OutputParams &params)
    {
      Ioss::Region &ioRegion = params.io_region();

      ioRegion.begin_mode( Ioss::STATE_MODEL );
      write_output_db_node_block(params);
      write_output_db_rest_of_mesh(params);
      ioRegion.end_mode( Ioss::STATE_MODEL );
    }

    //----------------------------------------------------------------------
    bool is_part_io_part(const stk::mesh::Part* part)
    {
      if(part == nullptr) { return false; }
      return is_part_io_part(*part);
    }

    bool is_part_io_part(const stk::mesh::Part &part)
    {
      return stk::mesh::impl::has_part_attribute<IossPartAttribute>(part);
    }

    // TODO: NOTE: The use of "FieldBase" here basically eliminates the use of the attribute
    // for any other fieldbase.  This is just being done now for a proof-of-concept for use
    // in the framework-based stk-mesh handoff to a stk-based stk-mesh...  If this looks promising,
    // then need to wrap it in a different classs...
    const stk::mesh::FieldBase *get_distribution_factor_field(const stk::mesh::Part &p)
    {
      return p.attribute<stk::mesh::FieldBase>();
    }

    void set_distribution_factor_field(stk::mesh::Part &p,
                                       const stk::mesh::FieldBase &dfField)
    {
      stk::mesh::MetaData &m = mesh::MetaData::get(p);
      if (const stk::mesh::FieldBase * existingDistFactField = p.attribute<stk::mesh::FieldBase>()) {
        m.remove_attribute(p, existingDistFactField);
      }

      m.declare_attribute_no_delete(p, &dfField);
    }

    const Ioss::Field::RoleType* get_field_role(const stk::mesh::FieldBase &f)
    {
      return f.attribute<Ioss::Field::RoleType>();
    }

    void set_field_role(stk::mesh::FieldBase &f, const Ioss::Field::RoleType &role)
    {
      Ioss::Field::RoleType *myRole = new Ioss::Field::RoleType(role);
      stk::mesh::MetaData &m = mesh::MetaData::get(f);
      const Ioss::Field::RoleType *check = m.declare_attribute_with_delete(f, myRole);
      if ( check != myRole ) {
        if (*check != *myRole) {
          std::ostringstream msg ;
          msg << " FAILED in IossBridge -- set_field_role:"
              << " The role type for field name= " << f.name()
              << " was already set to " << *check
              << ", so it is not possible to change it to " << *myRole;
          delete myRole;
          throw std::runtime_error( msg.str() );
        }
        delete myRole;
      }
    }

    namespace {
      void define_input_nodeblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::NodeBlockContainer& nodeBlocks = region.get_node_blocks();
        assert(nodeBlocks.size() == 1);

        Ioss::NodeBlock *nb = nodeBlocks[0];
        stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT,
                                  meta.universal_part(), stk::topology::NODE_RANK);
      }

      void define_input_elementblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::ElementBlockContainer& elemBlocks = region.get_element_blocks();
        for(size_t i=0; i < elemBlocks.size(); i++) {
          if (stk::io::include_entity(elemBlocks[i])) {
            stk::mesh::Part* const part = meta.get_part(elemBlocks[i]->name());
            assert(part != nullptr);
            stk::io::define_io_fields(elemBlocks[i], Ioss::Field::TRANSIENT,
                                      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_nodeset_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::NodeSetContainer& nodesets = region.get_nodesets();
        for(size_t i=0; i < nodesets.size(); i++) {
          if (stk::io::include_entity(nodesets[i])) {
            stk::mesh::Part* const part = meta.get_part(nodesets[i]->name());
            assert(part != nullptr);
            stk::io::define_io_fields(nodesets[i], Ioss::Field::TRANSIENT,
                                      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_sideset_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        unsigned sideRank = meta.side_rank();
        if (meta.spatial_dimension() <= sideRank) return;

        const Ioss::SideSetContainer& sideSets = region.get_sidesets();
        for(Ioss::SideSetContainer::const_iterator it = sideSets.begin();
            it != sideSets.end(); ++it) {
          Ioss::SideSet *entity = *it;
          if (stk::io::include_entity(entity)) {
            const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
            for(size_t i=0; i < blocks.size(); i++) {
              if (stk::io::include_entity(blocks[i])) {
                stk::mesh::Part* const part = meta.get_part(blocks[i]->name());
                assert(part != nullptr);
                stk::io::define_io_fields(blocks[i], Ioss::Field::TRANSIENT,
                                          *part, part_primary_entity_rank(*part));
              }
            }
          }
        }
      }

      void define_input_face_block_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::FaceBlockContainer& faceBlocks = region.get_face_blocks();
        for(size_t i=0; i < faceBlocks.size(); i++) {
          if (stk::io::include_entity(faceBlocks[i])) {
            stk::mesh::Part* const part = meta.get_part(faceBlocks[i]->name());
            assert(part != nullptr);
            stk::io::define_io_fields(faceBlocks[i], Ioss::Field::TRANSIENT,
                                      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_edge_block_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::EdgeBlockContainer& edgeBlocks = region.get_edge_blocks();
        for(size_t i=0; i < edgeBlocks.size(); i++) {
          if (stk::io::include_entity(edgeBlocks[i])) {
            stk::mesh::Part* const part = meta.get_part(edgeBlocks[i]->name());
            assert(part != nullptr);
            stk::io::define_io_fields(edgeBlocks[i], Ioss::Field::TRANSIENT,
                                      *part, part_primary_entity_rank(*part));
          }
        }
      }
    }

    // ========================================================================
    // Iterate over all Ioss entities in the input mesh Ioss region and
    // define a stk_field for all transient fields found.  The stk field
    // will have the same name as the field on the database.
    //
    // Note that all fields found on the database will have a
    // corresponding stk field defined.  If you want just a selected
    // subset of the defined fields, you will need to define the fields
    // manually.
    //
    // To populate the stk field with data from the database, call
    // process_input_request().
    void define_input_fields(Ioss::Region &region,  stk::mesh::MetaData &meta)
    {
      define_input_nodeblock_fields(region, meta);
      define_input_elementblock_fields(region, meta);
      define_input_nodeset_fields(region, meta);
      define_input_sideset_fields(region, meta);
      define_input_edge_block_fields(region, meta);
      define_input_face_block_fields(region, meta);
    }

    void insert_var_names_for_part(const Ioss::GroupingEntity* entity, stk::mesh::Part* part, FieldNameToPartVector& names)
    {
        Ioss::NameList geNames;
        entity->field_describe(Ioss::Field::TRANSIENT, &geNames);
        for(const std::string& geName : geNames)
            names.emplace_back(geName, part);
    }

    void insert_var_names(const Ioss::GroupingEntity *entity, FieldNameToPartVector &names, const stk::mesh::MetaData& meta)
    {
        if (stk::io::include_entity(entity)) {
            stk::mesh::Part* part = meta.get_part(entity->name());
            insert_var_names_for_part(entity, part, names);
        }
    }

    template <typename GroupingEntityVector>
    FieldNameToPartVector get_grouping_entity_var_names(const GroupingEntityVector &groupingEntities,
                                                        const stk::mesh::MetaData& meta)
    {
        FieldNameToPartVector names;
        for(size_t i=0; i < groupingEntities.size(); i++)
            insert_var_names(groupingEntities[i], names, meta);
        stk::util::sort_and_unique(names, FieldNameToPartLess());
        return names;
    }

    FieldNameToPartVector get_var_names(Ioss::Region &region, Ioss::EntityType type, const stk::mesh::MetaData& meta)
    {
        switch(type)
        {
            case Ioss::NODEBLOCK:
            {
                FieldNameToPartVector names;
                const Ioss::NodeBlockContainer nodeBlocks = region.get_node_blocks();
                STK_ThrowRequire(nodeBlocks.size() == 1);
                insert_var_names_for_part(nodeBlocks[0], &meta.universal_part(), names);
                return names;
            }
            case Ioss::ELEMENTBLOCK:
                return get_grouping_entity_var_names(region.get_element_blocks(), meta);
            case Ioss::NODESET:
                return get_grouping_entity_var_names(region.get_nodesets(), meta);
            case Ioss::SIDESET:
            {
                FieldNameToPartVector names;
                for(const Ioss::SideSet *sideset : region.get_sidesets())
                    for(const Ioss::SideBlock *sideBlock : sideset->get_side_blocks())
                        insert_var_names(sideBlock, names, meta);
                stk::util::sort_and_unique(names, FieldNameToPartLess());
                return names;
            }
            default:
                return FieldNameToPartVector();
        }
    }

    void put_field_data(stk::io::OutputParams& params,
                        stk::mesh::Part &part,
                        stk::mesh::EntityRank partType,
                        Ioss::GroupingEntity *ioEntity,
                        Ioss::Field::RoleType filterRole)
    {
      std::vector<stk::mesh::Entity> entities;
      stk::io::get_output_entity_list(ioEntity, partType, params, entities);

      const stk::mesh::BulkData &bulk = params.bulk_data();
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

      std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
      while (I != fields.end()) {
        const stk::mesh::FieldBase *f = *I; ++I;
        if (stk::io::is_valid_part_field(f, partType, part, filterRole)) {
          stk::io::field_data_to_ioss(bulk, f, entities, ioEntity, f->name(), filterRole);
        }
      }
    }

    void put_field_data(stk::mesh::BulkData &bulk,
                        stk::mesh::Part &part,
                        stk::mesh::EntityRank partType,
                        Ioss::GroupingEntity *ioEntity,
                        Ioss::Field::RoleType filterRole)
    {
      stk::io::OutputParams params(bulk);
      put_field_data(params, part, partType, ioEntity, filterRole);
    }

    struct DefineOutputFunctor
    {
      void operator()(stk::io::OutputParams& params, stk::mesh::Part &part, stk::mesh::EntityRank rank, Ioss::GroupingEntity *ge, Ioss::Field::RoleType role)
      {  stk::io::ioss_add_fields(part, rank, ge, role); }
    };

    struct ProcessOutputFunctor
    {
      void operator()(stk::io::OutputParams& params, stk::mesh::Part &part, stk::mesh::EntityRank rank, Ioss::GroupingEntity *ge, Ioss::Field::RoleType role)
      {  put_field_data(params, part, rank, ge, role); }
    };

    template <typename T>
    void process_field_loop(stk::io::OutputParams& params, T& callable)
    {
        Ioss::Region &region = params.io_region();
        const stk::mesh::BulkData &bulk = params.bulk_data();

        const stk::mesh::MetaData & meta = bulk.mesh_meta_data();

        Ioss::NodeBlock *nb = region.get_node_blocks()[0];
        callable(params, meta.universal_part(), stk::topology::NODE_RANK,
                 dynamic_cast<Ioss::GroupingEntity *>(nb), Ioss::Field::TRANSIENT);

        const stk::mesh::PartVector & allParts = meta.get_parts();
        for ( stk::mesh::PartVector::const_iterator
                ip = allParts.begin(); ip != allParts.end(); ++ip ) {

          stk::mesh::Part * const part = *ip;

          if (stk::io::is_part_io_part(*part)) {
            Ioss::GroupingEntity *entity = region.get_entity(part->name());
            STK_ThrowRequireMsg(entity != nullptr, "Could not find Ioss grouping entity with name: " + part->name() + ". Contact sierra-help.");

            if (entity->type() == Ioss::SIDESET) {
              Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
              int block_count = sset->block_count();

              for (int i=0; i < block_count; i++) {
                Ioss::SideBlock *fb = sset->get_block(i);
                callable(params, *part,
                         stk::mesh::EntityRank( part->primary_entity_rank() ),
                         dynamic_cast<Ioss::GroupingEntity *>(fb), Ioss::Field::TRANSIENT);
              }
            } else {
              callable(params, *part,
                       stk::mesh::EntityRank( part->primary_entity_rank() ),
                       entity, Ioss::Field::TRANSIENT);
            }
          }
        }
    }

    void process_output_request(stk::io::OutputParams& params, int step)
    {
      params.io_region().begin_state(step);
      ProcessOutputFunctor functor;
      process_field_loop(params, functor);
      params.io_region().end_state(step);
    }

    template <typename INT>
    void output_node_sharing_info( Ioss::CommSet* ioCs,  const EntitySharingInfo &nodeSharingInfo)
    {
        std::vector<INT> entity_proc(2*nodeSharingInfo.size());
        int counter = 0;
        for(auto &nodeProc : nodeSharingInfo)
        {
            entity_proc[counter]     = nodeProc.first;
            entity_proc[counter + 1] = nodeProc.second;
            counter += 2;
        }
        size_t size_field = entity_proc.size()*(sizeof(INT));
        ioCs->put_field_data("entity_processor", entity_proc.data(), size_field);
    }

    void write_node_sharing_info(Ioss::DatabaseIO *dbo, const EntitySharingInfo &nodeSharingInfo)
    {
        bool ints64bit =  db_api_int_size(dbo->get_region()) == 8;

        Ioss::CommSet* ioCs = dbo->get_region()->get_commset("commset_node");
        if(ioCs)
        {
          if (ints64bit)
            output_node_sharing_info<int64_t>(ioCs, nodeSharingInfo);
          else
            output_node_sharing_info<int>(ioCs, nodeSharingInfo);
        }
    }

    Ioss::DatabaseIO *create_database_for_subdomain(const std::string &baseFilename,
                                                    int indexSubdomain,
                                                    int numSubdomains)
    {
        std::string parallelFilename{construct_filename_for_serial_or_parallel(baseFilename, numSubdomains, indexSubdomain)};

        std::string dbtype("exodusII");
        Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, parallelFilename, Ioss::WRITE_RESULTS, MPI_COMM_SELF);

        return dbo;
    }

    void write_mesh_data_for_subdomain(stk::io::OutputParams& params, const EntitySharingInfo& nodeSharingInfo)
    {
      Ioss::Region &region = params.io_region();
      region.begin_mode(Ioss::STATE_DEFINE_MODEL);
      stk::io::define_output_db_within_state_define(params, {});
      Ioss::CommSet *commset = new Ioss::CommSet(region.get_database(), "commset_node", "node", nodeSharingInfo.size());
      commset->property_add(Ioss::Property("id", 1));
      region.add(commset);
      region.end_mode(Ioss::STATE_DEFINE_MODEL);

      region.begin_mode(Ioss::STATE_MODEL);
      stk::io::write_output_db_node_block(params);
      write_node_sharing_info(region.get_database(), nodeSharingInfo);
      stk::io::write_output_db_rest_of_mesh(params);
      region.end_mode(Ioss::STATE_MODEL);
    }

    int write_transient_data_for_subdomain(stk::io::OutputParams& params, double timeStep)
    {
      Ioss::Region &outRegion = params.io_region();

      if(!outRegion.transient_defined()) {
        outRegion.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        DefineOutputFunctor functor;
        process_field_loop(params, functor);
        outRegion.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      }

      outRegion.begin_mode(Ioss::STATE_TRANSIENT);
      int out_step = outRegion.add_state(timeStep);
      process_output_request(params, out_step);
      outRegion.end_mode(Ioss::STATE_TRANSIENT);

      return out_step;
    }

    void write_file_for_subdomain(stk::io::OutputParams& params,
                                  const EntitySharingInfo &nodeSharingInfo,
                                  int numSteps,
                                  double timeStep)
    {
        Ioss::Region &outRegion = params.io_region();

        Ioss::DatabaseIO *dbo = outRegion.get_database();
        STK_ThrowRequire(nullptr != dbo);

        write_mesh_data_for_subdomain(params, nodeSharingInfo);

        if(numSteps > 0) {
          write_transient_data_for_subdomain(params, timeStep);
        }
    }

    void add_properties_for_subdomain(stk::io::OutputParams& params,
                                      int indexSubdomain,
                                      int numSubdomains,
                                      int globalNumNodes,
                                      int globalNumElems)
    {
        Ioss::Region &outRegion = params.io_region();

        outRegion.property_add(Ioss::Property("processor_count", numSubdomains));
        outRegion.property_add(Ioss::Property("my_processor", indexSubdomain));
        outRegion.property_add(Ioss::Property("global_node_count", globalNumNodes));
        outRegion.property_add(Ioss::Property("global_element_count", globalNumElems));

        if(params.bulk_data().supports_large_ids()) {
            outRegion.property_add(Ioss::Property("INTEGER_SIZE_API" , 8));
            outRegion.property_add(Ioss::Property("INTEGER_SIZE_DB" , 8));

            Ioss::DatabaseIO *dbo = outRegion.get_database();
            dbo->set_int_byte_size_api(Ioss::USE_INT64_API);
        }
    }

    void write_file_for_subdomain(const std::string &baseFilename,
                                  int indexSubdomain,
                                  int numSubdomains,
                                  int globalNumNodes,
                                  int globalNumElems,
                                  stk::io::OutputParams& params,
                                  const EntitySharingInfo &nodeSharingInfo,
                                  int numSteps,
                                  double timeStep)
    {
        Ioss::DatabaseIO *dbo = create_database_for_subdomain(baseFilename, indexSubdomain, numSubdomains);
        Ioss::Region outRegion(dbo, "name");

        STK_ThrowRequireMsg(params.io_region_ptr() == nullptr, "OutputParams argument must have a NULL IORegion");
        params.set_io_region(&outRegion);
        add_properties_for_subdomain(params, indexSubdomain, numSubdomains, globalNumNodes, globalNumElems);

        write_file_for_subdomain(params, nodeSharingInfo, numSteps, timeStep);

        stk::io::delete_selector_property(outRegion);
        params.set_io_region(nullptr);
    }


    const stk::mesh::Part* get_parent_element_block_by_adjacency(const stk::mesh::BulkData& bulk,
                                                                 const std::string& name,
                                                                 const stk::mesh::Part* parentElementBlock)
    {
      const stk::mesh::Part* part = bulk.mesh_meta_data().get_part(name);
      if (part != nullptr) {
        std::vector<const stk::mesh::Part*> touchingParts = bulk.mesh_meta_data().get_blocks_touching_surface(part);
        if (touchingParts.size() == 1) {
          parentElementBlock = touchingParts[0];
        }
      }
      return parentElementBlock;
    }


    const stk::mesh::Part* get_parent_element_block(const stk::mesh::BulkData &bulk,
                                                    const Ioss::Region &ioRegion,
                                                    const std::string& name)
    {
        // Ugliness -- Handling of an embedded two-sided surface.  Currently
        // only handled in the case where the surfaces are "split by element
        // block".  Then, each side of the surface will have a different
        // parent element block and we can correctly assign the face to the
        // correct faceblock...

        std::string part_name = name + "_context";
        const stk::mesh::Part* parentElementBlock = bulk.mesh_meta_data().get_part(part_name);

        if(parentElementBlock == nullptr) {
            if(ioRegion.get_database()->get_surface_split_type() == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
                // If the surfaces were split by element block, then the surface
                // name will be of the form:  "name_block_id_facetopo_id" "name" is typically "surface".
                //                         or "name_blockname_facetopo_id"
                // See if we can extract a block name...
                std::vector<std::string> tokens;
                stk::util::tokenize(name, "_", tokens);
                if(tokens.size() >= 4) {
                    // Check whether the second-last token is a face topology
                    const Ioss::ElementTopology* faceTopo = Ioss::ElementTopology::factory(tokens[tokens.size() - 2], true);
                    if(faceTopo != nullptr) {
                        // Extract the blockname or "block"_id...
                        std::string ebName;
                        size_t lastToken = tokens.size() - 2;
                        for(size_t tok = 1; tok < lastToken; tok++) {
                            ebName += tokens[tok];
                            if(tok < lastToken - 1) ebName += "_";
                        }

                        stk::mesh::Part* elementBlock = bulk.mesh_meta_data().get_part(ebName);
                        if(elementBlock != nullptr && is_part_io_part(*elementBlock))
                            parentElementBlock = elementBlock;
                    }
                }
            }
        }

        if(parentElementBlock == nullptr) {
          parentElementBlock = get_parent_element_block_by_adjacency(bulk, name, parentElementBlock);
        }

        return parentElementBlock;
    }

    bool is_valid_for_output(stk::io::OutputParams &params, const stk::mesh::Part &part)
    {
        const stk::mesh::EntityRank rank = part.primary_entity_rank();
        const stk::mesh::Selector *outputSelector = params.get_output_selector(rank);

        bool isIoPart   = stk::io::is_part_io_part(part);
        bool isSelected = (outputSelector == nullptr) || (*outputSelector)(part);

        bool isEmptyElementBlock = false;

        if(rank == stk::topology::ELEM_RANK && params.get_filter_empty_entity_blocks()) {
          isEmptyElementBlock = is_empty_element_block(params, &part);
        }

        return (isIoPart && isSelected && !isEmptyElementBlock);
    }

    bool node_is_connected_to_local_element(const stk::mesh::BulkData &bulk, stk::mesh::Entity node, const stk::mesh::Selector *subsetSelector)
    {
        const stk::mesh::Entity *elems = bulk.begin_elements(node);
        const unsigned numElements = bulk.num_elements(node);
        bool isLocalElement = false;
        for(unsigned i = 0; i<numElements; ++i) {
            stk::mesh::Bucket &bucket = bulk.bucket(elems[i]);
            bool isSelected = (subsetSelector == nullptr) ? true : (*subsetSelector)(bucket);
            if(bucket.owned() && isSelected) {
                isLocalElement = true;
                break;
            }
        }
        return isLocalElement;
    }

    size_t count_selected_nodes(OutputParams &params, const stk::mesh::Selector &selector)
    {
        stk::mesh::EntityVector nodes;
        get_selected_nodes(params, selector, nodes);
        return nodes.size();
    }

    void filter_nodes_by_local_connectivity(const stk::mesh::BulkData& bulk,
                                            const stk::mesh::Selector* subsetSelector,
                                            stk::mesh::EntityVector& nodes)
    {
        unsigned nodesSize = nodes.size();
        for(unsigned i = 0; i < nodesSize; ++i)
        {
            unsigned index = nodesSize - i - 1;
            if(!node_is_connected_to_local_element(bulk, nodes[index], subsetSelector))
            {
                nodes.erase(nodes.begin() + index);
            }
        }
    }

    bool is_node_selected_for_output(OutputParams &params, stk::mesh::Entity node)
    {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        const stk::mesh::Selector* subsetSelector = params.get_subset_selector();

        bool hasGhosting = params.get_has_ghosting();
        bool hasAdaptivity = params.get_has_adaptivity();

        bool result = true;
        bool activeOnly = subsetSelector != nullptr;
        stk::mesh::Bucket &nodeBucket = bulk.bucket(node);
        if (hasAdaptivity) {
            result = activeOnly ? (*subsetSelector)(nodeBucket) : true;
        }
        if (hasGhosting && result) {
            // Now need to check whether this node is locally owned or is used by
            // a locally-owned element.
            if (nodeBucket.owned()) {
                result = true;
            } else {
                // Iterate elements that use this node and see if they are locally-owned.
                result = false;
                const stk::mesh::Entity* elements = bulk.begin_elements(node);

                for (unsigned i = 0, e = bulk.num_elements(node); i < e; ++i) {
                    stk::mesh::Entity elem = elements[i];
                    stk::mesh::Bucket &elemBucket = bulk.bucket(elem);
                    if (elemBucket.owned() && (!activeOnly || (activeOnly && (*subsetSelector)(elemBucket)))) {
                        result = true;
                        break;
                    }
                }
            }
        }
        return result;
    }

    void filter_nodes_by_ghosting(OutputParams &params, stk::mesh::EntityVector& nodes)
    {
        unsigned nodesSize = nodes.size();
        for(unsigned i = 0; i < nodesSize; ++i)
        {
            unsigned index = nodesSize - i - 1;
            if (!is_node_selected_for_output(params, nodes[index]))
            {
                nodes.erase(nodes.begin() + index);
            }
        }
    }

    void get_selected_nodes(OutputParams & params,
                            const stk::mesh::Selector &selector,
                            stk::mesh::EntityVector &nodes)
    {
        const stk::mesh::BulkData &bulk = params.bulk_data();
        Ioss::Region &ioRegion = params.io_region();
        nodes.clear();

        bool ignoreDisconnectedNodes = false;
        if(ioRegion.property_exists(stk::io::s_ignoreDisconnectedNodes)) {
            ignoreDisconnectedNodes = ioRegion.get_property(stk::io::s_ignoreDisconnectedNodes).get_int();
        }

        const bool sortById = true;
        stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, selector, nodes, sortById);
        filter_nodes_by_ghosting(params, nodes);
        if(!ignoreDisconnectedNodes) {
            return;
        }

        filter_nodes_by_local_connectivity(bulk, params.get_subset_selector(), nodes);
    }


  }//namespace io
}//namespace stk
