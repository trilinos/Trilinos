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

#include <stk_io/IossBridge.hpp>
#include <Ioss_NullEntity.h>            // for NullEntity
#include <assert.h>                     // for assert
#include <stdint.h>                     // for int64_t
#include <Shards_Array.hpp>             // for ArrayDimension
#include <algorithm>                    // for sort
#include <complex>                      // for complex
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <iomanip>
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for EntityLess, BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian, Matrix, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities, etc
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank, etc
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/tokenize.hpp>   // for tokenize
#include "Ioss_CommSet.h"               // for CommSet
#include "Ioss_DatabaseIO.h"            // for DatabaseIO
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_ElementTopology.h"       // for ElementTopology
#include "Ioss_EntityBlock.h"           // for EntityBlock
#include "Ioss_EntityType.h"            // for EntityType::ELEMENTBLOCK, etc
#include <Ioss_IOFactory.h>
#include "Ioss_Field.h"                 // for Field, Field::RoleType, etc
#include "Ioss_GroupingEntity.h"        // for GroupingEntity
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_NodeSet.h"               // for NodeSet
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, SideSetContainer, etc
#include "Ioss_SideBlock.h"             // for SideBlock
#include "Ioss_SideSet.h"               // for SideSet, SideBlockContainer
#include "Ioss_State.h"                 // for State::STATE_DEFINE_MODEL, etc
#include "Ioss_Utils.h"                 // for Utils
#include "Ioss_VariableType.h"          // for NameList, VariableType
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_topology/topology.hpp"    // for topology::num_nodes
#include "stk_util/util/PairIter.hpp"   // for PairIter

#include "SidesetTranslator.hpp"
#include "StkIoUtils.hpp"

namespace stk {
  namespace io {
    bool is_field_on_part(const stk::mesh::FieldBase *field,
                          const stk::mesh::EntityRank part_type,
                          const stk::mesh::Part &part);
  }
}

void STKIORequire(bool cond)
{
  if (!cond) throw std::runtime_error("");
}

namespace {

  const std::string internal_selector_name = "_stk_io_internal_selector";
  const std::string base_stk_part_name = "_base_stk_part_name";
  const std::string block_nodes_suffix = "_nodes";

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
        assert(sset != NULL);
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
        assert(sblk != NULL);
        int rank = sblk->topology()->parametric_dimension();
        if (rank == 2)
          return stk::topology::FACE_RANK;
        if (rank == 1)
          return stk::topology::EDGE_RANK;
        if (rank == 0)
          return stk::topology::NODE_RANK;
        else
          return stk::mesh::InvalidEntityRank;
      }
    default:
      return stk::mesh::InvalidEntityRank;
    }
  }


  template <typename T>
  const stk::mesh::FieldBase *declare_ioss_field_internal(stk::mesh::MetaData &meta,
                                                          stk::mesh::EntityRank type,
                                                          stk::mesh::Part &part,
                                                          const Ioss::Field &io_field,
                                                          bool use_cartesian_for_scalar)
  {
    std::string name = io_field.get_name();
    stk::mesh::FieldBase *field_ptr = meta.get_field(type, name);
    // If the field has already been declared, don't redeclare it.
    if (field_ptr != NULL && stk::io::is_field_on_part(field_ptr, type, part)) {
      return field_ptr;
    }

    std::string field_type = io_field.transformed_storage()->name();
    size_t num_components = io_field.transformed_storage()->component_count();
    stk::topology::rank_t entity_rank = static_cast<stk::topology::rank_t>(type);

    if (field_type == "scalar" || num_components == 1) {
      if (!use_cartesian_for_scalar) {
        stk::mesh::Field<double> & field = meta.declare_field<stk::mesh::Field<double> >(entity_rank, name);
        stk::mesh::put_field(field, part);
        field_ptr = &field;
      } else {
        stk::mesh::Field<double, stk::mesh::Cartesian> & field =
          meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(entity_rank, name);
        stk::mesh::put_field(field, part, 1);
        field_ptr = &field;
      }
    }
    else if (field_type == "vector_2d") {
      stk::mesh::Field<double, stk::mesh::Cartesian> & field =
        meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(entity_rank, name);
      stk::mesh::put_field(field, part, 2);
      field_ptr = &field;
    }
    else if (field_type == "vector_3d") {
      stk::mesh::Field<double, stk::mesh::Cartesian> & field =
        meta.declare_field<stk::mesh::Field<double,
        stk::mesh::Cartesian> >(entity_rank, name);
      stk::mesh::put_field(field, part, 3);
      field_ptr = &field;
    }
    else if (field_type == "sym_tensor_33") {
      stk::mesh::Field<double, stk::mesh::SymmetricTensor> & field =
        meta.declare_field<stk::mesh::Field<double,
        stk::mesh::SymmetricTensor> >(entity_rank, name);
      stk::mesh::put_field(field, part, 6);
      field_ptr = &field;
    }
    else if (field_type == "full_tensor_36") {
      stk::mesh::Field<double, stk::mesh::FullTensor> & field =
        meta.declare_field<stk::mesh::Field<double,
        stk::mesh::FullTensor> >(entity_rank, name);
      stk::mesh::put_field(field, part, 9);
      field_ptr = &field;
    }
    else if (field_type == "matrix_22") {
      stk::mesh::Field<double, stk::mesh::Matrix> & field =
        meta.declare_field<stk::mesh::Field<double,
        stk::mesh::Matrix> >(entity_rank, name);
      stk::mesh::put_field(field, part, 4);
      field_ptr = &field;
    }
    else if (field_type == "matrix_33") {
      stk::mesh::Field<double, stk::mesh::Matrix> & field =
        meta.declare_field<stk::mesh::Field<double,
        stk::mesh::Matrix> >(entity_rank, name);
      stk::mesh::put_field(field, part, 9);
      field_ptr = &field;
    }
    else {
      // Just create a field with the correct number of components...
      stk::mesh::Field<double,shards::ArrayDimension> & field =
        meta.declare_field<stk::mesh::Field<double,shards::ArrayDimension> >(entity_rank, name);
      stk::mesh::put_field(field, part, num_components);
      field_ptr = &field;
    }

    if (field_ptr != NULL) {
      stk::io::set_field_role(*field_ptr, io_field.get_role());
    }
    return field_ptr;
  }

  const stk::mesh::FieldBase *declare_ioss_field(stk::mesh::MetaData &meta,
                                                 stk::mesh::EntityRank type,
                                                 stk::mesh::Part &part,
                                                 const Ioss::Field &io_field,
                                                 bool use_cartesian_for_scalar)
  {
    const stk::mesh::FieldBase *field_ptr = NULL;
    if (io_field.get_type() == Ioss::Field::INTEGER) {
      field_ptr = declare_ioss_field_internal<int>(meta, type, part, io_field, use_cartesian_for_scalar);
    } else if (io_field.get_type() == Ioss::Field::INT64) {
      field_ptr = declare_ioss_field_internal<int64_t>(meta, type, part, io_field, use_cartesian_for_scalar);
    } else if (io_field.get_type() == Ioss::Field::REAL) {
      field_ptr = declare_ioss_field_internal<double>(meta, type, part, io_field, use_cartesian_for_scalar);
    } else if (io_field.get_type() == Ioss::Field::COMPLEX) {
      field_ptr = declare_ioss_field_internal<std::complex<double> >(meta, type, part, io_field, use_cartesian_for_scalar);
    } else {
      std::ostringstream errmsg;
      errmsg << "ERROR: Unrecognized field type for IO field '"
             << io_field.get_name() << "'.";
      throw std::runtime_error(errmsg.str());
    }
    return field_ptr;
  }

  template <typename T>
  void internal_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                     const Ioss::Field &io_field,
                                     const stk::mesh::FieldBase *field,
                                     std::vector<stk::mesh::Entity> &entities,
                                     Ioss::GroupingEntity *io_entity)
  {
    size_t field_component_count = io_field.transformed_storage()->component_count();

    std::vector<T> io_field_data;
    size_t io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
    assert(io_field_data.size() == entities.size() * field_component_count);

    size_t entity_count = entities.size();

    if (io_entity_count != entity_count) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << io_field.get_name()
             << "' on " << io_entity->type_string() << " " << io_entity->name()
             << ". The IO system has " << io_entity_count
             << " entries, but the stk:mesh system has " << entity_count
             << " entries. The two counts must match.";
      throw std::runtime_error(errmsg.str());
    }

    for (size_t i=0; i < entity_count; ++i) {
      if (mesh.is_valid(entities[i])) {
        T *fld_data = static_cast<T*>(stk::mesh::field_data(*field, entities[i]));
        if (fld_data !=NULL) {
          for(size_t j=0; j<field_component_count; ++j) {
            fld_data[j] = io_field_data[i*field_component_count+j];
          }
        }
      }
    }
  }

  template <typename T>
  void internal_subsetted_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                               const Ioss::Field &io_field,
                                               const stk::mesh::FieldBase *field,
                                               std::vector<stk::mesh::Entity> &entities,
                                               Ioss::GroupingEntity *io_entity,
                                               const stk::mesh::Part *stk_part)
  {
    size_t field_component_count = io_field.transformed_storage()->component_count();
    std::vector<T> io_field_data;
    size_t io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
    assert(io_field_data.size() == entities.size() * field_component_count);
    size_t entity_count = entities.size();
    if (io_entity_count != entity_count) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << io_field.get_name()
             << "' on " << io_entity->type_string() << " " << io_entity->name()
             << ". The IO system has " << io_entity_count
             << " entries, but the stk:mesh system has " << entity_count
             << " entries. The two counts must match.";
      throw std::runtime_error(errmsg.str());
    }

    stk::mesh::MetaData &meta = stk::mesh::MetaData::get(*stk_part);
    stk::mesh::Selector selector = (meta.globally_shared_part() | meta.locally_owned_part()) & *stk_part;

    for (size_t i=0; i < entity_count; ++i) {
      if (mesh.is_valid(entities[i])) {
        const stk::mesh::Bucket &bucket = mesh.bucket(entities[i]);
        if (selector(bucket)) {
          T *fld_data = static_cast<T*>(stk::mesh::field_data(*field, entities[i]));
          if (fld_data !=NULL) {
            for(size_t j=0; j<field_component_count; ++j) {
              fld_data[j] = io_field_data[i*field_component_count+j];
            }
          }
        }
      }
    }
  }

  template <typename T>
  void internal_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                                   const Ioss::Field &io_field,
                                   const stk::mesh::FieldBase *field,
                                   std::vector<stk::mesh::Entity> &entities,
                                   Ioss::GroupingEntity *io_entity)
  {
    size_t field_component_count = io_field.transformed_storage()->component_count();
    size_t entity_count = entities.size();

    std::vector<T> io_field_data(entity_count*field_component_count);

    for (size_t i=0; i < entity_count; ++i) {
      if (mesh.is_valid(entities[i]) && mesh.entity_rank(entities[i]) == field->entity_rank()) {
        T *fld_data = static_cast<T*>(stk::mesh::field_data(*field, entities[i]));
        if (fld_data != NULL) {
          for(size_t j=0; j<field_component_count; ++j) {
            io_field_data[i*field_component_count+j] = fld_data[j];
          }
        }
      }
    }

    size_t io_entity_count = io_entity->put_field_data(io_field.get_name(), io_field_data);
    assert(io_field_data.size() == entities.size() * field_component_count);

    if (io_entity_count != entity_count) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Field count mismatch for IO field '"
             << io_field.get_name()
             << "' on " << io_entity->type_string() << " " << io_entity->name()
             << ". The IO system has " << io_entity_count
             << " entries, but the stk:mesh system has " << entity_count
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
      if (stk::io::is_valid_part_field(f, rank, part, Ioss::Field::TRANSIENT)) {
        return true;
      }
    }
    return false;
  }

  void check_if_io_part_attribute_already_defined(const stk::mesh::Part& part)
  {
      if (part.attribute<Ioss::GroupingEntity>() != NULL) {
        std::string msg = "stk::io::put_io_part_attribute( ";
        msg += part.name();
        msg += " ) FAILED:";
        msg += " io_part_attribute is already defined";
        throw std::runtime_error( msg );
      }
  }

  void put_io_part_attribute(stk::mesh::Part & part, Ioss::GroupingEntity* entity)
  {
    check_if_io_part_attribute_already_defined(part);
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    meta.declare_attribute_no_delete(part, entity);
  }
}//namespace <empty>

namespace stk {
  namespace io {

    struct IossAlternatePartName
    {
      std::string value;
    };

    void set_alternate_part_name(stk::mesh::Part& part, const std::string& altPartName)
    {
      mesh::MetaData & meta = mesh::MetaData::get(part);
      stk::mesh::PartVector pv = meta.get_parts();
      for (unsigned ii = 0; ii < pv.size(); ++ii)
        {
          if (get_alternate_part_name(*pv[ii]) == altPartName && &part != pv[ii])
            {
              throw std::runtime_error(std::string("stk::io::set_alternate_part_name found another part with the same ")
                                       +"alternate part name attribute, part= "+part.name() +" altPartName= " +altPartName+" other part= " +pv[ii]->name());
            }
        }
      const IossAlternatePartName *altName = part.attribute<IossAlternatePartName>();
      if (!altName)
        {
          IossAlternatePartName *altName1 = new IossAlternatePartName();
          altName1->value = altPartName;
          meta.declare_attribute_with_delete(part, altName1);
        }
      else
        {
          if (altPartName != altName->value)
            {
              bool success = meta.remove_attribute(part, altName);
              if (!success)
                throw std::runtime_error(std::string("stk::io::set_alternate_part_name failed to remove")
                                         +"alternate part name attribute, part= "+part.name() +" altPartName= " +altPartName);
              IossAlternatePartName *altName1 = new IossAlternatePartName();
              altName1->value = altPartName;
              meta.declare_attribute_with_delete(part, altName1);
            }
        }
    }

    std::string get_alternate_part_name(stk::mesh::Part& part)
    {
      const IossAlternatePartName *altName = part.attribute<IossAlternatePartName>();
      if (altName)
        {
          return altName->value;
        }
      else
        return "";
    }

    std::string getPartName(stk::mesh::Part& part)
    {
      std::string apn = get_alternate_part_name(part);
      if (apn.length())
        {
          return apn;
        }
      else
        return part.name();
    }

    stk::mesh::Part *getPart(const stk::mesh::MetaData& meta_data, const std::string& name)
    {
      const mesh::PartVector & parts = meta_data.get_parts();
      for (unsigned ii=0; ii < parts.size(); ++ii)
        {
          stk::mesh::Part *pp = parts[ii];
          std::string altName = getPartName(*pp);
          if (altName == name)
            return pp;
        }
      return 0;
    }

    size_t db_api_int_size(const Ioss::GroupingEntity *entity)
    {
      return entity->get_database()->int_byte_size_api();
    }

    void initialize_spatial_dimension(stk::mesh::MetaData & meta, size_t spatial_dimension,
                                      const std::vector<std::string> &entity_rank_names)
    {
      if (!meta.is_initialized() ) {
        meta.initialize(spatial_dimension, entity_rank_names);
      }
    }

    bool is_field_on_part(const stk::mesh::FieldBase *field,
                          const stk::mesh::EntityRank part_type,
                          const stk::mesh::Part &part)
    {
      const stk::mesh::MetaData &meta = stk::mesh::MetaData::get(part);
      const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*field, part_type, part);
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
        if (part_type != stk::topology::NODE_RANK || part == meta.universal_part()) {
          return true;
        }

        const stk::mesh::FieldBase::Restriction &res_universe = stk::mesh::find_restriction(*field, part_type, meta.universal_part());
        if (res_universe.num_scalars_per_entity() <= 0) {
          // Field exists on current part, but not on the universal
          // set (and this part is not the universal part)
          return true;
        }
      }
      return false;
    }

    bool is_valid_part_field(const stk::mesh::FieldBase *field,
                             const stk::mesh::EntityRank part_type,
                             const stk::mesh::Part &part,
                             const Ioss::Field::RoleType filter_role)
    {
      const Ioss::Field::RoleType *role = stk::io::get_field_role(*field);

      if (role == NULL) {
        return false;
      }

      if (role != NULL && *role != filter_role)
        return false;

      return is_field_on_part(field, part_type, part);
    }

    void get_io_field_type(const stk::mesh::FieldBase *field,
                           const stk::mesh::FieldRestriction &res,
                           std::pair<std::string, Ioss::Field::BasicType> *result)
    {
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

#if 0
      // Not currently handled...
      static const std::string quaternion_2d("quaternion_2d");
      static const std::string quaternion_3d("quaternion_3d");
      static const std::string sym_tensor_13("sym_tensor_13");
      static const std::string sym_tensor_11("sym_tensor_11");
      static const std::string sym_tensor_10("sym_tensor_10");
      static const std::string asym_tensor_03("asym_tensor_03");
      static const std::string asym_tensor_02("asym_tensor_02");
      static const std::string asym_tensor_01("asym_tensor_01");
#endif

      const unsigned rank = field->field_array_rank();
      const shards::ArrayDimTag * const * const tags = field->dimension_tags();

      result->second = Ioss::Field::INVALID;

      if ( field->type_is<double>() ) {
        result->second = Ioss::Field::REAL;
      }
      else if ( field->type_is<int>() ) {
        result->second = Ioss::Field::INTEGER;
      }

      if ( 0 == rank ) {
        result->first = scalar ;
      }
      else if ( 1 == rank ) {
        size_t num_comp = res.num_scalars_per_entity();
        if ( tags[0] == & stk::mesh::Cartesian::tag() ) {
          if (1 == num_comp ) {
            result->first = scalar ;
          }
          else if ( 2 == num_comp ) {
            result->first = vector_2d ;
          }
          else if ( 3 == num_comp ) {
            result->first = vector_3d ;
          }
        }
        else if ( tags[0] == & stk::mesh::FullTensor::tag() ) {
          if ( 9 == num_comp ) {
            result->first = full_tensor_36 ;
          }
          else if ( 5 == num_comp ) {
            result->first = full_tensor_32 ;
          }
          else if ( 4 == num_comp ) {
            result->first = full_tensor_22 ;
          }
          else if ( 3 == num_comp ) {
            result->first = full_tensor_12 ;
          }
        }
        else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() ) {
          if ( 6 == num_comp ) {
            result->first = sym_tensor_33 ;
          }
          else if ( 4 == num_comp ) {
            result->first = sym_tensor_31 ;
          }
          else if ( 3 == num_comp ) {
            result->first = sym_tensor_21 ;
          }
        }
        else if ( tags[0] == & stk::mesh::Matrix::tag() ) {
          if (4 == num_comp ) {
            result->first =  matrix_22;
          }
          else if ( 9 == num_comp ) {
            result->first = matrix_33 ;
          }
        }
      }

      if ( result->first.empty() ) {
        size_t num_comp = res.num_scalars_per_entity();
        std::ostringstream tmp ;
        tmp << "Real[" << num_comp << "]" ;
        result->first = tmp.str();
      }
    }

    //----------------------------------------------------------------------
    const Ioss::GroupingEntity *get_associated_ioss_entity(const mesh::Part &part)
    {
      const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
      if (!entity || entity->type() == Ioss::INVALID_TYPE) {
        return NULL;
      } else {
        return entity;
      }
    }

    void put_io_part_attribute(mesh::Part & part)
    {
      check_if_io_part_attribute_already_defined(part);
      mesh::MetaData & meta = mesh::MetaData::get(part);
      Ioss::GroupingEntity *attr = new Ioss::NullEntity();
      meta.declare_attribute_with_delete(part, attr);
    }

    void remove_io_part_attribute(mesh::Part & part)
    {
      const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
      if (entity != NULL) {
        mesh::MetaData & meta = mesh::MetaData::get(part);
        bool success = meta.remove_attribute(part, entity);
        if (!success) {
          std::string msg = "stk::io::remove_io_part_attribute( ";
          msg += part.name();
          msg += " ) FAILED:";
          msg += " meta.remove_attribute(..) returned failure.";
          throw std::runtime_error( msg );
        }

        if (entity->type() == Ioss::INVALID_TYPE) {
          delete entity;
        }
      }
    }

    stk::topology get_start_topology(const Ioss::ElementTopology* topology, unsigned mesh_spatial_dimension)
    {
        if (topology->is_element() && topology->spatial_dimension() == (int)mesh_spatial_dimension)
        {
            return stk::topology::BEGIN_ELEMENT_RANK;
        }
        return stk::topology::BEGIN_TOPOLOGY;
    }

    stk::topology map_ioss_topology_to_stk(const Ioss::ElementTopology *topology,
                                           unsigned mesh_spatial_dimension)
    {
      stk::topology begin_topo = get_start_topology(topology, mesh_spatial_dimension);
      for (stk::topology topo=begin_topo; topo < stk::topology::END_TOPOLOGY; ++topo) {
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
      std::string tmpCopy = Ioss::Utils::lowercase(topology->name().substr(0,5));
      if (tmpCopy=="super")
        {
          return stk::create_superelement_topology(topology->number_nodes());
        }

      return stk::topology::INVALID_TOPOLOGY;
    }

    std::string map_stk_topology_to_ioss(stk::topology topo)
    {
      Ioss::ElementTopology *ioss_topo = Ioss::ElementTopology::factory(topo.name(), true);
      return ioss_topo != NULL ? ioss_topo->name() : "invalid";
    }

    void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta)
    {
      if (include_entity(entity)) {
        mesh::EntityRank type = get_entity_rank(entity, meta);
        stk::mesh::Part & part = meta.declare_part(entity->name(), type);
        if (entity->property_exists("id")) {
          meta.set_part_id(part, entity->get_property("id").get_int());
        }
        ::put_io_part_attribute(part, entity);
      }
    }

    void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta)
    {
      if (include_entity(entity)) {
        mesh::EntityRank type = get_entity_rank(entity, meta);
        stk::mesh::Part * part = NULL;
        part = &meta.declare_part(entity->name(), type);
        if (entity->property_exists("id")) {
            meta.set_part_id(*part, entity->get_property("id").get_int());
        }
        ::put_io_part_attribute(*part, entity);

        const Ioss::ElementTopology *topology = entity->topology();
        // Check spatial dimension of the element topology here so we can
        // issue a more meaningful error message.  If the dimension is bad
        // and we continue to the following calls, there is an exception
        // and we get unintelligible (to the user) error messages.  Could
        // also do a catch...

        if (entity->type() == Ioss::ELEMENTBLOCK) {
          assert(topology != NULL);
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
        }

        stk::topology stk_topology = map_ioss_topology_to_stk(topology, meta.spatial_dimension());
        stk::mesh::set_topology(*part, stk_topology);
        stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE, *part, type);
      }
    }

    //----------------------------------------------------------------------
    /** Add all stk::Fields on the entities of the specified part_type
     *  on the specified part of the specified role * to the specified
     *  Ioss::GroupingEntity
     */

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank part_type,
                         Ioss::GroupingEntity *entity,
                         const std::vector<FieldAndName> &namedFields,
                         const Ioss::Field::RoleType filter_role)
    {
      stk::mesh::EntityRank part_rank = part_primary_entity_rank(part);
      const stk::mesh::PartVector &blocks = part.subsets();
      bool check_subparts = (part_rank == 1 || part_rank == 2) && (blocks.size() > 0);

      for (size_t i=0;i<namedFields.size();i++)
        {
          const stk::mesh::FieldBase *f = namedFields[i].field();
          if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
            const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, part_type, part);
            std::pair<std::string, Ioss::Field::BasicType> field_type;
            get_io_field_type(f, res, &field_type);
            if (field_type.second != Ioss::Field::INVALID) {
              size_t entity_size = entity->get_property("entity_count").get_int();
              std::string name = namedFields[i].db_name();
              entity->field_add(Ioss::Field(name, field_type.second, field_type.first,
                                            filter_role, entity_size));
	      if (entity->type() == Ioss::NODEBLOCK) {
		namedFields[i].m_forceNodeblockOutput = true;
	      }
            }
          }
          // If this is a sideset, check the subset parts for the field also...
          if (check_subparts) {
            for (size_t j = 0; j < blocks.size(); j++) {
              mesh::Part & side_block_part = *blocks[j];
              if (stk::io::is_valid_part_field(f, part_type, side_block_part, filter_role)) {
                const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, part_type, side_block_part);
                std::pair<std::string, Ioss::Field::BasicType> field_type;
                get_io_field_type(f, res, &field_type);
                if (field_type.second != Ioss::Field::INVALID) {
                  size_t entity_size = entity->get_property("entity_count").get_int();
                  std::string name = namedFields[i].db_name();
                  entity->field_add(Ioss::Field(name, field_type.second, field_type.first,
                                                filter_role, entity_size));
		  if (entity->type() == Ioss::NODEBLOCK) {
		    namedFields[i].m_forceNodeblockOutput = true;
		  }
                }
              }
            }
          }
        }
    }

    void getNamedFields(const stk::mesh::MetaData &meta, Ioss::GroupingEntity *io_entity, std::vector<FieldAndName> &namedFields)
    {
      const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();
      namedFields.reserve(fields.size());
      std::vector<stk::mesh::FieldBase *>::const_iterator fieldIterator = fields.begin();
      for(;fieldIterator != fields.end();++fieldIterator)
        {
          std::string field_name = (*fieldIterator)->name();
          namedFields.push_back(FieldAndName(*fieldIterator, field_name));
        }
    }

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank part_type,
                         Ioss::GroupingEntity *entity,
                         const Ioss::Field::RoleType filter_role)
    {
      std::vector<FieldAndName> namedFields;
      stk::io::getNamedFields(mesh::MetaData::get(part), entity, namedFields);

      ioss_add_fields(part, part_type, entity, namedFields, filter_role);
    }

    void ioss_add_fields(const stk::mesh::Part &part,
                         const stk::mesh::EntityRank part_type,
                         Ioss::GroupingEntity *entity,
                         const std::vector<FieldAndName> &namedFields)
    {
      ioss_add_fields(part, part_type, entity, namedFields, Ioss::Field::Field::TRANSIENT);
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
                          stk::mesh::EntityRank part_type)
    {
      stk::mesh::MetaData &meta = mesh::MetaData::get(part);

      bool use_cartesian_for_scalar = false;
      if (role == Ioss::Field::ATTRIBUTE)
        use_cartesian_for_scalar = true;

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
        Ioss::Field io_field = entity->get_field(*I);
        declare_ioss_field(meta, part_type, part, io_field, use_cartesian_for_scalar);
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
      delete_selector_property(region.get_nodesets());
      delete_selector_property(region.get_commsets());

      const Ioss::SideSetContainer& side_sets = region.get_sidesets();
      for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
          it != side_sets.end(); ++it) {
        Ioss::SideSet *sset = *it;
        delete_selector_property(*it);
        delete_selector_property(sset->get_side_blocks());
      }
    }

    void delete_selector_property(Ioss::GroupingEntity *io_entity)
    {
      // If the Ioss::GroupingEntity has a property named 'selector' of
      // type 'pointer', delete the pointer and remove the property.
      if (io_entity->property_exists(internal_selector_name)) {
        mesh::Selector *select = reinterpret_cast<mesh::Selector*>(io_entity->get_property(internal_selector_name).get_pointer());
        delete select;
        io_entity->property_erase(internal_selector_name);
      }
    }



    template <typename INT>
    void get_entity_list(Ioss::GroupingEntity *io_entity,
                         stk::mesh::EntityRank part_type,
                         const stk::mesh::BulkData &bulk,
                         std::vector<stk::mesh::Entity> &entities)
    {
      if (io_entity->type() == Ioss::SIDEBLOCK) {
	std::vector<INT> elem_side ;
	io_entity->get_field_data("element_side", elem_side);
	size_t side_count = elem_side.size() / 2;
	for(size_t is=0; is<side_count; ++is)
	  entities.push_back(stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk, elem_side[is*2], elem_side[is*2+1]-1, part_type));
      }
      else {
	std::vector<INT> ids ;
	io_entity->get_field_data("ids", ids);

	size_t count = ids.size();
	entities.reserve(count);

	for(size_t i=0; i<count; ++i) {
	  entities.push_back(bulk.get_entity( part_type, ids[i] ));
	}
      }
    }
      

    void get_entity_list(Ioss::GroupingEntity *io_entity,
                         stk::mesh::EntityRank part_type,
                         const stk::mesh::BulkData &bulk,
                         std::vector<stk::mesh::Entity> &entities)
    {
      if (io_entity->get_database()->is_input()) {
        if (db_api_int_size(io_entity) == 4) {
          get_entity_list<int>(io_entity, part_type, bulk, entities);
        } else {
          get_entity_list<int64_t>(io_entity, part_type, bulk, entities);
        }
      } else {
        // Output database...
        assert(io_entity->property_exists(internal_selector_name));

        mesh::Selector *select = reinterpret_cast<mesh::Selector*>(io_entity->get_property(internal_selector_name).get_pointer());
        get_selected_entities(*select, bulk.buckets(part_type), entities);
      }
    }

    const std::string get_suffix_for_field_at_state(enum stk::mesh::FieldState field_state)
    {
      std::string suffix = "";
      switch(field_state)
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
          ThrowRequireMsg(false, "Internal Error: Unsupported stk::mesh::FieldState: " << field_state << ".\n");
        }
      return suffix;
    }

    std::string get_stated_field_name(const std::string &field_base_name, stk::mesh::FieldState state_identifier)
    {
      std::string field_name_with_suffix = field_base_name + get_suffix_for_field_at_state(state_identifier);
      return field_name_with_suffix;
    }

    bool field_state_exists_on_io_entity(const std::string& db_name, const stk::mesh::FieldBase* field, stk::mesh::FieldState state_identifier,
                                         Ioss::GroupingEntity *io_entity)
    {
        std::string field_name_with_suffix = get_stated_field_name(db_name, state_identifier);
        return io_entity->field_exists(field_name_with_suffix);
    }

    bool all_field_states_exist_on_io_entity(const std::string& db_name, const stk::mesh::FieldBase* field, Ioss::GroupingEntity *io_entity,
                                             std::vector<stk::mesh::FieldState> &missing_states)
    {
        bool all_states_exist = true;
        size_t state_count = field->number_of_states();
        for(size_t state = 0; state < state_count - 1; state++)
        {
            stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
            if (!field_state_exists_on_io_entity(db_name, field, state_identifier, io_entity))
            {
                all_states_exist = false;
                missing_states.push_back(state_identifier);
            }
        }

        return all_states_exist;
    }

    void multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                         const stk::mesh::FieldBase *field,
                                         std::vector<stk::mesh::Entity> &entity_list,
                                         Ioss::GroupingEntity *io_entity,
                                         const std::string &name,
                                         const size_t state_count,
                                         bool ignore_missing_fields)
    {
        for(size_t state = 0; state < state_count - 1; state++)
        {
            stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
            bool field_exists = field_state_exists_on_io_entity(name, field, state_identifier, io_entity);
            if (!field_exists && !ignore_missing_fields)
            {
                STKIORequire(field_exists);
            }
            if (field_exists)
            {
                stk::mesh::FieldBase *stated_field = field->field_state(state_identifier);
                std::string field_name_with_suffix = get_stated_field_name(name, state_identifier);
                stk::io::field_data_from_ioss(mesh, stated_field, entity_list, io_entity, field_name_with_suffix);
            }
        }
    }

    void subsetted_multistate_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                                   const stk::mesh::FieldBase *field,
                                                   std::vector<stk::mesh::Entity> &entity_list,
                                                   Ioss::GroupingEntity *io_entity,
                                                   const stk::mesh::Part *stk_part,
                                                   const std::string &name,
                                                   const size_t state_count,
                                                   bool ignore_missing_fields)
    {
        for(size_t state = 0; state < state_count - 1; state++)
        {
            stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
            bool field_exists = field_state_exists_on_io_entity(name, field, state_identifier, io_entity);
            if (!field_exists && !ignore_missing_fields)
            {
                STKIORequire(field_exists);
            }
            if (field_exists)
            {
                stk::mesh::FieldBase *stated_field = field->field_state(state_identifier);
                std::string field_name_with_suffix = get_stated_field_name(name, state_identifier);
                stk::io::subsetted_field_data_from_ioss(mesh, stated_field, entity_list,
                                                      io_entity, stk_part, field_name_with_suffix);
            }
        }
    }

    void field_data_from_ioss(const stk::mesh::BulkData& mesh,
                              const stk::mesh::FieldBase *field,
                              std::vector<stk::mesh::Entity> &entities,
                              Ioss::GroupingEntity *io_entity,
                              const std::string &io_fld_name)
    {
      /// \todo REFACTOR Need some additional compatibility checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != NULL && io_entity->field_exists(io_fld_name)) {
        const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
        if (field->type_is<double>()) {
          internal_field_data_from_ioss<double>(mesh, io_field, field, entities, io_entity);
        } else if (field->type_is<int>()) {
          // Make sure the IO field type matches the STK field type.
          // By default, all IO fields are created of type 'double'
          if (db_api_int_size(io_entity) == 4) {
            io_field.check_type(Ioss::Field::INTEGER);
            internal_field_data_from_ioss<int>(mesh, io_field, field, entities, io_entity);
          } else {
            io_field.check_type(Ioss::Field::INT64);
            internal_field_data_from_ioss<int64_t>(mesh, io_field, field, entities, io_entity);
          }
        }
      }
    }

    void subsetted_field_data_from_ioss(const stk::mesh::BulkData& mesh,
                                        const stk::mesh::FieldBase *field,
                                        std::vector<stk::mesh::Entity> &entities,
                                        Ioss::GroupingEntity *io_entity,
                                        const stk::mesh::Part *stk_part,
                                        const std::string &io_fld_name)
    {
      /// \todo REFACTOR Need some additional compatibility checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != NULL && io_entity->field_exists(io_fld_name)) {
        const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
        if (field->type_is<double>()) {
          internal_subsetted_field_data_from_ioss<double>(mesh, io_field, field, entities, io_entity, stk_part);
        } else if (field->type_is<int>()) {
          // Make sure the IO field type matches the STK field type.
          // By default, all IO fields are created of type 'double'
          if (db_api_int_size(io_entity) == 4) {
            io_field.check_type(Ioss::Field::INTEGER);
            internal_subsetted_field_data_from_ioss<int>(mesh, io_field, field, entities, io_entity, stk_part);
          } else {
            io_field.check_type(Ioss::Field::INT64);
            internal_subsetted_field_data_from_ioss<int64_t>(mesh, io_field, field, entities, io_entity,
                                                    stk_part);
          }
        }
      }
    }

    void multistate_field_data_to_ioss(const stk::mesh::BulkData& mesh,
                                       const stk::mesh::FieldBase *field,
                                       std::vector<stk::mesh::Entity> &entities,
                                       Ioss::GroupingEntity *io_entity,
                                       const std::string &io_fld_name,
                                       Ioss::Field::RoleType filter_role,
                                       const size_t state_count)
    {
      for(size_t state = 0; state < state_count - 1; state++)
        {
          stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
          std::string field_name_with_suffix = get_stated_field_name(io_fld_name, state_identifier);
          stk::mesh::FieldBase *stated_field = field->field_state(state_identifier);
          //STKIORequire(io_entity->field_exists(field_name_with_suffix));
          stk::io::field_data_to_ioss(mesh, stated_field, entities, io_entity, field_name_with_suffix, filter_role);
        }
    }

    void field_data_to_ioss(const stk::mesh::BulkData& mesh,
                            const stk::mesh::FieldBase *field,
                            std::vector<stk::mesh::Entity> &entities,
                            Ioss::GroupingEntity *io_entity,
                            const std::string &io_fld_name,
                            Ioss::Field::RoleType filter_role)
    {
      /// \todo REFACTOR Need some additional compatibility checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != NULL && io_entity->field_exists(io_fld_name)) {
        const Ioss::Field &io_field = io_entity->get_fieldref(io_fld_name);
        if (io_field.get_role() == filter_role) {
          if (field->type_is<double>()) {
            internal_field_data_to_ioss<double>(mesh, io_field, field, entities, io_entity);
          } else if (field->type_is<int>()) {
            io_field.check_type(Ioss::Field::INTEGER);
            internal_field_data_to_ioss<int>(mesh, io_field, field, entities, io_entity);
          } else if (field->type_is<int64_t>()) {
            io_field.check_type(Ioss::Field::INT64);
            internal_field_data_to_ioss<int64_t>(mesh, io_field, field, entities, io_entity);
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

    namespace {

      void define_side_block(const stk::mesh::BulkData &bulk,
                             stk::mesh::Selector selector,
                             const stk::mesh::Selector *subset_selector,
                             stk::mesh::Part &part,
                             Ioss::SideSet *sset,
                             int spatial_dimension)
      {
        stk::mesh::EntityRank type = part.primary_entity_rank();
        const stk::mesh::EntityRank siderank = stk::mesh::MetaData::get(part).side_rank();
        const stk::mesh::EntityRank edgerank = stk::topology::EDGE_RANK;
        STKIORequire(type == siderank || type == edgerank);

        stk::topology side_topology = part.topology();
        std::string io_topo = map_stk_topology_to_ioss(side_topology);
        std::string element_topo_name = "unknown";

        // Get sideblock parent element topology quantities...
        // Try to decode from part name...
        std::vector<std::string> tokens;
        stk::util::tokenize(getPartName(part), "_", tokens);
        const Ioss::ElementTopology *element_topo = nullptr;
        stk::topology stk_element_topology = stk::topology::INVALID_TOPOLOGY;
        if (tokens.size() >= 4) {
          // Name of form: "name_eltopo_sidetopo_id" or
          //               "name_block_id_sidetopo_id"
          // "name" is typically "surface".
          element_topo = Ioss::ElementTopology::factory(tokens[1], true);

          if (element_topo != NULL) {
            element_topo_name = element_topo->name();
            stk_element_topology = map_ioss_topology_to_stk(element_topo, bulk.mesh_meta_data().spatial_dimension());
          }
        }

        int sset_id = -1;
        if(sset->property_exists("id"))
        {
            sset_id = sset->get_property("id").get_int();
        }

        size_t side_count = get_number_sides_in_sideset(bulk, sset_id, selector, subset_selector, stk_element_topology, bulk.buckets(type));

        std::string name = getPartName(part);
        Ioss::SideBlock *side_block = sset->get_side_block(name);
        if(side_block == nullptr)
        {
            side_block = new Ioss::SideBlock(sset->get_database(), name, io_topo, element_topo_name, side_count);
            sset->add(side_block);
        }

        const mesh::FieldBase *df = get_distribution_factor_field(part);
        if (df != NULL) {
          int nodes_per_side = side_topology.num_nodes();
          std::string storage_type = "Real[";
          storage_type += Ioss::Utils::to_string(nodes_per_side);
          storage_type += "]";
          side_block->field_add(Ioss::Field("distribution_factors", Ioss::Field::REAL, storage_type,
                                            Ioss::Field::MESH, side_count));
        }

        selector &= bulk.mesh_meta_data().locally_owned_part();

        mesh::Selector *select = new mesh::Selector(selector);
        side_block->property_add(Ioss::Property(internal_selector_name, select, false));
        side_block->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), side_block, Ioss::Field::ATTRIBUTE);
      }

      void define_side_blocks(stk::mesh::Part &part,
                              const stk::mesh::BulkData &bulk_data,
                              Ioss::SideSet *sset,
                              stk::mesh::EntityRank type,
                              int spatial_dimension,
                              const stk::mesh::Selector *subset_selector)
      {
        STKIORequire(type == stk::topology::FACE_RANK || stk::topology::EDGE_RANK);

        const stk::mesh::PartVector &blocks = part.subsets();
        if (blocks.size() > 0) {
          for (size_t j = 0; j < blocks.size(); j++) {
            mesh::Part & side_block_part = *blocks[j];
            mesh::Selector selector = side_block_part;
            if (subset_selector) selector &= *subset_selector;

            define_side_block(bulk_data, selector, subset_selector, side_block_part,
                              sset, spatial_dimension);
          }
        } else {
          mesh::Selector selector = part;
          if (subset_selector) selector &= *subset_selector;
          define_side_block(bulk_data, selector, subset_selector, part, sset, spatial_dimension);
        }
      }

      //----------------------------------------------------------------------
      void define_node_block(stk::mesh::Part &part,
                             const stk::mesh::BulkData &bulk,
                             Ioss::Region &io_region,
                             const stk::mesh::Selector *subset_selector)
      {
        //--------------------------------
        // Set the spatial dimension:
        mesh::MetaData & meta = mesh::MetaData::get(part);

        //We now get spatial-dim from meta.spatial_dimension() rather than getting
        //it from the coordinate-field's restriction onto the universal part.
        //This is because some codes (sierra framework) don't put the coordinate
        //field on the universal part. (framework puts it on active and inactive parts)
        const int spatial_dim = meta.spatial_dimension();

        //--------------------------------
        // Create the special universal node block:

        mesh::Selector all_selector = meta.globally_shared_part() | meta.locally_owned_part();
        if (subset_selector) all_selector &= *subset_selector;

        mesh::Selector own_selector = meta.locally_owned_part();
        if (subset_selector) own_selector &= *subset_selector;

        int64_t all_nodes = count_selected_entities(all_selector, bulk.buckets(stk::topology::NODE_RANK));
        int64_t own_nodes = count_selected_entities(own_selector, bulk.buckets(stk::topology::NODE_RANK));

        const std::string name("nodeblock_1");

        Ioss::NodeBlock * nb = io_region.get_node_block(name);
        if(nb == nullptr)
        {
            nb = new Ioss::NodeBlock(io_region.get_database(),
                                                         name, all_nodes, spatial_dim);
            io_region.add( nb );
        }

        mesh::Selector *node_select = new mesh::Selector(all_selector);
        nb->property_add(Ioss::Property(internal_selector_name, node_select, false));
        nb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add locally-owned property...
        nb->property_add(Ioss::Property("locally_owned_count", own_nodes));
        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), nb, Ioss::Field::ATTRIBUTE);
      }

      void define_node_set(stk::mesh::Part &part,
                           const std::string &name,
                           const stk::mesh::BulkData &bulk,
                           Ioss::Region &io_region,
                           const stk::mesh::Selector *subset_selector)
      {
        mesh::MetaData & meta = mesh::MetaData::get(part);

        mesh::Selector all_selector = (meta.globally_shared_part() | meta.locally_owned_part()) & part;
        if (subset_selector) all_selector &= *subset_selector;

        mesh::Selector own_selector = meta.locally_owned_part() & part;
        if (subset_selector) own_selector &= *subset_selector;

        int64_t all_nodes = count_selected_entities(all_selector, bulk.buckets(stk::topology::NODE_RANK));
        int64_t own_nodes = count_selected_entities(own_selector, bulk.buckets(stk::topology::NODE_RANK));

        Ioss::NodeSet *ns = io_region.get_nodeset(name);
        if(ns == nullptr)
        {
            ns = new Ioss::NodeSet( io_region.get_database(), name, all_nodes);
            io_region.add(ns);
        }

        ns->property_add(Ioss::Property("locally_owned_count", own_nodes));

        mesh::Selector *select = new mesh::Selector(all_selector);
        ns->property_add(Ioss::Property(internal_selector_name, select, false));
        ns->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        if (part.id() != stk::mesh::Part::INVALID_ID)
        {
            ns->property_add(Ioss::Property("id", part.id()));
        }

        // Add the attribute fields.
        ioss_add_fields(part, stk::topology::NODE_RANK, ns, Ioss::Field::ATTRIBUTE);
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

      void define_element_block(stk::mesh::Part &part,
                                const stk::mesh::BulkData &bulk,
                                const std::vector<std::vector<int>> &attributeOrdering,
                                Ioss::Region &io_region,
                                const stk::mesh::Selector *subset_selector,
                                bool use_nodeset_for_nodal_fields,
                                bool order_blocks_by_creation_order)
      {
        mesh::MetaData & meta = mesh::MetaData::get(part);

        stk::topology topo = part.topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR: Part " << part.name() << " returned INVALID from get_topology()";
          throw std::runtime_error( msg.str() );
        }

        mesh::Selector selector = meta.locally_owned_part() & part;
        if (subset_selector) selector &= *subset_selector;

        const size_t num_elems = count_selected_entities( selector, bulk.buckets(stk::topology::ELEMENT_RANK));

        // Defer the counting of attributes until after we define the
        // element block so we can count them as we add them as fields to
        // the element block
        std::string name = getPartName(part);
        Ioss::ElementBlock *eb = io_region.get_element_block(name);
        if(eb == nullptr)
        {
            eb = new Ioss::ElementBlock(io_region.get_database() ,
                                                        name,
                                                        map_stk_topology_to_ioss(part.topology()),
                                                        num_elems);
            io_region.add(eb);
        }

        if (order_blocks_by_creation_order)
        {
            int ordinal = part.mesh_meta_data_ordinal();
            eb->property_add(Ioss::Property("original_block_order", ordinal));
        }

        if (part.id() != stk::mesh::Part::INVALID_ID)
        {
            eb->property_add(Ioss::Property("id", part.id()));
        }


        mesh::Selector *select = new mesh::Selector(selector);
        eb->property_add(Ioss::Property(internal_selector_name, select, false));
        eb->property_add(Ioss::Property(base_stk_part_name, getPartName(part)));

        // Add the attribute fields.
        ioss_add_fields(part, part_primary_entity_rank(part), eb, Ioss::Field::ATTRIBUTE);

        set_attribute_field_order(bulk, attributeOrdering, part, eb);

        // Check whether there are any transient fields defined on the nodes of this elementblock
        // that are to be output.  If so, create a nodeset named "part.name()"+block_nodes_suffix
        // and output the fields on that nodeset...
        if (use_nodeset_for_nodal_fields &&
            will_output_lower_rank_fields(part, stk::topology::NODE_RANK)) {
          std::string nodes_name = getPartName(part) + block_nodes_suffix;
          define_node_set(part, nodes_name, bulk, io_region, subset_selector);
        }
      }

      void define_communication_maps(const stk::mesh::BulkData &bulk,
                                     Ioss::Region &io_region,
                                     const stk::mesh::Selector *subset_selector)
      {
        if (bulk.parallel_size() > 1) {
          const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);
          const std::string cs_name("node_symm_comm_spec");

          mesh::Selector selector = meta.globally_shared_part();
          if (subset_selector) selector &= *subset_selector;

          std::vector<mesh::Entity> entities;
          get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), entities);

          std::vector<int> sharingProcs;
          size_t size = 0;
          for (size_t i=0; i < entities.size(); i++) {
            bulk.comm_shared_procs(bulk.entity_key(entities[i]), sharingProcs);
            size+=sharingProcs.size();
          }

          Ioss::DatabaseIO *dbo = io_region.get_database();
          Ioss::CommSet *io_cs = new Ioss::CommSet(dbo, cs_name, "node", size);
          io_region.add(io_cs);

          mesh::Selector *select = new mesh::Selector(selector);
          io_cs->property_add(Ioss::Property(internal_selector_name, select, false));

          // Update global node and element count...
          std::vector<size_t> entityCounts;
          stk::mesh::comm_mesh_counts(bulk, entityCounts);

          io_region.property_add(Ioss::Property("global_node_count",    static_cast<int64_t>(entityCounts[stk::topology::NODE_RANK])));
          io_region.property_add(Ioss::Property("global_element_count", static_cast<int64_t>(entityCounts[stk::topology::ELEMENT_RANK])));
        }
      }

      void define_side_set(stk::mesh::Part &part,
                           const stk::mesh::BulkData &bulk,
                           Ioss::Region &io_region,
                           const stk::mesh::Selector *subset_selector,
                           bool use_nodeset_for_nodal_fields)
      {
        const stk::mesh::EntityRank si_rank = mesh::MetaData::get(part).side_rank();

        bool create_sideset = true;
        if (part.subsets().empty()) {
          // Only define a sideset for this part if its superset part is
          // not a side-containing part..  (i.e., this part is not a subset part
          // in a surface...)
          const stk::mesh::PartVector &supersets = part.supersets();
          for (size_t i=0; i < supersets.size(); i++) {
            if (is_part_io_part(*supersets[i]) &&
                (supersets[i]->primary_entity_rank() == stk::topology::FACE_RANK ||
                 supersets[i]->primary_entity_rank() == stk::topology::EDGE_RANK)) {
              create_sideset = false;
              break;
            }
          }
        }

        if (create_sideset) {
          std::string name = getPartName(part);
          Ioss::SideSet *ss = io_region.get_sideset(name);
          if(ss == nullptr)
          {
              ss = new Ioss::SideSet(io_region.get_database(), name);
              io_region.add(ss);
          }

          if (part.id() != stk::mesh::Part::INVALID_ID)
          {
              ss->property_add(Ioss::Property("id",part.id()));
          }

          int spatial_dim = io_region.get_property("spatial_dimension").get_int();
          define_side_blocks(part, bulk, ss, si_rank, spatial_dim, subset_selector);

          if (use_nodeset_for_nodal_fields) {
            bool lower_rank_fields = will_output_lower_rank_fields(part, stk::topology::NODE_RANK);
            if (!lower_rank_fields) {
              // See if lower rank fields are defined on sideblock parts of this sideset...
              const stk::mesh::PartVector &blocks = part.subsets();
              for (size_t j = 0; j < blocks.size() && !lower_rank_fields; j++) {
                mesh::Part & side_block_part = *blocks[j];
                lower_rank_fields = will_output_lower_rank_fields(side_block_part, stk::topology::NODE_RANK);
              }
            }
            if (lower_rank_fields) {
              std::string nodes_name = getPartName(part) + block_nodes_suffix;
              define_node_set(part, nodes_name, bulk, io_region, subset_selector);
            }
          }
        }
      }

    } // namespace <blank>

    struct part_compare_by_name {
      bool operator() (stk::mesh::Part *i, stk::mesh::Part *j) { return (i->name() < j->name()); }
    };

    void define_output_db_within_state_define(Ioss::Region & io_region ,
                                              const mesh::BulkData &bulk_data,
                                              const std::vector<std::vector<int>> &attributeOrdering,
                                              const Ioss::Region *input_region = nullptr,
                                              const stk::mesh::Selector *subset_selector = nullptr,
                                              const bool sort_stk_parts_by_name = false,
                                              const bool use_nodeset_for_part_node_fields = true)
    {
       const mesh::MetaData & meta_data = mesh::MetaData::get(bulk_data);
       define_node_block(meta_data.universal_part(), bulk_data, io_region, subset_selector);

       // All parts of the meta data:
       const mesh::PartVector *parts = NULL;
       mesh::PartVector all_parts_sorted;

       const mesh::PartVector & all_parts = meta_data.get_parts();
       // sort parts so they go out the same on all processors (srk: this was induced by streaming refine)
       if (sort_stk_parts_by_name) {
         all_parts_sorted = all_parts;
         std::sort(all_parts_sorted.begin(), all_parts_sorted.end(), part_compare_by_name());
         parts = &all_parts_sorted;
       } else {
         parts = &all_parts;
       }

       const bool order_blocks_by_creation_order = (input_region == NULL) && !sort_stk_parts_by_name;

       for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {
         mesh::Part * const part = *i;

         if (is_part_io_part(*part)) {
           if (part->primary_entity_rank() == mesh::InvalidEntityRank)
             continue;
           else if (part->primary_entity_rank() == stk::topology::NODE_RANK)
             define_node_set(*part, getPartName(*part), bulk_data, io_region, subset_selector);
           else if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK) {
             define_element_block(*part, bulk_data, attributeOrdering, io_region, subset_selector,
                                  use_nodeset_for_part_node_fields, order_blocks_by_creation_order);
           }
           else if (part->primary_entity_rank() == stk::topology::FACE_RANK)
             define_side_set(*part, bulk_data, io_region, subset_selector,
                             use_nodeset_for_part_node_fields);
           else if (part->primary_entity_rank() == stk::topology::EDGE_RANK)
             define_side_set(*part, bulk_data, io_region, subset_selector,
                             use_nodeset_for_part_node_fields);
         }
       }

       define_communication_maps(bulk_data, io_region, subset_selector);

       if (input_region != NULL)
         io_region.synchronize_id_and_name(input_region, true);

       // for streaming refinement, each "pseudo-processor" doesn't know about others, so we pick a sort order
       //   and use it for all pseudo-procs - the original_block_order property is used to set the order
       //   on all procs.
       if (sort_stk_parts_by_name) {
         int offset=0;
         for (mesh::PartVector::const_iterator i = parts->begin(); i != parts->end(); ++i) {

           mesh::Part * const part = *i ;

           if (is_part_io_part(*part)) {
             if (part->primary_entity_rank() == mesh::InvalidEntityRank)
               continue;
             else if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK) {
               Ioss::GroupingEntity *element_block = io_region.get_entity(getPartName(*part));
               if (element_block) {
                 if (element_block->property_exists("original_block_order")) {
                   element_block->property_erase("original_block_order");
                 }
                 element_block->property_add(Ioss::Property("original_block_order", offset));
                 ++offset;
               }
             }
           }
         }
       }
    }

    void define_output_db(Ioss::Region & io_region ,
                          const mesh::BulkData &bulk_data,
                          const std::vector<std::vector<int>> &attributeOrdering,
                          const Ioss::Region *input_region,
                          const stk::mesh::Selector *subset_selector,
                          const bool sort_stk_parts_by_name,
                          const bool use_nodeset_for_part_node_fields)
    {
      io_region.begin_mode( Ioss::STATE_DEFINE_MODEL );
      define_output_db_within_state_define(io_region, bulk_data, attributeOrdering, input_region, subset_selector, sort_stk_parts_by_name, use_nodeset_for_part_node_fields);
      io_region.end_mode( Ioss::STATE_DEFINE_MODEL );
    }

    //----------------------------------------------------------------------

    namespace {



      template <typename INT>
      void write_side_data_to_ioss( Ioss::GroupingEntity & io ,
                                    mesh::Part * const part ,
                                    const mesh::BulkData & bulk_data,
                                    const Ioss::ElementTopology *element_topology,
                                    const stk::mesh::Selector *subset_selector )
      {
        std::vector<INT> elem_side_ids;
        stk::mesh::EntityVector sides;
        stk::topology stk_elem_topology = map_ioss_topology_to_stk(element_topology, bulk_data.mesh_meta_data().spatial_dimension());

        fill_element_and_side_ids(io, part, bulk_data, stk_elem_topology, subset_selector, sides, elem_side_ids);
        size_t num_sides = sides.size();

        const size_t num_side_written = io.put_field_data("element_side",elem_side_ids);

        if ( num_sides != num_side_written ) {
          std::ostringstream msg ;

          msg << "stk::io::write_side_data_to_ioss FAILED for " ;
          msg << io.name();
          msg << " in Ioss::GroupingEntity::put_field_data:" ;
          msg << " num_sides = " << num_sides ;
          msg << " , num_side_written = " << num_side_written ;
          throw std::runtime_error( msg.str() );
        }

        const mesh::FieldBase *df = get_distribution_factor_field(*part);
        if (df != NULL) {
          field_data_to_ioss(bulk_data, df, sides, &io, "distribution_factors", Ioss::Field::MESH);
        }

        const mesh::MetaData & meta_data = mesh::MetaData::get(*part);

        const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
            stk::io::field_data_to_ioss(bulk_data, f, sides, &io, f->name(), Ioss::Field::ATTRIBUTE);
          }
        }
      }

      //----------------------------------------------------------------------
      template <typename INT>
      void output_node_block(Ioss::NodeBlock &nb,
                             stk::mesh::Part &part,
                             const stk::mesh::BulkData &bulk,
                             const stk::mesh::Selector *subset_selector)
      {
        //----------------------------------
        // Exactly one node block to obtain the nodal coordinates and ids:
        // Note that the "ids" field of the nodes needs to be written
        // before any other bulk data that uses node ids since it sets up
        // the global->local mapping of nodes for the output database.
        // Similarly for the element "ids" field related to bulk data
        // using element ids.
        std::vector<mesh::Entity> nodes ;
        size_t num_nodes = get_entities(part, stk::topology::NODE_RANK,
                                        bulk, nodes, true, subset_selector);

        std::vector<INT> node_ids; node_ids.reserve(num_nodes);
        for(size_t i=0; i<num_nodes; ++i) {
          const mesh::Entity node = nodes[i] ;
          node_ids.push_back(bulk.identifier(node));
        }

        size_t num_ids_written = nb.put_field_data("ids", node_ids);
        if ( num_nodes != num_ids_written) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::NodeBlock::put_field_data:" ;
          msg << " num_nodes = " << num_nodes ;
          msg << " , num_ids_written = " << num_ids_written ;
          throw std::runtime_error( msg.str() );
        }

        if (nb.get_database()->needs_shared_node_information()) {
          std::vector<int> owning_processor; owning_processor.reserve(num_nodes);
          for(size_t i=0; i<num_nodes; ++i) {
            owning_processor.push_back(bulk.parallel_owner_rank(nodes[i]));
          }
          nb.put_field_data("owning_processor", owning_processor);
        }

        const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
        const mesh::FieldBase *coord_field = meta_data.coordinate_field();
        assert(coord_field != NULL);
        field_data_to_ioss(bulk, coord_field, nodes, &nb, "mesh_model_coordinates", Ioss::Field::MESH);

        const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          if (stk::io::is_valid_part_field(f, part_primary_entity_rank(part), part,
                                           Ioss::Field::ATTRIBUTE)) {
            stk::io::field_data_to_ioss(bulk, f, nodes, &nb, f->name(), Ioss::Field::ATTRIBUTE);
          }
        }
      }

      template <typename INT>
      void output_element_block(Ioss::ElementBlock *block,
                                const stk::mesh::BulkData &bulk,
                                const stk::mesh::Selector *subset_selector)
      {
        const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
        const std::string& name = block->name();
        mesh::Part* part = getPart( meta_data, name);

        assert(part != NULL);
        std::vector<mesh::Entity> elements;
        stk::mesh::EntityRank type = part_primary_entity_rank(*part);
        size_t num_elems = get_entities(*part, type, bulk, elements, false, subset_selector);

        stk::topology topo = part->topology();
        if (topo == stk::topology::INVALID_TOPOLOGY) {
          std::ostringstream msg ;
          msg << " INTERNAL_ERROR: Part " << part->name() << " returned INVALID from get_topology()";
          throw std::runtime_error( msg.str() );
        }
        size_t nodes_per_elem = topo.num_nodes();

        std::vector<INT> elem_ids; elem_ids.reserve(num_elems);
        std::vector<INT> connectivity; connectivity.reserve(num_elems*nodes_per_elem);

        for (size_t i = 0; i < num_elems; ++i) {

          elem_ids.push_back(bulk.identifier(elements[i]));

          stk::mesh::Entity const * elem_nodes = bulk.begin_nodes(elements[i]);

          for (size_t j = 0; j < nodes_per_elem; ++j) {
            connectivity.push_back(bulk.identifier(elem_nodes[j]));
          }
        }

        const size_t num_ids_written = block->put_field_data("ids", elem_ids);
        const size_t num_con_written = block->put_field_data("connectivity", connectivity);

        if ( num_elems != num_ids_written || num_elems != num_con_written ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::ElementBlock::put_field_data:" << std::endl ;
          msg << "  num_elems = " << num_elems << std::endl ;
          msg << "  num_ids_written = " << num_ids_written << std::endl ;
          msg << "  num_connectivity_written = " << num_con_written << std::endl ;
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
        const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, elem_rank, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, elements, block, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }
      }

      template <typename INT>
      void output_node_set(Ioss::NodeSet *ns, const stk::mesh::BulkData &bulk,
                           const stk::mesh::Selector *subset_selector)
      {
        const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
        const std::string& name = ns->name();
        mesh::Part* part = getPart( meta_data, name);

        // If part is null, then it is possible that this nodeset is a "viz nodeset" which
        // means that it is a nodeset containing the nodes of an element block.
        // See if there is a property base_stk_part_name and if so, get the part with
        // that name.
        if (part == NULL) {
          if (ns->property_exists(base_stk_part_name)) {
            std::string base_name = ns->get_property(base_stk_part_name).get_string();
            part = getPart( meta_data, base_name);
          }
          if (part == NULL) {
            std::ostringstream msg ;
            msg << " FAILED in Ioss::NodeSet::output_node_set:"
                << " Could not find stk part corresponding to nodeset named '"
                << name << "'";
            throw std::runtime_error( msg.str() );
          }
        }

        std::vector<stk::mesh::Entity> nodes ;
        size_t num_nodes = get_entities(*part, stk::topology::NODE_RANK, bulk, nodes, true, subset_selector);

        std::vector<INT> node_ids; node_ids.reserve(num_nodes);
        for(size_t i=0; i<num_nodes; ++i) {
          const stk::mesh::Entity node = nodes[i] ;
          node_ids.push_back(bulk.identifier(node));
        }

        size_t num_ids_written = ns->put_field_data("ids", node_ids);
        if ( num_nodes != num_ids_written ) {
          std::ostringstream msg ;
          msg << " FAILED in Ioss::NodeSet::output_node_set:"
              << " num_nodes = " << num_nodes
              << ", num_ids_written = " << num_ids_written;
          throw std::runtime_error( msg.str() );
        }

        stk::mesh::Field<double> *df_field = meta_data.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors_" + name);
        if (df_field != NULL) {
          const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*df_field, stk::topology::NODE_RANK, *part);
          if (res.num_scalars_per_entity() > 0) {
            stk::io::field_data_to_ioss(bulk, df_field, nodes, ns, "distribution_factors", Ioss::Field::MESH);
          }
        }

        const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
        std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
        while (I != fields.end()) {
          const mesh::FieldBase *f = *I ; ++I ;
          const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
          if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
            const mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*f, stk::topology::NODE_RANK, *part);
            if (res.num_scalars_per_entity() > 0) {
              stk::io::field_data_to_ioss(bulk, f, nodes, ns, f->name(), Ioss::Field::ATTRIBUTE);
            }
          }
        }
      }

      template <typename INT>
      void output_communication_maps(Ioss::Region &io_region,
                                     const stk::mesh::BulkData &bulk,
                                     const stk::mesh::Selector *subset_selector)
      {
        if (bulk.parallel_size() > 1) {
          const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);
          mesh::Selector selector = meta.globally_shared_part();
          if (subset_selector) selector &= *subset_selector;

          std::vector<mesh::Entity> entities;
          get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), entities);

          const std::string cs_name("node_symm_comm_spec");
          Ioss::CommSet * io_cs = io_region.get_commset(cs_name);
          STKIORequire(io_cs != NULL);

          // Allocate data space to store <id, processor> pair
          assert(io_cs->field_exists("entity_processor"));
          size_t size = io_cs->get_field("entity_processor").raw_count();

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
          io_cs->put_field_data("entity_processor", ep);
        }
      }

      template <typename INT>
      void output_side_set(Ioss::SideSet *ss,
                           const stk::mesh::BulkData &bulk,
                           const stk::mesh::Selector *subset_selector)
      {
        const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
        size_t block_count = ss->block_count();
        for (size_t i=0; i < block_count; i++) {
          Ioss::SideBlock *block = ss->get_block(i);
          if (stk::io::include_entity(block)) {
            stk::mesh::Part * part = getPart(meta_data, block->name());
            const Ioss::ElementTopology *parent_topology = block->parent_element_topology();
            stk::io::write_side_data_to_ioss<INT>(*block, part, bulk, parent_topology, subset_selector);
          }
        }
      }

    } // namespace <blank>

    void write_output_db_node_block(Ioss::Region& io_region,
                                    const stk::mesh::BulkData& bulk,
                                    const stk::mesh::Selector *subset_selector = nullptr)
    {
        const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);

        bool ints64bit = db_api_int_size(&io_region) == 8;

        Ioss::NodeBlock & nb = *io_region.get_node_blocks()[0];

        if (ints64bit)
          output_node_block<int64_t>(nb, meta.universal_part(), bulk, subset_selector);
        else
          output_node_block<int>(nb, meta.universal_part(), bulk, subset_selector);
    }

    void write_output_db_rest_of_mesh(Ioss::Region& io_region,
                                      const stk::mesh::BulkData& bulk,
                                      const stk::mesh::Selector *subset_selector = nullptr)
    {
      bool ints64bit = db_api_int_size(&io_region) == 8;

      //----------------------------------
      const Ioss::ElementBlockContainer& elem_blocks = io_region.get_element_blocks();
      for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
          it != elem_blocks.end(); ++it) {
        if (ints64bit)
          output_element_block<int64_t>(*it, bulk, subset_selector);
        else
          output_element_block<int>(*it, bulk, subset_selector);
      }

      //----------------------------------
      const Ioss::NodeSetContainer& node_sets = io_region.get_nodesets();
      for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
          it != node_sets.end(); ++it) {
        if (ints64bit)
          output_node_set<int64_t>(*it, bulk, subset_selector);
        else
          output_node_set<int>(*it, bulk, subset_selector);
      }

      //----------------------------------
      const Ioss::SideSetContainer& side_sets = io_region.get_sidesets();
      for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
          it != side_sets.end(); ++it) {
        if (ints64bit)
          output_side_set<int64_t>(*it, bulk, subset_selector);
        else
          output_side_set<int>(*it, bulk, subset_selector);
      }

      if (ints64bit)
        output_communication_maps<int64_t>(io_region, bulk, subset_selector);
      else
        output_communication_maps<int>(io_region, bulk, subset_selector);
    }

    void write_output_db(Ioss::Region& io_region,
                         const stk::mesh::BulkData& bulk,
                         const stk::mesh::Selector *subset_selector)
    {
      io_region.begin_mode( Ioss::STATE_MODEL );
      write_output_db_node_block(io_region, bulk, subset_selector);
      write_output_db_rest_of_mesh(io_region, bulk, subset_selector);
      io_region.end_mode( Ioss::STATE_MODEL );
    }

    //----------------------------------------------------------------------
    bool is_part_io_part(const stk::mesh::Part &part)
    {
      return NULL != part.attribute<Ioss::GroupingEntity>();
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
                                       const stk::mesh::FieldBase &df_field)
    {
      stk::mesh::MetaData &m = mesh::MetaData::get(p);
      m.declare_attribute_no_delete(p,&df_field);
    }

    const Ioss::Field::RoleType* get_field_role(const stk::mesh::FieldBase &f)
    {
      return f.attribute<Ioss::Field::RoleType>();
    }

    void set_field_role(stk::mesh::FieldBase &f, const Ioss::Field::RoleType &role)
    {
      Ioss::Field::RoleType *my_role = new Ioss::Field::RoleType(role);
      stk::mesh::MetaData &m = mesh::MetaData::get(f);
      const Ioss::Field::RoleType *check = m.declare_attribute_with_delete(f, my_role);
      if ( check != my_role ) {
        if (*check != *my_role) {
          std::ostringstream msg ;
          msg << " FAILED in IossBridge -- set_field_role:"
              << " The role type for field name= " << f.name() 
              << " was already set to " << *check
              << ", so it is not possible to change it to " << *my_role;
          delete my_role;
          throw std::runtime_error( msg.str() );
        }
        delete my_role;
      }
    }

    namespace {
      void define_input_nodeblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
        assert(node_blocks.size() == 1);

        Ioss::NodeBlock *nb = node_blocks[0];
        stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT,
                                  meta.universal_part(), stk::topology::NODE_RANK);
      }

      void define_input_elementblock_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
        for(size_t i=0; i < elem_blocks.size(); i++) {
          if (stk::io::include_entity(elem_blocks[i])) {
            stk::mesh::Part* const part = meta.get_part(elem_blocks[i]->name());
            assert(part != NULL);
            stk::io::define_io_fields(elem_blocks[i], Ioss::Field::TRANSIENT,
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
            assert(part != NULL);
            stk::io::define_io_fields(nodesets[i], Ioss::Field::TRANSIENT,
                                      *part, part_primary_entity_rank(*part));
          }
        }
      }

      void define_input_sideset_fields(Ioss::Region &region, stk::mesh::MetaData &meta)
      {
        if (meta.spatial_dimension() <= meta.side_rank())
          return;

        const Ioss::SideSetContainer& side_sets = region.get_sidesets();
        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {
          Ioss::SideSet *entity = *it;
          if (stk::io::include_entity(entity)) {
            const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
            for(size_t i=0; i < blocks.size(); i++) {
              if (stk::io::include_entity(blocks[i])) {
                stk::mesh::Part* const part = meta.get_part(blocks[i]->name());
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
    }

    void insert_var_names_for_part(const Ioss::GroupingEntity* entity, stk::mesh::Part* part, FieldNameToPartVector& names)
    {
        Ioss::NameList geNames;
        entity->field_describe(Ioss::Field::TRANSIENT, &geNames);
        for(const std::string& geName : geNames)
            names.push_back(std::make_pair(geName, part));
    }

    void insert_var_names(const Ioss::GroupingEntity *entity, FieldNameToPartVector &names, stk::mesh::MetaData& meta)
    {
        if (stk::io::include_entity(entity)) {
            stk::mesh::Part* part = meta.get_part(entity->name());
            insert_var_names_for_part(entity, part, names);
        }
    }

    template <typename GroupingEntityVector>
    FieldNameToPartVector get_grouping_entity_var_names(const GroupingEntityVector &groupingEntities,stk::mesh::MetaData& meta)
    {
        FieldNameToPartVector names;
        for(size_t i=0; i < groupingEntities.size(); i++)
            insert_var_names(groupingEntities[i], names, meta);
        stk::util::sort_and_unique(names, FieldNameToPartLess());
        return names;
    }

    FieldNameToPartVector get_var_names(Ioss::Region &region, Ioss::EntityType type, stk::mesh::MetaData& meta)
    {
        switch(type)
        {
            case Ioss::NODEBLOCK:
            {
                FieldNameToPartVector names;
                const Ioss::NodeBlockContainer nodeBlocks = region.get_node_blocks();
                ThrowRequire(nodeBlocks.size() == 1);
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

    void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                        stk::mesh::EntityRank part_type,
                        Ioss::GroupingEntity *io_entity,
                        Ioss::Field::RoleType filter_role)
    {
      std::vector<stk::mesh::Entity> entities;
      stk::io::get_entity_list(io_entity, part_type, bulk, entities);

      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

      std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
      while (I != fields.end()) {
        const stk::mesh::FieldBase *f = *I; ++I;
        if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
          stk::io::field_data_to_ioss(bulk, f, entities, io_entity, f->name(), filter_role);
        }
      }
    }

    struct DefineOutputFunctor
    {
      void operator()(stk::mesh::BulkData &bulk, stk::mesh::Part &part, stk::mesh::EntityRank rank, Ioss::GroupingEntity *ge, Ioss::Field::RoleType role)
      {  stk::io::ioss_add_fields(part, rank, ge, role); }
    };

    struct ProcessOutputFunctor
    {
      void operator()(stk::mesh::BulkData &bulk, stk::mesh::Part &part, stk::mesh::EntityRank rank, Ioss::GroupingEntity *ge, Ioss::Field::RoleType role)
      {  put_field_data(bulk, part, rank, ge, role); }
    };

    template <typename T>
    void process_field_loop(Ioss::Region &region,
                            stk::mesh::BulkData &bulk, T& callable)
    {
        stk::mesh::MetaData & meta = bulk.mesh_meta_data();

        Ioss::NodeBlock *nb = region.get_node_blocks()[0];
        callable(bulk, meta.universal_part(), stk::topology::NODE_RANK,
                 dynamic_cast<Ioss::GroupingEntity *>(nb), Ioss::Field::TRANSIENT);

        const stk::mesh::PartVector & all_parts = meta.get_parts();
        for ( stk::mesh::PartVector::const_iterator
                ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

          stk::mesh::Part * const part = *ip;

          if (stk::io::is_part_io_part(*part)) {
            Ioss::GroupingEntity *entity = region.get_entity(part->name());
            ThrowRequireMsg(entity != nullptr, "Could not find Ioss grouping entity with name: " + part->name() + ". Contact sierra-help.");

            if (entity->type() == Ioss::SIDESET) {
              Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
              int block_count = sset->block_count();

              for (int i=0; i < block_count; i++) {
                Ioss::SideBlock *fb = sset->get_block(i);
                callable(bulk, *part,
                         stk::mesh::EntityRank( part->primary_entity_rank() ),
                         dynamic_cast<Ioss::GroupingEntity *>(fb), Ioss::Field::TRANSIENT);
              }
            } else {
              callable(bulk, *part,
                       stk::mesh::EntityRank( part->primary_entity_rank() ),
                       entity, Ioss::Field::TRANSIENT);
            }
          }
        }
    }

    void process_output_request(Ioss::Region &region,
                                stk::mesh::BulkData &bulk,
                                int step)
    {
      region.begin_state(step);
      ProcessOutputFunctor functor;
      process_field_loop(region, bulk, functor);
      region.end_state(step);
    }

    void write_node_sharing_info(Ioss::DatabaseIO *dbo, const EntitySharingInfo &nodeSharingInfo)
    {
        Ioss::CommSet* io_cs = dbo->get_region()->get_commset("commset_node");
        if(io_cs)
        {
            std::vector<int> entity_proc(2*nodeSharingInfo.size());
            int counter = 0;
            for(auto &nodeProc : nodeSharingInfo)
            {
                entity_proc[counter]     = nodeProc.first;
                entity_proc[counter + 1] = nodeProc.second;
                counter += 2;
            }
            size_t size_field = entity_proc.size()*(sizeof(int));
            io_cs->put_field_data("entity_processor", entity_proc.data(), size_field);
        }
    }

    void write_file_for_subdomain(const std::string &baseFilename,
                                  int index_subdomain,
                                  int num_subdomains,
                                  int global_num_nodes,
                                  int global_num_elems,
                                  stk::mesh::BulkData& bulkData,
                                  const EntitySharingInfo &nodeSharingInfo,
                                  int numSteps,
                                  double timeStep)
    {
        int width = std::log10(static_cast<double>(num_subdomains - 1))+1;
        std::ostringstream os;
        os << baseFilename << "." << num_subdomains << "." << std::setfill('0') << std::setw(width) << index_subdomain;

        std::string dbtype("exodusII");
        Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, os.str(), Ioss::WRITE_RESULTS, MPI_COMM_SELF);

        Ioss::Region out_region(dbo, "name");
        dbo->get_region()->property_add(Ioss::Property("processor_count", num_subdomains));
        dbo->get_region()->property_add(Ioss::Property("my_processor", index_subdomain));
        dbo->get_region()->property_add(Ioss::Property("global_node_count", global_num_nodes));
        dbo->get_region()->property_add(Ioss::Property("global_element_count", global_num_elems));

        //////////////////////////

        out_region.begin_mode(Ioss::STATE_DEFINE_MODEL);
        stk::io::define_output_db_within_state_define(out_region, bulkData, {});
        Ioss::CommSet *commset = new Ioss::CommSet(dbo, "commset_node", "node", nodeSharingInfo.size());
        commset->property_add(Ioss::Property("id", 1));
        dbo->get_region()->add(commset);
        out_region.end_mode(Ioss::STATE_DEFINE_MODEL);

        //////////////////////////

        out_region.begin_mode(Ioss::STATE_MODEL);
        stk::io::write_output_db_node_block(out_region, bulkData);
        write_node_sharing_info(dbo, nodeSharingInfo);
        stk::io::write_output_db_rest_of_mesh(out_region, bulkData);
        out_region.end_mode(Ioss::STATE_MODEL);

        //////////////////////////

        if(numSteps > 0)
        {
          out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
          DefineOutputFunctor functor;
          process_field_loop(out_region, bulkData, functor);
          out_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

          //////////////////////////

          out_region.begin_mode(Ioss::STATE_TRANSIENT);
          int out_step = out_region.add_state(timeStep);
          process_output_request(out_region, bulkData, out_step);
          out_region.end_mode(Ioss::STATE_TRANSIENT);
        }

        stk::io::delete_selector_property(out_region);
    }

  }//namespace io
}//namespace stk

