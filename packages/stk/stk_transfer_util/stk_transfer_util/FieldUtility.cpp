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
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_transfer_util/FieldUtility.hpp"
#include "stk_util/util/string_case_compare.hpp"
#include <vector>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

FieldSpec::FieldSpec(const stk::mesh::FieldBase* fieldPtr)
{
  if(nullptr != fieldPtr) {
    field = fieldPtr;
    name  = fieldPtr->name();
    state = fieldPtr->state();
  }
}

FieldSpec::FieldSpec(const stk::mesh::FieldBase* fieldPtr, unsigned int fieldIndex)
{
  if(nullptr != fieldPtr) {
    field = fieldPtr;
    name  = fieldPtr->name();
    state = fieldPtr->state();
    index = fieldIndex;
  }
}

stk::mesh::FieldState state_name(const char* name)
{
  using StateNameValue = std::pair<const std::string, const stk::mesh::FieldState>;

  static const std::vector<StateNameValue> table{
    StateNameValue("NONE", stk::mesh::StateNone),     StateNameValue("STATE_NONE", stk::mesh::StateNone),
    StateNameValue("NEW", stk::mesh::StateNew),       StateNameValue("OLD", stk::mesh::StateOld),
    StateNameValue("NP1", stk::mesh::StateNP1),       StateNameValue("N", stk::mesh::StateN),
    StateNameValue("NM1", stk::mesh::StateNM1),       StateNameValue("NM2", stk::mesh::StateNM2),
    StateNameValue("NM3", stk::mesh::StateNM3),       StateNameValue("NM4", stk::mesh::StateNM4),
    StateNameValue("STATE_NEW", stk::mesh::StateNew), StateNameValue("STATE_OLD", stk::mesh::StateOld),
    StateNameValue("STATE_NP1", stk::mesh::StateNP1), StateNameValue("STATE_N", stk::mesh::StateNP1),
    StateNameValue("STATE_NM1", stk::mesh::StateNM1), StateNameValue("STATE_NM2", stk::mesh::StateNM2),
    StateNameValue("STATE_NM3", stk::mesh::StateNM3), StateNameValue("STATE_NM4", stk::mesh::StateNM4),
  };

  for(const StateNameValue& entry : table) {
    if(stk::equal_case(entry.first, name)) return entry.second;
  }

  return stk::mesh::StateInvalid;
}

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData* metaData, const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  return (metaData != nullptr) ? get_fields(*metaData, fieldSpecs) : std::vector<IndexedField>{};
}

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData& metaData, const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<IndexedField> indexedFieldVec;
  indexedFieldVec.resize(fieldSpecs.size());

  // provide field names
  unsigned key = 0;
  for(const auto& spec : fieldSpecs) {
    const stk::mesh::FieldBase* field = (nullptr != spec.field) ? spec.field : stk::mesh::get_field_by_name(spec.name, metaData);
    STK_ThrowRequireMsg(field != nullptr, "field not found in get_fields, name= " + spec.name);
    indexedFieldVec[key] = { field->field_state(spec.state), spec.index, stk::transfer::get_size_of_field(field) };
    key++;
  }
  return indexedFieldVec;
}

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData* metaData, const std::vector<stk::transfer::FieldSpec>& fieldSpecs, const stk::mesh::EntityRank rank)
{
  return (metaData != nullptr) ? get_fields(*metaData, fieldSpecs, rank) : std::vector<IndexedField>{};
}

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData& metaData, const std::vector<stk::transfer::FieldSpec>& fieldSpecs, const stk::mesh::EntityRank rank)
{
  std::vector<IndexedField> indexedFieldVec;
  indexedFieldVec.resize(fieldSpecs.size());

  // provide field names
  unsigned key = 0;
  for(const auto& spec : fieldSpecs) {
    const stk::mesh::FieldBase* field = (nullptr != spec.field) ? spec.field : metaData.get_field(rank, spec.name);
    STK_ThrowRequireMsg(field != nullptr, "field not found in get_fields, name= " + spec.name);
    indexedFieldVec[key] = { field->field_state(spec.state), spec.index, stk::transfer::get_size_of_field(field) };
    key++;
  }
  return indexedFieldVec;
}

std::vector<const stk::mesh::FieldBase*> extract_field_pointers(const std::vector<IndexedField>& indexedFields)
{
  std::vector<const stk::mesh::FieldBase*> fields;

  for(auto indexedField : indexedFields) {
    fields.push_back(indexedField.field);
  }

  return fields;
}

std::vector<FieldTransform> get_pre_transforms(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<FieldTransform> preTransforms(fieldSpecs.size());

  for(unsigned i = 0; i < fieldSpecs.size(); ++i) {
    preTransforms[i] = fieldSpecs[i].preTransform;
  }
  return preTransforms;
}

std::vector<FieldTransform> get_post_transforms(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<FieldTransform> postTransforms(fieldSpecs.size());

  for(unsigned i = 0; i < fieldSpecs.size(); ++i) {
    postTransforms[i] = fieldSpecs[i].postTransform;
  }
  return postTransforms;
}

std::vector<double> get_upper_bounds(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<double> upperBounds(fieldSpecs.size());

  for(unsigned i = 0; i < fieldSpecs.size(); ++i) {
    upperBounds[i] = fieldSpecs[i].upperBound;
  }
  return upperBounds;
}

std::vector<double> get_lower_bounds(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<double> lowerBounds(fieldSpecs.size());

  for(unsigned i = 0; i < fieldSpecs.size(); ++i) {
    lowerBounds[i] = fieldSpecs[i].lowerBound;
  }
  return lowerBounds;
}

std::vector<double> get_default_field_values(const std::vector<stk::transfer::FieldSpec>& fieldSpecs)
{
  std::vector<double> defaultFieldValues(fieldSpecs.size());

  for(unsigned i = 0; i < fieldSpecs.size(); ++i) {
    defaultFieldValues[i] = fieldSpecs[i].defaultValue;
  }
  return defaultFieldValues;
}

void apply_bounds(const unsigned length, double* fieldData, const double lowerBound, const double upperBound)
{
  static constexpr double doubleMax = std::numeric_limits<double>::max();
  static constexpr double doubleMin = std::numeric_limits<double>::lowest();

  if(doubleMin != lowerBound) {
    for(unsigned i(0); i < length; ++i) {
      if(fieldData[i] < lowerBound) {
        fieldData[i] = lowerBound;
      }
    }
  }
  if(doubleMax != upperBound) {
    for(unsigned i(0); i < length; ++i) {
      if(upperBound < fieldData[i]) {
        fieldData[i] = upperBound;
      }
    }
  }
}

void apply_bounds(std::vector<double>& fieldData, const double lowerBound, const double upperBound)
{
  apply_bounds(fieldData.size(), fieldData.data(), lowerBound, upperBound);
}

unsigned get_size_of_field(const stk::mesh::FieldBase* field)
{
  if(!field->restrictions().empty()) {
    auto size = field->restrictions().begin()->num_scalars_per_entity();
    return size;
  }


  return 0;
}

} // namespace transfer
} // namespace stk

