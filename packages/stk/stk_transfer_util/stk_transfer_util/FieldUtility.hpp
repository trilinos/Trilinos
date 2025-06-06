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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_FIELDUTILITY_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_FIELDUTILITY_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/FieldState.hpp"  // for FieldState, StateNP1, StateNM1
#include "stk_mesh/base/Types.hpp"       // for EntityRank
#include "stk_transfer/TransferTypes.hpp"
#include "stk_util/util/string_case_compare.hpp"

#include <functional>                    // for function
#include <limits>                        // for numeric_limits
#include <string>                        // for string, operator==, basic_st...
#include <utility>                       // for pair
#include <vector>                        // for vector
namespace stk::mesh { class BulkData; }
namespace stk::mesh { class FieldBase; }
namespace stk::mesh { class MetaData; }
namespace stk::mesh { class Selector; }
namespace stk::mesh { struct Entity; }
namespace stk::mesh { class Part; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

struct FieldSpec {
  const stk::mesh::FieldBase* field{nullptr};
  std::string name;
  stk::mesh::FieldState state{stk::mesh::FieldState::StateNone};
  unsigned int index{0};
  double upperBound{ std::numeric_limits<double>::max()};
  double lowerBound{-std::numeric_limits<double>::max()};
  FieldTransform  preTransform{impl::default_transform};
  FieldTransform postTransform{impl::default_transform};
  double defaultValue{0.0};

  FieldSpec() {}
  FieldSpec(const std::string& fieldName)
    : name(fieldName)
  {
  }
  FieldSpec(const std::string& fieldName, stk::mesh::FieldState fieldState)
    : name(fieldName)
    , state(fieldState)
  {
  }
  FieldSpec(const std::string& fieldName, stk::mesh::FieldState fieldState, unsigned int fieldIndex)
    : name(fieldName)
    , state(fieldState)
    , index(fieldIndex)
  {
  }
  FieldSpec(const stk::mesh::FieldBase* fieldPtr);
  FieldSpec(const stk::mesh::FieldBase* fieldPtr, unsigned int fieldIndex);

  bool operator < (const FieldSpec& rhs) const { return stk::less_case(name, rhs.name); }
  bool operator != (const FieldSpec& rhs) const { return stk::equal_case(name, rhs.name); }
  bool operator == (const FieldSpec& rhs) const { return stk::not_equal_case(name, rhs.name); }

  operator const stk::mesh::FieldBase*() { return field; }
};

using FieldSpecVector = std::vector<FieldSpec>;

struct IndexedField {
  const stk::mesh::FieldBase* field{nullptr};
  unsigned int index{0};
  unsigned fieldSize{0};
};

stk::mesh::FieldState state_name(const char* name);

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData& meta, const std::vector<stk::transfer::FieldSpec>& fieldnames);

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData& meta, const std::vector<stk::transfer::FieldSpec>& fieldSpecs, const stk::mesh::EntityRank rank);

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData* meta, const std::vector<stk::transfer::FieldSpec>& fieldnames);

std::vector<IndexedField>
get_fields(const stk::mesh::MetaData* meta, const std::vector<stk::transfer::FieldSpec>& fieldSpecs, const stk::mesh::EntityRank rank);

std::vector<const stk::mesh::FieldBase*> extract_field_pointers(const std::vector<IndexedField>& indexedFields);

std::vector<double> get_upper_bounds(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);
std::vector<double> get_lower_bounds(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);

std::vector<double> get_default_field_values(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);

void apply_bounds(const unsigned length, double* fieldData, const double lowerBound, const double upperBound);
void apply_bounds(std::vector<double>& fieldData, const double lowerBound, const double upperBound);

std::vector<FieldTransform> get_pre_transforms(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);
std::vector<FieldTransform> get_post_transforms(const std::vector<stk::transfer::FieldSpec>& fieldSpecs);

unsigned get_size_of_field(const stk::mesh::FieldBase* field);

} // namespace transfer
} // namespace stk



#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_FIELDUTILITY_HPP_ */
