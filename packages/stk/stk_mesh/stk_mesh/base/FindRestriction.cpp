// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stk_mesh/base/FindRestriction.hpp>
#include <stddef.h>                     // for size_t
#include <ostream>                      // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Part.hpp>       // for Part
#include <string>                       // for operator<<, char_traits
#include <vector>                       // for vector, etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase::Restriction, etc
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsg

namespace stk {
namespace mesh {

namespace {

void throw_conflicting_restrictions(const FieldBase& field, const PartVector& parts, const FieldBase::Restriction& first, const FieldBase::Restriction& second)
{
  std::ostringstream oss;
  for ( size_t i = 0, ie = parts.size(); i < ie; ++i) {
    oss << parts[i]->name() << ", ";
  }
  ThrowErrorMsg("find_restriction calculation failed with different field-restriction selectors giving incompatible sizes.  Field name = "
                << field.name() << "on bucket with parts = " << oss.str() << "with one selector = " << first.selector() << ", num_scalars = "
                << first.num_scalars_per_entity()
                << " and another selector = " << second.selector() << ", dimension = " << second.num_scalars_per_entity() << "!");

}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
inline
const FieldBase::Restriction & empty_field_restriction()
{
  //NOT THREAD SAFE
  static const FieldBase::Restriction empty ;
  return empty ;
}
#ifdef __clang__
#pragma clang diagnostic pop
#endif

}

//Given a field and a range of parts, determine whether the field has a restriction
//for the part range. (Common usage is to provide the part-range from a bucket; i.e.,
//determine whether the field should be allocated in the bucket.)
const FieldBase::Restriction& find_restriction(const FieldBase& field,
                                               EntityRank erank,
                                               const PartVector& parts)
{
  if(static_cast<unsigned>(field.entity_rank()) == erank)
  {
      const std::vector<FieldBase::Restriction> & restrictions = field.restrictions();
      for(std::vector<FieldBase::Restriction>::const_iterator it=restrictions.begin(), it_end=restrictions.end(); it != it_end; ++it) {
        const Selector& selector = it->selector();
        if (selector(parts)) {
          return *it;
        }
      }
  }
  return empty_field_restriction();
}

const FieldBase::Restriction& find_restriction(const FieldBase& field,
                                               EntityRank erank,
                                               const Part & part)
{
  if(static_cast<unsigned>(field.entity_rank()) == erank)
  {
      const std::vector<FieldBase::Restriction> & restrictions = field.restrictions();
      for(std::vector<FieldBase::Restriction>::const_iterator it=restrictions.begin(), it_end=restrictions.end(); it != it_end; ++it) {
        const Selector& selector = it->selector();
        if (selector(part)) {
          return *it;
        }
      }
  }

  return empty_field_restriction();
}

const FieldBase::Restriction& find_and_check_restriction(const FieldBase& field,
                                                         EntityRank erank,
                                                         const PartVector& parts)
{
  const FieldBase::Restriction & empty = empty_field_restriction();
  const FieldBase::Restriction * restriction = & empty;

  if(static_cast<unsigned>(field.entity_rank()) == erank)
  {
      const std::vector<FieldBase::Restriction> & restrictions = field.restrictions();
      for(std::vector<FieldBase::Restriction>::const_iterator it=restrictions.begin(), it_end=restrictions.end(); it != it_end; ++it) {
        const Selector& selector = it->selector();
        if (selector(parts)) {
          if (restriction == &empty) {
            restriction = &*it;
          }
          else if (it->num_scalars_per_entity() != restriction->num_scalars_per_entity()) {
            throw_conflicting_restrictions(field, parts, *restriction, *it);
          }
        }
      }
  }
  return *restriction;
}

} // namespace mesh
} // namespace stk
