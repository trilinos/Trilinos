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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SEARCHFIELD_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SEARCHFIELD_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_mesh/base/FieldBase.hpp>     // for FieldBase
#include "stk_topology/topology.hpp"       // for topology
#include <map>                             // for map, map<>::value_compare
#include <vector>                          // for vector
#include <string>

namespace stk { namespace mesh { class BulkData; } }

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

class SearchField
{
public:
  using Transform = double (*)(double);
  
  SearchField(const stk::mesh::FieldBase* field = nullptr,
              double (*transformFunc)(double) = default_transform,
              double defaultValue = 0.0)
  : m_field(field)
  , m_transform(transformFunc)
  , m_defaultValue(defaultValue)
  {
    if(field && !field->restrictions().empty()) {
      m_fieldSize = field->restrictions().begin()->num_scalars_per_entity();
    }
  }

  const stk::mesh::FieldBase* get_field() const { return m_field; }

  double default_value() const { return m_defaultValue; }

  void set_default_value(double defaultValue)
  {
    m_defaultValue = defaultValue;
  }

  unsigned get_field_size() const { return m_fieldSize; }

  const stk::mesh::BulkData& get_mesh() const
  {
    return m_field->get_mesh();
  }

  const std::string& name() const
  {
    return m_field->name();
  }

  Transform get_transform() const { return m_transform; }

  void set_transform(Transform transformFunc)
  {
    m_transform = transformFunc;
  }
  
  double transform(double x) const
  {
    return m_transform(x);
  }
  
protected:
  const stk::mesh::FieldBase* m_field{nullptr};
  Transform m_transform = default_transform;
  double m_defaultValue{0.0};
  unsigned m_fieldSize{0};

  static double default_transform(double x)
  {
    return x;
  }
  
  SearchField() = delete;
};

} // namespace stk
} // namespace search

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SEARCHFIELD_HPP_ */
