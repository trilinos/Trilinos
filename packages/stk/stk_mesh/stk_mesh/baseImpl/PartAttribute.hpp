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

#ifndef STK_STK_MESH_STK_MESH_BASEIMPL_PARTATTRIBUTE_HPP_
#define STK_STK_MESH_STK_MESH_BASEIMPL_PARTATTRIBUTE_HPP_

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <string>
#include <vector>
#include <utility>

namespace stk
{
namespace mesh
{
namespace impl
{
template <typename T>
bool has_part_attribute(const stk::mesh::Part& part)
{
  return part.attribute<T>() != nullptr;
}

template <typename T>
const decltype(T::value) get_part_attribute(const stk::mesh::Part& part)
{
  const T* altAttr = part.attribute<T>();

  if (altAttr == nullptr)
    throw std::runtime_error(std::string("stk::mesh::impl::get_part_attribute called with no ") + typeid(T).name() +
                             " attribute on part= " + part.name());

  return altAttr->value;
}

template <typename T>
void set_part_attribute(stk::mesh::Part& part, const decltype(T::value)& value)
{
  mesh::MetaData& meta = mesh::MetaData::get(part);

  const T* altAttr = part.attribute<T>();
  if (!altAttr) {
    T* altAttr1 = new T();
    altAttr1->value = value;
    meta.declare_attribute_with_delete(part, altAttr1);
  } else {
    if (value != altAttr->value) {
      bool success = meta.remove_attribute(part, altAttr);
      if (!success)
        throw std::runtime_error(std::string("stk::mesh::impl::set_part_attribute failed to remove") +
                                 typeid(T).name() + " attribute, part= " + part.name());
      T* altAttr1 = new T();
      altAttr1->value = value;
      meta.declare_attribute_with_delete(part, altAttr1);
    }
  }
}

template <typename T>
std::pair<bool, stk::mesh::Part*> is_unique_part_attribute(stk::mesh::Part& part, const decltype(T::value)& value)
{
  mesh::MetaData& meta = mesh::MetaData::get(part);
  stk::mesh::PartVector pv = meta.get_parts();
  for (unsigned ii = 0; ii < pv.size(); ++ii) {
    if (has_part_attribute<T>(*pv[ii]) && get_part_attribute<T>(*pv[ii]) == value && &part != pv[ii]) {
      return std::make_pair(false, pv[ii]);
    }
  }

  return std::make_pair(true, nullptr);
}

template <typename T>
void set_unique_part_attribute(stk::mesh::Part& part, const decltype(T::value)& value)
{
  auto isUnique = is_unique_part_attribute<T>(part, value);
  if(!isUnique.first) {
    throw std::runtime_error(
        std::string("stk::mesh::impl::set_unique_part_attribute found another part with the same ") +
        typeid(T).name() + " attribute, part= " + part.name() + " other part= " + isUnique.second->name());
  }

  set_part_attribute<T>(part, value);
}

}  // namespace impl
}  // namespace mesh
}  // namespace stk

#endif /* STK_STK_MESH_STK_MESH_BASEIMPL_PARTATTRIBUTE_HPP_ */
