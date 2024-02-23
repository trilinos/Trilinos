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

#ifndef stk_mesh_baseImpl_FieldRestriction_hpp
#define stk_mesh_baseImpl_FieldRestriction_hpp

#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/Types.hpp>      // for FieldArrayRank
#include <string>                       // for string
#include <vector>                       // for vector


namespace stk {
namespace mesh {

/**
 * A field restrictions is one of the three fundamental components of a field specification
 * (see Field.hpp for a full discusssion); it defines a set of entities that have a field.
 *
 * This class encapsulates a minimal set of data for a field restriction. The API for
 * declaring field restrictions is in MetaData.hpp.
 */
class FieldRestriction {
  public:

  typedef int size_type ;

  FieldRestriction()
    : m_part(nullptr),
      m_selector(),
      m_num_scalars_per_entity(0),
      m_dimension(0)
  {
  }

  FieldRestriction( const FieldRestriction & rhs )
    : m_part(rhs.m_part),
      m_selector( rhs.m_selector ),
      m_num_scalars_per_entity( rhs.m_num_scalars_per_entity ),
      m_dimension(rhs.m_dimension)
  {
  }

  FieldRestriction & operator = ( const FieldRestriction & rhs )
  {
    m_part = rhs.m_part;
    m_selector = rhs.m_selector;
    m_num_scalars_per_entity = rhs.m_num_scalars_per_entity;
    m_dimension = rhs.m_dimension;
    return *this ;
  }

  void add_union(const Selector& otherSelector)
  {
    if (m_part != nullptr) {
      m_selector = *m_part;
      m_part = nullptr;
    }
    m_selector |= otherSelector;
  }

  explicit FieldRestriction( const Part& input_part)
   : m_part(&input_part),
     m_selector(input_part),
     m_num_scalars_per_entity(0),
     m_dimension(0)
  {
  }

  explicit FieldRestriction( const Selector& input_selector)
   : m_part(nullptr),
     m_selector(input_selector),
     m_num_scalars_per_entity(0),
     m_dimension(0)
  {
  }

  const Selector& selector() const { return m_selector; }

  bool selects(const Part& part) const;

  void set_num_scalars_per_entity(size_type value) { m_num_scalars_per_entity = value; }
  size_type num_scalars_per_entity() const { return m_num_scalars_per_entity; }

  void set_dimension(size_type dim) { m_dimension = dim; }
  size_type dimension() const { return m_dimension; }

  bool operator < ( const FieldRestriction & rhs ) const
  {
    return m_part!=nullptr&&rhs.m_part!=nullptr ? m_part->mesh_meta_data_ordinal() < rhs.m_part->mesh_meta_data_ordinal() : m_selector < rhs.m_selector;
  }
  bool operator == ( const FieldRestriction & rhs ) const
  {
    return m_part!=nullptr&&rhs.m_part!=nullptr ? m_part->mesh_meta_data_ordinal() == rhs.m_part->mesh_meta_data_ordinal() : this->m_selector == rhs.m_selector;
  }
  bool operator != ( const FieldRestriction & rhs ) const
  {
    return this->m_selector != rhs.m_selector;
  }

  void print(std::ostream & os, const Selector & selector) const;

  private:
  const Part* m_part;
  Selector m_selector;
  size_type m_num_scalars_per_entity;
  size_type m_dimension;
};

typedef std::vector<FieldRestriction> FieldRestrictionVector;

std::string print_restriction(const FieldRestriction & restr, const Selector& selector);

} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_FieldRestriction_hpp
