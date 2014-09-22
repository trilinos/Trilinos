/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_baseImpl_FieldRestriction_hpp
#define stk_mesh_baseImpl_FieldRestriction_hpp

#include <Shards_Array.hpp>             // for int_t
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

  typedef shards::array_traits::int_t size_type ;

  FieldRestriction()
    : m_selector(),
      m_num_scalars_per_entity(0),
      m_dimension(0)
  {
  }

  FieldRestriction( const FieldRestriction & rhs )
    : m_selector( rhs.m_selector ),
      m_num_scalars_per_entity( rhs.m_num_scalars_per_entity ),
      m_dimension(rhs.m_dimension)
  {
  }

  FieldRestriction & operator = ( const FieldRestriction & rhs )
  {
    m_selector = rhs.m_selector;
    m_num_scalars_per_entity = rhs.m_num_scalars_per_entity;
    m_dimension = rhs.m_dimension;
    return *this ;
  }

  explicit FieldRestriction( const Selector& input_selector)
   : m_selector(input_selector),
     m_num_scalars_per_entity(0),
     m_dimension(0)
  {
  }

  const Selector& selector() const { return m_selector; }

  void set_num_scalars_per_entity(size_type value) { m_num_scalars_per_entity = value; }
  const size_type num_scalars_per_entity() const { return m_num_scalars_per_entity; }

  void set_dimension(size_type dim) { m_dimension = dim; }
  size_type dimension() const { return m_dimension; }

  bool operator < ( const FieldRestriction & rhs ) const
  {
    return m_selector < rhs.m_selector;
  }
  bool operator == ( const FieldRestriction & rhs ) const
  {
    return this->m_selector == rhs.m_selector;
  }
  bool operator != ( const FieldRestriction & rhs ) const
  {
    return this->m_selector != rhs.m_selector;
  }

  void print(
      std::ostream & os,
      const Selector & selector,
      FieldArrayRank field_rank
      ) const;

  private:
  Selector m_selector;
  size_type m_num_scalars_per_entity;
  size_type m_dimension;
};

typedef std::vector<FieldRestriction> FieldRestrictionVector;

std::string print_restriction(
    const FieldRestriction & restr,
    const Selector& selector,
    FieldArrayRank field_rank
    );

} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_FieldRestriction_hpp
