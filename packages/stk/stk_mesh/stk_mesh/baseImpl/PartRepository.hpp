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

#ifndef stk_mesh_PartRepository_hpp
#define stk_mesh_PartRepository_hpp

#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank
#include <string>                       // for string
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace mesh {


namespace impl {

class PartRepository {
public:
  explicit PartRepository(MetaData * meta);
  ~PartRepository();

  Part * universal_part() const;

  const PartVector & get_all_parts()  const;  // returns all parts
  const PartVector   get_mesh_parts() const; // returns the non-internal parts

  Part * declare_part( const std::string & arg_name , EntityRank arg_rank, bool force_no_induce=false );
  void declare_subset( Part & superset, Part & subset );

  template<class T>
  const T * declare_attribute_with_delete( Part & , const T *);
  template<class T>
  const T * declare_attribute_no_delete( Part & , const T *);
  template<class T>
  bool remove_attribute( Part & , const T *);

private:
  PartRepository();
  PartRepository(const PartRepository & );
  PartRepository & operator = ( const PartRepository & );

  Part * declare_part_impl( const std::string & name, EntityRank rank, bool force_no_induce );
  void declare_subset_impl( Part & superset, Part & subset );

  MetaData * m_meta_data;
  Part * m_universal_part;
  PartVector m_all_parts;
};

template<class T>
inline
const T *
PartRepository::declare_attribute_with_delete( Part & p, const T * a )
{
  return p.m_partImpl.declare_attribute_with_delete<T>( a );
}

template<class T>
inline
const T *
PartRepository::declare_attribute_no_delete( Part & p, const T * a )
{
  return p.m_partImpl.declare_attribute_no_delete<T>( a );
}

template<class T>
inline
bool
PartRepository::remove_attribute( Part & p, const T * a )
{
  return p.m_partImpl.remove_attribute<T>( a );
}

bool is_internal_part(const Part& part);
std::string convert_to_internal_name(const std::string& part_name);

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_PartRepository_hpp
