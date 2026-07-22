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

#ifndef stk_mesh_PartRepository_hpp
#define stk_mesh_PartRepository_hpp

#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank
#include <stk_util/util/string_case_compare.hpp>
#include <string>                       // for string
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace mesh {


namespace impl {

struct StringLessCase
{
    bool operator()(const std::string &a, const std::string &b) const
    {
        return stk::less_case(a, b);
    }
};

class PartRepository {
public:
  explicit PartRepository(MetaData * meta);
  ~PartRepository();

  Part * universal_part() const;

  Part * get_part_by_name(const std::string &name) const;
  void rename(Part* part, const std::string& newName);
  const PartVector & get_all_parts()  const;  // returns all parts
  const PartVector   get_mesh_parts() const; // returns the non-internal parts

  size_t size() const;

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
  void add_part(Part* part);

  MetaData * m_meta_data;
  Part * m_universal_part;
  PartVector m_all_parts;
  std::map<std::string, stk::mesh::Part*, StringLessCase> m_name_to_parts_map;
};

template<class T>
inline
const T *
PartRepository::declare_attribute_with_delete( Part & p, const T * a )
{
  return p.declare_attribute_with_delete<T>( a );
}

template<class T>
inline
const T *
PartRepository::declare_attribute_no_delete( Part & p, const T * a )
{
  return p.declare_attribute_no_delete<T>( a );
}

template<class T>
inline
bool
PartRepository::remove_attribute( Part & p, const T * a )
{
  return p.remove_attribute<T>( a );
}

inline
size_t PartRepository::size() const 
{
  return m_all_parts.size();
}

bool is_internal_part_name(const std::string& partName);
bool is_internal_part(const Part& part);
std::string convert_to_internal_name(const std::string& part_name);

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_PartRepository_hpp
