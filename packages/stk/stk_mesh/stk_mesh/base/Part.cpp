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

#include <stk_mesh/base/Part.hpp>
#include <algorithm>                    // for lower_bound, sort, unique
#include <ostream>                      // for operator<<, basic_ostream, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_util/util/string_case_compare.hpp>  // for not_equal_case


namespace stk {
namespace mesh {

Part * find( const PartVector & parts , const std::string & name )
{
  PartVector::const_iterator i = parts.begin();

  while ( i != parts.end() && not_equal_case((*i)->name(),name) ) { ++i ; }

  return i != parts.end() ? *i : NULL ;
}

//----------------------------------------------------------------------

std::ostream &
print( std::ostream & os , const char * const lead , const Part & p )
{
  const PartVector & supersets = p.supersets();
  const PartVector & subsets   = p.subsets();

  std::vector<Part*>::const_iterator i;

  if ( lead != NULL ) { os << lead ; }
  os << "Part[ ";
  os << "name: \"";
  os << p.name();
  os << "\" , ord: ";
  os << p.mesh_meta_data_ordinal();
  os << " , rank: ";
  if (p.primary_entity_rank() == stk::topology::INVALID_RANK) {
    os << "INVALID_RANK";
  }
  else {
    os << p.primary_entity_rank();
  }
  os << " ]";
  os << std::endl;

  if ( lead != NULL ) { os << lead ; }
  os << "  Supersets {";
  for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
    const std::string & n = (*i)->name();
    os << " \"" << n << "\"";
  }
  os << "  }" << std::endl;

  if ( lead != NULL ) { os << lead; }
  os << "  Subsets {";
  if (&p == &MetaData::get(p).universal_part() ) {
    os << " *all_parts*";
  }
  else {
    for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
      const std::string & n = (*i)->name();
      os << " \"" << n << "\"";
    }
  }
  os << "  }" << std::endl;

  return os;
}

//----------------------------------------------------------------------


bool insert( ConstPartVector & v , const Part & part )
{
  const ConstPartVector::iterator e = v.end();
        ConstPartVector::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool new_member = i == e || *i != & part ;

  if ( new_member ) { v.insert( i , &part ); }
  return new_member ;
}

bool insert( PartVector & v , Part & part )
{
  const PartVector::iterator e = v.end();
        PartVector::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool new_member = i == e || *i != & part ;

  if ( new_member ) { v.insert( i , &part ); }
  return new_member ;
}

void get_part_and_all_subsets(const Part& part, ConstPartVector& part_and_all_subsets)
{
  insert(part_and_all_subsets, part);
  
  const PartVector& subsets = part.subsets();
  for(size_t i=0; i<subsets.size(); ++i) {
    get_part_and_all_subsets(*subsets[i], part_and_all_subsets);
  }
}

void remove( PartVector & v , Part & part )
{
  const PartVector::iterator e = v.end();
        PartVector::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  if ( i != e && *i == & part ) { v.erase( i ); }
}


size_t intersect( const PartVector & v , const PartVector & p , PartVector & r )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartVector::const_iterator ev = v.end();
        PartVector::const_iterator iv = v.begin();

  const PartVector::const_iterator ep = p.end();
        PartVector::const_iterator ip = p.begin();

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { r.push_back( q ); }
  }

  return r.size() ;
}

bool intersect( const Part & a , const Part & b )
{
  const PartVector & a_sub = a.subsets();
  const PartVector & b_sub = b.subsets();
  return contain( a_sub , b ) ||
         contain( b_sub , a ) ||
         intersect( b_sub , a_sub );
}

bool contains( const PartVector & v , const Part & part )
{
  const PartVector::const_iterator e = v.end();
        PartVector::const_iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool match = (i != e) && (*i == &part);
  return match;
}

BulkData & Part::mesh_bulk_data() const
{
    return mesh_meta_data().mesh_bulk_data();
}

bool Part::contains(const Part& part) const
{
  if (this == &part) { // same part
    return true;
  }
  return m_subsetsEmpty ? false : stk::mesh::contains(subsets(), part);
}

void Part::set_name(const std::string& newName)
{
  m_name = newName;
}

void Part::set_primary_entity_rank( EntityRank entity_rank )
{
  if ( entity_rank == InvalidEntityRank || entity_rank == m_entity_rank ) return;

  const bool rank_already_set = m_entity_rank != InvalidEntityRank && entity_rank != m_entity_rank;

  if (rank_already_set) {
    std::stringstream str;
    str << "A part with name '" << m_name << "' "
        << "is defined in two different contexts " << m_mesh_meta_data->entity_rank_name(entity_rank) << " and "
        << m_mesh_meta_data->entity_rank_name(m_entity_rank) << ".";
    STK_ThrowErrorMsg(str.str());
  }

  m_entity_rank = entity_rank;
}

bool Part::add_part_to_subset(Part & part)
{
  m_subsetsEmpty = false;
  return insert(m_subsets, part);
}

bool Part::add_part_to_superset(Part & part)
{
  return insert(m_supersets, part);
}

namespace impl {

stk::CSet & get_attributes(stk::mesh::Part & part) {
  return part.get_attributes();
}

}

} // namespace mesh
} // namespace stk
