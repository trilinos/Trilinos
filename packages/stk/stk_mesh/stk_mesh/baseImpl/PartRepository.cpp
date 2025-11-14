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

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stddef.h>                     // for NULL, size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Part.hpp>       // for Part, contain, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank
#include <stk_util/util/ReportHandler.hpp>  // for ThrowErrorMsgIf
#include <vector>                       // for vector, etc
#include "stk_topology/topology.hpp"    // for topology, etc

namespace stk {
namespace mesh {
namespace impl {

namespace {

inline
std::string universal_part_name()
{
  static const std::string name = convert_to_internal_name("UNIVERSAL");
  return name;
}

void assert_same_universe( const Part & superset ,
                           const Part & subset,
                           const char * method )
{
  const PartVector & a = superset.supersets();
  const PartVector & b = subset.supersets();

  STK_ThrowErrorMsgIf( a.empty() || b.empty() || a[0] != b[0],
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << superset.name() << "] and " <<
                   "Part[" << subset.name() << "] are in the same " <<
                   universal_part_name() );
}

void assert_not_same( const Part & part1 ,
                      const Part & part2 ,
                      const char * method )
{
  STK_ThrowErrorMsgIf( & part1 == & part2,
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << part1.name() << "] and " <<
                   "Part[" << part2.name() << "] are not the same" );
}

void assert_not_superset( const Part & superset ,
                          const Part & subset ,
                          const char * method )
{
  STK_ThrowErrorMsgIf( contain( subset.supersets() , superset ),
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << superset.name() << "] " <<
                   "is not a superset of " <<
                   "Part[" << subset.name() << "]" );
}

void assert_rank_ordering( const Part & superset ,
                           const Part & subset ,
                           const char * method )
{
  STK_ThrowErrorMsgIf( superset.primary_entity_rank() < subset.primary_entity_rank(),
                   method << " Part '" << subset.name() << "' (rank="
                   << subset.primary_entity_rank() << ", topology=" 
                   << subset.topology() << ") can't be a subset of part '" << superset.name()
                   << "' (rank="<<superset.primary_entity_rank()<<", topology="
                   << superset.topology() << "), needs to have rank <= "
                   << superset.primary_entity_rank());
}

} // namespace

Part * PartRepository::universal_part() const
{
  return m_universal_part;
}

Part * PartRepository::get_part_by_name(const std::string &name) const
{
    auto iter = m_name_to_parts_map.find(name);
    if(iter != m_name_to_parts_map.end())
        return iter->second;
    return nullptr;
}

void PartRepository::rename(Part* part, const std::string& newName)
{
  const std::string& oldName = part->name();
  auto iter = m_name_to_parts_map.find(oldName);
  STK_ThrowRequireMsg(iter != m_name_to_parts_map.end(), "PartRepository::rename: Failed to find part with oldName="<<oldName);
  m_name_to_parts_map.erase(oldName);
  m_name_to_parts_map[newName] = part;
  part->set_name(newName);
}

void PartRepository::add_part(Part* part)
{
    m_all_parts.push_back(part);
    m_name_to_parts_map[part->name()] = part;
}

const PartVector & PartRepository::get_all_parts() const
{
  return m_all_parts;
}

const PartVector PartRepository::get_mesh_parts() const
{
  PartVector non_internal_parts;
  for (size_t i=0 ; i<m_all_parts.size() ; ++i) {
    if (!is_internal_part(*m_all_parts[i])) {
      non_internal_parts.push_back(m_all_parts[i]);
    }
  }
  return non_internal_parts;
}

Part * PartRepository::declare_part( const std::string & arg_name , EntityRank arg_rank, bool force_no_induce )
{
  Part * p = get_part_by_name( arg_name );

  if ( p == NULL ) {
    p = declare_part_impl( arg_name, arg_rank, force_no_induce );
  }
  else {
    p->set_primary_entity_rank(arg_rank);
    p->set_force_no_induce(force_no_induce);
  }

  return p;
}

Part * PartRepository::declare_part_impl( const std::string & name, EntityRank rank, bool force_no_induce )
{
  size_t ordinal = get_all_parts().size();
  Part * part = new Part(m_meta_data, name, rank, ordinal, force_no_induce);
  declare_subset_impl(*m_universal_part, *part);
  add_part(part);
  return part;
}

void PartRepository::declare_subset_impl( Part & superset_part, Part & subset_part )
{
  superset_part.add_part_to_subset( subset_part );
  subset_part.add_part_to_superset( superset_part );
}

void PartRepository::declare_subset( Part & superset, Part & subset )
{
  static const char method[] = "stk::mesh::impl::PartRepository::declare_subset" ;

  if ( ! contain( subset.supersets() , superset ) ) {

    assert_not_same(      superset , subset , method );
    assert_not_superset(  superset , subset , method );
    assert_same_universe( superset , subset , method );
    assert_rank_ordering( superset , subset , method );

    // Insert this symmetric relationship first
    // so that it does not get revisited.

    declare_subset_impl( superset, subset );

    // Transitive:

    const PartVector & subset_subsets = subset.subsets();
    for ( PartVector::const_iterator
          i =  subset_subsets.begin() ; i != subset_subsets.end() ; ++i ) {
      declare_subset( superset, **i );
    }

    const PartVector & superset_supersets = superset.supersets();
    for ( PartVector::const_iterator
          i =  superset_supersets.begin() ; i != superset_supersets.end() ; ++i ) {
      declare_subset( **i, subset );
    }

  }
}

PartRepository::PartRepository(MetaData * meta)
  : m_meta_data(meta),
    m_universal_part(NULL),
    m_all_parts(),
    m_name_to_parts_map()
{
  m_universal_part = new Part( m_meta_data, universal_part_name(), stk::topology::INVALID_RANK, 0 /*ordinal*/);
  add_part(m_universal_part);
}

PartRepository::~PartRepository()
{
  try {
    for ( PartVector::const_iterator i = m_all_parts.begin() ;
          i != m_all_parts.end() ; ++i) {
      Part * part = *i ;
      try { delete part ; } catch(...) {}
    }
  } catch(...){}
}

static const char INTERNAL_PART_PREFIX  = '{';
static const char INTERNAL_PART_POSTFIX = '}';

bool is_internal_part_name(const std::string& partName)
{
    return partName.size() > 2 && partName.front() == INTERNAL_PART_PREFIX && partName.back() == INTERNAL_PART_POSTFIX;
}

bool is_internal_part(const Part& part)
{
    return is_internal_part_name(part.name());
}
std::string convert_to_internal_name(const std::string& part_name)
{
  return std::string(INTERNAL_PART_PREFIX+part_name+INTERNAL_PART_POSTFIX);
}

} // namespace impl
} // namespace mesh
} // namespace stk

