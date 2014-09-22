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

#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stddef.h>                     // for NULL, size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Part.hpp>       // for Part, contain, etc
#include <stk_mesh/base/Trace.hpp>      // for TraceIf
#include <stk_mesh/base/Types.hpp>      // for PartVector, EntityRank
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowErrorMsgIf
#include <vector>                       // for vector, etc
#include "stk_mesh/baseImpl/PartImpl.hpp"  // for PartImpl
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

  ThrowErrorMsgIf( a.empty() || b.empty() || a[0] != b[0],
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << superset.name() << "] and " <<
                   "Part[" << subset.name() << "] are in the same " <<
                   universal_part_name() );
}

void assert_not_same( const Part & part1 ,
                      const Part & part2 ,
                      const char * method )
{
  ThrowErrorMsgIf( & part1 == & part2,
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << part1.name() << "] and " <<
                   "Part[" << part2.name() << "] are not the same" );
}

void assert_not_superset( const Part & superset ,
                          const Part & subset ,
                          const char * method )
{
  ThrowErrorMsgIf( contain( subset.supersets() , superset ),
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << superset.name() << "] " <<
                   "is not a superset of " <<
                   "Part[" << subset.name() << "]" );
}

void assert_rank_ordering( const Part & superset ,
                           const Part & subset ,
                           const char * method )
{
  ThrowErrorMsgIf( superset.primary_entity_rank() < subset.primary_entity_rank(),
                   method << "(...) FAILED Requirement that " <<
                   "Part[ " << superset.name() <<
                   " , rank(" << superset.primary_entity_rank() <<
                   ") ] has greater rank than " <<
                   "Part[ " << subset.name() <<
                   " , rank(" << subset.primary_entity_rank() << ") ]");
}

} // namespace

Part * PartRepository::universal_part() const
{
  return m_universal_part;
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
  TraceIf("stk::mesh::impl::PartRepository::declare_part", LOG_PART);

  const PartVector & all_parts = get_all_parts();
  Part * p = find( all_parts, arg_name );

  if ( p == NULL ) {
    p = declare_part_impl( arg_name, arg_rank, force_no_induce );
  }
  else {
    p->m_partImpl.set_primary_entity_rank(arg_rank);
    p->m_partImpl.set_force_no_induce(force_no_induce);
  }

  return p;
}

Part * PartRepository::declare_part_impl( const std::string & name, EntityRank rank, bool force_no_induce )
{
  size_t ordinal = get_all_parts().size();
  Part * part = new Part(m_meta_data, name, rank, ordinal, force_no_induce);
  declare_subset_impl(*m_universal_part, *part);
  m_all_parts.push_back(part);
  return part;
}

void PartRepository::declare_subset_impl( Part & superset_part, Part & subset_part )
{
  superset_part.m_partImpl.add_part_to_subset( subset_part );
  subset_part.m_partImpl.add_part_to_superset( superset_part );
}

void PartRepository::declare_subset( Part & superset, Part & subset )
{
  static const char method[] = "stk::mesh::impl::PartRepository::declare_subset" ;
  TraceIf(method, LOG_PART);

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
    m_all_parts()
{
  m_universal_part = new Part( m_meta_data, universal_part_name(), stk::topology::INVALID_RANK, 0 /*ordinal*/);
  m_all_parts.push_back(m_universal_part);
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
bool is_internal_part(const Part& part)
{
    return part.name().size() > 2 && *part.name().begin() == INTERNAL_PART_PREFIX && *part.name().rbegin() == INTERNAL_PART_POSTFIX;
}
std::string convert_to_internal_name(const std::string& part_name)
{
  std::ostringstream out;
  out << INTERNAL_PART_PREFIX << part_name << INTERNAL_PART_POSTFIX;
  std::string out_str = out.str();
  return out_str;
}

} // namespace impl
} // namespace mesh
} // namespace stk

