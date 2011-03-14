/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>

#include <stdlib.h>

#include <iostream>

namespace stk {
namespace mesh {
namespace impl {

namespace {

const char * universal_part_name()
{
  static const char name[] = "{UNIVERSAL}" ;
  return name ;
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

void assert_same( const Part & part1 ,
                  const Part & part2,
                  const char * method )
{
  ThrowErrorMsgIf( & part1 != & part2,
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << part1.name() << "] and " <<
                   "Part[" << part2.name() << "] are the same" );
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

void assert_superset( Part & superset ,
                      Part & subset ,
                      const char * method )
{
  ThrowErrorMsgIf( ! contain( subset.supersets() , superset ),
                   method << "(...) FAILED Requirement that " <<
                   "Part[" << superset.name() << "] " <<
                   "is a superset of " <<
                   "Part[" << subset.name() << "]" );
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

Part * PartRepository::declare_part( const std::string & arg_name , EntityRank arg_rank )
{
  const PartVector & all_parts = get_all_parts();
  Part * p = find( all_parts, arg_name );

  if ( p == NULL ) {
    p = declare_part_impl( arg_name, arg_rank );
  }
  else {
    p->m_partImpl.set_primary_entity_rank(arg_rank);
  }

  return p;
}

Part * PartRepository::declare_part( const PartVector & part_intersect )
{
  static const char method[] = "stk::mesh::PartRepository::declare_part" ;

  PartVector pset_clean ;

  for ( PartVector::const_iterator
        i = part_intersect.begin() ; i != part_intersect.end() ; ++i ) {
    Part * const p = *i ;
    assert_superset( *m_universal_part, *p , method );

    // If 'p' is a superset of another member
    // then it is redundant in this intersection.
    // Only keep non-redundant intersections.

    PartVector::const_iterator j = part_intersect.begin();
    for ( ; j != part_intersect.end() &&
            ! contain( (*j)->supersets() , *p ) ; ++j );
    if ( j == part_intersect.end() ) {
      pset_clean.push_back( p );
    }
  }

  // Sort and unique the intersection
  order( pset_clean );

  Part * p = NULL ;
  if ( 1 == pset_clean.size() ) {
    // Only one remaining part, it is the subset.
    p = pset_clean[0] ;
  }
  else {
    const char separator[] = "^" ;
    // Generate a name and rank reflecting the intersection.
    // Rank is the minimum rank of the intersection members.

    std::string p_name ;
    EntityRank p_rank = InvalidEntityRank;

    p_name.assign("{");
    for ( PartVector::iterator
          i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
      if ( i != pset_clean.begin() ) { p_name.append( separator ); }
      p_name.append( (*i)->name() );
      if ( (*i)->primary_entity_rank() < p_rank ) {
        p_rank = (*i)->primary_entity_rank();
      }
    }
    p_name.append("}");

    const PartVector & all_parts = get_all_parts();
    p = find( all_parts, p_name );
    if ( p == NULL ) {
      // Create the part:

      p = declare_part_impl( p_name , p_rank );

      // Define the part to be an intersection of the given parts:

      p->m_partImpl.set_intersection_of( pset_clean );

      for ( PartVector::iterator
            i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
        declare_subset( **i, *p );
      }
    }
    else {
      // This error is "inconceivable" and is
      // only possible by heroic malicious abuse.
      ThrowInvalidArgMsgIf( pset_clean != p->intersection_of(),
                            p_name << " FAILED FROM MALICIOUS ABUSE" );
    }
  }

  return p ;
}


Part * PartRepository::declare_part_impl( const std::string & name, EntityRank rank)
{
  size_t ordinal = get_all_parts().size();
  Part * part = new Part(m_meta_data,name,rank,ordinal);
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
  static const char method[] = "stk::mesh::PartRepository::declare_subset" ;

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

    // Induced intersection-part membership:

    const PartVector & superset_subsets = superset.subsets();
    for ( PartVector::const_iterator
          i =  superset_subsets.begin() ;
          i != superset_subsets.end() ; ++i ) {

      Part & pint = **i ;

      if ( ! pint.intersection_of().empty() && ( & pint != & subset ) ) {

        // If 'subset' is a subset of every member of 'pint.intersection_of()'
        // then it is by definition a subset of 'pint.
        if ( contain( subset.supersets() , pint.intersection_of() ) ) {
          declare_subset( pint, subset );
        }
      }
    }
  }
}


void PartRepository::declare_part_relation( Part & root_part, PartRelation relation, Part & target_part )
{
  static const char method[] = "stk::mesh::PartRepository::declare_part_relation" ;

  assert_not_same(      root_part   , target_part        , method );
  assert_same_universe( root_part   , target_part        , method );
  assert_same(          root_part   , *relation.m_root   , method );
  assert_same(          target_part , *relation.m_target , method );

  root_part.m_partImpl.add_relation( relation );
  target_part.m_partImpl.add_relation( relation );
}


PartRepository::PartRepository(MetaData * meta)
  : m_meta_data(meta),
    m_universal_part(NULL),
    m_all_parts()
{
  m_universal_part = new Part( m_meta_data, universal_part_name(), ~0u, 0 );
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

} // namespace impl
} // namespace mesh
} // namespace stk

