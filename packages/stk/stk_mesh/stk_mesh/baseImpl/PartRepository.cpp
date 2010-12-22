/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>

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

void append_part_method(
  std::string & msg , const Part & part , const char * method )
{
  msg.append( "stk::mesh::Part[" );
  msg.append( part.name() );
  msg.append( "]." );
  msg.append( method );
}


void assert_contain( Part & superset , Part & subset , const char * method )
{
  if ( ! contain( subset.supersets() , superset ) ) {
    std::string msg ;
    append_part_method( msg , superset , method );
    msg.append( "(...) FAILED Requirement that " );
    msg.append( "Part[" );
    msg.append( subset.name() );
    msg.append( "] is a subset" );
    throw std::runtime_error(msg);
  }
}

void assert_same_universe( const Part & part_superset ,
                           const char * method ,
                           const Part & part_subset )
{
  const PartVector & a = part_superset.supersets();
  const PartVector & b = part_subset.supersets();

  if ( a.empty() || b.empty() || a[0] != b[0] ) {
    std::string msg ;
    append_part_method( msg , part_superset , method );
    msg.append( "(...) FAILED Requirement that Part[" );
    msg.append( part_subset.name() );
    msg.append( "] are in the same " );
    msg.append( universal_part_name() );
    throw std::runtime_error(msg);
  }
}

void assert_not_same( const Part & part ,
                      const char * method ,
                      const Part & arg_part )
{
  if ( & part == & arg_part ) {
    std::string msg ;
    append_part_method( msg , part , method );
    msg.append( "(...) FAILED Requirement that Part[" );
    msg.append( arg_part.name() );
    msg.append( "] is not the same" );
    throw std::runtime_error(msg);
  }
}

void assert_same( const Part & part ,
                      const char * method ,
                      const Part & arg_part )
{
  if ( & part != & arg_part ) {
    std::string msg ;
    append_part_method( msg , part , method );
    msg.append( "(...) FAILED Requirement that Part[" );
    msg.append( arg_part.name() );
    msg.append( "] is the same" );
    throw std::runtime_error(msg);
  }
}

void assert_not_superset( const Part & part ,
                          const char * method ,
                          const Part & arg_part )
{
  if ( contain( part.supersets() , arg_part ) ) {
    std::string msg ;
    append_part_method( msg , part , method );
    msg.append( "(...) FAILED Requirement that Part[" );
    msg.append( arg_part.name() );
    msg.append( "] is not a superset" );
    throw std::runtime_error(msg);
  }
}

void assert_rank_ordering( const Part & superset ,
                           const char * method ,
                           const Part & subset )
{
  if ( superset.primary_entity_rank() < subset.primary_entity_rank() ) {
    std::ostringstream msg ;
    msg << "stk::mesh::Part[ " << superset.name();
    msg << " , rank(" << superset.primary_entity_rank();
    msg << ") ]." << method ;
    msg << "( Part[ " << subset.name();
    msg << " , rank(" << subset.primary_entity_rank();
    msg << ") ] ) FAILED Rank ordering requirement" ;
    throw std::runtime_error( msg.str() );
  }
}

} // namespace


Part * PartRepository::universal_part() const
{
  return m_universal_part;
}

const PartVector & PartRepository::all_parts() const
{
  return m_universal_part->subsets();
}

Part * PartRepository::declare_part( const std::string & arg_name , EntityRank arg_rank )
{
  const PartVector & all_parts = m_universal_part->subsets();
  Part * p = find( all_parts, arg_name );

  if ( p == NULL ) {
    p = declare_part_impl( arg_name, arg_rank );
  }

  if ( p->primary_entity_rank() != arg_rank ) {
    std::ostringstream msg;
    msg << "stk::mesh::Part[ " << arg_name ;
    msg << ",rank(" << p->primary_entity_rank() << ")" ;
    msg << "] : FAILED to declare part; " ;
    msg << "Part of name '" << arg_name ;
    msg << "' of rank " << p->primary_entity_rank() ;
    msg << " already exists";
    msg << " User cannot redeclare " << arg_name ;
    msg << " with different rank, " << arg_rank ;
    throw std::runtime_error ( msg.str() );
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
    assert_contain( *m_universal_part, *p , method );

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

    const PartVector & all_parts = m_universal_part->subsets();
    p = find( all_parts, p_name );
    if ( p == NULL ) {
      // Create the part:

      p = new Part( m_meta_data , p_name , p_rank , all_parts.size() );
      declare_subset_impl( *m_universal_part, *p );

      // Define the part to be an intersection of the given parts:

      p->m_partImpl.set_intersection_of( pset_clean );

      for ( PartVector::iterator
            i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
        declare_subset( **i, *p );
      }
    }
    else if ( pset_clean != p->intersection_of()) {
      // This error is "inconceivable" and is
      // only possible by heroic malicious abuse.
      std::string msg ;
      msg.append(method).append(p_name).append(" FAILED FROM MALICIOUS ABUSE");
      throw std::invalid_argument(msg);
    }
  }

  return p ;
}


Part * PartRepository::declare_part_impl( const std::string & name, EntityRank rank)
{
  size_t ordinal = m_universal_part->subsets().size();
  Part * part = new Part(m_meta_data,name,rank,ordinal);
  declare_subset_impl(*m_universal_part, *part);
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

    assert_not_same( superset , method , subset );
    assert_not_superset(  superset , method , subset );
    assert_same_universe( superset , method , subset );
    assert_rank_ordering( superset , method , subset );

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

  assert_not_same( root_part , method , target_part );
  assert_same_universe( root_part, method, target_part );
  assert_same( root_part, method, *relation.m_root );
  assert_same( target_part, method, *relation.m_target );

  root_part.m_partImpl.add_relation( relation );
  target_part.m_partImpl.add_relation( relation );
}


PartRepository::PartRepository(MetaData * meta)
  : m_meta_data(meta)
{
  m_universal_part = new Part( m_meta_data, universal_part_name(), ~0u, 0 );
  m_universal_part->m_partImpl.add_part_to_subset(*m_universal_part);
}

PartRepository::~PartRepository()
{
  // The universal part is the 0^th entry in the subset vector.
  // Delete all but the universal part in the loop, deleting
  // the universal part will invalidate the universal part subset vector.
  // Thus delete the universal part outside of the loop.

  try {
    for ( PartVector::const_iterator
          i = m_universal_part->subsets().end() ;
          --i != m_universal_part->subsets().begin() ; ) {
      Part * part = *i ;
      try { delete part ; } catch(...) {}
    }
    try { delete m_universal_part ; } catch(...) {}
    m_universal_part = NULL ;
  } catch(...){}
}

} // namespace impl 
} // namespace mesh 
} // namespace stk 


