
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_mesh/base/Part.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

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
  const PartVector & intersection = p.intersection_of();

  std::vector<Part*>::const_iterator i ;

  if ( lead != NULL ) { os << lead ; }
  os << "Part[ " ;
  os << p.name() ;
  os << " , " ;
  os << p.mesh_meta_data_ordinal() ;
  os << " ] {" ;
  os << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Supersets {" ;
  for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Intersection_Of {" ;
  for ( i = intersection.begin() ; i != intersection.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " } }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Subsets {" ;
  for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

void order( PartVector & v )
{
  PartVector::iterator ev = v.end();
  PartVector::iterator iv = v.begin();
  std::sort( iv , ev , PartLess() );
  iv = std::unique( iv , ev );
  v.erase( iv , ev );
}

bool insert( PartVector & v , Part & part )
{
  const PartVector::iterator e = v.end();
        PartVector::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool new_member = i == e || *i != & part ;

  if ( new_member ) { Part * const tmp = & part ; v.insert( i , tmp ); }
  return new_member ;
}

void remove( PartVector & v , Part & part )
{
  const PartVector::iterator e = v.end();
        PartVector::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  if ( i != e && *i == & part ) { v.erase( i ); }
}

bool contain( const PartVector & v , const Part & part )
{
  const PartVector::const_iterator e = v.end();
        PartVector::const_iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  return i != e && *i == & part ;
}

bool contain( const PartVector & super , const PartVector & sub )
{
  bool result = ( ! sub.empty() ) && ( sub.size() <= super.size() );

  if ( result ) {
    PartLess comp ;

    const PartVector::const_iterator ev = super.end();
          PartVector::const_iterator iv = super.begin();

    const PartVector::const_iterator ep = sub.end();
          PartVector::const_iterator ip = sub.begin();

    while ( result && ip != ep ) {
      Part * const q = *ip ; ++ip ;
      iv = std::lower_bound( iv , ev , q , comp );
      result = iv != ev && *iv == q ; 
    }
  }

  return result ;
}

size_t intersect( const PartVector & v , const PartVector & p )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartVector::const_iterator ev = v.end();
        PartVector::const_iterator iv = v.begin();

  const PartVector::const_iterator ep = p.end();
        PartVector::const_iterator ip = p.begin();

  size_t count = 0 ;

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { ++count ; }
  }

  return count ;
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

//----------------------------------------------------------------------
//----------------------------------------------------------------------

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

void assert_universal_part( Part & part , const char * method )
{
  if ( ! part.supersets().empty() ) {
    std::string msg ;
    append_part_method( msg , part , method );
    msg.append( "(...) FAILED: Only valid for Part[" );
    msg.append( universal_part_name() );
    msg.append( "]" );
    throw std::runtime_error( msg );
  }
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

void assert_same_universe( const Part & part ,
                           const char * method ,
                           const Part & arg_part )
{
  const PartVector & a = part.supersets();
  const PartVector & b = arg_part.supersets();

  if ( a.empty() || b.empty() || a[0] != b[0] ) {
    std::string msg ;
    append_part_method( msg , part , method );
    msg.append( "(...) FAILED Requirement that Part[" );
    msg.append( arg_part.name() );
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
  if ( superset.primary_entity_type() < subset.primary_entity_type() ) {
    std::ostringstream msg ;
    msg << "stk::mesh::Part[ " << superset.name();
    msg << " , rank(" << superset.primary_entity_type();
    msg << ") ]." << method ;
    msg << "( Part[ " << subset.name();
    msg << " , rank(" << subset.primary_entity_type();
    msg << ") ] ) FAILED Rank ordering requirement" ;
    throw std::runtime_error( msg.str() ); 
  }
}

}

Part::~Part()
{
  // If I am the universal part then I need to delete my subsets:
  try {
    if ( m_supersets.empty() && this == m_subsets[0] ) {
      PartVector::iterator i = m_subsets.begin();
      for ( ++i ; i != m_subsets.end() ; ++i ) {
        try { delete *i ; } catch(...) {}
      }
    }
  } catch(...){}
}

// Universal part constructor:
Part::Part( MetaData * arg_meta_data )
  : m_name( universal_part_name() ),
    m_attribute(),
    m_subsets() , m_supersets() , m_intersect() , m_relations() ,
    m_mesh_meta_data( arg_meta_data ) ,
    m_universe_ordinal( 0 ),
    m_entity_rank( ~0u )
{
  m_subsets.push_back( this );
}

// Subset part constructor:
Part::Part( MetaData          * arg_meta_data ,
            const std::string & arg_name ,
            EntityType            arg_rank ,
            size_t              arg_ordinal )
  : m_name( arg_name ),
    m_attribute(),
    m_subsets() , m_supersets() , m_intersect() , m_relations() ,
    m_mesh_meta_data( arg_meta_data ),
    m_universe_ordinal( arg_ordinal ),
    m_entity_rank( arg_rank )
{}

//----------------------------------------------------------------------
// A universal part may declare named subsets

Part & Part::declare_part( const std::string & arg_name , EntityType arg_rank )
{
  assert_universal_part( *this , "declare_part" );

  Part * p = find( m_subsets , arg_name );

  if ( p == NULL ) {
    p = new Part( m_mesh_meta_data , arg_name , arg_rank , m_subsets.size() );
    m_subsets.push_back( p );
    p->m_supersets.push_back( this );
  }

  return *p ;
}

//----------------------------------------------------------------------
// A universal part may declare intersection subsets

Part & Part::declare_part( const PartVector & part_intersect )
{
  static const char method[] = "declare_part" ;

  assert_universal_part( *this , method );

  PartVector pset_clean ;

  for ( PartVector::const_iterator
        i = part_intersect.begin() ; i != part_intersect.end() ; ++i ) {
    Part * const p = *i ;
    assert_contain( *this , *p , method );

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
    EntityType p_rank = std::numeric_limits<EntityType>::max();

    p_name.assign("{");
    for ( PartVector::iterator
          i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
      if ( i != pset_clean.begin() ) { p_name.append( separator ); }
      p_name.append( (*i)->name() );
      if ( (*i)->primary_entity_type() < p_rank ) {
        p_rank = (*i)->primary_entity_type();
      }
    }
    p_name.append("}");

    p = find( m_subsets , p_name );

    if ( p == NULL ) {

      // Create the part:

      p = new Part( m_mesh_meta_data , p_name , p_rank , m_subsets.size() );
      m_subsets.push_back( p );
      p->m_supersets.push_back( this );

      // Define the part to be an intersection of the given parts:

      p->m_intersect = pset_clean ; // Copy

      for ( PartVector::iterator
            i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
        (*i)->declare_subset( *p );
      }
    }
    else if ( pset_clean != p->m_intersect ) {
      // This error is "inconceivable" and is
      // only possible by heroic malicious abuse.
      std::string msg ;
      msg.append(method).append(p_name).append(" FAILED FROM MALICIOUS ABUSE");
      throw std::invalid_argument(msg);
    }
  }

  return *p ;
}

//----------------------------------------------------------------------

void Part::declare_subset( Part & subset )
{
  static const char method[] = "declare_subset" ;

  if ( ! contain( subset.supersets() , *this ) ) {

    assert_not_same( *this , method , subset );
    assert_not_superset(  *this , method , subset );
    assert_same_universe( *this , method , subset );
    assert_rank_ordering( *this , method , subset );

    // Insert this symmetric relationship first
    // so that it does not get revisited.

    insert( this->m_subsets , subset );
    insert( subset.m_supersets , *this );

    // Transitive:

    for ( PartVector::iterator
          i =  subset.m_subsets.begin() ; i != subset.m_subsets.end() ; ++i ) {
      declare_subset( **i );
    }

    for ( PartVector::iterator
          i =  m_supersets.begin() ; i != m_supersets.end() ; ++i ) {
      (*i)->declare_subset( subset );
    }

    // Induced intersection-part membership:

    for ( PartVector::iterator
          i =  m_subsets.begin() ;
          i != m_subsets.end() ; ++i ) {

      Part & pint = **i ;

      if ( ! pint.m_intersect.empty() && ( & pint != & subset ) ) {

        // If 'subset' is a subset of every member of 'pint.m_intersect'
        // then it is by definition a subset of 'pint.
        if ( contain( subset.m_supersets , pint.m_intersect ) ) {
          pint.declare_subset( subset );
        }
      }
    }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

