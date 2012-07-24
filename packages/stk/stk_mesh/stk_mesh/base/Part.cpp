/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <algorithm>
#include <ostream>
#include <sstream>

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

  if ( new_member ) { v.insert( i , &part ); }
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

std::string convert_to_internal_name(const std::string& part_name)
{
  std::ostringstream out;
  out << INTERNAL_PART_PREFIX << part_name << INTERNAL_PART_POSTFIX;
  return out.str();
}


//----------------------------------------------------------------------
//----------------------------------------------------------------------


} // namespace mesh
} // namespace stk

