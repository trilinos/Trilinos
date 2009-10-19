//----------------------------------------------------------------------

#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

SelectorInterface::~SelectorInterface() {}

//----------------------------------------------------------------------

namespace {

/** \brief  Does [ibeg,iend) contain [jbeg,jend) */
template< typename iType , typename jType >
inline
bool contains( iType i , const iType iend ,
               jType j , const jType jend )
{
  while ( j < jend && i < iend ) {
    if      ( *j < *i ) { break ; } // [i,iend) cannot contain *j
    else if ( *i < *j ) { ++i ; }
    else { ++i , ++j ; }
  }
  return j == jend ; // [i,iend) contains every [j,jend)
}

/** \brief  Does [ibeg,iend) contain value */
template< typename iType , typename vType >
inline
bool contains( iType i , const iType iend , const vType & v )
{
  while ( i < iend && *i < v ) { ++i ; }
  return i < iend && *i == v ;
}

}

//----------------------------------------------------------------------

bool Selector::select_parts(
  const Bucket & candidate ,
  PartVector & selected_parts_from_candidate ) const
{
  selected_parts_from_candidate.clear();

  const std::pair<const unsigned *,const unsigned *>
    part_ordinals = candidate.superset_part_ordinals();

  bool result = contains( part_ordinals.first , part_ordinals.second ,
                          m_intersection.begin(), m_intersection.end() );

  if ( result && m_all_parts ) { // Union check requested
    std::vector<unsigned>::const_iterator i = m_union.begin();
    const unsigned * j = part_ordinals.first ;

    // Determine intersection of union with the candidate

    result = false ;

    while ( i != m_union.end() && j != part_ordinals.second ) {
      if      ( *i < *j ) { ++i ; }
      else if ( *j < *i ) { ++j ; }
      else {
        result = true ;
        // Find every match:
        Part * const part = (*m_all_parts)[ *i ];
        selected_parts_from_candidate.push_back( part );
        ++i ; ++j ;
      }
    }
  }

  return result ;
}

bool Selector::select( const Bucket & candidate ) const
{
  const std::pair<const unsigned *,const unsigned *>
    part_ordinals = candidate.superset_part_ordinals();

  bool result = contains( part_ordinals.first , part_ordinals.second ,
                          m_intersection.begin(), m_intersection.end() );

  if ( result && m_all_parts ) { // Union check requested

    std::vector<unsigned>::const_iterator i = m_union.begin();
    const unsigned * j = part_ordinals.first ;

    // Need at least one member of union

    result = false ;

    while ( i != m_union.end() && j != part_ordinals.second ) {
      if      ( *i < *j ) { ++i ; }
      else if ( *j < *i ) { ++j ; }
      else {
        result = true ;
        break ; // Only need one match
      }
    }
  }
  return result ;
}

//----------------------------------------------------------------------

namespace {

void verify_meta_data( const Part & p , const PartVector & v )
{
  const MetaData * const md  = & p.mesh_meta_data();

  PartVector::const_iterator i = v.begin() , j = v.end();
  for ( ; i != j ; ++i ) {
    if ( md != & (*i)->mesh_meta_data() ) {
      std::ostringstream msg ;
      msg << "Parts \"" << p.name() << "\" and \"" << (*i)->name()
          << "\" are not members of the same meta-data" ;
      throw std::runtime_error( msg.str() );
    }
  }
}

const PartVector * all_parts( const PartVector & v )
{
  const PartVector * all = ! v.empty() ?
                           & v[0]->mesh_meta_data().get_parts() : NULL ;
  return all ;
}

void copy_ids( std::vector<unsigned> & v , const Part & p )
{
  v.resize( 1 );
  v[0] = p.mesh_meta_data_ordinal();
}

void copy_ids( std::vector<unsigned> & v , const PartVector & p )
{
  {
    const size_t n = p.size();
    v.resize( n );
    for ( size_t k = 0 ; k < n ; ++k ) {
      v[k] = p[k]->mesh_meta_data_ordinal();
    }
  }

  {
    std::vector<unsigned>::iterator i = v.begin() , j = v.end();
    std::sort( i , j );
    i = std::unique( i , j );
    v.erase( i , j );
  }
}

}

Selector::Selector( const Part & required_part )
  : SelectorInterface(),
    m_intersection(),
    m_union(),
    m_all_parts( NULL ) // Only need if union exists
{
  copy_ids( m_intersection , required_part );
}

Selector::Selector( const PartVector & part_intersection )
  : SelectorInterface(),
    m_intersection(),
    m_union(),
    m_all_parts( NULL ) // Only need if union exists
{
  if ( ! part_intersection.empty() ) {
    verify_meta_data( * part_intersection[0] , part_intersection );
  }
  copy_ids( m_intersection , part_intersection );
}

Selector::Selector( const Part & required_part ,
                    const PartVector & part_union )
  : SelectorInterface(),
    m_intersection(),
    m_union(),
    m_all_parts( all_parts( part_union ) )
{
  verify_meta_data( required_part , part_union );
  copy_ids( m_intersection , required_part );
  copy_ids( m_union , part_union );
}

Selector::Selector( const PartVector & part_intersection ,
                    const PartVector & part_union )
  : SelectorInterface(),
    m_intersection(),
    m_union(),
    m_all_parts( all_parts( part_union ) )
{
  if ( ! part_intersection.empty() ) {
    verify_meta_data( * part_intersection[0] , part_intersection );
    verify_meta_data( * part_intersection[0] , part_union );
  }
  else if ( ! part_union.empty() ) {
    verify_meta_data( * part_union[0] , part_union );
  }
  copy_ids( m_intersection , part_intersection );
  copy_ids( m_union , part_union );
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


