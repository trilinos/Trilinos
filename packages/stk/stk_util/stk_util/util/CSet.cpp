/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <assert.h>

#include <stk_util/util/CSet.hpp>

namespace stk {

namespace {

typedef void (* DeleteFunction )( void * );

typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

// Comparison for sorted vector

struct less_cset {
  bool operator()( const Manager        & lhs ,
                   const std::type_info & rhs ) const ;
  bool operator()( const std::type_info & lhs ,
                   const Manager        & rhs ) const ;
};

// On some systems, namely AIX, std::type_info::before(...)
// has a bug where it returns true instead of false for equality.
// Thus we pay a small price on all systems to specifically
// test for and eliminate equality.

bool less_cset::operator()( const Manager        & lhs ,
                            const std::type_info & rhs ) const
{ return lhs.first->before( rhs ) && * lhs.first != rhs ; }

bool less_cset::operator()( const std::type_info & lhs ,
                            const Manager        & rhs ) const
{ return lhs.before( *rhs.first ) && lhs != *rhs.first ; }


std::vector< Manager >::const_iterator
lower_bound( const std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::const_iterator i = v.begin();
  std::vector< Manager >::const_iterator j = v.end();

  return std::lower_bound( i , j , t , less_cset() );
}

std::vector< Manager >::iterator
lower_bound( std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::iterator i = v.begin();
  std::vector< Manager >::iterator j = v.end();

  return std::lower_bound( i , j , t , less_cset() );
}

}

//----------------------------------------------------------------------

const void * CSet::p_get( const std::type_info & t ) const
{
  const void * result = NULL ;

  const std::vector< Manager >::const_iterator im = lower_bound(m_manager,t);

  if ( im < m_manager.end() && t == * im->first ) {
    const size_t offset = im - m_manager.begin();
    result = m_value[ offset ];
  }

  return result ;
}

const void *
CSet::p_insert( const Manager & m , const void * v )
{
  std::vector< Manager >::iterator im = lower_bound( m_manager , * m.first );

  const size_t offset = im - m_manager.begin();

  assert(m_value.size() == m_manager.size());
  std::vector<const void *>::iterator iv = m_value.begin();
  std::advance( iv , offset );

  if ( im == m_manager.end() || * m.first != * im->first ) {
    im = m_manager.insert( im , m );
    iv = m_value  .insert( iv , v );
  }

  assert(iv != m_value.end());
  return *iv ;
}

bool CSet::p_remove( const std::type_info & t , const void * v )
{
  bool result = false;
  const std::vector< Manager >::iterator im = lower_bound( m_manager , t );

  if (im != m_manager.end()) {
    const size_t offset = im - m_manager.begin();

    if (offset <= m_value.size()) {
      std::vector<const void *>::iterator iv = m_value.begin();
      std::advance( iv , offset );

      result = t == * im->first && v == * iv ;

      if ( result ) {
	m_manager.erase( im );
	m_value  .erase( iv );
      }
    }
  }
  return result ;
}

//----------------------------------------------------------------------

CSet::~CSet()
{
  try {
    const size_t n = m_manager.size();
    for ( size_t i = 0 ; i < n ; ++i ) {
      try {
        if ( m_manager[i].second ) {
          (*m_manager[i].second)( const_cast<void*>( m_value[i] ) );
        }
      } catch(...) {}
    }
  } catch(...) {}
}

CSet::CSet() : m_manager(), m_value() {}

} // namespace stk


