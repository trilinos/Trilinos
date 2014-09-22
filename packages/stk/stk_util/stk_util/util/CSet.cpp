// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_util/util/CSet.hpp>
#include <assert.h>                     // for assert
#include <stddef.h>                     // for size_t, NULL
#include <algorithm>                    // for lower_bound
#include <iterator>                     // for advance


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


std::vector< Manager >::iterator
lower_bound( std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::iterator i = v.begin();
  std::vector< Manager >::iterator j = v.end();

  return std::lower_bound( i , j , t , less_cset() );
}

}

struct equal_cset {
  bool operator()(const Manager& lhs, const std::type_info& rhs) const
  { return *lhs.first == rhs; }
  bool operator()(const std::type_info& lhs, const Manager& rhs) const
  { return lhs == *rhs.first; }
};

//----------------------------------------------------------------------

const void * CSet::p_get( const std::type_info & t ) const
{
  for(std::vector<Manager>::const_iterator it=m_manager.begin(), end=m_manager.end(); it!=end; ++it) {
    if (*it->first == t) return m_value[it-m_manager.begin()];
  }

  return NULL ;
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


