/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <ostream>
#include <Kokkos_BaseMappedArray.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

void BaseMappedArray::require_not_owned_view() const
{
  if ( m_map_on_host ) {
    std::string msg("Kokkos::BasedMappedArray::require_not_owned_view FAILED");
    throw std::logic_error( msg );
  }
}

void BaseMappedArray::require_not_allocated() const
{
  if ( m_ptr_on_device ) {
    std::string msg("Kokkos::BasedMappedArray::require_not_allocated FAILED");
    throw std::runtime_error( msg );
  }
}

//----------------------------------------------------------------------------

BaseMappedArray::BaseMappedArray()
  : m_map_on_host( NULL )
  , m_next_on_host( this )
  , m_ptr_on_device( NULL )
  , m_rank( 0 )
{
  m_dimension[0] = 0 ;
  m_dimension[1] = 0 ;
  m_dimension[2] = 0 ;
  m_dimension[3] = 0 ;
  m_dimension[4] = 0 ;
  m_dimension[5] = 0 ;
  m_dimension[6] = 0 ;
  m_dimension[7] = 0 ;
}

  /** \brief  Construct a view of the array */
BaseMappedArray::BaseMappedArray( const BaseMappedArray & rhs )
  : m_map_on_host( NULL )
  , m_next_on_host(  rhs.m_next_on_host )
  , m_ptr_on_device( rhs.m_ptr_on_device )
  , m_rank(          rhs.m_rank )
{
  m_dimension[0] = rhs.m_dimension[0] ;
  m_dimension[1] = rhs.m_dimension[1] ;
  m_dimension[2] = rhs.m_dimension[2] ;
  m_dimension[3] = rhs.m_dimension[3] ;
  m_dimension[4] = rhs.m_dimension[4] ;
  m_dimension[5] = rhs.m_dimension[5] ;
  m_dimension[6] = rhs.m_dimension[6] ;
  m_dimension[7] = rhs.m_dimension[7] ;

  const_cast<BaseMappedArray&>(rhs).m_next_on_host = this ;
}

BaseMappedArray &
BaseMappedArray::operator = ( const BaseMappedArray & rhs )
{
  require_not_owned_view();

  prev_on_host()->m_next_on_host = m_next_on_host ;

  m_next_on_host  = rhs.m_next_on_host ;
  m_ptr_on_device = rhs.m_ptr_on_device ;
  m_rank          = rhs.m_rank ;
  m_dimension[0]  = rhs.m_dimension[0] ;
  m_dimension[1]  = rhs.m_dimension[1] ;
  m_dimension[2]  = rhs.m_dimension[2] ;
  m_dimension[3]  = rhs.m_dimension[3] ;
  m_dimension[4]  = rhs.m_dimension[4] ;
  m_dimension[5]  = rhs.m_dimension[5] ;
  m_dimension[6]  = rhs.m_dimension[6] ;
  m_dimension[7]  = rhs.m_dimension[7] ;

  const_cast<BaseMappedArray&>(rhs).m_next_on_host = this ;
  
  return *this ;
}

BaseMappedArray::~BaseMappedArray()
{
  prev_on_host()->m_next_on_host = m_next_on_host ;

  m_map_on_host   = NULL ;
  m_next_on_host  = this ;
  m_ptr_on_device = NULL ;
  m_rank          = 0 ;
  m_dimension[0]  = 0 ;
  m_dimension[1]  = 0 ;
  m_dimension[2]  = 0 ;
  m_dimension[3]  = 0 ;
  m_dimension[4]  = 0 ;
  m_dimension[5]  = 0 ; 
  m_dimension[6]  = 0 ;
  m_dimension[7]  = 0 ;
}

//----------------------------------------------------------------------------

BaseMapInterface *
BaseMappedArray::base_map() const
{
  require_not_owned_view();

  for ( BaseMappedArray *
        ptr = m_next_on_host ; this != ptr ; ptr = ptr->m_next_on_host ) {
    if ( ptr->m_map_on_host ) { return ptr->m_map_on_host ; }
  }

  return NULL ;
}

BaseMappedArray *
BaseMappedArray::prev_on_host() const
{
  BaseMappedArray * ptr = m_next_on_host ;
  for ( ; this != ptr->m_next_on_host ; ptr = ptr->m_next_on_host );
  return ptr ;
}

//----------------------------------------------------------------------------

void * BaseMappedArray::clear( std::list< BaseMappedArray > & mapped_arrays )
{
  void * const ptr = m_ptr_on_device ;

  if ( ptr ) {

    m_ptr_on_device = NULL ;
    m_rank          = 0 ;
    m_dimension[0]  = 0 ;
    m_dimension[1]  = 0 ;
    m_dimension[2]  = 0 ;
    m_dimension[3]  = 0 ;
    m_dimension[4]  = 0 ;
    m_dimension[5]  = 0 ;
    m_dimension[6]  = 0 ;
    m_dimension[7]  = 0 ;

    for ( BaseMappedArray *
          a = m_next_on_host ; a != this ; a = a->m_next_on_host ) {
      a->m_ptr_on_device = NULL ;
      a->m_rank          = 0 ;
      a->m_dimension[0]  = 0 ;
      a->m_dimension[1]  = 0 ;
      a->m_dimension[2]  = 0 ;
      a->m_dimension[3]  = 0 ;
      a->m_dimension[4]  = 0 ;
      a->m_dimension[5]  = 0 ;
      a->m_dimension[6]  = 0 ;
      a->m_dimension[7]  = 0 ;
    }

    { // Destroy the owned view, any other views are now not-owned

      BaseMappedArray * owned = this ;

      if ( NULL == owned->m_map_on_host ) {
        for ( owned = m_next_on_host ;
              owned != this && NULL == owned->m_map_on_host ;
              owned = owned->m_next_on_host );
      }

      std::list< BaseMappedArray >::iterator i = mapped_arrays.begin();

      while ( i != mapped_arrays.end() && owned != & *i ) { ++i ; }

      if ( i == mapped_arrays.end() ) {
        std::string msg("Kokkos::BaseMappedArray::clear FAILED");
        throw std::logic_error( msg );
      }

      mapped_arrays.erase( i );
    }
  }

  return ptr ;
}

void BaseMappedArray::assign(
  std::list< BaseMappedArray > & mapped_arrays ,
  BaseMapInterface * const map ,
  void * const ptr_on_device , 
  const BaseMappedArray::size_type rank ,
  const BaseMappedArray::size_type dimension[] )
{
  require_not_allocated();

  // Create the owned view of the input array
  mapped_arrays.push_back( *this );

  BaseMappedArray & owned = mapped_arrays.back();

  // Assigned the owned view to this host map,

  owned.m_map_on_host   = map ;
  owned.m_ptr_on_device = ptr_on_device ;
  owned.m_rank          = rank ;
  owned.m_dimension[0]  = dimension[0] ;
  owned.m_dimension[1]  = dimension[1] ;
  owned.m_dimension[2]  = dimension[2] ;
  owned.m_dimension[3]  = dimension[3] ;
  owned.m_dimension[4]  = dimension[4] ;
  owned.m_dimension[5]  = dimension[5] ;
  owned.m_dimension[6]  = dimension[6] ;
  owned.m_dimension[7]  = dimension[7] ;

  for ( BaseMappedArray *
        a = owned.m_next_on_host ; a != & owned ; a = a->m_next_on_host ) {
    a->m_ptr_on_device = ptr_on_device ;
    a->m_rank          = rank ;
    a->m_dimension[0]  = dimension[0] ;
    a->m_dimension[1]  = dimension[1] ;
    a->m_dimension[2]  = dimension[2] ;
    a->m_dimension[3]  = dimension[3] ;
    a->m_dimension[4]  = dimension[4] ;
    a->m_dimension[5]  = dimension[5] ;
    a->m_dimension[6]  = dimension[6] ;
    a->m_dimension[7]  = dimension[7] ;
  }
}

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


