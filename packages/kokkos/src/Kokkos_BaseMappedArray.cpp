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

inline
void BaseMappedArray::require_not_owned_view() const
{
  if ( m_map_on_host ) {
    std::string msg("Kokkos::BasedMappedArray::require_not_owned_view FAILED");
    throw std::logic_error( msg );
  }
}

inline
void BaseMappedArray::require_next_exists() const
{
  if ( NULL == m_next_on_host ) {
    std::string msg("Kokkos::BasedMappedArray::require_next_exists FAILED");
    throw std::logic_error( msg );
  }
}

//----------------------------------------------------------------------------

inline
BaseMappedArray *
BaseMappedArray::next_on_host() const
{
  require_next_exists();
  return m_next_on_host ;
}

// Only non-owned views are ever queried for the owned view.
// Only the owned view has a non-null 'm_map'.
BaseMappedArray *
BaseMappedArray::owned_view()
{
  require_not_owned_view();

  for ( BaseMappedArray *
        ptr = next_on_host() ; this != ptr ; ptr = ptr->next_on_host() ) {
    if ( ptr->m_map_on_host ) { return ptr ; }
  }

  return NULL ;
}

// Only non-owned views are ever queried for the map.
// Only the owned view has a non-null 'm_map'.
BaseMappedArrayRepository *
BaseMappedArray::repository() const
{
  require_not_owned_view();

  for ( BaseMappedArray *
        ptr = next_on_host() ; this != ptr ; ptr = ptr->next_on_host() ) {
    if ( ptr->m_map_on_host ) { return ptr->m_map_on_host ; }
  }

  return NULL ;
}

BaseMappedArray *
BaseMappedArray::prev_on_host() const
{
  BaseMappedArray * ptr = next_on_host() ;
  for ( ; this != ptr->m_next_on_host ; ptr = ptr->next_on_host() );
  return ptr ;
}

void BaseMappedArray::clear()
{
  m_map_on_host   = NULL ;
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
        a = next_on_host() ; a != this ; a = a->next_on_host() ) {
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
}

void BaseMappedArray::assign(
  BaseMappedArrayRepository * const map , void * const ptr_on_device , 
  const BaseMappedArray::size_type rank ,
  const BaseMappedArray::size_type dimension[] )
{
  m_map_on_host   = map ;
  m_ptr_on_device = ptr_on_device ;
  m_rank          = rank ;
  m_dimension[0]  = dimension[0] ;
  m_dimension[1]  = dimension[1] ;
  m_dimension[2]  = dimension[2] ;
  m_dimension[3]  = dimension[3] ;
  m_dimension[4]  = dimension[4] ;
  m_dimension[5]  = dimension[5] ;
  m_dimension[6]  = dimension[6] ;
  m_dimension[7]  = dimension[7] ;

  for ( BaseMappedArray *
        a = next_on_host() ; a != this ; a = a->next_on_host() ) {
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

BaseMappedArrayRepository::BaseMappedArrayRepository( BaseDeviceMemory & device ,
                              BaseMappedArrayRepository::size_type parallel_work_count )
  : m_base_device( device )
  , m_allocated_arrays()
  , m_parallel_work_count( parallel_work_count )
{}

BaseMappedArrayRepository::~BaseMappedArrayRepository()
{
  // Deallocate all member arrays.
  for ( std::list< BaseMappedArray >::iterator
        i  = m_allocated_arrays.begin() ;
        i != m_allocated_arrays.end() ; ++i ) {

    void * const ptr = i->pointer_on_device<void>();

    // Assign the pointer and chunk size to zero
    i->clear();

    if ( ptr ) { m_base_device.deallocate( ptr ); }
  }

  m_allocated_arrays.clear();
}

void BaseMappedArrayRepository::deallocate( BaseMappedArray & array )
{
  std::list< BaseMappedArray >::iterator
    i = require_array_is_member( array );

  void * const ptr = i->pointer_on_device<void>();

  // Assign the pointer and chunk size to zero
  i->clear();

  if ( ptr ) { m_base_device.deallocate( ptr ); }

  // Destroy the owned view, any other views are now not-owned
  m_allocated_arrays.erase( i );
}

void BaseMappedArrayRepository::allocate( BaseMappedArray              & array ,
                              BaseMappedArrayRepository::size_type       sizeof_value ,
                              BaseMappedArrayRepository::size_type       rank ,
                              const BaseMappedArrayRepository::size_type dimension[] )
{
  require_array_is_not_owned( array );

  size_type dim[8] ;

  size_type n = 1 ;

  for ( size_type i = 0 ; i < rank - 1 ; ++i ) {
    n *= ( dim[i] = dimension[i] );
  }

  dim[ rank - 1 ] = m_parallel_work_count ;

  // Allocation and initialization

  void * const pointer =
    n ? m_base_device.allocate( sizeof_value , n , m_parallel_work_count )
      : NULL ;

  // Create an owned view of the input array
  m_allocated_arrays.push_back( array );

  // Assigned the owned view to this host map,
  // assignment operator pushes the pointer and chunking to
  // the other views.
  m_allocated_arrays.back().assign( this , pointer , rank , dim );
}

void BaseMappedArrayRepository::require_array_is_not_owned( BaseMappedArray & array )
{
  if ( NULL != array.owned_view() ) {
    std::string msg( "Kokkos::BaseMappedArrayRepository::allocate FAILED required condition that array is not already allocated" );
    throw std::runtime_error( msg );
  }
}

std::list< BaseMappedArray >::iterator
BaseMappedArrayRepository::require_array_is_member( BaseMappedArray & array )
{
  BaseMappedArray * const owned = array.owned_view();

  std::list< BaseMappedArray >::iterator i = m_allocated_arrays.begin();

  while ( i != m_allocated_arrays.end() && owned != & *i ) { ++i ; }

  if ( i == m_allocated_arrays.end() ) {
    std::string msg("Kokkos::BaseMappedArrayRepository::deallocate FAILED require_array_is_member");
     throw std::logic_error( msg );
  }

  return i ;
}

} // namespace Kokkos


