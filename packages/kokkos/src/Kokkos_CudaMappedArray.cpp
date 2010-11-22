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
#include <Kokkos_CudaMappedArray.hpp>

namespace Kokkos {

CudaMappedArray *
CudaMappedArray::owned_view()
{
  for ( CudaMappedArray *
        ptr = m_next_on_host ; this != ptr ; ptr = ptr->m_next_on_host ) {
    if ( ptr->m_map_on_host ) { return ptr ; }
  }
  return NULL ;
}

CudaMappedArray::map_type *
CudaMappedArray::map() const
{
  for ( const CudaMappedArray *
        ptr = m_next_on_host ; this != ptr ; ptr = ptr->m_next_on_host ) {
    if ( ptr->m_map_on_host ) { return ptr->m_map_on_host ; }
  }
  return NULL ;
}

void CudaMappedArray::assign(
  map_type * map , void * pointer ,
  CudaMappedArray::size_type work_count ,
  CudaMappedArray::size_type chunk_count )
{
  m_map_on_host   = map ;
  m_ptr_on_device = pointer ;
  m_work_count    = work_count ;
  m_chunk_count   = chunk_count ;

  for ( CudaMappedArray *
        a = m_next_on_host ; a != this ; a = a->m_next_on_host ) {
    a->m_ptr_on_device = pointer ;
    a->m_work_count    = work_count ;
    a->m_chunk_count   = chunk_count ;
  }
}

/*--------------------------------------------------------------------------*/

CudaMap::CudaMap(
  CudaMap::device_type & device ,
  CudaMap::size_type parallel_work_length )
  : m_device( device )
  , m_allocated_arrays()
  , m_parallel_work_length( parallel_work_length )
{}

CudaMap::~CudaMap()
{
  // Deallocate all member arrays.
  for ( std::list< mapped_array_type >::iterator
        i  = m_allocated_arrays.begin() ;
        i != m_allocated_arrays.end() ; ++i ) {

    void * const ptr = i->m_ptr_on_device ;

    // Assign the pointer and chunk size to zero
    i->assign( NULL , NULL , 0 , 0 );
    
    if ( ptr ) { cudaFree( ptr ); }
  }

  m_allocated_arrays.clear();
}

void CudaMap::deallocate( CudaMap::mapped_array_type & array )
{
  std::list< mapped_array_type >::iterator
    i = require_array_is_member( array );

  void * const ptr = i->m_ptr_on_device ;

  // Assign the pointer and chunk size to zero
  i->assign( NULL , NULL , 0 , 0 );

  if ( ptr ) { cudaFree( ptr ); }

  // Destroy the owned view, any other views are now not-owned
  m_allocated_arrays.erase( i );
}

void CudaMap::allocate( CudaMap::mapped_array_type & array ,
                        CudaMap::size_type values_per_chunk ,
                        CudaMap::size_type sizeof_value )
{
  require_array_is_not_owned( array );

  // Allocation and initialization
  const size_type n = m_parallel_work_length * values_per_chunk ;

  void * ptr_on_device = NULL ;

  if ( n ) {
    cudaMalloc( & ptr_on_device , sizeof_value * n );
    // Initialize the memory to zero
 
  }

  // Create an owned view of the input array
  m_allocated_arrays.push_back( array );

  // Assigned the owned view to this host map,
  // assignment operator pushes the pointer and chunking to
  // the other views.
  m_allocated_arrays.back().assign( this , ptr_on_device ,
                                    m_parallel_work_length ,
                                    values_per_chunk );
}

void CudaMap::require_array_is_not_owned( CudaMappedArray & array )
{
  if ( NULL != array.owned_view() ) {
    std::string msg( "Kokkos::CudaMap::allocate FAILED required condition that array is not already allocated" );
    throw std::runtime_error( msg );
  }
}

std::list< CudaMap::mapped_array_type >::iterator
CudaMap::require_array_is_member( CudaMappedArray & array )
{
  CudaMappedArray * const owned = array.owned_view();

  std::list< mapped_array_type >::iterator i = m_allocated_arrays.begin();

  while ( i != m_allocated_arrays.end() && owned != & *i ) { ++i ; }

  if ( i == m_allocated_arrays.end() ) {
    std::string msg(" Kokkos::CudaMap::deallocate FAILED required condition that array is member of this host map" );
     throw std::runtime_error( msg );
  }

  return i ;
}

} // namespace Kokkos


