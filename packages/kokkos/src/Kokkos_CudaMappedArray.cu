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
namespace {

void device_free( void * pointer_on_device )
{
  cudaFree( pointer_on_device );
}

void * device_allocate( int sizeof_value ,
                        int chunk_count ,
                        int work_count )
{
  void * pointer_on_device = NULL ;

  cudaMalloc( & pointer_on_device , sizeof_value * work_count * chunk_count );

  return pointer_on_device ;
}

}

CudaMap::CudaMap( CudaMap::size_type parallel_work_count )
  : m_allocated_arrays()
  , m_parallel_work_count( parallel_work_count )
{}

CudaMap::~CudaMap()
{
  while ( ! m_allocated_arrays.empty() ) {
    void * ptr = m_allocated_arrays.back().clear( m_allocated_arrays );
    device_free( ptr );
  }
}

void CudaMap::deallocate( BaseMappedArray & array )
{
  // Clear all views and destroy the owned view
  void * ptr = array.clear( m_allocated_arrays );

  device_free( ptr );
}

void CudaMap::allocate( BaseMappedArray                 & array ,
                        BaseMapInterface::size_type       sizeof_value ,
                        BaseMapInterface::size_type       rank ,
                        const BaseMapInterface::size_type * const dimension )
{
  array.require_not_allocated();

  size_type dim[8] ;

  size_type n = 1 ;

  for ( size_type i = 0 ; i < rank - 1 ; ++i ) {
    n *= ( dim[i] = dimension[i] );
  }

  dim[ rank - 1 ] = m_parallel_work_count ;

  void * const pointer =
    device_allocate( sizeof_value , n , m_parallel_work_count );

  if ( pointer ) {
    array.assign( m_allocated_arrays , this , pointer , rank , dim );
  }
}

} // namespace Kokkos


