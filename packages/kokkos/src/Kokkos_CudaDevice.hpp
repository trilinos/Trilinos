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

#ifndef KOKKOS_CUDADEVICE_HPP
#define KOKKOS_CUDADEVICE_HPP

#include <string>
#include <iosfwd>
#include <Kokkos_MDArrayViewRawData.hpp>

#undef  KOKKOS_DEVICE_FUNCTION
#define KOKKOS_DEVICE_FUNCTION inline __device__

#undef  KOKKOS_DEVICE_AND_HOST_FUNCTION
#define KOKKOS_DEVICE_AND_HOST_FUNCTION inline __device__ __host__

namespace Kokkos {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class CudaDeviceImpl ;

class CudaDevice {
public:
  typedef size_t size_type ;

  static CudaDevice & singleton();

  template< typename ValueType >
  void allocate( MDArrayViewRawData<ValueType,CudaDevice> & array ,
                 const std::string & label )
    {
      if ( array.m_rank ) {
        size_type count = 1 ;
        for ( size_type r = 0 ; r < array.m_rank ; ++r ) {
          count *= array.m_dimension[r] ;
        }
        array.m_ptr_on_device = reinterpret_cast<ValueType*>(
          allocate_memory( sizeof(ValueType) , count , label ) );
      }
    }

  template< typename ValueType >
  void deallocate( MDArrayViewRawData<ValueType,CudaDevice> & array )
    { deallocate_memory( array.m_ptr_on_device ); }

  void print_allocations( std::ostream & ) const ;

private:

  void * allocate_memory( size_type member_size ,
                          size_type member_count ,
                          const std::string & label );

  void deallocate_memory( void * );

  CudaDeviceImpl * m_impl ;

  CudaDevice();
  CudaDevice( const CudaDevice & );
  CudaDevice & operator = ( const CudaDevice & );
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* KOKKOS_CUDADEVICE_HPP */


