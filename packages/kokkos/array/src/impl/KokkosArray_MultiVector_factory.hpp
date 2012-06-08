/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_MULTIVECTOR_FACTORY_HPP
#define KOKKOS_IMPL_MULTIVECTOR_FACTORY_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< typename ValueType , class Device >
struct Factory< MultiVector< ValueType , Device > , void >
{
  typedef MultiVector< ValueType , Device > output_type ;

  static inline
  output_type create( const std::string & label ,
                      size_t length ,
                      size_t count ,
                      size_t stride = 0 )
  {
    typedef Impl::MemoryManager<typename Device::memory_space> memory_manager ;

    if ( 0 == stride ) {
      stride = length ;
      if ( 1 < count ) {
        stride = memory_manager::template preferred_alignment<ValueType>(length);
      }
    }

    output_type output ;

    output.m_length = length ;
    output.m_count  = count ;
    output.m_stride = stride ;
    output.m_memory.allocate( output.m_count * output.m_stride , label );
    output.m_ptr_on_device = output.m_memory.ptr_on_device();

    return output ;
  }
};

template< typename ValueType , class DeviceOutput >
struct Factory< MultiVector< ValueType , DeviceOutput > , MirrorUseView >
{
  typedef MultiVector< ValueType , DeviceOutput > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const MultiVector<ValueType,DeviceInput> & input )
  {
    typedef MultiVector<ValueType,DeviceInput> input_type ;
    return Factory<output_type,input_type>::create( input );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_IMPL_MULTIVECTOR_FACTORY_HPP */

