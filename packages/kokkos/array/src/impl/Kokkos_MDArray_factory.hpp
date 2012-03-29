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

#ifndef KOKKOS_IMPL_MDARRAY_FACTORY_HPP
#define KOKKOS_IMPL_MDARRAY_FACTORY_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename ValueType , class Device >
struct Factory< MDArray< ValueType , Device > , void >
{
  typedef MDArray< ValueType , Device > output_type ;

  static output_type create( const std::string & label ,
                             size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                             size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    typedef DeepCopyKernelMDArray< output_type, ValueType, 1 > kernel_1 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 2 > kernel_2 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 3 > kernel_3 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 4 > kernel_4 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 5 > kernel_5 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 6 > kernel_6 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 7 > kernel_7 ;
    typedef DeepCopyKernelMDArray< output_type, ValueType, 8 > kernel_8 ;

    output_type array ;

    array.m_map.template assign< ValueType >(nP,n1,n2,n3,n4,n5,n6,n7);
    array.m_memory.allocate( array.m_map.allocation_size() , label );

    switch( array.m_map.rank() ) {
    case 1 : parallel_for( nP , kernel_1( array , 0 ) ); break ;
    case 2 : parallel_for( nP , kernel_2( array , 0 ) ); break ;
    case 3 : parallel_for( nP , kernel_3( array , 0 ) ); break ;
    case 4 : parallel_for( nP , kernel_4( array , 0 ) ); break ;
    case 5 : parallel_for( nP , kernel_5( array , 0 ) ); break ;
    case 6 : parallel_for( nP , kernel_6( array , 0 ) ); break ;
    case 7 : parallel_for( nP , kernel_7( array , 0 ) ); break ;
    case 8 : parallel_for( nP , kernel_8( array , 0 ) ); break ;
    }

    return array ;
  }
};

//----------------------------------------------------------------------------

/** \brief Mirror with view optimization */
template< typename ValueType , class Device >
struct Factory< MDArray< ValueType , Device > , Impl::MirrorUseView >
{
  typedef MDArray< ValueType , Device > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const MDArray< ValueType , DeviceInput > & input )
  {
    typedef MDArray< ValueType , DeviceInput > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_IMPL_MDARRAY_FACTORY_HPP */


