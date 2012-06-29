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
  output_type create( const std::string & label , size_t length , size_t count )
  {
    output_type output ;

    output.m_memory = KokkosArray::create< typename output_type::view_type >( label , length , count );

    return output ;
  }
};

template< typename ValueType , class DeviceDst , class DeviceSrc >
struct Factory< MultiVector< ValueType , DeviceDst > ,
                MultiVector< ValueType , DeviceSrc > >
{
  typedef MultiVector< ValueType , DeviceDst > output_type ;
  typedef MultiVector< ValueType , DeviceSrc > input_type ;

  inline
  static void deep_copy( const output_type & output ,
                         const input_type  & input )
  {
    KokkosArray::deep_copy( output.m_memory , input.m_memory );
  }

  inline static 
  output_type create( const input_type & input )
  {
    output_type output ;
    output.m_memory = Factory< typename output_type::view_type ,
                               typename input_type::view_type >
                       ::create( input.m_memory );
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

