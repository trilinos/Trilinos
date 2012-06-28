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

#ifndef KOKKOS_IMPL_PREFIXSUM_FACTORY_HPP
#define KOKKOS_IMPL_PREFIXSUM_FACTORY_HPP

#include <vector>
#include <impl/KokkosArray_MemoryView.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename IntType , class DeviceDst , class DeviceSrc >
struct Factory< PrefixSum< IntType , DeviceDst > ,
                PrefixSum< IntType , DeviceSrc > >
{
  typedef PrefixSum< IntType, DeviceDst > output_type ;
  typedef PrefixSum< IntType, DeviceSrc > input_type ;

  typedef typename output_type::view_type output_view_type ;
  typedef typename input_type::view_type  input_view_type ;

  typedef Factory< output_view_type , input_view_type > view_factory ;

  static inline
  void deep_copy( output_type & output , const input_type & input )
  {
    view_factory::deep_copy( output.m_data , input.m_data );
    output.m_sum = input.m_sum ;
  }

  static inline
  output_type create( const input_type & input )
  {
    output_type output ;

    output.m_data = view_factory::create( input.m_data );

    view_factory::deep_copy( output.m_data , input.m_data );

    output.m_sum = input.m_sum ;

    return output ;
  }
};

//----------------------------------------------------------------------------

template< typename IntType , class DeviceOutput >
struct Factory< PrefixSum< IntType , DeviceOutput > , MirrorUseView >
{
  typedef PrefixSum< IntType , DeviceOutput > output_type ;

  static inline
  const output_type & create( const output_type & input ) { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const PrefixSum< IntType , DeviceInput > & input )
  {
    typedef PrefixSum< IntType , DeviceInput > input_type ;
    return Factory< output_type , input_type >::create( input );
  }
};

//----------------------------------------------------------------------------

template< typename IntTypeOutput , class DeviceOutput ,
          typename IntTypeInput >
struct Factory< PrefixSum< IntTypeOutput , DeviceOutput > ,
                std::vector< IntTypeInput > >
{
  typedef PrefixSum< IntTypeOutput , DeviceOutput > output_type ;
  typedef std::vector< IntTypeInput > input_type ;

  static
  output_type create( const std::string & label , const input_type & input )
  {
    typedef typename output_type::view_type output_view_type ;
    typedef typename output_view_type::HostMirror  mirror_view_type ;

    typedef Factory< output_view_type , unsigned_<1> > output_view_factory ;

    typedef MemoryView< IntTypeOutput ,
                        typename DeviceOutput::memory_space > memory_output ;

    typedef typename memory_output::HostMirror  memory_mirror ;

    const size_t count = input.size();

    output_type output ;

    output.m_data = output_view_factory::create( label , count + 1 );

    // If same memory space then a view:
    typename output_view_type::HostMirror tmp =
      create_mirror( output.m_data , MirrorUseView() );

    // Could be made parallel
    output.m_sum = tmp[0] = 0 ;
    for ( size_t i = 0 ; i < count ; ++i ) {
      tmp[i+1] = output.m_sum += input[i] ;
    }

    deep_copy( output.m_data , tmp );

    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_PREFIXSUM_FACTORY_HPP */

