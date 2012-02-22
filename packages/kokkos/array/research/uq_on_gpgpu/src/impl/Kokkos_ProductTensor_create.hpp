/*
//@HEADER
// ************************************************************************
// 
//                         Kokkos Array
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_PRODUCTTENSOR_CREATE_HPP
#define KOKKOS_PRODUCTTENSOR_CREATE_HPP

#include <map>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_MDArray.hpp>
#include <Kokkos_CrsMap.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Create a sparse product tensor on the device
 *          from a map of tensor indices to values.
 *
 *  The std::map input guarantees uniqueness and proper sorting of
 *  the product tensor's symmetric entries.
 */
template< unsigned Rank , typename ValueType , class Device , class D >
class CreateSparseProductTensor<
  SparseProductTensor< Rank , ValueType , Device > ,
  std::map< ProductTensorIndex<Rank,D> , ValueType > >
{
public:
  typedef SparseProductTensor<Rank,ValueType,Device> type ;
  typedef std::map< ProductTensorIndex< Rank , D > , ValueType > input_type ;

  static
  type create( const input_type & input )
  {
    typedef ValueType                   value_type ;
    typedef Device                      device_type ;
    typedef typename Device::size_type  size_type ;

    typedef MDArray< size_type, device_type>  coord_array_type ;
    typedef MultiVector< value_type, device_type> value_array_type ;

    enum { is_host_memory =
             SameType< typename device_type::memory_space , Host >::value };

    type tensor ;

    tensor.m_coord = create_mdarray< size_type , device_type >( input.size(), 3 );
    tensor.m_value = create_multivector< value_type, device_type >( input.size() );

    // Create mirror, is a view if is host memory

    typename coord_array_type::HostMirror
      host_coord = Kokkos::Impl::CreateMirror< coord_array_type , is_host_memory >::create( tensor.m_coord );

    typename value_array_type::HostMirror
      host_value = Kokkos::Impl::CreateMirror< value_array_type , is_host_memory >::create( tensor.m_value );

    size_type n = 0 ;

    tensor.m_dimen = 0 ;

    for ( typename input_type::const_iterator
          iter = input.begin() ; iter != input.end() ; ++iter , ++n ) {

      host_value(n) = (*iter).second ;

      for ( size_type c = 0 ; c < Rank ; ++c ) {

        const size_type i = (*iter).first.coord(c);

        host_coord(n,c) = i ;

        tensor.m_dimen = std::max( tensor.m_dimen , i + 1 );
      }
    }

    Kokkos::deep_copy( tensor.m_coord , host_coord );
    Kokkos::deep_copy( tensor.m_value , host_value );

    return tensor ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_PRODUCTTENSOR_CREATE_HPP */


