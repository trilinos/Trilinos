/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_FEINT_HPP
#define KOKKOS_EXAMPLE_FEINT_HPP

#include <iostream>
#include <BoxElemFixture.hpp>
#include <ElemFunctor.hpp>
#include <feint_fwd.hpp>

namespace Kokkos {
namespace Example {

struct MyFunctionType {

  enum { value_count = 8 };

  // Integral on a unit cube:
  // { 1 , 1/2 , 1/2 , 1/2 , 1/4 , 1/4 , 1/4 , 1/8 }

  // Evaluate function at coordinate.
  template< typename CoordType , typename ValueType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const CoordType point[] , ValueType value[] ) const
    {
      value[0] = 1 ;
      value[1] = point[0] ;
      value[2] = point[1] ;
      value[3] = point[2] ;
      value[4] = point[0] * point[1] ;
      value[5] = point[1] * point[2] ;
      value[6] = point[2] * point[0] ;
      value[7] = point[0] * point[1] * point[2] ;
    }
};

template < class Device , bool UseAtomic >
void feint(
  const unsigned global_elem_nx ,
  const unsigned global_elem_ny ,
  const unsigned global_elem_nz )
{
  const unsigned mpi_rank = 0 ;
  const unsigned mpi_size = 1 ;

  static const Kokkos::Example::BoxElemPart::Decompose
    decompose = Kokkos::Example::BoxElemPart:: DecomposeElem ; // DecomposeElem | DecomposeNode ;

  // typedef Kokkos::Example::BoxElemFixture< Device , Kokkos::Example::BoxElemPart:: ElemLinear > BoxFixtureType ;
  typedef Kokkos::Example::BoxElemFixture< Device , Kokkos::Example::BoxElemPart:: ElemQuadratic > BoxFixtureType ;

  typedef Kokkos::Example::FiniteElementIntegration< BoxFixtureType , MyFunctionType , UseAtomic > FeintType ;

  typedef Kokkos::Example::LumpElemToNode< typename FeintType::NodeValueType ,
                                           typename FeintType::ElemValueType ,
                                           UseAtomic > LumpType ;

  const BoxFixtureType fixture( decompose , mpi_size , mpi_rank ,
                                global_elem_nx , global_elem_ny , global_elem_nz );

  const FeintType feint( fixture , MyFunctionType() );

  typename FeintType::value_type elem_integral ;

  Kokkos::parallel_reduce( fixture.elem_count() , feint , elem_integral );

  if ( elem_integral.error ) {
    std::cout << "Spatial jacobian error" << std::endl ;
    return ;
  }

  std::cout << "Elem integral =" ;
  for ( int i = 0 ; i < MyFunctionType::value_count ; ++i ) {
    std::cout << " " << elem_integral.value[i] ;
  }
  std::cout << std::endl ;
 
  const LumpType lump( feint.m_node_lumped ,
                       feint.m_elem_integral ,
                       fixture.elem_node() );

  typename LumpType ::value_type node_sum ;

  Kokkos::parallel_reduce( fixture.node_count() , lump , node_sum );

  std::cout << "Node lumped sum =" ;
  for ( int i = 0 ; i < MyFunctionType::value_count ; ++i ) {
    std::cout << " " << node_sum.value[i] ;
  }
  std::cout << std::endl ;
}

} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FEINT_HPP */

