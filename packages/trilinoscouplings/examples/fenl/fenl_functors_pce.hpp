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

#ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_PCE_HPP
#define KOKKOS_EXAMPLE_FENLFUNCTORS_PCE_HPP

#include <fenl_functors.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

// Specialization of ResponseComputation for PCE scalar types.  We currently
// have to do this because parallel_reduce doesn't support non-POD types
template< typename FixtureType ,
          typename S , typename L , typename D , typename M >
class ResponseComputation< FixtureType,
                           Kokkos::View<S,L,D,M,Kokkos::Impl::ViewPCEContiguous>
                         >
{
public:

  typedef FixtureType fixture_type ;
  typedef Kokkos::View<S,L,D,M,Kokkos::Impl::ViewPCEContiguous> vector_type ;
  typedef typename vector_type::device_type device_type ;
  typedef typename vector_type::value_type scalar_type ;

  // Hack to get parallel_reduce to work
  typedef typename vector_type::intrinsic_scalar_type value_type[] ;
  typedef typename vector_type::cijk_type cijk_type ;

  typedef Kokkos::Example::HexElement_Data< fixture_type::ElemNode > element_data_type ;
  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;
  static const unsigned ElemNodeCount    = element_data_type::element_node_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;

  //------------------------------------
  // Computational data:

  const element_data_type    elem_data ;
  const fixture_type         fixture ;
  const vector_type          solution ;
  const unsigned             value_count ;
  const cijk_type            cijk ;

  ResponseComputation( const ResponseComputation & rhs )
    : elem_data()
    , fixture( rhs.fixture )
    , solution( rhs.solution )
    , value_count( rhs.value_count )
    , cijk( rhs.cijk )
    {}

  ResponseComputation( const fixture_type& arg_fixture ,
                       const vector_type & arg_solution )
    : elem_data()
    , fixture( arg_fixture )
    , solution( arg_solution )
    , value_count( solution.sacado_size() )
    , cijk( solution.cijk() )
    {}

  //------------------------------------

  scalar_type apply() const
  {
    scalar_type response(cijk, value_count);
    //Kokkos::parallel_reduce( fixture.elem_count() , *this , response.coeff() );
    Kokkos::parallel_reduce( solution.dimension_0() , *this , response.coeff() );
    return response;
  }

  //------------------------------------

   KOKKOS_INLINE_FUNCTION
  float compute_detJ(
    const float grad[][ ElemNodeCount ] , // Gradient of bases master element
    const double x[] ,
    const double y[] ,
    const double z[] ) const
  {
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Jacobian accumulation:

    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( unsigned i = 0; i < ElemNodeCount ; ++i ) {
      const double x1 = x[i] ;
      const double x2 = y[i] ;
      const double x3 = z[i] ;

      const float g1 = grad[0][i] ;
      const float g2 = grad[1][i] ;
      const float g3 = grad[2][i] ;

      J[j11] += g1 * x1 ;
      J[j12] += g1 * x2 ;
      J[j13] += g1 * x3 ;

      J[j21] += g2 * x1 ;
      J[j22] += g2 * x2 ;
      J[j23] += g2 * x3 ;

      J[j31] += g3 * x1 ;
      J[j32] += g3 * x2 ;
      J[j33] += g3 * x3 ;
    }

    // Inverse jacobian:

    float invJ[ TensorDim ] = {
      static_cast<float>( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      static_cast<float>( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      static_cast<float>( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      static_cast<float>( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      static_cast<float>( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      static_cast<float>( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      static_cast<float>( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      static_cast<float>( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      static_cast<float>( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const float detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    return detJ ;
  }

  KOKKOS_INLINE_FUNCTION
  scalar_type contributeResponse(
    const scalar_type dof_values[] ,
    const float  detJ ,
    const float  integ_weight ,
    const float  bases_vals[] ) const
  {
    // $$ g_i = \int_{\Omega} T^2 d \Omega $$

    scalar_type value_at_pt = 0 ;
    for ( unsigned m = 0 ; m < ElemNodeCount ; m++ ) {
      value_at_pt += dof_values[m] * bases_vals[m] ;
    }

    scalar_type elem_response =
      value_at_pt * value_at_pt * detJ * integ_weight ;

    return elem_response;
  }

  /*
  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem , value_type response ) const
  {
    // Gather nodal coordinates and solution vector:

    double x[ ElemNodeCount ] ;
    double y[ ElemNodeCount ] ;
    double z[ ElemNodeCount ] ;
    scalar_type val[ ElemNodeCount ] ;

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = fixture.elem_node( ielem , i );

      x[i] = fixture.node_coord( ni , 0 );
      y[i] = fixture.node_coord( ni , 1 );
      z[i] = fixture.node_coord( ni , 2 );

      val[i] = solution( ni );
    }

    scalar_type response_pce( cijk, value_count, response, false );

    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const float detJ = compute_detJ( elem_data.gradients[i] , x , y , z );

      response_pce += contributeResponse( val , detJ , elem_data.weights[i] ,
                                          elem_data.values[i] );
    }
  }
  */

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned i , value_type response ) const
  {
    const scalar_type& u = solution(i);
    scalar_type response_pce( cijk, value_count, response, false );
    response_pce += (u * u) / fixture.node_count_global();
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type response ) const
  {
    for (unsigned i=0; i<value_count; ++i)
      response[i] = 0 ;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type  response ,
             volatile const value_type input ) const
  {
    for (unsigned i=0; i<value_count; ++i)
      response[i] += input[i] ;
  }

}; /* ResponseComputation */

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template <typename Storage>
inline
Sacado::UQ::PCE<Storage>
all_reduce( const Sacado::UQ::PCE<Storage>& local ,
            const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
  const int sz = local.size();
  Sacado::UQ::PCE<Storage> global(local.cijk(), sz) ;
  Teuchos::reduceAll(
    *comm , Teuchos::REDUCE_SUM , sz , local.coeff() , global.coeff() );
  return global ;
}

} // namespace Example
} // namespace Kokkos

#endif /* #ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_PCE_HPP */
