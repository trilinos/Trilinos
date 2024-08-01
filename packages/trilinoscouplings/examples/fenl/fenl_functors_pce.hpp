//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_PCE_HPP
#define KOKKOS_EXAMPLE_FENLFUNCTORS_PCE_HPP

#include <fenl_functors.hpp>
#include "Kokkos_Parallel_MP_Vector.hpp"
#if defined( KOKKOS_ENABLE_CUDA )
#include "Stokhos_Cuda_WarpShuffle.hpp"
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

// Specialization of ResponseComputation for PCE scalar types.  We currently
// have to do this because parallel_reduce doesn't support non-POD types
template< typename FixtureType ,
          typename S , typename ... P , typename ... Q >
class ResponseComputation< FixtureType,
                           Kokkos::View<Sacado::UQ::PCE<S>*,P...>,
                           Kokkos::View<Sacado::UQ::PCE<S>**,Q...>
                         >
{
public:

  typedef FixtureType fixture_type ;
  typedef Kokkos::View<Sacado::UQ::PCE<S>*,P...> vector_type ;
  typedef Kokkos::View<Sacado::UQ::PCE<S>**,Q...> multi_vector_type ;
  typedef typename vector_type::execution_space execution_space ;
  typedef typename vector_type::value_type scalar_type ;

  // Hack to get parallel_reduce to work
  typedef typename Kokkos::IntrinsicScalarType<vector_type>::type reduce_scalar_type;
  typedef reduce_scalar_type value_type[] ;
  typedef typename Kokkos::CijkType<vector_type>::type cijk_type ;

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
  const multi_vector_type    solution_dp ;
  const unsigned             pce_size ;
  const unsigned             num_param ;
  const unsigned             value_count ;
  const cijk_type            cijk ;

  ResponseComputation( const ResponseComputation & rhs )
    : elem_data()
    , fixture( rhs.fixture )
    , solution( rhs.solution )
    , solution_dp( rhs.solution_dp )
    , pce_size( rhs.pce_size )
    , num_param( rhs.num_param )
    , value_count( rhs.value_count )
    , cijk( rhs.cijk )
    {}

  ResponseComputation( const fixture_type& arg_fixture ,
                       const vector_type & arg_solution ,
                       const multi_vector_type& arg_solution_dp =
                         multi_vector_type() )
    : elem_data()
    , fixture( arg_fixture )
    , solution( arg_solution )
    , solution_dp( arg_solution_dp )
    , pce_size( Kokkos::dimension_scalar(solution) )
    , num_param( solution_dp.extent(1) )
    , value_count( pce_size*(num_param+1) )
    , cijk( Kokkos::cijk(solution) )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        num_param > 0, std::logic_error, "PCE gradients not implemented!");
    }

  //------------------------------------

  scalar_type apply() const
  {
    scalar_type response(cijk, value_count);
    //Kokkos::parallel_reduce( fixture.elem_count() , *this , response.coeff() );
    Kokkos::View<reduce_scalar_type*, Kokkos::HostSpace> resp_view(response.coeff(), value_count);
    Kokkos::parallel_reduce( solution.extent(0) , *this , resp_view );
    return response;
  }

  Teuchos::Array<scalar_type> apply_gradient() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "PCE gradients not implemented!");
  }

  //------------------------------------

  /*
  KOKKOS_INLINE_FUNCTION
  double compute_detJ(
    const double grad[][ ElemNodeCount ] , // Gradient of bases master element
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

      const double g1 = grad[0][i] ;
      const double g2 = grad[1][i] ;
      const double g3 = grad[2][i] ;

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

    double invJ[ TensorDim ] = {
      static_cast<double>( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      static_cast<double>( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      static_cast<double>( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      static_cast<double>( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      static_cast<double>( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      static_cast<double>( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      static_cast<double>( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      static_cast<double>( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      static_cast<double>( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const double detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    return detJ ;
  }

  KOKKOS_INLINE_FUNCTION
  scalar_type contributeResponse(
    const scalar_type dof_values[] ,
    const double  detJ ,
    const double  integ_weight ,
    const double  bases_vals[] ) const
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

      const double detJ = compute_detJ( elem_data.gradients[i] , x , y , z );

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
  void join( value_type  response ,
             const value_type input ) const
  {
    for (unsigned i=0; i<value_count; ++i)
      response[i] += input[i] ;
  }

}; /* ResponseComputation */

//------------------------------------

template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename ensemble_scalar_type,
           int EnsembleSize,
           typename Device = typename pce_view_type::execution_space >
struct EvaluatePCE {
  typedef Device execution_space;
  typedef typename pce_view_type::array_type pce_array_type;

  const pce_array_type   pce_view;
  const scalar_view_type scalar_view;
  const quad_values_type quad_values;
  unsigned               qp;
  const unsigned         num_pce;

  EvaluatePCE( const pce_view_type&    arg_pce_view,
               const scalar_view_type& arg_scalar_view,
               const quad_values_type& arg_quad_values) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ) {}

  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    Kokkos::parallel_for( pce_view.extent(1), *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const unsigned row ) const {
    ensemble_scalar_type s = 0.0;
    for (unsigned pce=0; pce<num_pce; ++pce)
      for (int i=0; i<EnsembleSize; ++i)
        s.fastAccessCoeff(i) += pce_view(pce,row)*quad_values(qp+i,pce);
    scalar_view(row) = s;
  }
};

template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename quad_weights_type,
           int EnsembleSize,
           typename Device = typename pce_view_type::execution_space >
struct AssemblePCE {
  typedef Device execution_space;
  typedef typename pce_view_type::array_type  pce_array_type;

  const pce_array_type    pce_view;
  const scalar_view_type  scalar_view;
  const quad_values_type  quad_values;
  const quad_weights_type quad_weights;
  unsigned                qp;
  const unsigned          num_pce;

  AssemblePCE( const pce_view_type&     arg_pce_view,
               const scalar_view_type&  arg_scalar_view,
               const quad_values_type&  arg_quad_values,
               const quad_weights_type& arg_quad_weights ) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    quad_weights( arg_quad_weights ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ) {}

  void apply( const unsigned arg_qp ) {
    qp = arg_qp;
    Kokkos::parallel_for( pce_view.extent(1), *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const unsigned row ) const {
    for (unsigned pce=0; pce<num_pce; ++pce) {
      for (int i=0; i<EnsembleSize; ++i)
        pce_view(pce,row) +=
          quad_weights(qp+i)*scalar_view(row).fastAccessCoeff(i)*quad_values(qp+i,pce);
    }
  }
};

template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename quad_weights_type,
           int EnsembleSize,
           typename Device = typename pce_view_type::execution_space >
struct AssembleRightPCE {
  typedef Device execution_space;
  typedef typename pce_view_type::array_type  pce_array_type;

  const pce_array_type    pce_view;
  const scalar_view_type  scalar_view;
  const quad_values_type  quad_values;
  const quad_weights_type quad_weights;
  unsigned                qp;
  const unsigned          num_pce;

  AssembleRightPCE( const pce_view_type&     arg_pce_view,
                    const scalar_view_type&  arg_scalar_view,
                    const quad_values_type&  arg_quad_values,
                    const quad_weights_type& arg_quad_weights ) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    quad_weights( arg_quad_weights ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ) {}

  void apply( const unsigned arg_qp ) {
    qp = arg_qp;
    Kokkos::parallel_for( pce_view.extent(0), *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const unsigned row ) const {
    for (unsigned pce=0; pce<num_pce; ++pce) {
      for (int i=0; i<EnsembleSize; ++i)
        pce_view(row,pce) +=
          quad_weights(qp+i)*scalar_view(row).fastAccessCoeff(i)*quad_values(qp+i,pce);
    }
  }
};

#if defined( KOKKOS_ENABLE_CUDA ) && defined( __CUDACC__ )
template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename ensemble_scalar_type,
           int EnsembleSize >
struct EvaluatePCE< pce_view_type,
                    scalar_view_type,
                    quad_values_type,
                    ensemble_scalar_type,
                    EnsembleSize,
                    Kokkos::Cuda > {
  typedef Kokkos::Cuda execution_space;
  typedef typename pce_view_type::array_type pce_array_type;
  typedef typename pce_array_type::value_type scalar_type;
  typedef typename quad_values_type::value_type quad_scalar_type;

  const pce_array_type   pce_view;
  const scalar_view_type scalar_view;
  const quad_values_type quad_values;
  unsigned               qp;
  const unsigned         num_pce;
  const unsigned         row_count;

  EvaluatePCE( const pce_view_type&    arg_pce_view,
               const scalar_view_type& arg_scalar_view,
               const quad_values_type& arg_quad_values) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ),
    row_count( pce_view.extent(1) ) {}

  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    Kokkos::parallel_for( Kokkos::MPVectorWorkConfig<execution_space>( row_count,
                                                      EnsembleSize ),
                          *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const unsigned row, const unsigned tid ) const {
    scalar_type s = 0.0;
    for (unsigned pce=0; pce<num_pce; ++pce) {
      s += pce_view(pce,row)*quad_values(qp+tid,pce);
    }
    scalar_view(row).fastAccessCoeff(tid) = s;
  }

/*
  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    const unsigned threads_per_vector = EnsembleSize;
    const unsigned rows_per_block = 6;
    const unsigned num_blocks = (row_count+rows_per_block-1)/rows_per_block;
    const unsigned shared = EnsembleSize*EnsembleSize*sizeof(quad_scalar_type);
    const dim3 block( threads_per_vector, rows_per_block, 1 );
    const dim3 grid( num_blocks, 1 );

    Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( *this );
    Kokkos::parallel_for( Kokkos::MPVectorWorkConfig<execution_space>( row_count,
                                                      EnsembleSize ),
                          *this );
  }

  __device__ inline
  void operator() () const {
    __shared__ volatile quad_scalar_type sh_qv[EnsembleSize][EnsembleSize];

    const unsigned row = blockIdx.x * blockDim.y + threadIdx.y;
    if ( row < row_count ) {
      scalar_type s = 0.0;
      for (unsigned pce_block=0; pce_block<num_pce; pce_block+=EnsembleSize) {
        const unsigned block_size = pce_block+EnsembleSize<num_pce ?
          EnsembleSize : num_pce-pce_block;
        for (unsigned j=threadIdx.y; j<block_size; j+=blockDim.y) {
          sh_qv[j][threadIdx.x] = quad_values(qp+threadIdx.x,pce_block+j);
        }
        __syncthreads();
        for (unsigned j=0; j<block_size; ++j) {
          s += pce_view(pce_block+j,row)*sh_qv[j][threadIdx.x];
        }
      }
      scalar_view(row).fastAccessCoeff(threadIdx.x) = s;
    }
  }
*/
};

template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename quad_weights_type,
           int EnsembleSize >
struct AssemblePCE< pce_view_type,
                    scalar_view_type,
                    quad_values_type,
                    quad_weights_type,
                    EnsembleSize,
                    Kokkos::Cuda > {
  typedef Kokkos::Cuda execution_space;
  typedef typename pce_view_type::array_type  pce_array_type;
  typedef typename pce_array_type::value_type scalar_type;
  typedef typename quad_weights_type::value_type weights_scalar_type;

  const pce_array_type    pce_view;
  const scalar_view_type  scalar_view;
  const quad_values_type  quad_values;
  const quad_weights_type quad_weights;
  unsigned                qp;
  const unsigned          num_pce;

  AssemblePCE( const pce_view_type&     arg_pce_view,
               const scalar_view_type&  arg_scalar_view,
               const quad_values_type&  arg_quad_values,
               const quad_weights_type& arg_quad_weights ) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    quad_weights( arg_quad_weights ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ) {}

  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    Kokkos::parallel_for( Kokkos::MPVectorWorkConfig<execution_space>( pce_view.extent(1),
                                                      EnsembleSize ),
                          *this );
  }

  __device__ inline
  void operator() ( const unsigned row, const unsigned tid ) const {
    const weights_scalar_type w = quad_weights(qp+tid);
    const scalar_type v = w*scalar_view(row).fastAccessCoeff(tid);
    for (unsigned pce=0; pce<num_pce; ++pce) {
      scalar_type s = v*quad_values(qp+tid,pce);
      if (EnsembleSize >= 2) s += Stokhos::shfl_down(s, 1, EnsembleSize);
      if (EnsembleSize >= 4) s += Stokhos::shfl_down(s, 2, EnsembleSize);
      if (EnsembleSize >= 8) s += Stokhos::shfl_down(s, 4, EnsembleSize);
      if (EnsembleSize >= 16) s += Stokhos::shfl_down(s, 8, EnsembleSize);
      if (EnsembleSize >= 32) s += Stokhos::shfl_down(s, 16, EnsembleSize);
      if (tid == 0) pce_view(pce,row) += s;
    }
  }
};

template < typename pce_view_type,
           typename scalar_view_type,
           typename quad_values_type,
           typename quad_weights_type,
           int EnsembleSize >
struct AssembleRightPCE< pce_view_type,
                         scalar_view_type,
                         quad_values_type,
                         quad_weights_type,
                         EnsembleSize,
                         Kokkos::Cuda > {
  typedef Kokkos::Cuda execution_space;
  typedef typename pce_view_type::array_type  pce_array_type;
  typedef typename pce_array_type::value_type scalar_type;
  typedef typename quad_weights_type::value_type weights_scalar_type;
  typedef typename quad_values_type::value_type quad_scalar_type;

  static const unsigned threads_per_vector = EnsembleSize;
  static const unsigned rows_per_block = 256 / EnsembleSize;

  const pce_array_type    pce_view;
  const scalar_view_type  scalar_view;
  const quad_values_type  quad_values;
  const quad_weights_type quad_weights;
  unsigned                qp;
  const unsigned          num_pce;
  const unsigned          row_count;

  AssembleRightPCE( const pce_view_type&     arg_pce_view,
                    const scalar_view_type&  arg_scalar_view,
                    const quad_values_type&  arg_quad_values,
                    const quad_weights_type& arg_quad_weights ) :
    pce_view( arg_pce_view ),
    scalar_view( arg_scalar_view ),
    quad_values( arg_quad_values ),
    quad_weights( arg_quad_weights ),
    qp( 0 ),
    num_pce( Kokkos::dimension_scalar(arg_pce_view) ),
    row_count( pce_view.extent(0) ) {}

  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    Kokkos::parallel_for( Kokkos::MPVectorWorkConfig<execution_space>( row_count,
                                                      EnsembleSize ),
                          *this );
  }

  __device__ inline
  void operator() ( const unsigned row, const unsigned tid ) const {
    const weights_scalar_type w = quad_weights(qp+tid);
    const scalar_type v = w*scalar_view(row).fastAccessCoeff(tid);
    for (unsigned pce=0; pce<num_pce; ++pce) {
      scalar_type s = v*quad_values(qp+tid,pce);
      if (EnsembleSize >= 2) s += Stokhos::shfl_down(s, 1, EnsembleSize);
      if (EnsembleSize >= 4) s += Stokhos::shfl_down(s, 2, EnsembleSize);
      if (EnsembleSize >= 8) s += Stokhos::shfl_down(s, 4, EnsembleSize);
      if (EnsembleSize >= 16) s += Stokhos::shfl_down(s, 8, EnsembleSize);
      if (EnsembleSize >= 32) s += Stokhos::shfl_down(s, 16, EnsembleSize);
      if (tid == 0) pce_view(row,pce) += s;
    }
  }

  /*
  void apply(const unsigned arg_qp) {
    qp = arg_qp;
    const unsigned num_blocks = (row_count+rows_per_block-1)/rows_per_block;
    const unsigned shared = 0;
    const dim3 block( threads_per_vector, rows_per_block, 1 );
    const dim3 grid( num_blocks, 1 );

    Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid, block, shared >>>
        ( *this );
  }

  __device__ inline
  void operator() () const {
    __shared__ volatile quad_scalar_type sh_qv[rows_per_block][EnsembleSize];

    const unsigned row = blockIdx.x * blockDim.y + threadIdx.y;
    if ( row < row_count ) {
      const weights_scalar_type w = quad_weights(qp+threadIdx.x);
      const scalar_type v = w*scalar_view(row).fastAccessCoeff(threadIdx.x);
      for (unsigned pce_block=0; pce_block<num_pce; pce_block+=rows_per_block) {
        const unsigned block_size = pce_block+rows_per_block<num_pce ?
          rows_per_block : num_pce-pce_block;
        if (threadIdx.y < block_size)
          sh_qv[threadIdx.y][threadIdx.x] =
            quad_values(qp+threadIdx.x,pce_block+threadIdx.y);
        __syncthreads();
        for (unsigned j=0; j<block_size; ++j) {
          scalar_type s = v*sh_qv[j][threadIdx.x];
          if (EnsembleSize >= 2) s += Stokhos::shfl_down(s, 1, EnsembleSize);
          if (EnsembleSize >= 4) s += Stokhos::shfl_down(s, 2, EnsembleSize);
          if (EnsembleSize >= 8) s += Stokhos::shfl_down(s, 4, EnsembleSize);
          if (EnsembleSize >= 16) s += Stokhos::shfl_down(s, 8, EnsembleSize);
          if (EnsembleSize >= 32) s += Stokhos::shfl_down(s, 16, EnsembleSize);
          if (threadIdx.x == 0) pce_view(row,pce_block+j) += s;
        }
      }
    }
  }
  */
};
#endif

template< typename ExecutionSpace ,
          BoxElemPart::ElemOrder Order ,
          typename CoordinateMap ,
          typename StorageType ,
          typename OrdinalType ,
          typename MemoryTraits ,
          typename SizeType ,
          class CoeffFunctionType>
class ElementComputation<
  Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap >,
  KokkosSparse::CrsMatrix< Sacado::UQ::PCE<StorageType> , OrdinalType , Kokkos::Device<ExecutionSpace, typename ExecutionSpace::memory_space> , MemoryTraits , SizeType >,
  CoeffFunctionType >
{
public:

  typedef Kokkos::Example::BoxElemFixture< ExecutionSpace, Order, CoordinateMap >  mesh_type ;
  typedef Kokkos::Example::HexElement_Data< mesh_type::ElemNode >              element_data_type ;
  typedef Sacado::UQ::PCE<StorageType> ScalarType;

  //------------------------------------

  typedef ExecutionSpace   execution_space ;
  typedef Kokkos::Device<execution_space, typename execution_space::memory_space> DeviceType;
  typedef ScalarType   scalar_type ;

  typedef KokkosSparse::CrsMatrix< ScalarType , OrdinalType , DeviceType , MemoryTraits , SizeType >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType sparse_graph_type ;
  typedef typename sparse_matrix_type::values_type matrix_values_type ;
  typedef Kokkos::View< scalar_type* , Kokkos::LayoutLeft, execution_space > vector_type ;

  //------------------------------------

  typedef typename scalar_type::value_type scalar_value_type;
  typedef typename scalar_type::ordinal_type ordinal_type;
  static const int EnsembleSize = 32;
  typedef Stokhos::StaticFixedStorage<ordinal_type,scalar_value_type,EnsembleSize,ExecutionSpace> ensemble_storage_type;
  typedef Sacado::MP::Vector<ensemble_storage_type> ensemble_scalar_type;
  typedef typename Sacado::mpl::apply<CoeffFunctionType, ensemble_scalar_type>::type scalar_coeff_function_type;
  typedef ElementComputation<
    Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap >,
    KokkosSparse::CrsMatrix< ensemble_scalar_type , OrdinalType , DeviceType , MemoryTraits , SizeType >,
    scalar_coeff_function_type > scalar_element_computation_type;
  typedef typename scalar_element_computation_type::sparse_matrix_type scalar_sparse_matrix_type ;
  typedef typename scalar_sparse_matrix_type::values_type scalar_matrix_values_type ;
  typedef typename scalar_element_computation_type::vector_type scalar_vector_type ;
  typedef typename scalar_element_computation_type::elem_graph_type elem_graph_type;

  //------------------------------------
  // Computational data:

  const vector_type         solution ;
  const vector_type         residual ;
  const sparse_matrix_type  jacobian ;

  const scalar_vector_type              scalar_solution ;
  const scalar_vector_type              scalar_residual ;
  const scalar_sparse_matrix_type       scalar_jacobian ;
  const scalar_coeff_function_type      scalar_diffusion_coefficient;
  const scalar_element_computation_type scalar_element_computation ;

  typedef QuadratureData<ExecutionSpace> QD;
  typename QD::quad_weights_type quad_weights;
  typename QD::quad_values_type quad_points;
  typename QD::quad_values_type quad_values;

  ElementComputation( const ElementComputation & rhs )
    : solution( rhs.solution )
    , residual( rhs.residual )
    , jacobian( rhs.jacobian )
    , scalar_solution( rhs.scalar_solution )
    , scalar_residual( rhs.scalar_residual )
    , scalar_jacobian( rhs.scalar_jacobian )
    , scalar_diffusion_coefficient( rhs.scalar_diffusion_coefficient )
    , scalar_element_computation( rhs.scalar_element_computation )
    , quad_weights( rhs.quad_weights )
    , quad_points( rhs.quad_points )
    , quad_values( rhs.quad_values )
    {}

  // If the element->sparse_matrix graph is provided then perform atomic updates
  // Otherwise fill per-element contributions for subequent gather-add into a residual and jacobian.
  ElementComputation( const mesh_type          & arg_mesh ,
                      const CoeffFunctionType  & arg_coeff_function ,
                      const bool                 arg_isotropic ,
                      const double             & arg_coeff_source ,
                      const double             & arg_coeff_advection ,
                      const vector_type        & arg_solution ,
                      const elem_graph_type    & arg_elem_graph ,
                      const sparse_matrix_type & arg_jacobian ,
                      const vector_type        & arg_residual ,
                      const KokkosSparse::DeviceConfig arg_dev_config ,
                      const QD& qd )
    : solution( arg_solution )
    , residual( arg_residual )
    , jacobian( arg_jacobian )
    , scalar_solution( "scalar_solution", solution.extent(0) )
    , scalar_residual( "scalar_residual", residual.extent(0) )
    , scalar_jacobian( "scalar_jacobian", jacobian.graph, maximum_entry(jacobian.graph) + 1 )
    , scalar_diffusion_coefficient( arg_coeff_function.m_mean,
                                    arg_coeff_function.m_variance,
                                    arg_coeff_function.m_corr_len,
                                    arg_coeff_function.m_num_rv,
                                    arg_coeff_function.m_use_exp,
                                    arg_coeff_function.m_exp_shift,
                                    arg_coeff_function.m_exp_scale,
                                    arg_coeff_function.m_use_disc_exp_scale )
    , scalar_element_computation( arg_mesh,
                                  scalar_diffusion_coefficient,
                                  arg_isotropic,
                                  arg_coeff_source,
                                  arg_coeff_advection,
                                  scalar_solution,
                                  arg_elem_graph,
                                  scalar_jacobian,
                                  scalar_residual,
                                  arg_dev_config )
    , quad_weights( qd.weights_view )
    , quad_points( qd.points_view )
    , quad_values( qd.values_view )
    {
      // Set global vector size -- this is mandatory
      Kokkos::global_sacado_mp_vector_size = EnsembleSize;
    }

  //------------------------------------

  void apply() const
  {
    typedef EvaluatePCE< vector_type, scalar_vector_type, typename QD::quad_values_type, ensemble_scalar_type, EnsembleSize> evaluate_solution_type;
    typedef AssemblePCE< vector_type, scalar_vector_type, typename QD::quad_values_type, typename QD::quad_weights_type, EnsembleSize> assemble_residual_type;
    typedef AssembleRightPCE< matrix_values_type, scalar_matrix_values_type, typename QD::quad_values_type, typename QD::quad_weights_type, EnsembleSize> assemble_jacobian_type;

    typedef scalar_coeff_function_type KL;
    typedef typename KL::RandomVariableView RV;
    typedef typename RV::HostMirror HRV;
    RV rv = scalar_diffusion_coefficient.getRandomVariables();
    HRV hrv = Kokkos::create_mirror_view(rv);
    auto hqp = Kokkos::create_mirror_view(quad_points);
    Kokkos::deep_copy(hqp, quad_points);

    // Note:  num_quad_points is aligned to the ensemble size to make
    // things easier
    const unsigned num_quad_points = quad_points.extent(0);
    const unsigned dim = quad_points.extent(1);

    evaluate_solution_type evaluate_pce(solution, scalar_solution, quad_values);
    assemble_residual_type assemble_res(residual, scalar_residual,
                                        quad_values, quad_weights);
    assemble_jacobian_type assemble_jac(jacobian.values, scalar_jacobian.values,
                                        quad_values, quad_weights);

    for (unsigned qp=0; qp<num_quad_points; qp+=EnsembleSize) {
      // Zero out residual, and Jacobian
      Kokkos::deep_copy( scalar_residual, 0.0 );
      Kokkos::deep_copy( scalar_jacobian.values, 0.0 );

      // Evaluate PCE solution at quadrature point
      evaluate_pce.apply(qp);

      // Set quadrature point in diffusion coefficient
      for (unsigned i=0; i<dim; ++i)
        for (unsigned j=0; j< unsigned(EnsembleSize); ++j)
          hrv(i).fastAccessCoeff(j) = hqp(qp+j,i);
      Kokkos::deep_copy( rv, hrv );

      // Compute element residual/Jacobian at quadrature point
      scalar_element_computation.apply();

      // Assemble element residual/Jacobian into PCE residual/Jacobian
      assemble_res.apply(qp);
      assemble_jac.apply(qp);
    }
  }

  //------------------------------------
}; /* ElementComputation */

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
