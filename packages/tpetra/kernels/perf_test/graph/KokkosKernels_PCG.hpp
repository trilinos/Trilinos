/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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

#ifndef KOKKOS_EXAMPLE_CG_SOLVE
#define KOKKOS_EXAMPLE_CG_SOLVE

#include <cmath>
#include <limits>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <Kokkos_Atomic.hpp>
#include <Kokkos_MemoryTraits.hpp>

//#include <WrapMPI.hpp>

#include <iostream>
#include "KokkosKernels_GaussSeidel.hpp"
#include "KokkosKernels_Handle.hpp"
//#include "Kokkos_Sparse_MV.hpp"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include <Kokkos_Sparse.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


//#define PRINTRES
namespace KokkosKernels {

namespace Experimental{

namespace Example {

template< typename ValueType ,typename Idx_Type, class Space >
struct CrsMatrix {
  struct StaticCrsGraphType{
    Idx_Type nv, ne;
    const typename Kokkos::View<Idx_Type *, Space > row_map;
    const typename Kokkos::View<Idx_Type *, Space > entries;
    StaticCrsGraphType():row_map(), entries(){}
    StaticCrsGraphType(
        Idx_Type nv_, Idx_Type ne_,
        const Kokkos::View<Idx_Type *, Space > row_map_,
        const Kokkos::View<Idx_Type *, Space > entries_):
          nv(nv_), ne(ne_),
          row_map(row_map_), entries(entries_){}
  } graph;

  typedef typename Kokkos::View< ValueType * , Space > coeff_type ;

  coeff_type coeff ;

  CrsMatrix() : graph(), coeff() {}

  CrsMatrix(
            Idx_Type nv_, Idx_Type ne_,
            const Kokkos::View<Idx_Type *, Space > row_map,
            const Kokkos::View<Idx_Type *, Space > entries,
            coeff_type coeff_)
    : graph( nv_,ne_, row_map,  entries)
    , coeff(coeff_)
    {}
};


template< typename view_y, typename crsMat_t
        , typename view_x>
struct Multiply {

  typedef typename crsMat_t::ordinal_type Idx_Type;
  const crsMat_t    m_A ;
  const view_x m_x ;
  const view_y m_y ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int iRow ) const
    {
      const Idx_Type iEntryBegin = m_A.graph.row_map[iRow];
      const Idx_Type iEntryEnd   = m_A.graph.row_map[iRow+1];

      double sum = 0 ;

      for ( Idx_Type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.values(iEntry) * m_x( m_A.graph.entries(iEntry) );
      }

      m_y(iRow) = sum ;
    }

  Multiply( const view_y & y
          , const crsMat_t    & A
          , const view_x & x
          )
  : m_A( A ), m_x( x ), m_y( y )
  {}
};

template< typename view_y, typename crsMat_t
        , typename view_x
        , class Space >
inline
void multiply( const int nrow
             , const view_y    & y
             , const crsMat_t & A
             , const view_x    & x
             )
{
  Kokkos::parallel_for( Kokkos::RangePolicy<Space>(0,nrow), Multiply<view_y,crsMat_t, view_x>( y , A , x ) );
}

template< typename view_w, typename view_x,typename view_y >
struct WAXPBY {
  const view_w  m_x ;
  const view_x  m_y ;
  const view_y  m_w ;
  const double m_alpha ;
  const double m_beta ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i ) const
    { m_w(i) = m_alpha * m_x(i) + m_beta * m_y(i); }

  WAXPBY( const view_w  & arg_w
        , const double arg_alpha
        , const view_x  & arg_x
        , const double arg_beta
        , const view_y  & arg_y
        )
    : m_x( arg_x )
    , m_y( arg_y )
    , m_w( arg_w )
    , m_alpha( arg_alpha )
    , m_beta( arg_beta )
    {}
};

template< typename view_w, typename view_x, typename view_y , typename Space >
void waxpby( const int n
           , const view_w & arg_w
           , const double arg_alpha
           , const view_x & arg_x
           , const double arg_beta
           , const view_y & arg_y
           )
{
  Kokkos::parallel_for( Kokkos::RangePolicy<Space>(0,n), WAXPBY<view_w,view_x, view_y>(arg_w,arg_alpha,arg_x,arg_beta,arg_y) );
}


template< typename view_x , typename view_y>
struct Dot {
  typedef double value_type ;

  const view_x  m_x ;
  const view_y  m_y ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i , value_type & update ) const
    { update += m_x(i) * m_y(i); }

  Dot( const view_x  & arg_x
     , const view_y & arg_y
     )
    : m_x(arg_x), m_y(arg_y) {}
};

template< typename view_x , typename view_y, typename Space >
double dot( const int n
          , const view_x & arg_x
          , const view_y & arg_y
          )
{
  double result = 0 ;
  Kokkos::parallel_reduce( Kokkos::RangePolicy<Space>(0,n) , Dot<view_x,view_y>( arg_x , arg_y ) , result );
  return result ;
}


template <typename view_type>
struct InitZeroView{
  view_type myview;
  InitZeroView(view_type myview_): myview(myview_){}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i ) const {
    myview(i) = 0;
  }
};

} // namespace Example
} // namespace Kokkos
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosKernels {
namespace Experimental{
namespace Example {


struct CGSolveResult {
  size_t  iteration ;
  double  iter_time ;
  double  matvec_time ;
  double  norm_res ;
  double precond_time;
  double precond_init_time;
};

//#define PRECOND_NORM
//#define PRE
template< typename KernelHandle_t,
          typename crsMatrix_t,
          typename y_vector_t,
          typename x_vector_t
          >
void pcgsolve(
               KernelHandle_t &kh
            ,  const crsMatrix_t &crsMat
            ,  const y_vector_t &y_vector
            ,  x_vector_t x_vector
            ,  const size_t  maximum_iteration = 200
            ,  const double  tolerance = std::numeric_limits<double>::epsilon()
            ,  CGSolveResult * result = 0
            ,  bool use_sgs = true) {

  typedef typename KernelHandle_t::HandleExecSpace Space;
  typedef y_vector_t  VectorType ;

  const size_t count_total = crsMat.numRows();

  size_t  iteration = 0 ;
  double  iter_time = 0 ;
  double  matvec_time = 0 ;
  double  norm_res = 0 ;
  double precond_time = 0;
  double precond_init_time = 0;

  Kokkos::Impl::Timer wall_clock ;
  Kokkos::Impl::Timer timer;

  // Need input vector to matvec to be owned + received
  VectorType pAll ( "cg::p" , count_total );

  VectorType p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_total) );
  VectorType r ( "cg::r" , count_total );
  VectorType Ap( "cg::Ap", count_total );

  /* r = b - A * x ; */

  /* p  = x       */  Kokkos::deep_copy( p , x_vector );

  ///* Ap = A * p   */  multiply<VectorType,crsMatrix_t,VectorType, Space>  ( count_total , Ap , crsMat , pAll );
  KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);

  /* r = b - Ap   */  waxpby<VectorType, y_vector_t, VectorType, Space>( count_total , r , 1.0 , y_vector , -1.0 , Ap );
  /* p  = r       */  Kokkos::deep_copy( p , r );

  //double old_rdot = Kokkos::Example::all_reduce( dot( count_owned , r , r ) , import.comm );
  double old_rdot = dot<VectorType,VectorType, Space>( count_total , r , r );

  norm_res  = sqrt( old_rdot );



  int apply_count = 1;
  VectorType z;
  //double precond_old_rdot = Kokkos::Example::all_reduce( dot( count_owned , r , z ) , import.comm );
  double precond_old_rdot = 1;
#ifdef PRECOND_NORM
  double precond_norm_res  = 1;
#endif
  Kokkos::deep_copy( p , z );

  //typename KernelHandle::GaussSeidelHandleType *gsHandler;
  bool owner_handle = false;
  if (use_sgs){
    if (kh.get_gs_handle() == NULL){

      owner_handle = true;
      kh.create_gs_handle();
    }
    //gsHandler = kh.get_gs_handle();
    timer.reset();

    KokkosKernels::Experimental::Graph::gauss_seidel_numeric
      (&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values);

    Space::fence();
    precond_init_time += timer.seconds();

    z = VectorType( "pcg::z" , count_total );
    Space::fence();
    timer.reset();

    KokkosKernels::Experimental::Graph::symmetric_gauss_seidel_apply
        (&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values, z, r, true, apply_count);

    Space::fence();
    precond_time += timer.seconds();
    //double precond_old_rdot = Kokkos::Example::all_reduce( dot( count_owned , r , z ) , import.comm );
    precond_old_rdot = dot<VectorType,VectorType, Space>( count_total , r , z );
#ifdef PRECOND_NORM
    precond_norm_res  = sqrt( precond_old_rdot );
#endif

    Kokkos::deep_copy( p , z );
  }

  iteration = 0 ;

#ifdef PRINTRES

  std::cout << "norm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;
#ifdef PRECOND_NORM
  if (use_sgs)
  std::cout << "precond_norm_res:" << precond_norm_res << " precond_old_rdot:" << precond_old_rdot<<  std::endl;
#endif

#endif
  while ( tolerance < norm_res && iteration < maximum_iteration ) {

    /* pAp_dot = dot( p , Ap = A * p ) */

    timer.reset();
    ///* import p    */  import( pAll );
    ///* Ap = A * p   */  multiply<VectorType,crsMatrix_t,VectorType, Space>( count_total , Ap , crsMat , pAll );
    KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);


    Space::fence();
    matvec_time += timer.seconds();

    //const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p , Ap ) , import.comm );
    const double pAp_dot = dot<VectorType,VectorType, Space>( count_total , p , Ap ) ;

    double alpha  = 0;
    if (use_sgs){
      alpha = precond_old_rdot / pAp_dot ;
    }
    else {
      alpha = old_rdot / pAp_dot ;
    }

    /* x +=  alpha * p ;  */ waxpby<VectorType, y_vector_t, VectorType, Space>( count_total , x_vector ,  alpha, p  , 1.0 , x_vector );
    /* r += -alpha * Ap ; */ waxpby<VectorType, y_vector_t, VectorType, Space>( count_total , r , -alpha, Ap , 1.0 , r );

    //const double r_dot = Kokkos::Example::all_reduce( dot( count_owned , r , r ) , import.comm );
    const double r_dot = dot<VectorType,VectorType, Space>( count_total , r , r );
    const double beta_original  = r_dot / old_rdot ;

    double precond_r_dot = 1;
    double precond_beta = 1;
    if (use_sgs){
      Space::fence();
      timer.reset();
      KokkosKernels::Experimental::Graph::symmetric_gauss_seidel_apply(
          &kh,
          count_total, count_total,
          crsMat.graph.row_map,
          crsMat.graph.entries,
          crsMat.values, z, r, true,
          apply_count);

      Space::fence();
      precond_time += timer.seconds();
      //const double precond_r_dot = Kokkos::Example::all_reduce( dot( count_owned , r , z ) , import.comm );
      precond_r_dot = dot<VectorType,VectorType, Space>( count_total , r , z );
      precond_beta  = precond_r_dot / precond_old_rdot ;
    }

    double beta  = 1;
    if (!use_sgs){
      beta = beta_original;
      /* p = r + beta * p ; */ waxpby<VectorType, y_vector_t, VectorType, Space>( count_total , p , 1.0 , r , beta , p );
    }
    else {
      beta = precond_beta;
      waxpby<VectorType, y_vector_t, VectorType, Space>( count_total , p , 1.0 , z , beta , p );
    }

#ifdef PRINTRES
    std::cout << "\tbeta_original:" << beta_original <<  std::endl;

    if (use_sgs)
    std::cout << "\tprecond_beta:" << precond_beta <<  std::endl;

#endif


    norm_res = sqrt( old_rdot = r_dot );
#ifdef PRECOND_NORM
    if (use_sgs){
      precond_norm_res = sqrt( precond_old_rdot = precond_r_dot );
    }
#else
    precond_old_rdot = precond_r_dot;
#endif

#ifdef PRINTRES
    std::cout << "\tnorm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;
#ifdef PRECOND_NORM

    if (use_sgs)
    std::cout << "\tprecond_norm_res:" << precond_norm_res << " precond_old_rdot:" << precond_old_rdot<<  std::endl;
#endif
#endif
    ++iteration ;
  }

  Space::fence();
  iter_time = wall_clock.seconds();

  if ( 0 != result ) {
    result->iteration   = iteration ;
    result->iter_time   = iter_time ;
    result->matvec_time = matvec_time ;
    result->norm_res    = norm_res ;
    result->precond_time = precond_time;
    result->precond_init_time = precond_init_time;
  }

  if (use_sgs & owner_handle ){

    kh.destroy_gs_handle();
  }
}

} // namespace Example
} // namespace Kokkos
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_CG_SOLVE */


