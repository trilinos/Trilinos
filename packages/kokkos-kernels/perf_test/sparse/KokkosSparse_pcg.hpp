/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
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

#include <iostream>
#include "KokkosKernels_Handle.hpp"
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


//#define KK_TICTOCPRINT
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

template< typename KernelHandle_t,
          typename crsMatrix_t,
          typename y_vector_t,
          typename x_vector_t
          >
void block_pcgsolve(
               KernelHandle_t &kh
            ,  const crsMatrix_t &point_crsMat
			,  const crsMatrix_t &_block_crsMat, int block_size
            ,  const y_vector_t &y_vector
            ,  x_vector_t x_vector
            ,  const size_t  maximum_iteration = 200
            ,  const double  tolerance = std::numeric_limits<double>::epsilon()
            ,  CGSolveResult * result = 0
            ,  bool use_sgs = true) {

  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  typedef typename KernelHandle_t::HandleExecSpace Space;


  const size_t count_total = point_crsMat.numRows();

  size_t  iteration = 0 ;
  double  iter_time = 0 ;
  double  matvec_time = 0 ;
  double  norm_res = 0 ;
  double precond_time = 0;
  double precond_init_time = 0;

  Kokkos::Impl::Timer wall_clock ;
  Kokkos::Impl::Timer timer;

  // Need input vector to matvec to be owned + received
  y_vector_t pAll ( "cg::p" , count_total );

  y_vector_t p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_total) );
  y_vector_t r ( "cg::r" , count_total );
  y_vector_t Ap( "cg::Ap", count_total );

  /* r = b - A * x ; */
  /* p  = x       */  Kokkos::deep_copy( p , x_vector );

  /* Ap = A * p   */  KokkosSparse::spmv("N", 1, point_crsMat, pAll, 0, Ap);

  /* r  = Ap       */  Kokkos::deep_copy( r , Ap );

  /* r = b - r   */  KokkosBlas::axpby(1.0, y_vector, -1.0, r);

  /* p  = r       */  Kokkos::deep_copy( p , r );
;
  double old_rdot = KokkosBlas::dot( r , r );
  norm_res  = sqrt( old_rdot );

  int apply_count = 1;
  y_vector_t z;

  double precond_old_rdot = 1;
  //Kokkos::deep_copy( p , z );

  bool owner_handle = false;

  KernelHandle_t block_kh; block_kh.create_gs_handle(); block_kh.get_gs_handle()->set_block_size(block_size);
    //block_kh.set_shmem_size(8032);
  if (use_sgs){
    if (kh.get_gs_handle() == NULL){
      owner_handle = true;
      kh.create_gs_handle();
    }

    timer.reset();
/*
    gauss_seidel_numeric
      (&kh, count_total, count_total, point_crsMat.graph.row_map, point_crsMat.graph.entries, point_crsMat.values);

    Space::fence();
    timer.reset();
*/

    //block_kh.set_verbose(true);
    block_gauss_seidel_numeric
          (&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(), block_size, _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values);

    precond_init_time += timer.seconds();

    z = y_vector_t( "pcg::z" , count_total );
    Space::fence();
    timer.reset();
    symmetric_block_gauss_seidel_apply
            (&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(),block_size,  _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values,
            		z, r, true, true, apply_count);
/*
    symmetric_gauss_seidel_apply
        (&kh, count_total, count_total, point_crsMat.graph.row_map, point_crsMat.graph.entries, point_crsMat.values, z, r, true, true, apply_count);
*/
    Space::fence();
    precond_time += timer.seconds();
    precond_old_rdot = KokkosBlas::dot( r , z );
    Kokkos::deep_copy( p , z );
  }

  iteration = 0 ;

#ifdef KK_TICTOCPRINT

  std::cout << "norm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;

#endif
  while ( tolerance < norm_res && iteration < maximum_iteration ) {


    timer.reset();
    /* Ap = A * p   */  KokkosSparse::spmv("N", 1, point_crsMat, pAll, 0, Ap);


    Space::fence();
    matvec_time += timer.seconds();

    //const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p , Ap ) , import.comm );
    //const double pAp_dot = dot<y_vector_t,y_vector_t, Space>( count_total , p , Ap ) ;

    /* pAp_dot = dot(Ap , p ) */ const double pAp_dot = KokkosBlas::dot( p , Ap ) ;


    double alpha  = 0;
    if (use_sgs){
      alpha = precond_old_rdot / pAp_dot ;
    }
    else {
      alpha = old_rdot / pAp_dot ;
    }

    /* x +=  alpha * p ;  */  KokkosBlas::axpby(alpha, p, 1.0, x_vector);

    /* r += -alpha * Ap ; */  KokkosBlas::axpby(-alpha, Ap, 1.0, r);

    const double r_dot = KokkosBlas::dot( r , r );

    const double beta_original  = r_dot / old_rdot ;
    double precond_r_dot = 1;
    double precond_beta = 1;
    if (use_sgs){
      Space::fence();
      timer.reset();
      symmetric_block_gauss_seidel_apply
                  (&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(),block_size, _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values,
                  		z, r, true, true, apply_count);
      /*

      symmetric_gauss_seidel_apply(
          &kh,
          count_total, count_total,
          point_crsMat.graph.row_map,
          point_crsMat.graph.entries,
          point_crsMat.values, z, r, true,
          apply_count);
          */

      Space::fence();
      precond_time += timer.seconds();
      precond_r_dot = KokkosBlas::dot(r , z );
      precond_beta  = precond_r_dot / precond_old_rdot ;
    }

    double beta  = 1;
    if (!use_sgs){
      beta = beta_original;
      /* p = r + beta * p ; */  KokkosBlas::axpby(1.0, r, beta, p);
    }
    else {
      beta = precond_beta;
      KokkosBlas::axpby(1.0, z, beta, p);
    }

#ifdef KK_TICTOCPRINT
    std::cout << "\tbeta_original:" << beta_original <<  std::endl;
    if (use_sgs)
    std::cout << "\tprecond_beta:" << precond_beta <<  std::endl;

#endif


    norm_res = sqrt( old_rdot = r_dot );
    precond_old_rdot = precond_r_dot;

#ifdef KK_TICTOCPRINT
    std::cout << "\tnorm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;
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

  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  typedef typename KernelHandle_t::HandleExecSpace Space;


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
  y_vector_t pAll ( "cg::p" , count_total );

  y_vector_t p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_total) );
  y_vector_t r ( "cg::r" , count_total );
  y_vector_t Ap( "cg::Ap", count_total );

  /* r = b - A * x ; */
  /* p  = x       */  Kokkos::deep_copy( p , x_vector );

  /* Ap = A * p   */  KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);

  /* r  = Ap       */  Kokkos::deep_copy( r , Ap );

  /* r = b - r   */  KokkosBlas::axpby(1.0, y_vector, -1.0, r);

  /* p  = r       */  Kokkos::deep_copy( p , r );
;
  double old_rdot = KokkosBlas::dot( r , r );
  norm_res  = sqrt( old_rdot );

  int apply_count = 1;
  y_vector_t z;

  double precond_old_rdot = 1;
  //Kokkos::deep_copy( p , z );

  bool owner_handle = false;
  if (use_sgs){
    if (kh.get_gs_handle() == NULL){
      owner_handle = true;
      kh.create_gs_handle();
    }

    timer.reset();
    //kh.set_verbose(true);

    gauss_seidel_numeric
      (&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values);

    Space::fence();

    precond_init_time += timer.seconds();
    z = y_vector_t( "pcg::z" , count_total );
    Space::fence();
    timer.reset();

    symmetric_gauss_seidel_apply
        (&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values, z, r, true, true, apply_count);

    Space::fence();
    precond_time += timer.seconds();
    precond_old_rdot = KokkosBlas::dot( r , z );
    Kokkos::deep_copy( p , z );
  }

  iteration = 0 ;

#ifdef KK_TICTOCPRINT

  std::cout << "norm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;

#endif
  while (tolerance < norm_res && iteration < maximum_iteration ) {


    timer.reset();
    /* Ap = A * p   */  KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);


    Space::fence();
    matvec_time += timer.seconds();

    //const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p , Ap ) , import.comm );
    //const double pAp_dot = dot<y_vector_t,y_vector_t, Space>( count_total , p , Ap ) ;

    /* pAp_dot = dot(Ap , p ) */ const double pAp_dot = KokkosBlas::dot( p , Ap ) ;


    double alpha  = 0;
    if (use_sgs){
      alpha = precond_old_rdot / pAp_dot ;
    }
    else {
      alpha = old_rdot / pAp_dot ;
    }

    /* x +=  alpha * p ;  */  KokkosBlas::axpby(alpha, p, 1.0, x_vector);

    /* r += -alpha * Ap ; */  KokkosBlas::axpby(-alpha, Ap, 1.0, r);

    const double r_dot = KokkosBlas::dot( r , r );

    const double beta_original  = r_dot / old_rdot ;
    double precond_r_dot = 1;
    double precond_beta = 1;
    if (use_sgs){
      Space::fence();
      timer.reset();
      symmetric_gauss_seidel_apply(
          &kh,
          count_total, count_total,
          crsMat.graph.row_map,
          crsMat.graph.entries,
          crsMat.values, z, r, true,
          apply_count);

      Space::fence();
      precond_time += timer.seconds();
      precond_r_dot = KokkosBlas::dot(r , z );
      precond_beta  = precond_r_dot / precond_old_rdot ;
    }

    double beta  = 1;
    if (!use_sgs){
      beta = beta_original;
      /* p = r + beta * p ; */  KokkosBlas::axpby(1.0, r, beta, p);
    }
    else {
      beta = precond_beta;
      KokkosBlas::axpby(1.0, z, beta, p);
    }

#ifdef KK_TICTOCPRINT
    std::cout << "\tbeta_original:" << beta_original <<  std::endl;
    if (use_sgs)
    std::cout << "\tprecond_beta:" << precond_beta <<  std::endl;

#endif


    norm_res = sqrt( old_rdot = r_dot );
    precond_old_rdot = precond_r_dot;

#ifdef KK_TICTOCPRINT
    std::cout << "\tnorm_res:" << norm_res << " old_rdot:" << old_rdot<<  std::endl;
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


