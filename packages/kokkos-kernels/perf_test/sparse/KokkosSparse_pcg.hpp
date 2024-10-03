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

#ifndef KOKKOS_EXAMPLE_CG_SOLVE
#define KOKKOS_EXAMPLE_CG_SOLVE

#include <cmath>
#include <limits>
#include <Kokkos_Core.hpp>

#include <iostream>
#include "KokkosKernels_Handle.hpp"
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosSparse_sor_sequential_impl.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// #define KK_TICTOCPRINT
namespace KokkosKernels {
namespace Experimental {
namespace Example {

struct CGSolveResult {
  size_t iteration;
  double iter_time;
  double matvec_time;
  double norm_res;
  double precond_time;
  double precond_init_time;
};

template <typename KernelHandle_t, typename crsMatrix_t, typename y_vector_t, typename x_vector_t>
void block_pcgsolve(KernelHandle_t &kh, const crsMatrix_t &point_crsMat, const crsMatrix_t &_block_crsMat,
                    int block_size, const y_vector_t &y_vector, x_vector_t x_vector,
                    const size_t maximum_iteration = 200,
                    const double tolerance = std::numeric_limits<double>::epsilon(), CGSolveResult *result = 0,
                    bool use_sgs = true) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  typedef typename KernelHandle_t::HandleExecSpace Space;

  const size_t count_total = point_crsMat.numRows();

  size_t iteration         = 0;
  double iter_time         = 0;
  double matvec_time       = 0;
  double norm_res          = 0;
  double precond_time      = 0;
  double precond_init_time = 0;

  Kokkos::Timer wall_clock;
  Kokkos::Timer timer;

  // Need input vector to matvec to be owned + received
  y_vector_t pAll("cg::p", count_total);

  y_vector_t p = Kokkos::subview(pAll, std::pair<size_t, size_t>(0, count_total));
  y_vector_t r("cg::r", count_total);
  y_vector_t Ap("cg::Ap", count_total);

  // r = b - A * x ;
  // p  = x
  Kokkos::deep_copy(p, x_vector);

  // Ap = A * p
  KokkosSparse::spmv("N", 1, point_crsMat, pAll, 0, Ap);

  // r  = Ap
  Kokkos::deep_copy(r, Ap);

  // r = b - r
  KokkosBlas::axpby(1.0, y_vector, -1.0, r);

  // p  = r
  Kokkos::deep_copy(p, r);
  ;
  double old_rdot = KokkosBlas::dot(r, r);
  norm_res        = sqrt(old_rdot);

  int apply_count = 1;
  y_vector_t z;

  double precond_old_rdot = 1;
  // Kokkos::deep_copy( p , z );

  bool owner_handle = false;

  KernelHandle_t block_kh;
  block_kh.create_gs_handle();
  block_kh.get_point_gs_handle()->set_block_size(block_size);
  // block_kh.set_shmem_size(8032);
  if (use_sgs) {
    if (kh.get_gs_handle() == NULL) {
      owner_handle = true;
      kh.create_gs_handle();
    }

    timer.reset();

    // gauss_seidel_numeric
    //  (&kh, count_total, count_total, point_crsMat.graph.row_map,
    //  point_crsMat.graph.entries, point_crsMat.values);

    // Space().fence();
    // timer.reset();

    // block_kh.set_verbose(true);
    block_gauss_seidel_numeric(&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(), block_size,
                               _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values);

    precond_init_time += timer.seconds();

    z = y_vector_t("pcg::z", count_total);
    Space().fence();
    timer.reset();
    symmetric_block_gauss_seidel_apply(&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(), block_size,
                                       _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values,
                                       z, r, true, true, 1.0, apply_count);

    // symmetric_gauss_seidel_apply
    //    (&kh, count_total, count_total, point_crsMat.graph.row_map,
    //    point_crsMat.graph.entries, point_crsMat.values, z, r, true, true,
    //    apply_count);
    Space().fence();
    precond_time += timer.seconds();
    precond_old_rdot = KokkosBlas::dot(r, z);
    Kokkos::deep_copy(p, z);
  }

  iteration = 0;

#ifdef KK_TICTOCPRINT

  std::cout << "norm_res:" << norm_res << " old_rdot:" << old_rdot << std::endl;

#endif
  while (tolerance < norm_res && iteration < maximum_iteration) {
    timer.reset();
    // Ap = A * p
    KokkosSparse::spmv("N", 1, point_crsMat, pAll, 0, Ap);

    Space().fence();
    matvec_time += timer.seconds();

    // const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p
    // , Ap ) , import.comm ); const double pAp_dot = dot<y_vector_t,y_vector_t,
    // Space>( count_total , p , Ap ) ;

    // pAp_dot = dot(Ap , p);
    const double pAp_dot = KokkosBlas::dot(p, Ap);

    double alpha = 0;
    if (use_sgs) {
      alpha = precond_old_rdot / pAp_dot;
    } else {
      alpha = old_rdot / pAp_dot;
    }

    // x +=  alpha * p ;
    KokkosBlas::axpby(alpha, p, 1.0, x_vector);

    // r += -alpha * Ap ;
    KokkosBlas::axpby(-alpha, Ap, 1.0, r);

    const double r_dot = KokkosBlas::dot(r, r);

    const double beta_original = r_dot / old_rdot;
    double precond_r_dot       = 1;
    double precond_beta        = 1;
    if (use_sgs) {
      Space().fence();
      timer.reset();
      symmetric_block_gauss_seidel_apply(&block_kh, _block_crsMat.numRows(), _block_crsMat.numCols(), block_size,
                                         _block_crsMat.graph.row_map, _block_crsMat.graph.entries, _block_crsMat.values,
                                         z, r, true, true, 1.0, apply_count);

      // symmetric_gauss_seidel_apply(
      //    &kh,
      //    count_total, count_total,
      //    point_crsMat.graph.row_map,
      //    point_crsMat.graph.entries,
      //    point_crsMat.values, z, r, true,
      //    apply_count);

      Space().fence();
      precond_time += timer.seconds();
      precond_r_dot = KokkosBlas::dot(r, z);
      precond_beta  = precond_r_dot / precond_old_rdot;
    }

    double beta = 1;
    if (!use_sgs) {
      beta = beta_original;
      // p = r + beta * p ;
      KokkosBlas::axpby(1.0, r, beta, p);
    } else {
      beta = precond_beta;
      KokkosBlas::axpby(1.0, z, beta, p);
    }

#ifdef KK_TICTOCPRINT
    std::cout << "\tbeta_original:" << beta_original << std::endl;
    if (use_sgs) std::cout << "\tprecond_beta:" << precond_beta << std::endl;

#endif

    norm_res         = sqrt(old_rdot = r_dot);
    precond_old_rdot = precond_r_dot;

#ifdef KK_TICTOCPRINT
    std::cout << "\tnorm_res:" << norm_res << " old_rdot:" << old_rdot << std::endl;
#endif
    ++iteration;
  }

  Space().fence();
  iter_time = wall_clock.seconds();

  if (0 != result) {
    result->iteration         = iteration;
    result->iter_time         = iter_time;
    result->matvec_time       = matvec_time;
    result->norm_res          = norm_res;
    result->precond_time      = precond_time;
    result->precond_init_time = precond_init_time;
  }

  if (use_sgs & owner_handle) {
    kh.destroy_gs_handle();
  }
}

template <typename KernelHandle_t, typename crsMatrix_t, typename y_vector_t, typename x_vector_t>
void pcgsolve(KernelHandle_t &kh, const crsMatrix_t &crsMat, const y_vector_t &y_vector, x_vector_t x_vector,
              const size_t maximum_iteration = 200, const double tolerance = std::numeric_limits<double>::epsilon(),
              CGSolveResult *result = 0, bool use_sgs = true, int /*clusterSize*/ = 1,
              bool use_sequential_sgs = false) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  using size_type = typename KernelHandle_t::size_type;
  using nnz_lno_t = typename KernelHandle_t::nnz_lno_t;
  using Space     = typename KernelHandle_t::HandleExecSpace;
  static_assert(std::is_same<double, typename KernelHandle_t::nnz_scalar_t>::value,
                "The PCG performance test only works with scalar = double.");

  const nnz_lno_t count_total = crsMat.numRows();

  size_t iteration         = 0;
  double iter_time         = 0;
  double matvec_time       = 0;
  double norm_res          = 0;
  double precond_time      = 0;
  double precond_init_time = 0;

  Kokkos::Timer wall_clock;
  Kokkos::Timer timer;

  // Need input vector to matvec to be owned + received
  y_vector_t pAll("cg::p", count_total);

  y_vector_t p = Kokkos::subview(pAll, std::pair<size_t, size_t>(0, count_total));
  y_vector_t r("cg::r", count_total);
  y_vector_t Ap("cg::Ap", count_total);

  /* r = b - A * x ; */
  /* p  = x       */ Kokkos::deep_copy(p, x_vector);

  /* Ap = A * p   */ KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);

  /* r  = Ap       */ Kokkos::deep_copy(r, Ap);

  /* r = b - r   */ KokkosBlas::axpby(1.0, y_vector, -1.0, r);

  /* p  = r       */ Kokkos::deep_copy(p, r);

  double old_rdot = KokkosBlas::dot(r, r);
  norm_res        = sqrt(old_rdot);

  int apply_count = 1;
  y_vector_t z;

  double precond_old_rdot = 1;
  // Kokkos::deep_copy( p , z );

  bool use_par_sgs = use_sgs && !use_sequential_sgs;

  auto ptrHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crsMat.graph.row_map);
  auto indHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crsMat.graph.entries);
  auto valHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crsMat.values);
  Kokkos::View<double *, Kokkos::HostSpace> diagHost;
  if (use_sequential_sgs) {
    diagHost = Kokkos::View<double *, Kokkos::HostSpace>("Diag for Seq SOR", count_total);
    for (int i = 0; i < count_total; i++) {
      for (size_type j = ptrHost(i); j < ptrHost(i + 1); j++) {
        if (indHost(j) == i) diagHost(i) = 1.0 / valHost(j);
      }
    }
  }
  auto xHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x_vector);
  auto yHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_vector);

  if (use_sgs) {
    timer.reset();
    z = y_vector_t("pcg::z", count_total);
    if (use_par_sgs) {
      gauss_seidel_numeric(&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values);

      Space().fence();

      precond_init_time += timer.seconds();
      Space().fence();
      timer.reset();

      symmetric_gauss_seidel_apply(&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries,
                                   crsMat.values, z, r, true, true, 1.0, apply_count);

      Space().fence();
    } else if (use_sequential_sgs) {
      // z = LHS (aka x), r RHS (aka y or b)
      Kokkos::deep_copy(z, 0.0);
      auto zhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), z);
      auto rhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
      // as with par_sgs, init unknown to 0
      timer.reset();
      for (int sweep = 0; sweep < apply_count; sweep++) {
        KokkosSparse::Impl::Sequential::gaussSeidel<nnz_lno_t, size_type, double, double, double>(
            count_total,  // rows = cols of the matrix
            1,            // number of vectors in X and B
            ptrHost.data(), indHost.data(), valHost.data(), rhost.data(),
            count_total,  // raw ptr to B vector, and B column stride (for when
                          // multiple RHS gets added to MTSGS)
            zhost.data(),
            count_total,  // raw ptr to X vector, and X column stride
            diagHost.data(), 1.0, "F");
        KokkosSparse::Impl::Sequential::gaussSeidel<nnz_lno_t, size_type, double, double, double>(
            count_total, 1, ptrHost.data(), indHost.data(), valHost.data(), rhost.data(), count_total, zhost.data(),
            count_total, diagHost.data(), 1.0, "B");
      }
      // result is in z (but r doesn't change)
      Kokkos::deep_copy(z, zhost);
      Kokkos::deep_copy(r, rhost);
    }
    precond_time += timer.seconds();
    precond_old_rdot = KokkosBlas::dot(r, z);
    Kokkos::deep_copy(p, z);
  }

  iteration = 0;

#ifdef KK_TICTOCPRINT

  std::cout << "norm_res:" << norm_res << " old_rdot:" << old_rdot << std::endl;

#endif
  while (tolerance < norm_res && iteration < maximum_iteration) {
    std::cout << "Running CG iteration " << iteration << ", current resnorm = " << norm_res << '\n';

    timer.reset();
    /* Ap = A * p   */ KokkosSparse::spmv("N", 1, crsMat, pAll, 0, Ap);

    Space().fence();
    matvec_time += timer.seconds();

    // const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p
    // , Ap ) , import.comm ); const double pAp_dot = dot<y_vector_t,y_vector_t,
    // Space>( count_total , p , Ap ) ;

    /* pAp_dot = dot(Ap , p ) */ const double pAp_dot = KokkosBlas::dot(p, Ap);

    double alpha = 0;
    if (use_sgs) {
      alpha = precond_old_rdot / pAp_dot;
    } else {
      alpha = old_rdot / pAp_dot;
    }

    /* x +=  alpha * p ;  */ KokkosBlas::axpby(alpha, p, 1.0, x_vector);

    /* r += -alpha * Ap ; */ KokkosBlas::axpby(-alpha, Ap, 1.0, r);

    const double r_dot = KokkosBlas::dot(r, r);

    const double beta_original = r_dot / old_rdot;
    double precond_r_dot       = 1;
    double precond_beta        = 1;
    if (use_sgs) {
      Space().fence();
      timer.reset();
      if (use_par_sgs) {
        symmetric_gauss_seidel_apply(&kh, count_total, count_total, crsMat.graph.row_map, crsMat.graph.entries,
                                     crsMat.values, z, r, true, true, 1.0, apply_count);
      } else if (use_sequential_sgs) {
        // z = LHS (aka x), r RHS (aka y or b)
        Kokkos::deep_copy(z, 0.0);
        auto zhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), z);
        auto rhost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), r);
        // as with the par_sgs version, init unknown (here, zhost) to 0
        for (int sweep = 0; sweep < apply_count; sweep++) {
          KokkosSparse::Impl::Sequential::gaussSeidel<nnz_lno_t, size_type, double, double, double>(
              count_total, 1, ptrHost.data(), indHost.data(), valHost.data(), rhost.data(), count_total, zhost.data(),
              count_total, diagHost.data(), 1.0, "F");
          KokkosSparse::Impl::Sequential::gaussSeidel<nnz_lno_t, size_type, double, double, double>(
              count_total, 1, ptrHost.data(), indHost.data(), valHost.data(), rhost.data(), count_total, zhost.data(),
              count_total, diagHost.data(), 1.0, "B");
        }
        Kokkos::deep_copy(z, zhost);
        Kokkos::deep_copy(r, rhost);
      }
      precond_time += timer.seconds();
      precond_r_dot = KokkosBlas::dot(r, z);
      precond_beta  = precond_r_dot / precond_old_rdot;
    }
    double beta = 1;
    if (!use_sgs) {
      beta = beta_original;
      /* p = r + beta * p ; */ KokkosBlas::axpby(1.0, r, beta, p);
    } else {
      beta = precond_beta;
      KokkosBlas::axpby(1.0, z, beta, p);
    }

#ifdef KK_TICTOCPRINT
    std::cout << "\tbeta_original:" << beta_original << std::endl;
    if (use_sgs) std::cout << "\tprecond_beta:" << precond_beta << std::endl;

#endif

    norm_res         = sqrt(old_rdot = r_dot);
    precond_old_rdot = precond_r_dot;

#ifdef KK_TICTOCPRINT
    std::cout << "\tnorm_res:" << norm_res << " old_rdot:" << old_rdot << std::endl;
#endif
    ++iteration;
  }

  Space().fence();
  iter_time = wall_clock.seconds();

  if (0 != result) {
    result->iteration         = iteration;
    result->iter_time         = iter_time;
    result->matvec_time       = matvec_time;
    result->norm_res          = norm_res;
    result->precond_time      = precond_time;
    result->precond_init_time = precond_init_time;
  }
}

}  // namespace Example
}  // namespace Experimental
}  // namespace KokkosKernels
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_CG_SOLVE */
