// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

//#define SACADO_DISABLE_FAD_VIEW_SPEC

#include "Sacado.hpp"

#include "mat_vec.hpp"

#include "Kokkos_Timer.hpp"

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c) {
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;

  const int m = A.extent(0);
  const int n = A.extent(1);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>( 0,m ),
    KOKKOS_LAMBDA (const int i) {
      scalar_type t = 0.0;
      for (int j=0; j<n; ++j)
        t += A(i,j)*b(j);
      c(i) = t;
    }
  );
}

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_scratch(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> Policy;
  typedef typename Policy::member_type team_member;
  typedef Kokkos::View<scalar_type*,Kokkos::LayoutLeft, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged> TmpScratchSpace;

  const int m = A.extent(0);
  const int n = A.extent(1);
  const int p = dimension_scalar(A);

#ifdef KOKKOS_ENABLE_CUDA
  const bool is_cuda = std::is_same<execution_space,Kokkos::Cuda>::value;
#else
  const bool is_cuda = false;
#endif
  const int TeamSize = is_cuda ? 128 : 1;
  const int N = (m+TeamSize-1)/TeamSize;
  Policy policy(N, TeamSize, 1);
  const size_t bytes = TmpScratchSpace::shmem_size(TeamSize,p);
  Kokkos::parallel_for(
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes)),
    KOKKOS_LAMBDA (const team_member& team) {
      const int team_rank = team.team_rank();
      const int team_size = team.team_size();
      TmpScratchSpace t(team.team_scratch(0), team_size, p);
      const int i = team.league_rank()*team_size + team_rank;
      if (i < m) {
        t(team_rank) = 0.0;
        for (int j=0; j<n; ++j)
          t(team_rank) += A(i,j)*b(j);
        c(i) = t(team_rank);
      }
    }
  );
}

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  typedef typename ViewTypeC::execution_space execution_space;

  const int m = A.extent(0);
  const int n = A.extent(1);
  const int p = A.extent(2)-1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>( 0,m ),
    KOKKOS_LAMBDA (const int i) {
      c(i,p) = 0.0;
      for (int k=0; k<p; ++k)
        c(i,k) = 0.0;
      for (int j=0; j<n; ++j) {
        c(i,p) += A(i,j,p)*b(j,p);
        for (int k=0; k<p; ++k) {
          c(i,k) += A(i,j,k)*b(j,p) + A(i,j,p)*b(j,k);
        }
      }
    }
  );
}

template <int MaxP, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv_sl(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;

  const int m = A.extent(0);
  const int n = A.extent(1);
  const int p = A.extent(2)-1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>( 0,m ),
    KOKKOS_LAMBDA (const int i) {
      scalar_type cv = 0.0;
      scalar_type t[MaxP];
      for (int k=0; k<p; ++k)
        t[k] = 0.0;

      for (int j=0; j<n; ++j) {
        scalar_type av = A(i,j,p);
        scalar_type bv = b(j,p);
        cv += av*bv;
        for (int k=0; k<p; ++k) {
          t[k] += A(i,j,k)*bv + av*b(j,k);
        }
      }

      for (int k=0; k<p; ++k)
        c(i,k) = t[k];
      c(i,p) = cv;
    }
  );
}

template <int p, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv_s(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;

  const int m = A.extent(0);
  const int n = A.extent(1);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>( 0,m ),
    KOKKOS_LAMBDA (const int i) {
      scalar_type cv = 0.0;
      scalar_type t[p];
      for (int k=0; k<p; ++k)
        t[k] = 0.0;

      for (int j=0; j<n; ++j) {
        const scalar_type av = A(i,j,p);
        const scalar_type bv = b(j,p);
        cv += av*bv;

// Using simd here results in much better performance.  Othewise the compiler
// appears to try and vectorize the j loop with gather instructions, which
// doesn't work very well.
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma simd
#endif
        for (int k=0; k<p; ++k) {
          t[k] += A(i,j,k)*bv + av*b(j,k);
        }
      }

      for (int k=0; k<p; ++k)
        c(i,k) = t[k];
      c(i,p) = cv;
    }
  );
}

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
check_val(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  const double tol = 1.0e-14;
  typedef typename ViewTypeC::value_type value_type;
  typename ViewTypeC::HostMirror h_c = Kokkos::create_mirror_view(c);
  Kokkos::deep_copy(h_c, c);
  const size_t m = A.extent(0);
  const size_t n = A.extent(1);
  for (size_t i=0; i<m; ++i) {
    value_type t = n;
    if (std::abs(h_c(i)- t) > tol) {
      std::cout << "Comparison failed!  " << i << " : " << h_c(i) << " , " << t
                << std::endl;
    }
  }
}

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
check_deriv(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  const double tol = 1.0e-14;
  typedef typename ViewTypeC::value_type value_type;
  typename ViewTypeC::HostMirror h_c = Kokkos::create_mirror_view(c);
  Kokkos::deep_copy(h_c, c);
  const size_t m = A.extent(0);
  const size_t n = A.extent(1);
  const size_t p = A.extent(2);
  for (size_t i=0; i<m; ++i) {
    for (size_t j=0; j<p; ++j) {
      value_type t = (j == p-1 ? n : 2*n);
      if (std::abs(h_c(i,j)- t) > tol) {
        std::cout << "Comparison failed!  " << i << "," << j << " : "
                  << h_c(i,j) << " , " << t << std::endl;
      }
    }
  }
}

template <typename ... ViewArgs>
Perf
do_time_val(const size_t m, const size_t n, const size_t nloop,
            const bool check)
{
  typedef Kokkos::View<double**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  ViewTypeA A("A",m,n);
  ViewTypeB b("B",n);
  ViewTypeC c("c",m);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*2;
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_val(A,b,c);

  return perf;
}

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad(const size_t m, const size_t n, const size_t p, const size_t nloop,
            const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  // Set amount of memory available for dynamic memory allocation on GPU
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<execution_space,Kokkos::Cuda>::value &&
      std::is_same<FadType,Sacado::Fad::DFad<double> >::value) {
    const size_t concurrency = execution_space().concurrency();
    const size_t mem = std::min(m,concurrency) * p * sizeof(double);
    //std::cout << "mem = " << mem / (1024*1024) << " MB" << std::endl;
    cudaDeviceSetLimit(cudaLimitMallocHeapSize, mem);
  }
#endif

#ifndef SACADO_DISABLE_FAD_VIEW_SPEC
  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);
#else
  ViewTypeA A("A",m,n);
  ViewTypeB b("B",n);
  ViewTypeC c("c",m);
#endif

  // FadType a(p, 1.0);
  // for (size_t k=0; k<p; ++k)
  //   a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(typename ViewTypeA::array_type(A), 1.0);
  Kokkos::deep_copy(typename ViewTypeB::array_type(b), 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

#ifndef SACADO_DISABLE_FAD_VIEW_SPEC
  if (check) {
    typename ViewTypeA::array_type A_flat = A;
    typename ViewTypeB::array_type b_flat = b;
    typename ViewTypeC::array_type c_flat = c;
    check_deriv(A_flat, b_flat, c_flat);
  }
#endif

  return perf;
}

template <typename FadType, typename ... ViewArgs>
Perf
do_time_scratch(const size_t m, const size_t n, const size_t p, const size_t nloop,
                const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

#ifndef SACADO_DISABLE_FAD_VIEW_SPEC
  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);
#else
  ViewTypeA A("A",m,n);
  ViewTypeB b("B",n);
  ViewTypeC c("c",m);
#endif

  // FadType a(p, 1.0);
  // for (size_t k=0; k<p; ++k)
  //   a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(typename ViewTypeA::array_type(A), 1.0);
  Kokkos::deep_copy(typename ViewTypeB::array_type(b), 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_scratch( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_scratch( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

#ifndef SACADO_DISABLE_FAD_VIEW_SPEC
  if (check) {
    typename ViewTypeA::array_type A_flat = A;
    typename ViewTypeB::array_type b_flat = b;
    typename ViewTypeC::array_type c_flat = c;
    check_deriv(A_flat, b_flat, c_flat);
  }
#endif

  return perf;
}

template <typename ... ViewArgs>
Perf
do_time_analytic(const size_t m, const size_t n, const size_t p,
                 const size_t nloop, const bool check)
{
  typedef Kokkos::View<double***, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv( A, b, c );
  execution_space().fence();

  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_deriv(A,b,c);

  return perf;
}

template <int MaxP, typename ... ViewArgs>
Perf
do_time_analytic_sl(const size_t m, const size_t n, const size_t p,
                    const size_t nloop, const bool check)
{
  typedef Kokkos::View<double***, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv_sl<MaxP>( A, b, c );
  execution_space().fence();

  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv_sl<MaxP>( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_deriv(A,b,c);

  return perf;
}

template <int p, typename ... ViewArgs>
Perf
do_time_analytic_s(const size_t m, const size_t n,
                   const size_t nloop, const bool check)
{
  typedef Kokkos::View<double**[p+1], ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  ViewTypeA A("A",m,n);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv_s<p>( A, b, c );
  execution_space().fence();

  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv_s<p>( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_deriv(A,b,c);

  return perf;
}

typedef Sacado::Fad::SFad<double,SFadSize> SFad_type;
typedef Sacado::Fad::SLFad<double,SLFadSize> SLFad_type;
typedef Sacado::Fad::DFad<double> DFad_type;

#define INST_FUNC_VAL_DEV(DEV) \
  template Perf do_time_val< Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check ); \
  template Perf do_time_val< Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check ); \
  template Perf do_time_val< DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check ); \
  template Perf do_time_analytic< Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic< Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic< DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic_sl< SLFadSize, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic_sl< SLFadSize, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic_sl< SLFadSize, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check); \
  template Perf do_time_analytic_s< SFadSize, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check); \
  template Perf do_time_analytic_s< SFadSize, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check); \
  template Perf do_time_analytic_s< SFadSize, DEV > ( const size_t m, const size_t n, const size_t nloop, const bool check);

#define INST_FUNC_FAD_DEV(FAD,DEV)                                      \
  template Perf do_time_fad< FAD, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad< FAD, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad< FAD, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_scratch< FAD, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_scratch< FAD, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_scratch< FAD, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check );

#define INST_FUNC_DEV(DEV)                                       \
  INST_FUNC_VAL_DEV( DEV )                                       \
  INST_FUNC_FAD_DEV( SFad_type, DEV )   \
  INST_FUNC_FAD_DEV( SLFad_type, DEV ) \
  INST_FUNC_FAD_DEV( DFad_type, DEV )

#ifdef KOKKOS_ENABLE_SERIAL
INST_FUNC_DEV(Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
INST_FUNC_DEV(Kokkos::OpenMP)
#endif

#ifdef KOKKOS_ENABLE_THREADS
INST_FUNC_DEV(Kokkos::Threads)
#endif

#ifdef KOKKOS_ENABLE_CUDA
INST_FUNC_DEV(Kokkos::Cuda)
#endif

#ifdef KOKKOS_ENABLE_HIP
INST_FUNC_DEV(Kokkos::HIP)
#endif
