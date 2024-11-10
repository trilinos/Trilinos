// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#define SACADO_VIEW_CUDA_HIERARCHICAL 1
#define SACADO_ALIGN_SFAD 1

#include "Sacado.hpp"

#include "mat_vec_hierarchical.hpp"

#include "Kokkos_Timer.hpp"

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec_hierarchical(const ViewTypeA& A, const ViewTypeB& b,
                              const ViewTypeC& c) {
  typedef typename Kokkos::ThreadLocalScalarType<ViewTypeC>::type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;

#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
  const unsigned vector_size = is_cuda ? 32 : 1;
  const unsigned team_size = is_cuda ? 128 / vector_size : 1;
#elif defined (KOKKOS_ENABLE_HIP)
  const bool is_hip = std::is_same<execution_space, Kokkos::HIP>::value;
  const unsigned vector_size = is_hip ? 64 : 1;
  const unsigned team_size = is_hip ? 128 / vector_size : 1;
#else
  const unsigned vector_size = 1;
  const unsigned team_size = 1;
#endif

  const int m = A.extent(0);
  const int n = A.extent(1);
  const int range = (m+team_size-1)/team_size;

  typedef Kokkos::TeamPolicy<execution_space> Policy;
  Kokkos::parallel_for(
    Policy( range,team_size,vector_size ),
    KOKKOS_LAMBDA (const typename Policy::member_type& team) {
      const int i = team.league_rank()*team.team_size() + team.team_rank();
      if (i >= m)
        return;

      scalar_type t = 0.0;
      for (int j=0; j<n; ++j)
        t += A(i,j)*b(j);
      c(i) = t;
    }
  );
}

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
check_deriv_hierarchical(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  const double tol = 1.0e-14;
  typedef typename ViewTypeC::value_type value_type;
  typename ViewTypeC::HostMirror h_c = Kokkos::create_mirror_view(c);
  Kokkos::deep_copy(h_c, c);
  const size_t m = A.extent(0);
  const size_t n = A.extent(1);
  const size_t p = Kokkos::dimension_scalar(A);
  for (size_t i=0; i<m; ++i) {
    for (size_t j=0; j<p; ++j) {
      value_type t = (j == p-1 ? n : 2*n);
      if (std::abs(h_c(i).fastAccessDx(j)- t) > tol) {
        std::cout << "Comparison failed!  " << i << "," << j << " : "
                  << h_c(i).fastAccessDx(j) << " , " << t << std::endl;
      }
    }
  }
}

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad_hierarchical(const size_t m, const size_t n, const size_t p,
                         const size_t nloop, const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
  const int FadStride = is_cuda ? 32 : 1;
#elif defined (KOKKOS_ENABLE_HIP)
  const bool is_hip = std::is_same<execution_space, Kokkos::HIP>::value;
  const int FadStride = is_hip ? 64 : 1;
#else
  const int FadStride = 1;
#endif

#if defined(SACADO_ALIGN_SFAD)
  const int N = Sacado::StaticSize<FadType>::value;
  const int Nalign = ((N+FadStride-1)/FadStride)*FadStride;
  const size_t pa = N > 0 ? ((p+FadStride-1)/FadStride)*FadStride : p;
  typedef typename FadType::template apply_N<Nalign>::type AlignedFadType;
#else
  typedef FadType AlignedFadType;
  const size_t pa = p;
#endif

  typedef Kokkos::LayoutContiguous<typename ViewTypeA::array_layout,FadStride> ConLayoutA;
  typedef Kokkos::LayoutContiguous<typename ViewTypeB::array_layout,FadStride> ConLayoutB;
  typedef Kokkos::LayoutContiguous<typename ViewTypeC::array_layout,FadStride> ConLayoutC;

  typedef Kokkos::View<AlignedFadType**, ConLayoutA, execution_space> ConViewTypeA;
  typedef Kokkos::View<AlignedFadType*,  ConLayoutB, execution_space> ConViewTypeB;
  typedef Kokkos::View<AlignedFadType*,  ConLayoutC, execution_space> ConViewTypeC;

  ConViewTypeA A("A",m,n,pa+1);
  ConViewTypeB b("B",n,pa+1);
  ConViewTypeC c("c",m,pa+1);

  // AlignedFadType a(pa, 1.0);
  // for (size_t k=0; k<pa; ++k)
  //   a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(typename ConViewTypeA::array_type(A), 1.0);
  Kokkos::deep_copy(typename ConViewTypeB::array_type(b), 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_hierarchical( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_hierarchical( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check) {
    check_deriv_hierarchical(A, b, c);
  }

  return perf;
}

typedef Sacado::Fad::SFad<double,HierSFadSize> SFad_type;
typedef Sacado::Fad::SLFad<double,HierSLFadSize> SLFad_type;

#define INST_FUNC_FAD_DEV(FAD,DEV)                                      \
  template Perf do_time_fad_hierarchical< FAD, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical< FAD, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical< FAD, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check );

#define INST_FUNC_DEV(DEV)              \
  INST_FUNC_FAD_DEV( SFad_type, DEV )   \
  INST_FUNC_FAD_DEV( SLFad_type, DEV )

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
