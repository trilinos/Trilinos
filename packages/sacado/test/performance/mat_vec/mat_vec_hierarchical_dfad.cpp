// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD 1
#define SACADO_KOKKOS_USE_MEMORY_POOL 1

#include "Sacado.hpp"

#include "mat_vec_hierarchical_dfad.hpp"

#include "Kokkos_Timer.hpp"

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec_hierarchical_dfad(const ViewTypeA& A, const ViewTypeB& b,
                                   const ViewTypeC& c) {
  typedef typename ViewTypeC::value_type scalar_type;
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
void run_mat_vec_hierarchical_dfad_scratch(
  const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c) {
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;
  typedef Kokkos::TeamPolicy<execution_space> Policy;
  typedef typename Policy::member_type team_member;
  typedef Kokkos::View<scalar_type*,Kokkos::LayoutLeft, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged> TmpScratchSpace;

#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
  const unsigned VectorSize = is_cuda ? 32 : 1;
  const unsigned TeamSize = is_cuda ? 128 / VectorSize : 1;
#elif defined (KOKKOS_ENABLE_HIP)
  const bool is_hip = std::is_same<execution_space, Kokkos::HIP>::value;
  const unsigned VectorSize = is_hip ? 64 : 1;
  const unsigned TeamSize = is_hip ? 128 / VectorSize : 1;
#else
  const unsigned VectorSize = 1;
  const unsigned TeamSize = 1;
#endif

  const int m = A.extent(0);
  const int n = A.extent(1);
  const int p = dimension_scalar(A);
  const int N = (m+TeamSize-1)/TeamSize;

  Policy policy(N, TeamSize, VectorSize);
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
check_deriv_hierarchical_dfad(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
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
do_time_fad_hierarchical_dfad(const size_t m, const size_t n, const size_t p,
                              const size_t nloop, const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;
  typedef Kokkos::LayoutContiguous<typename ViewTypeA::array_layout> ConLayoutA;
  typedef Kokkos::LayoutContiguous<typename ViewTypeB::array_layout> ConLayoutB;
  typedef Kokkos::LayoutContiguous<typename ViewTypeC::array_layout> ConLayoutC;
  typedef Kokkos::View<FadType**, ConLayoutA, execution_space> ConViewTypeA;
  typedef Kokkos::View<FadType*,  ConLayoutB, execution_space> ConViewTypeB;
  typedef Kokkos::View<FadType*,  ConLayoutC, execution_space> ConViewTypeC;

  ConViewTypeA A("A",m,n,p+1);
  ConViewTypeB b("B",n,p+1);
  ConViewTypeC c("c",m,p+1);

  // FadType a(p, 1.0);
  // for (size_t k=0; k<p; ++k)
  //   a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(typename ConViewTypeA::array_type(A), 1.0);
  Kokkos::deep_copy(typename ConViewTypeB::array_type(b), 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
  const size_t warp_dim = is_cuda ? 32 : 1;
#elif defined (KOKKOS_ENABLE_HIP)
  const bool is_hip = std::is_same<execution_space, Kokkos::HIP>::value;
  const size_t warp_dim = is_hip ? 64 : 1;
#else
  const size_t warp_dim = 1;
#endif

  const size_t concurrency = execution_space().concurrency();
  const size_t block_size = p*sizeof(double);
  const size_t nkernels = concurrency / warp_dim;
  const size_t mem_pool_size =
    static_cast<size_t>(1.2*nkernels*block_size);
  const size_t superblock_size = std::max<size_t>(nkernels / 100, 1) * block_size;
  execution_space space;
  Sacado::createGlobalMemoryPool(space, mem_pool_size,
      block_size,
      block_size,
      superblock_size
      );

  // Execute the kernel once to warm up
  run_mat_vec_hierarchical_dfad( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_hierarchical_dfad( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check) {
    check_deriv_hierarchical_dfad(A, b, c);
  }

  Sacado::destroyGlobalMemoryPool(space);

  return perf;
}

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad_hierarchical_dfad_scratch(
  const size_t m, const size_t n, const size_t p, const size_t nloop,
  const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;
  typedef Kokkos::LayoutContiguous<typename ViewTypeA::array_layout> ConLayoutA;
  typedef Kokkos::LayoutContiguous<typename ViewTypeB::array_layout> ConLayoutB;
  typedef Kokkos::LayoutContiguous<typename ViewTypeC::array_layout> ConLayoutC;
  typedef Kokkos::View<FadType**, ConLayoutA, execution_space> ConViewTypeA;
  typedef Kokkos::View<FadType*,  ConLayoutB, execution_space> ConViewTypeB;
  typedef Kokkos::View<FadType*,  ConLayoutC, execution_space> ConViewTypeC;

  ConViewTypeA A("A",m,n,p+1);
  ConViewTypeB b("B",n,p+1);
  ConViewTypeC c("c",m,p+1);

  // FadType a(p, 1.0);
  // for (size_t k=0; k<p; ++k)
  //   a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(typename ConViewTypeA::array_type(A), 1.0);
  Kokkos::deep_copy(typename ConViewTypeB::array_type(b), 1.0);

  Kokkos::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_hierarchical_dfad_scratch( A, b, c );
  execution_space().fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_hierarchical_dfad_scratch( A, b, c );
  }
  execution_space().fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check) {
    check_deriv_hierarchical_dfad(A, b, c);
  }

  return perf;
}

typedef Sacado::Fad::DFad<double> DFad_type;

#define INST_FUNC_FAD_DEV(FAD,DEV)                                      \
  template Perf do_time_fad_hierarchical_dfad< FAD, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical_dfad< FAD, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical_dfad< FAD, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical_dfad_scratch< FAD, Kokkos::LayoutLeft, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical_dfad_scratch< FAD, Kokkos::LayoutRight, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check ); \
  template Perf do_time_fad_hierarchical_dfad_scratch< FAD, DEV > ( const size_t m, const size_t n, const size_t p, const size_t nloop, const bool check );

#define INST_FUNC_DEV(DEV) \
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
