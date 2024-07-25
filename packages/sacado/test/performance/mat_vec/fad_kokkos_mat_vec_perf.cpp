// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#define SACADO_VIEW_CUDA_HIERARCHICAL 1
//#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD 1
//#define SACADO_KOKKOS_USE_MEMORY_POOL 1
#define SACADO_ALIGN_SFAD 1

//#define SACADO_DISABLE_FAD_VIEW_SPEC
#include "Sacado.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Time.hpp"

#include "Kokkos_Timer.hpp"

// For vtune
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>

// A performance test that computes the derivative of a simple Kokkos kernel
// using various Fad classes

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

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec_hierarchical(const ViewTypeA& A, const ViewTypeB& b,
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

#elif defined(SACADO_VIEW_CUDA_HIERARCHICAL)

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

#else

template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec_hierarchical(const ViewTypeA& A, const ViewTypeB& b,
                              const ViewTypeC& c) {
  typedef typename ViewTypeC::value_type scalar_type;
  typedef typename ViewTypeC::execution_space execution_space;

#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
#else
  const bool is_cuda = false;
#endif
  const unsigned vector_size = 1;
  const unsigned team_size = is_cuda ? 128 / vector_size : 1;

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

#endif

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

struct Perf {
  double time;
  double flops;
  double throughput;
};

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad(const size_t m, const size_t n, const size_t p, const size_t nloop,
            const bool check)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;
  typedef typename ViewTypeA::execution_space execution_space;

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  FadType a(p, 1.0);
  for (size_t k=0; k<p; ++k)
    a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(A, a);
  Kokkos::deep_copy(b, a);

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

  if (check) {
    check_deriv(A, b, c);
  }

  return perf;
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

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL)
#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
  const int FadStride = is_cuda ? 32 : 1;
#elif defined(KOKKOS_ENABLE_HIP)
  const bool is_hip = std::is_same<execution_space, Kokkos::HIP>::value;
  const int FadStride = is_hip ? 64 : 1;
#else
  const int FadStride 1;
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
#else
  const int FadStride = 1;
  typedef FadType AlignedFadType;
  const size_t pa = p;
#endif

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  typedef Kokkos::LayoutContiguous<typename ViewTypeA::array_layout,FadStride> ConLayoutA;
  typedef Kokkos::LayoutContiguous<typename ViewTypeB::array_layout,FadStride> ConLayoutB;
  typedef Kokkos::LayoutContiguous<typename ViewTypeC::array_layout,FadStride> ConLayoutC;
#else
  typedef typename ViewTypeA::array_layout ConLayoutA;
  typedef typename ViewTypeB::array_layout ConLayoutB;
  typedef typename ViewTypeC::array_layout ConLayoutC;
  (void) FadStride;
#endif


  typedef Kokkos::View<AlignedFadType**, ConLayoutA, execution_space> ConViewTypeA;
  typedef Kokkos::View<AlignedFadType*,  ConLayoutB, execution_space> ConViewTypeB;
  typedef Kokkos::View<AlignedFadType*,  ConLayoutC, execution_space> ConViewTypeC;

  ConViewTypeA A("A",m,n,pa+1);
  ConViewTypeB b("B",n,pa+1);
  ConViewTypeC c("c",m,pa+1);

  AlignedFadType a(pa, 1.0);
  for (size_t k=0; k<pa; ++k)
    a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(A, a);
  Kokkos::deep_copy(b, a);

  Kokkos::Timer wall_clock;
  Perf perf;

#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  const size_t concurrency = execution_space().concurrency();

#if defined (KOKKOS_ENABLE_CUDA)
  const size_t warp_dim = is_cuda ? 32 : 1;
#elif defined (KOKKOS_ENABLE_HIP)
  const size_t warp_dim = is_hip ? 64 : 1;
#else
  const size_t warp_dim = 1;
#endif

  const size_t block_size = pa*sizeof(double);
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
#endif

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
    check_deriv(A, b, c);
  }

#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  Sacado::destroyGlobalMemoryPool(space);
#endif

  return perf;
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

void
print_perf(const Perf& perf, const Perf& perf_base, const size_t p,
           const std::string& name)
{
  std::cout << name << "\t "
            << perf.time << "\t "
            << perf.throughput << "\t "
            << perf.time / (perf_base.time*p)
            << std::endl;
}

template <int SFadSize, int SLFadSize, int HierSFadSize, int HierSLFadSize,
          typename ... ViewArgs>
void
do_times(const size_t m,
         const size_t n,
         const size_t p,
         const size_t ph,
         const size_t nloop,
         const bool value,
         const bool sfad,
         const bool slfad,
         const bool dfad,
         const bool hierarchical,
         const bool check)
{
  Perf perf_value;
  perf_value.time = 1.0;

  // Run value
  if (value) {
    Perf perf = do_time_val<ViewArgs...>(m,n,nloop,check);
    perf_value = perf;
    print_perf(perf, perf_value, p, "Value     ");
  }

  // Run SFad
  if (sfad && p == SFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SFad<double,SFadSize>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_value, p, "SFad      ");
  }

  // Run SLFad
  if (slfad && p <= SLFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SLFad<double,SLFadSize>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_value, p, "SLFad     ");
  }

  // Run DFad
  if (dfad) {
    Perf perf =
      do_time_fad<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_value, p, "DFad      ");
  }

  // Run hierarchical
  if (hierarchical) {
    if (sfad && ph == HierSFadSize) {
      Perf perf =
        do_time_fad_hierarchical<Sacado::Fad::SFad<double,HierSFadSize>, ViewArgs...>(m,n,ph,nloop,check);
      print_perf(perf, perf_value, ph, "Hier SFad ");
    }
    if (slfad && ph <= HierSLFadSize) {
      Perf perf =
        do_time_fad_hierarchical<Sacado::Fad::SLFad<double,HierSLFadSize>, ViewArgs...>(m,n,ph,nloop,check);
      print_perf(perf, perf_value, ph, "Hier SLFad");
    }
    if (dfad) {
      Perf perf =
        do_time_fad_hierarchical<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,ph,nloop,check);
      print_perf(perf, perf_value, ph, "Hier DFad ");
    }
  }

}

enum LayoutType {
  LAYOUT_LEFT=0,
  LAYOUT_RIGHT,
  LAYOUT_DEFAULT
};
const int num_layout_types = 3;
const LayoutType layout_values[] = {
  LAYOUT_LEFT, LAYOUT_RIGHT, LAYOUT_DEFAULT };
const char *layout_names[] = { "left", "right", "default" };

template <int SFadSize, int SLFadSize, int HierSFadSize, int HierSLFadSize,
          typename Device>
void
do_times_layout(const size_t m,
                const size_t n,
                const size_t p,
                const size_t ph,
                const size_t nloop,
                const bool value,
                const bool sfad,
                const bool slfad,
                const bool dfad,
                const bool hierarchical,
                const bool check,
                const LayoutType& layout,
                const std::string& device)
{
  int prec = 2;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(prec);
  std::cout << std::endl
            << device
            << " performance for layout "
            << layout_names[layout]
            << " m = " << m << " n = " << n << " p = " << p << " ph = " << ph
            << std::endl << std::endl;
  std::cout << "Computation \t Time     \t Throughput \t Ratio" << std::endl;

  if (layout == LAYOUT_LEFT)
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::LayoutLeft,Device>(
      m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check);
  else if (layout == LAYOUT_RIGHT)
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::LayoutRight,Device>(
      m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check);
  else
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Device>
      (m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check);
}

// Connect executable to vtune for profiling
void connect_vtune() {
  std::stringstream cmd;
  pid_t my_os_pid=getpid();
  const std::string vtune_loc =
    "amplxe-cl";
  const std::string output_dir = "./vtune";
  cmd << vtune_loc
      << " -collect hotspots -result-dir " << output_dir
      << " -target-pid " << my_os_pid << " &";
  std::cout << cmd.str() << std::endl;
  system(cmd.str().c_str());
  system("sleep 10");
}

//const int SFadSize  = 8;
const int SFadSize  = 32;
const int SLFadSize = SFadSize;
//const int HierSFadSize  = 50;
const int HierSFadSize  = 32;
const int HierSLFadSize = HierSFadSize;

int main(int argc, char* argv[]) {
  bool success = true;
  try {

    // Set up command line options
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("This program tests the speed of various forward mode AD implementations for simple Kokkos kernel");
    int m = 100000;
    clp.setOption("m", &m, "Number of matrix rows");
    int n = 100;
    clp.setOption("n", &n, "Number of matrix columns");
    int p = SFadSize;
    clp.setOption("p", &p, "Number of derivative components");
    int ph = HierSFadSize;
    clp.setOption("ph", &ph, "Number of derivative components for hierarchical");
    int nloop = 10;
    clp.setOption("nloop", &nloop, "Number of loops");
#ifdef KOKKOS_ENABLE_SERIAL
    bool serial = 0;
    clp.setOption("serial", "no-serial", &serial, "Whether to run Serial");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    int openmp = 0;
    clp.setOption("openmp", &openmp, "Number of OpenMP threads");
#endif
#ifdef KOKKOS_ENABLE_THREADS
    int threads = 0;
    clp.setOption("threads", &threads, "Number of pThreads threads");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = 0;
    clp.setOption("cuda", "no-cuda", &cuda, "Whether to run CUDA");
#endif
#ifdef KOKKOS_ENABLE_HIP
    bool hip = 0;
    clp.setOption("hip", "no-hip", &cuda, "Whether to run HIP");
#endif
    int numa = 0;
    clp.setOption("numa", &numa,
                  "Number of NUMA domains to use (set to 0 to use all NUMAs");
    int cores_per_numa = 0;
    clp.setOption("cores-per-numa", &cores_per_numa,
                  "Number of CPU cores per NUMA to use (set to 0 to use all cores)");
    bool print_config = false;
    clp.setOption("print-config", "no-print-config", &print_config,
                  "Whether to print Kokkos device configuration");
    LayoutType layout = LAYOUT_DEFAULT;
    clp.setOption("layout", &layout, num_layout_types, layout_values,
                  layout_names, "View layout");
    bool vtune = false;
    clp.setOption("vtune", "no-vtune", &vtune, "Profile with vtune");
    bool value = true;
    clp.setOption("value", "no-value", &value, "Run value calculation");
    bool sfad = true;
    clp.setOption("sfad", "no-sfad", &sfad, "Run SFad derivative calculation");
    bool slfad = true;
    clp.setOption("slfad", "no-slfad", &slfad, "Run SLFad derivative calculation");
    bool dfad = true;
    clp.setOption("dfad", "no-dfad", &dfad, "Run DFad derivative calculation");
    bool hierarchical = true;
    clp.setOption("hierarchical", "no-hierarchical", &hierarchical, "Run hierarchical Fad derivative calculation");
    bool check = false;
    clp.setOption("check", "no-check", &check, "Check calculations are correct");

    // Parse options
    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
        return 0;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
        return 1;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
        break;
    }

    if (vtune)
      connect_vtune();

    Kokkos::InitializationSettings init_args;
    init_args.set_num_threads(cores_per_numa);

    Kokkos::initialize(init_args);

    if (print_config)
      Kokkos::print_configuration(std::cout, true);

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Serial>(
        m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check,layout,"Serial");
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::OpenMP>(
        m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check,layout,"OpenMP");
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Threads>(
        m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check,layout,"Threads");
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Cuda>(
        m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check,layout,"Cuda");
    }
#endif

#ifdef KOKKOS_ENABLE_HIP
    if (hip) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::HIP>(
        m,n,p,ph,nloop,value,sfad,slfad,dfad,hierarchical,check,layout,"HIP");
    }
#endif

    Kokkos::finalize();

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return !success;
}
