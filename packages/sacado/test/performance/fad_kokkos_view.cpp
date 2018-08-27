// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

//#define SACADO_DISABLE_FAD_VIEW_SPEC
#include "Sacado.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Time.hpp"

#include "impl/Kokkos_Timer.hpp"

// For vtune
#include <sys/types.h>
#include <unistd.h>

// A performance test that computes the derivative of a simple Kokkos kernel
// using various Fad classes

// Our Kokkos kernel: computes c = A*b for A mxn and b nx1
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
struct MatVecFunctor {

  // The scalar type used in this calculation, e.g., double
  typedef typename ViewTypeC::value_type scalar_type;

  // The best ordinal type for the architecture we are running on,
  // e.g., int or size_t
  //typedef typename ViewTypeC::size_type size_type;
  typedef int size_type;

  // The execution space where this functor will run
  typedef typename ViewTypeC::execution_space execution_space;

  // Data needed by functor
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_type n;

  // Constructor
  MatVecFunctor(const ViewTypeA& A_arg,
                const ViewTypeB& b_arg,
                const ViewTypeC& c_arg) :
    A(A_arg), b(b_arg), c(c_arg), n(A.extent(1))
  {}

  // Function to compute matrix-vector product for a given row i
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    scalar_type t = 0.0;
    for (size_type j=0; j<n; ++j)
      t += A(i,j)*b(j);
    c(i) = t;
  }

};

// Computes the derivative of c = A*b for A mxnx(p+1) and b nx1x(p+1)
// where p is the number of derivatives
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
struct MatVecDerivFunctor {

  // The scalar type used in this calculation, e.g., double
  typedef typename ViewTypeC::value_type scalar_type;

  // The best ordinal type for the architecture we are running on,
  // e.g., int or size_t
  //typedef typename ViewTypeC::size_type size_type;
  typedef int size_type;

  // The execution space where this functor will run
  typedef typename ViewTypeC::execution_space execution_space;

  // Data needed by functor
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_type n;
  const size_type p;

  // Constructor
  MatVecDerivFunctor(const ViewTypeA& A_arg,
                     const ViewTypeB& b_arg,
                     const ViewTypeC& c_arg) :
    A(A_arg), b(b_arg), c(c_arg), n(A.extent(1)), p(A.extent(2)-1)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    c(i,p) = 0.0;
    for (size_type k=0; k<p; ++k)
      c(i,k) = 0.0;
    for (size_type j=0; j<n; ++j) {
      c(i,p) += A(i,j,p)*b(j,p);
      for (size_type k=0; k<p; ++k) {
        c(i,k) += A(i,j,k)*b(j,p) + A(i,j,p)*b(j,k);
      }
    }
  }

};

// Computes the derivative of c = A*b for A mxnx(p+1) and b nx1x(p+1)
// where p is the number of derivatives
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC,
          int MaxP>
struct SLMatVecDerivFunctor {

  // The scalar type used in this calculation, e.g., double
  typedef typename ViewTypeC::value_type scalar_type;

  // The best ordinal type for the architecture we are running on,
  // e.g., int or size_t
  //typedef typename ViewTypeC::size_type size_type;
  typedef int size_type;

  // The execution space where this functor will run
  typedef typename ViewTypeC::execution_space execution_space;

  // Data needed by functor
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_type n;
  const size_type p;

  // Constructor
  SLMatVecDerivFunctor(const ViewTypeA& A_arg,
                       const ViewTypeB& b_arg,
                       const ViewTypeC& c_arg) :
    A(A_arg), b(b_arg), c(c_arg), n(A.extent(1)), p(A.extent(2)-1)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    scalar_type cv = 0.0;
    scalar_type t[MaxP];
    for (size_type k=0; k<p; ++k)
      t[k] = 0.0;

    for (size_type j=0; j<n; ++j) {
      scalar_type av = A(i,j,p);
      scalar_type bv = b(j,p);
      cv += av*bv;
      for (size_type k=0; k<p; ++k) {
        t[k] += A(i,j,k)*bv + av*b(j,k);
      }
    }

    for (size_type k=0; k<p; ++k)
      c(i,k) = t[k];
    c(i,p) = cv;
  }

};

// Computes the derivative of c = A*b for A mxnx(p+1) and b nx1x(p+1)
// where p is the number of derivatives
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC,
          int p>
struct SMatVecDerivFunctor {

  // The scalar type used in this calculation, e.g., double
  typedef typename ViewTypeC::value_type scalar_type;

  // The best ordinal type for the architecture we are running on,
  // e.g., int or size_t
  //typedef typename ViewTypeC::size_type size_type;
  typedef int size_type;

  // The execution space where this functor will run
  typedef typename ViewTypeC::execution_space execution_space;

  // Data needed by functor
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_type n;

  // Constructor
  SMatVecDerivFunctor(const ViewTypeA& A_arg,
                      const ViewTypeB& b_arg,
                      const ViewTypeC& c_arg) :
    A(A_arg), b(b_arg), c(c_arg), n(A.extent(1))
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    scalar_type cv = 0.0;
    scalar_type t[p];
    for (size_type k=0; k<p; ++k)
      t[k] = 0.0;

    for (size_type j=0; j<n; ++j) {
      const scalar_type av = A(i,j,p);
      const scalar_type bv = b(j,p);
      cv += av*bv;

// Using simd here results in much better performance.  Othewise the compiler
// appears to try and vectorize the j loop with gather instructions, which
// doesn't work very well.
#if defined(__INTEL_COMPILER) && ! defined(__CUDA_ARCH__)
#pragma simd
#endif
      for (size_type k=0; k<p; ++k) {
        t[k] += A(i,j,k)*bv + av*b(j,k);
      }
    }

    for (size_type k=0; k<p; ++k)
      c(i,k) = t[k];
    c(i,p) = cv;
  }

};

// Create a mat-vec functor from given A, b, c
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  MatVecFunctor<ViewTypeA, ViewTypeB, ViewTypeC> f( A, b, c );
  Kokkos::parallel_for( A.extent(0), f );
}

// Create a mat-vec derivative functor from given A, b, c
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  MatVecDerivFunctor<ViewTypeA, ViewTypeB, ViewTypeC> f( A, b, c );
  Kokkos::parallel_for( A.extent(0), f );
}

// Create a mat-vec derivative functor from given A, b, c
template <int MaxP, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv_sl(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  SLMatVecDerivFunctor<ViewTypeA, ViewTypeB, ViewTypeC, MaxP> f( A, b, c );
  Kokkos::parallel_for( A.extent(0), f );
}

// Create a mat-vec derivative functor from given A, b, c
template <int p, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv_s(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  SMatVecDerivFunctor<ViewTypeA, ViewTypeB, ViewTypeC, p> f( A, b, c );
  Kokkos::parallel_for( A.extent(0), f );
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

  Kokkos::Impl::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );
  execution_space::fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  execution_space::fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check) {
    typename ViewTypeA::array_type A_flat = A;
    typename ViewTypeB::array_type b_flat = b;
    typename ViewTypeC::array_type c_flat = c;
    check_deriv(A_flat, b_flat, c_flat);
  }

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

  Kokkos::Impl::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv( A, b, c );
  execution_space::fence();

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv( A, b, c );
  }
  execution_space::fence();
  timer.stop();

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

  Kokkos::Impl::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv_sl<MaxP>( A, b, c );
  execution_space::fence();

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv_sl<MaxP>( A, b, c );
  }
  execution_space::fence();
  timer.stop();

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

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  Kokkos::Impl::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec_deriv_s<p>( A, b, c );
  execution_space::fence();

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv_s<p>( A, b, c );
  }
  execution_space::fence();
  timer.stop();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_deriv(A,b,c);

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

  Kokkos::Impl::Timer wall_clock;
  Perf perf;

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );
  execution_space::fence();

  wall_clock.reset();
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  execution_space::fence();

  perf.time = wall_clock.seconds() / nloop;
  perf.flops = m*n*2;
  perf.throughput = perf.flops / perf.time / 1.0e9;

  if (check)
    check_val(A,b,c);

  return perf;
}

void
print_perf(const Perf& perf, const Perf& perf_base, const std::string& name)
{
  std::cout << name << "\t "
            << perf.time << "\t "
            << perf.throughput << "\t "
            << perf.time / perf_base.time
            << std::endl;
}

template <int SFadSize, int SLFadSize, typename ... ViewArgs>
void
do_times(const size_t m,
         const size_t n,
         const size_t p,
         const size_t nloop,
         const bool value,
         const bool analytic,
         const bool sfad,
         const bool slfad,
         const bool dfad,
         const bool check)
{
  Perf perf_analytic;
  perf_analytic.time = 1.0;

  // Run analytic
  if (analytic) {
    perf_analytic = do_time_analytic<ViewArgs...>(m,n,p,nloop,check);
  }

  // Run value
  if (value) {
    Perf perf = do_time_val<ViewArgs...>(m,n,nloop,check);
    print_perf(perf, perf_analytic, "Value     ");
  }

  if (analytic) {
    print_perf(perf_analytic, perf_analytic, "Analytic  ");
  }

  if(analytic && p == SFadSize) {
    Perf perf =
      do_time_analytic_s<SFadSize, ViewArgs...>(m,n,nloop,check);
    print_perf(perf, perf_analytic, "Analytic-s");
  }

  if(analytic && p <= SLFadSize) {
    Perf perf =
      do_time_analytic_sl<SLFadSize, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_analytic, "Analytic-sl");
  }

  // Run SFad
  if (sfad && p == SFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SFad<double,SFadSize>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_analytic, "SFad      ");
  }

  // Run SLFad
  if (slfad && p <= SLFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SLFad<double,SLFadSize>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_analytic, "SLFad     ");
  }

  // Run DFad
  if (dfad) {
    Perf perf =
      do_time_fad<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
    print_perf(perf, perf_analytic, "DFad      ");
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

template <int SFadSize, int SLFadSize, typename Device>
void
do_times_layout(const size_t m,
                const size_t n,
                const size_t p,
                const size_t nloop,
                const bool value,
                const bool analytic,
                const bool sfad,
                const bool slfad,
                const bool dfad,
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
            << " m = " << m << " n = " << n << " p = " << p
            << std::endl << std::endl;
  std::cout << "Computation \t Time     \t Throughput \t Ratio" << std::endl;

  if (layout == LAYOUT_LEFT)
    do_times<SFadSize,SLFadSize,Kokkos::LayoutLeft,Device>(
      m,n,p,nloop,value,analytic,sfad,slfad,dfad,check);
  else if (layout == LAYOUT_RIGHT)
    do_times<SFadSize,SLFadSize,Kokkos::LayoutRight,Device>(
      m,n,p,nloop,value,analytic,sfad,slfad,dfad,check);
  else
    do_times<SFadSize,SLFadSize,Device>
      (m,n,p,nloop,value,analytic,sfad,slfad,dfad,check);
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

const int SFadSize  = 8;
const int SLFadSize = SFadSize;

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
    bool analytic = true;
    clp.setOption("analytic", "no-analytic", &analytic,
                  "Run analytic derivative calculation");
    bool sfad = true;
    clp.setOption("sfad", "no-sfad", &sfad, "Run SFad derivative calculation");
    bool slfad = true;
    clp.setOption("slfad", "no-slfad", &slfad, "Run SLFad derivative calculation");
#if defined(KOKKOS_ENABLE_CUDA_UVM)
    bool dfad = true;
    clp.setOption("dfad", "no-dfad", &dfad, "Run DFad derivative calculation");
#else
    bool dfad = false;
#endif
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

    Kokkos::InitArguments init_args;
    init_args.num_threads = cores_per_numa;
    init_args.num_numa = numa;

    Kokkos::initialize(init_args);

    if (print_config)
      Kokkos::print_configuration(std::cout, true);

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      do_times_layout<SFadSize,SLFadSize,Kokkos::Serial>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,check,layout,"Serial");
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      do_times_layout<SFadSize,SLFadSize,Kokkos::OpenMP>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,check,layout,"OpenMP");
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      do_times_layout<SFadSize,SLFadSize,Kokkos::Threads>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,check,layout,"Threads");
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      do_times_layout<SFadSize,SLFadSize,Kokkos::Cuda>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,check,layout,"Cuda");
    }
#endif
    Kokkos::finalize();

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return !success;
}
