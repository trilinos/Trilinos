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

#include "Sacado.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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
  typedef typename ViewTypeC::size_type size_type;

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
    A(A_arg), b(b_arg), c(c_arg), n(A.dimension_1())
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
  typedef typename ViewTypeC::size_type size_type;

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
    A(A_arg), b(b_arg), c(c_arg), n(A.dimension_1()), p(A.dimension_2()-1)
  {}

  // Function to compute matrix-vector product for a given row i
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    scalar_type t = 0.0;
    for (size_type j=0; j<n; ++j)
      t += A(i,j,0)*b(j,0);
    c(i,0) = t;
    for (size_type k=0; k<p; ++k) {
      t = 0.0;
      for (size_type j=0; j<n; ++j)
        t += A(i,j,k+1)*b(j,0) + A(i,j,0)*b(j,k+1);
      c(i,k+1) = t;
    }
  }

};

// Create a mat-vec functor from given A, b, c
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  MatVecFunctor<ViewTypeA, ViewTypeB, ViewTypeC> f( A, b, c );
  Kokkos::parallel_for( A.dimension_0(), f );
  Kokkos::fence();
}

// Create a mat-vec derivative functor from given A, b, c
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void
run_mat_vec_deriv(const ViewTypeA& A, const ViewTypeB& b, const ViewTypeC& c)
{
  MatVecDerivFunctor<ViewTypeA, ViewTypeB, ViewTypeC> f( A, b, c );
  Kokkos::parallel_for( A.dimension_0(), f );
  Kokkos::fence();
}

struct Perf {
  double time;
  double flops;
  double throughput;
};

template <typename FadType, typename ... ViewArgs>
Perf
do_time_fad(const size_t m, const size_t n, const size_t p, const size_t nloop)
{
  typedef Kokkos::View<FadType**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<FadType*,  ViewArgs...> ViewTypeC;

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  FadType a(p, 1.0);
  for (size_t k=0; k<p; ++k)
    a.fastAccessDx(k) = 1.0;
  Kokkos::deep_copy(A, a);
  Kokkos::deep_copy(b, a);

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  timer.stop();

  Perf perf;
  perf.time = timer.totalElapsedTime() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  return perf;
}

template <typename ... ViewArgs>
Perf
do_time_analytic(const size_t m, const size_t n, const size_t p, const size_t nloop)
{
  typedef Kokkos::View<double***, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double**,  ViewArgs...> ViewTypeC;

  ViewTypeA A("A",m,n,p+1);
  ViewTypeB b("B",n,p+1);
  ViewTypeC c("c",m,p+1);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  // Execute the kernel once to warm up
  run_mat_vec_deriv( A, b, c );

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec_deriv( A, b, c );
  }
  timer.stop();

  Perf perf;
  perf.time = timer.totalElapsedTime() / nloop;
  perf.flops = m*n*(2+4*p);
  perf.throughput = perf.flops / perf.time / 1.0e9;

  return perf;
}

template <typename ... ViewArgs>
Perf
do_time_val(const size_t m, const size_t n, const size_t nloop)
{
  typedef Kokkos::View<double**, ViewArgs...> ViewTypeA;
  typedef Kokkos::View<double*,  ViewArgs...> ViewTypeB;
  typedef Kokkos::View<double*,  ViewArgs...> ViewTypeC;

  ViewTypeA A("A",m,n);
  ViewTypeB b("B",n);
  ViewTypeC c("c",m);

  Kokkos::deep_copy(A, 1.0);
  Kokkos::deep_copy(b, 1.0);

  // Execute the kernel once to warm up
  run_mat_vec( A, b, c );

  Teuchos::Time timer("mult", false);
  timer.start(true);
  for (size_t l=0; l<nloop; l++) {
    run_mat_vec( A, b, c );
  }
  timer.stop();

  Perf perf;
  perf.time = timer.totalElapsedTime() / nloop;
  perf.flops = m*n*2;
  perf.throughput = perf.flops / perf.time / 1.0e9;

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
do_times(const size_t m, const size_t n, const size_t p, const size_t nloop)
{
  // Run analytic
  Perf perf_analytic = do_time_analytic<ViewArgs...>(m,n,p,nloop);

  // Run value
  {
    Perf perf = do_time_val<ViewArgs...>(m,n,nloop);
    print_perf(perf, perf_analytic, "Value     ");
  }

  print_perf(perf_analytic, perf_analytic, "Analytic  ");

  // Run SFad
  if (p == SFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SFad<double,SFadSize>, ViewArgs...>(m,n,p,nloop);
    print_perf(perf, perf_analytic, "SFad      ");
  }

  // Run SLFad
  if (p <= SLFadSize) {
    Perf perf =
      do_time_fad<Sacado::Fad::SLFad<double,SLFadSize>, ViewArgs...>(m,n,p,nloop);
    print_perf(perf, perf_analytic, "SLFad     ");
  }

  // Run DFad
  {
    Perf perf =
      do_time_fad<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop);
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
do_times_layout(const size_t m, const size_t n, const size_t p,
                const size_t nloop, const LayoutType& layout,
                const std::string& device)
{
  int prec = 2;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(prec);
  std::cout << std::endl
            << device
            << " performance for layout "
            << layout_names[layout]
            << std::endl << std::endl;
  std::cout << "Computation \t Time     \t Throughput \t Ratio" << std::endl;

  if (layout == LAYOUT_LEFT)
    do_times<SFadSize,SLFadSize,Kokkos::LayoutLeft,Device>(m,n,p,nloop);
  else if (layout == LAYOUT_RIGHT)
    do_times<SFadSize,SLFadSize,Kokkos::LayoutRight,Device>(m,n,p,nloop);
  else
    do_times<SFadSize,SLFadSize,Device>(m,n,p,nloop);
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
const int SLFadSize = 8;

int main(int argc, char* argv[]) {
  bool success = true;
  try {

    // Set up command line options
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("This program tests the speed of various forward mode AD implementations for simple Kokkos kernel");
    int m = 1000000;
    clp.setOption("m", &m, "Number of matrix rows");
    int n = 10;
    clp.setOption("n", &m, "Number of matrix columns");
    int p = 8;
    clp.setOption("p", &p, "Number of derivative components");
    int nloop = 10;
    clp.setOption("nloop", &nloop, "Number of loops");
#ifdef KOKKOS_HAVE_SERIAL
    bool serial = 0;
    clp.setOption("serial", "no-serial", &serial, "Whether to run Serial");
#endif
#ifdef KOKKOS_HAVE_OPENMP
    int openmp = 0;
    clp.setOption("openmp", &openmp, "Number of OpenMP threads");
#endif
#ifdef KOKKOS_HAVE_PTHREAD
    int threads = 0;
    clp.setOption("threads", &threads, "Number of pThreads threads");
#endif
#ifdef KOKKOS_HAVE_CUDA
    bool cuda = 0;
    clp.setOption("cuda", "no-cuda", &cuda, "Whether to run CUDA");
#endif
    LayoutType layout = LAYOUT_DEFAULT;
    clp.setOption("layout", &layout, num_layout_types, layout_values,
                  layout_names, "View layout");
    bool vtune = false;
    clp.setOption("vtune", "no-vtune", &vtune,  "Profile with vtune");

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

#ifdef KOKKOS_HAVE_SERIAL
    if (serial) {
      Kokkos::Serial::initialize();
      do_times_layout<SFadSize,SLFadSize,Kokkos::Serial>(m,n,p,nloop,layout,
                                                         "Serial");
      Kokkos::Serial::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    if (openmp) {
      Kokkos::OpenMP::initialize(openmp);
      do_times_layout<SFadSize,SLFadSize,Kokkos::OpenMP>(m,n,p,nloop,layout,
                                                         "OpenMP");
      Kokkos::OpenMP::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_THREADS
    if (threads) {
      Kokkos::Threads::initialize(threads);
      do_times_layout<SFadSize,SLFadSize,Kokkos::Threads>(m,n,p,nloop,layout,
                                                          "Threads");
      Kokkos::Threads::finalize();
    }
#endif

#ifdef KOKKOS_HAVE_CUDA
    if (cuda) {
      Kokkos::HostSpace::execution_space::initialize();
      Kokkos::Cuda::initialize();
      do_times_layout<SFadSize,SLFadSize,Kokkos::Cuda>(m,n,p,nloop,layout,
                                                       "Cuda");
      Kokkos::HostSpace::execution_space::finalize();
      Kokkos::Cuda::finalize();
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return success;
}
