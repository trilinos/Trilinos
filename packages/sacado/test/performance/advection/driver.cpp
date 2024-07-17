// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// A performance test that computes the derivative of a simple Kokkos kernel
// using various Fad classes

#include "Sacado.hpp"

#include "advection.hpp"
#include "advection_hierarchical.hpp"
#include "advection_hierarchical_dfad.hpp"
#include "common.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

template <typename ExecSpace>
void run(const int cell_begin, const int cell_end, const int cell_step,
         const int nbasis, const int npoint, const int ntrial, const bool check)
{
  const int ndim = 3;
  printf("ncell %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "flat sfad", "flat slfad", "flat dfad", "dfad sc", "analytic", "const", "team", "hier sfad", "hier slfad", "hier dfad", "h dfad sc");
  for(int i=cell_begin; i<=cell_end; i+=cell_step) {
    double sfad_flat = time_fad_flat<SFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double slfad_flat = time_fad_flat<SLFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double dfad_flat = time_fad_flat<DFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double dfad_scratch = time_fad_scratch<DFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double analytic = time_analytic_flat<fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double analytic_const = time_analytic_const<fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double analytic_team = time_analytic_team<fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double sfad_hierarchical = time_fad_hierarchical_team<SFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double slfad_hierarchical = time_fad_hierarchical_team<SLFadType,fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double dfad_hierarchical = time_dfad_hierarchical_team<fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    double dfad_hierarchical_scratch =
      time_dfad_hierarchical_team_scratch<fad_dim,ExecSpace>(
      i,nbasis,npoint,ndim,ntrial,check);
    printf("%5d %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n",i,sfad_flat,slfad_flat,dfad_flat,dfad_scratch,analytic,analytic_const,analytic_team,sfad_hierarchical,slfad_hierarchical,dfad_hierarchical,dfad_hierarchical_scratch);
  }
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc,argv);

  bool success = true;
  try {

    // Set up command line options
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("This program tests the speed of various forward mode AD implementations for simple Kokkos kernel");
#ifdef KOKKOS_ENABLE_SERIAL
    bool serial = 0;
    clp.setOption("serial", "no-serial", &serial, "Whether to run Serial");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    bool openmp = 0;
    clp.setOption("openmp", "no-openmp", &openmp, "Whether to run OpenMP");
#endif
#ifdef KOKKOS_ENABLE_THREADS
    bool threads = 0;
    clp.setOption("threads", "no-threads", &threads, "Whether to run Threads");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = 0;
    clp.setOption("cuda", "no-cuda", &cuda, "Whether to run CUDA");
#endif
    bool print_config = false;
    clp.setOption("print-config", "no-print-config", &print_config,
                  "Whether to print Kokkos device configuration");
    int cell_begin = 100;
    clp.setOption("begin", &cell_begin, "Starting number of cells");
    int cell_end = 8000;
    clp.setOption("end", &cell_end, "Ending number of cells");
    int cell_step = 100;
    clp.setOption("step", &cell_step, "Cell increment");
    int nbasis = 8;
    clp.setOption("basis", &nbasis, "Number of basis functions");
    int npoint = 8;
    clp.setOption("point", &npoint, "Number of integration points");
    int ntrial = 5;
    clp.setOption("trial", &ntrial, "Number of trials");
    bool check = false;
    clp.setOption("check", "no-check", &check,
                  "Check correctness of results");

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

    if (print_config)
      Kokkos::print_configuration(std::cout, true);

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      using Kokkos::Serial;
      run<Serial>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      using Kokkos::OpenMP;
      run<OpenMP>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      using Kokkos::Threads;
      run<Threads>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      using Kokkos::Cuda;
      run<Cuda>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  Kokkos::finalize();

  return !success;
}
