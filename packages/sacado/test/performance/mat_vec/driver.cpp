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

#include "mat_vec.hpp"
#include "mat_vec_hierarchical.hpp"
#include "mat_vec_hierarchical_dfad.hpp"

#include "Sacado.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// For vtune
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>

void
print_perf(const Perf& perf, const Perf& perf_base, const size_t p,
           const std::string& name)
{
  std::cout << name << "\t "
            << perf.time << "\t "
            << perf.throughput << "\t "
            << perf.time / perf_base.time
            << std::endl;
}

template <int SFadSize, int SLFadSize, int HierSFadSize, int HierSLFadSize,
          typename ... ViewArgs>
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
         const bool flat,
         const bool hierarchical,
         const bool check)
{
  Perf perf_value;
  perf_value.time = 1.0;

  // Run value
  if (value) {
    try {
      Perf perf = do_time_val<ViewArgs...>(m,n,nloop,check);
      perf_value = perf;
      print_perf(perf, perf_value, p, "Value     ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run analytic
  if (analytic) {
    try {
      Perf perf =
        do_time_analytic<ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "Analytic  ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }
  if(analytic && p == SFadSize) {
    try {
      Perf perf =
        do_time_analytic_s<SFadSize, ViewArgs...>(m,n,nloop,check);
      print_perf(perf, perf_value, p, "Analytic-s");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }
  if(analytic && p <= SLFadSize) {
    try {
      Perf perf =
        do_time_analytic_sl<SLFadSize, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "Analytic-sl");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run flat SFad
  if (flat && sfad && p == SFadSize) {
    try {
      Perf perf =
        do_time_fad<Sacado::Fad::SFad<double,SFadSize>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "SFad      ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run flat SLFad
  if (flat && slfad && p <= SLFadSize) {
    try {
      Perf perf =
        do_time_fad<Sacado::Fad::SLFad<double,SLFadSize>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "SLFad     ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run flat DFad
  if (flat && dfad) {
    try {
      Perf perf =
        do_time_fad<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "DFad      ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    try {
      Perf perf_scratch =
        do_time_scratch<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf_scratch, perf_value, p, "DFad Scratch");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run hierarchical SFad
  if (hierarchical && sfad && p == HierSFadSize) {
    try {
      Perf perf =
        do_time_fad_hierarchical<Sacado::Fad::SFad<double,HierSFadSize>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "H. SFad   ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run hierarchical SLFad
  if (hierarchical && slfad && p <= HierSLFadSize) {
    try {
      Perf perf =
        do_time_fad_hierarchical<Sacado::Fad::SLFad<double,HierSLFadSize>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "H. SLFad  ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  // Run hierarchical DFad
  if (hierarchical && dfad) {
    try {
      Perf perf =
        do_time_fad_hierarchical_dfad<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf, perf_value, p, "H. DFad   ");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    try {
      Perf perf_scratch =
        do_time_fad_hierarchical_dfad_scratch<Sacado::Fad::DFad<double>, ViewArgs...>(m,n,p,nloop,check);
      print_perf(perf_scratch, perf_value, p, "H. DFad Scratch");
    }
    catch(std::exception& e) {
      std::cout << e.what() << std::endl;
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
                const size_t nloop,
                const bool value,
                const bool analytic,
                const bool sfad,
                const bool slfad,
                const bool dfad,
                const bool flat,
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
            << " m = " << m << " n = " << n << " p = " << p
            << std::endl << std::endl;
  std::cout << "Computation \t Time     \t Throughput \t Ratio" << std::endl;

  if (layout == LAYOUT_LEFT)
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::LayoutLeft,Device>(
      m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check);
  else if (layout == LAYOUT_RIGHT)
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::LayoutRight,Device>(
      m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check);
  else
    do_times<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Device>
      (m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check);
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

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc,argv);

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
#ifdef KOKKOS_ENABLE_HIP
    bool hip = 0;
    clp.setOption("hip", "no-hip", &hip, "Whether to run HIP");
#endif
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
    bool dfad = true;
    clp.setOption("dfad", "no-dfad", &dfad, "Run DFad derivative calculation");
    bool flat = true;
    clp.setOption("flat", "no-flat", &flat, "Run flat Fad derivative calculation");
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

    if (print_config)
      Kokkos::print_configuration(std::cout, true);

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Serial>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check,layout,"Serial");
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::OpenMP>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check,layout,"OpenMP");
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Threads>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check,layout,"Threads");
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::Cuda>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check,layout,"Cuda");
    }
#endif

#ifdef KOKKOS_ENABLE_HIP
    if (hip) {
      do_times_layout<SFadSize,SLFadSize,HierSFadSize,HierSLFadSize,Kokkos::HIP>(
        m,n,p,nloop,value,analytic,sfad,slfad,dfad,flat,hierarchical,check,layout,"HIP");
    }
#endif

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  Kokkos::finalize();

  return !success;
}
