// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "fe_jac_fill_funcs.hpp"

std::vector<double>
do_times(int work_count, int num_eqns_begin, int num_eqns_end,
         int num_eqns_delta,
         double (*func)(unsigned int,unsigned int,double)) {
  std::vector<double> times;
  for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
       num_eqns += num_eqns_delta) {
    int num_nodes = work_count / num_eqns;
    double mesh_spacing = 1.0 / (num_nodes - 1);
    times.push_back(func(num_nodes, num_eqns, mesh_spacing));
  }
  return times;
}

template <template <typename> class FadType>
std::vector<double>
do_times_fad(int work_count, int num_eqns_begin, int num_eqns_end,
             int num_eqns_delta) {
  std::vector<double> times;
  for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
       num_eqns += num_eqns_delta) {
    int num_nodes = work_count / num_eqns;
    double mesh_spacing = 1.0 / (num_nodes - 1);
    times.push_back(fad_jac_fill<FadType<double> >(num_nodes, num_eqns, mesh_spacing));
  }
  return times;
}

template <template <typename,int> class FadType>
std::vector<double>
do_times_slfad(int work_count, int num_eqns_begin, int num_eqns_end,
             int num_eqns_delta) {
  const int slfad_max = 130; // Maximum number of derivative components for SLFad
  std::vector<double> times;
  for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
       num_eqns += num_eqns_delta) {
    int num_nodes = work_count / num_eqns;
    double mesh_spacing = 1.0 / (num_nodes - 1);
    if (num_eqns*2 < slfad_max)
      times.push_back(fad_jac_fill<FadType<double,slfad_max> >(num_nodes, num_eqns, mesh_spacing));
  }
  return times;
}

template <template <typename,int> class FadType>
std::vector<double>
do_times_sfad(int work_count, int num_eqns_begin, int num_eqns_end,
             int num_eqns_delta) {
  std::vector<double> times;
  for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
       num_eqns += num_eqns_delta) {
    int num_nodes = work_count / num_eqns;
    double mesh_spacing = 1.0 / (num_nodes - 1);
    if (num_eqns*2 == 10)
      times.push_back(fad_jac_fill<FadType<double,10> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 30)
      times.push_back(fad_jac_fill<FadType<double,30> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 50)
      times.push_back(fad_jac_fill<FadType<double,50> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 70)
      times.push_back(fad_jac_fill<FadType<double,70> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 90)
      times.push_back(fad_jac_fill<FadType<double,90> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 110)
      times.push_back(fad_jac_fill<FadType<double,110> >(num_nodes, num_eqns, mesh_spacing));
    else if (num_eqns*2 == 130)
      times.push_back(fad_jac_fill<FadType<double,130> >(num_nodes, num_eqns, mesh_spacing));
  }
  return times;
}

void print_times(const std::vector<double>& times,
                 const std::vector<double>& base,
                 const std::string& name, int p, int w, int w_name) {
  std::cout.setf(std::ios::scientific);
  std::cout.precision(p);
  std::cout.setf(std::ios::right);
  std::cout << std::setw(w_name) << name << " ";
  std::cout.setf(std::ios::right);
  for (unsigned int i=0; i<times.size(); i++)
    std::cout << std::setw(w) << times[i]/base[i] << " ";
  std::cout << std::endl;
}

int main(int argc, char* argv[]) {
  int ierr = 0;
  int p = 1;
  int w = p+7;
  int w_name = 13;

  try {

    // Set up command line options
    Teuchos::CommandLineProcessor clp;
    clp.setDocString("This program tests the speed of various forward mode AD implementations for a finite-element-like Jacobian fill");
    int work_count = 200000;
    int num_eqns_begin = 5;
    int num_eqns_end = 65;
    int num_eqns_delta = 10;
    int rt = 0;
    clp.setOption("wc", &work_count, "Work count = num_nodes*num_eqns");
    clp.setOption("p_begin", &num_eqns_begin, "Intitial number of equations");
    clp.setOption("p_end", &num_eqns_end, "Final number of equations");
    clp.setOption("p_delta", &num_eqns_delta, "Step in number of equations");
    clp.setOption("rt", &rt, "Include ADOL-C retaping test");

    // Parse options
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn= clp.parse(argc, argv);
    if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    // Print header
    std::cout.setf(std::ios::right);
    std::cout << std::setw(w_name) << "Name" << " ";
    for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
         num_eqns += num_eqns_delta)
      std::cout << std::setw(w) << num_eqns << " ";
    std::cout << std::endl;
    for (int j=0; j<w_name; j++)
      std::cout << '=';
    std::cout << " ";
    for (int num_eqns = num_eqns_begin; num_eqns <= num_eqns_end;
         num_eqns += num_eqns_delta) {
      for (int j=0; j<w; j++)
        std::cout << '=';
      std::cout << " ";
    }
    std::cout << std::endl;

    // Analytic
    std::vector<double> times_analytic =
      do_times(work_count, num_eqns_begin, num_eqns_end, num_eqns_delta,
               analytic_jac_fill);
    print_times(times_analytic, times_analytic, "Analytic", p, w, w_name);

#ifdef HAVE_ADIC
    // Note there seems to be a bug in ADIC where doing more than one num_eqns
    // value results in incorrect timings after the first.  Doing one value
    // at a time seems to give correct values though.
    std::vector<double> times_adic =
      do_times(work_count, num_eqns_begin, num_eqns_end, num_eqns_delta,
               adic_jac_fill);
    print_times(times_adic, times_analytic, "ADIC", p, w, w_name);
#endif

    // Original Fad
    std::vector<double> times_sfad =
      do_times_sfad<Sacado::Fad::SFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_sfad, times_analytic, "SFAD", p, w, w_name);

    std::vector<double> times_slfad =
      do_times_sfad<Sacado::Fad::SLFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_slfad, times_analytic, "SLFAD", p, w, w_name);

    std::vector<double> times_dfad =
      do_times_fad<Sacado::Fad::DFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_dfad, times_analytic, "DFAD", p, w, w_name);


    // ELR Fad
    std::vector<double> times_elr_sfad =
      do_times_sfad<Sacado::ELRFad::SFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_elr_sfad, times_analytic, "ELRSFAD", p, w, w_name);

    std::vector<double> times_elr_slfad =
      do_times_sfad<Sacado::ELRFad::SLFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_elr_slfad, times_analytic, "ELRSLFAD", p, w, w_name);

    std::vector<double> times_elr_dfad =
      do_times_fad<Sacado::ELRFad::DFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_elr_dfad, times_analytic, "ELRDFAD", p, w, w_name);


    // Cache Fad
    std::vector<double> times_cache_sfad =
      do_times_sfad<Sacado::CacheFad::SFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_sfad, times_analytic, "CacheSFAD", p, w, w_name);

    std::vector<double> times_cache_slfad =
      do_times_sfad<Sacado::CacheFad::SLFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_slfad, times_analytic, "CacheSLFAD", p, w, w_name);

    std::vector<double> times_cache_dfad =
      do_times_fad<Sacado::CacheFad::DFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_dfad, times_analytic, "CacheDFAD", p, w, w_name);

    // ELR Cache Fad
    std::vector<double> times_cache_elr_sfad =
      do_times_sfad<Sacado::ELRCacheFad::SFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_elr_sfad, times_analytic, "ELRCacheSFAD", p, w, w_name);

    std::vector<double> times_cache_elr_slfad =
      do_times_sfad<Sacado::ELRCacheFad::SLFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_elr_slfad, times_analytic, "ELRCacheSLFAD", p, w, w_name);

    std::vector<double> times_cache_elr_dfad =
      do_times_fad<Sacado::ELRCacheFad::DFad>(
        work_count, num_eqns_begin, num_eqns_end, num_eqns_delta);
    print_times(times_cache_elr_dfad, times_analytic, "ELRCacheDFAD", p, w, w_name);

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  return ierr;
}
