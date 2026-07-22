// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <cstdlib>

#include "Kokkos_Core.hpp"

#include "Stokhos_ConfigDefs.h"

template <typename scalar, typename device>
int mainHost(bool test_flat, bool test_orig, bool test_deg, bool test_lin,
             bool test_block, bool symmetric, bool mkl);

template <typename scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_lin,
             bool test_block, bool symmetric, int device_id);

int main(int argc, char *argv[])
{
  // Defaults
  bool test_host = true;
#ifdef KOKKOS_ENABLE_CUDA
  bool test_cuda = true;
  int device = 0;
#endif
  bool test_block = true;
  bool test_flat = true;
  bool test_orig = true;
  bool test_deg = false;
  bool test_lin = false;
  bool symmetric = true;
  bool single = false;
  bool mkl = false;
#ifdef KOKKOS_ENABLE_SERIAL
  bool serial = true;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  bool omp = true;
#endif
#ifdef KOKKOS_ENABLE_THREADS
  bool threads = true;
#endif

  // Parse command line arguments
  bool print_usage = false;
  int i=1;
  while (i<argc) {
    std::string s(argv[i]);
    if (s == "host")
      test_host = true;
    else if (s == "no-host")
      test_host = false;
#ifdef KOKKOS_ENABLE_CUDA
    else if (s == "cuda")
      test_cuda = true;
    else if (s == "no-cuda")
      test_cuda = false;
    else if (s == "device") {
      ++i;
      device = std::atoi(argv[i]);
    }
#endif
    else if (s == "block")
      test_block = true;
    else if (s == "no-block")
      test_block = false;
    else if (s == "flat")
      test_flat = true;
    else if (s == "no-flat")
      test_flat = false;
    else if (s == "orig")
      test_orig = true;
    else if (s == "no-orig")
      test_orig = false;
    else if (s == "deg")
      test_deg = true;
    else if (s == "no-deg")
      test_deg = false;
    else if (s == "linear")
      test_lin = true;
    else if (s == "no-linear")
      test_lin = false;
    else if (s == "symmetric")
      symmetric = true;
    else if (s == "no-symmetric")
      symmetric = false;
    else if (s == "mkl")
      mkl = true;
    else if (s == "no-mkl")
      mkl = false;
    else if (s == "single")
      single = true;
    else if (s == "double")
      single = false;
#ifdef KOKKOS_ENABLE_SERIAL
    else if (s == "serial")
      serial = true;
    else if (s == "no-serial")
      serial = false;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    else if (s == "omp")
      omp = true;
    else if (s == "no-omp")
      omp = false;
#endif
#ifdef KOKKOS_ENABLE_THREADS
    else if (s == "threads")
      threads = true;
    else if (s == "no-threads")
      threads = false;
#endif
    else if (s == "-h" || s == "--help")
      print_usage = true;
    else {
      std::cout << "Invalid argument:  " << s << std::endl;
      print_usage = true;
    }
    ++i;
  }
  if (print_usage) {
    std::cout << "Usage:" << std::endl
              << "\t" << argv[0]
              << " [no-][cuda|host|serial|omp|threads|block|flat|orig|deg|linear|symmetric] [single|double] [device device_id]"
              << std::endl << "Defaults are all enabled." << std::endl;
    return -1;
  }

  if (test_host) {

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      if (single)
        mainHost<float,Kokkos::Serial>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
      else
        mainHost<double,Kokkos::Serial>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      if (single)
        mainHost<float,Kokkos::Threads>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
      else
        mainHost<double,Kokkos::Threads>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (omp) {
      if (single)
        mainHost<float,Kokkos::OpenMP>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
      else
        mainHost<double,Kokkos::OpenMP>(
          test_flat, test_orig, test_deg, test_lin, test_block, symmetric, mkl);
    }
#endif

  }

#ifdef KOKKOS_ENABLE_CUDA
  if (test_cuda) {
    if (single)
      mainCuda<float>(test_flat, test_orig, test_lin, test_block, symmetric, device);
    else
      mainCuda<double>(test_flat, test_orig, test_lin, test_block, symmetric, device);
  }
#endif

  return 0 ;
}
