//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
//
// Created by Berger-Vergiat, Luc on 2/6/23.
//

#ifndef KOKKOSKERNELS_PERF_TEST_UTILITIES_HPP
#define KOKKOSKERNELS_PERF_TEST_UTILITIES_HPP

#include "KokkosKernels_TestUtils.hpp"  // for string_compare_no_case

// Namepsace that defines common utilities
// for performance tests
namespace perf_test {

struct CommonInputParams {
  int use_cuda    = 0;
  int use_hip     = 0;
  int use_sycl    = 0;
  int use_openmp  = 0;
  int use_threads = 0;

  int repeat      = 0;
  bool print_help = false;
};

std::string list_common_options() {
  std::ostringstream common_options;
  common_options << "\t[Required] Backend: the available backends are:\n"
#ifdef KOKKOS_ENABLE_THREADS
                 << "\t\t'--threads [numThreads]'\n"
#endif
#ifdef KOKKOS_ENABLE_OPENMP
                 << "\t\t'--openmp [numThreads]'\n"
#endif
#ifdef KOKKOS_ENABLE_CUDA
                 << "\t\t'--cuda [deviceIndex]'\n"
#endif
#ifdef KOKKOS_ENABLE_HIP
                 << "\t\t'--hip [deviceIndex]'\n"
#endif
#ifdef KOKKOS_ENABLE_SYCL
                 << "\t\t'--sycl [deviceIndex]'\n"
#endif
#ifdef KOKKOS_ENABLE_SERIAL
                 << "\t\tIf no parallel backend is requested, Serial will be used.\n"
#endif
                 << "\n"
                 << "\t The following backends are not available because Kokkos was not "
                    "configured with them:\n"
#ifndef KOKKOS_ENABLE_THREADS
                 << "\t\t'--threads [numThreads]'\n"
#endif
#ifndef KOKKOS_ENABLE_OPENMP
                 << "\t\t'--openmp [numThreads]'\n"
#endif
#ifndef KOKKOS_ENABLE_CUDA
                 << "\t\t'--cuda [deviceIndex]'\n"
#endif
#ifndef KOKKOS_ENABLE_HIP
                 << "\t\t'--hip [deviceIndex]'\n"
#endif
#ifndef KOKKOS_ENABLE_SYCL
                 << "\t\t'--sycl [deviceIndex]'\n"
#endif
#ifndef KOKKOS_ENABLE_SERIAL
                 << "\t\tSerial is not enabled so a parallel backend must be selected.\n"
#endif
                 << "\n"
                 << "\t[Optional]:\n"
                 << "\t\t'-h', '--help': show available options\n\n";

  return common_options.str();
}

void process_arg_int(char const* str_val, int& val) {
  errno = 0;
  char* ptr_end;
  val = std::strtol(str_val, &ptr_end, 10);

  if (str_val == ptr_end) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << str_val << "' to an integer.\n";
    throw std::invalid_argument(ss.str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for command line argument '" << str_val << "' falls out of range.\n";
    throw std::invalid_argument(ss.str());
  }
}

void process_arg_double(char const* str_val, double& val) {
  errno = 0;
  char* ptr_end;
  val = std::strtod(str_val, &ptr_end);

  if (str_val == ptr_end) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << str_val << "' to a double.\n";
    throw std::invalid_argument(ss.str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for command line argument '" << str_val << "' falls out of range.\n";
    throw std::invalid_argument(ss.str());
  }
}

bool check_arg_int(int const i, int const argc, char** argv, char const* name, int& val) {
  if (0 != Test::string_compare_no_case(argv[i], name)) {
    return false;
  }

  if (i < argc - 1) {
    process_arg_int(argv[i + 1], val);
  } else {
    std::stringstream msg;
    msg << name << " input argument needs to be followed by an int";
    throw std::invalid_argument(msg.str());
  }
  return true;
}

bool check_arg_double(int const i, int const argc, char** argv, char const* name, double& val) {
  if (0 != Test::string_compare_no_case(argv[i], name)) {
    return false;
  }

  if (i < argc - 1) {
    process_arg_double(argv[i + 1], val);
  } else {
    std::stringstream msg;
    msg << name << " input argument needs to be followed by a real number";
    throw std::invalid_argument(msg.str());
  }
  return true;
}

bool check_arg_bool(int const i, int const /*argc*/, char** argv, char const* name, bool& val) {
  if (0 != Test::string_compare_no_case(argv[i], name)) {
    return false;
  }
  val = true;
  return true;
}

bool check_arg_str(int const i, int const argc, char** argv, char const* name, std::string& val) {
  if (0 != Test::string_compare_no_case(argv[i], name)) {
    return false;
  }

  if (i < argc - 1) {
    val = std::string(argv[i + 1]);
  } else {
    std::stringstream msg;
    msg << name << " input argument needs to be followed by a string";
    throw std::invalid_argument(msg.str());
  }
  return true;
}

void parse_common_options(int& argc, char** argv, CommonInputParams& params) {
  // Skip the program name, start with argIdx=1
  int argIdx = 1;
  // Note: after parsing a GPU device ID, always add 1 to it.
  // If e.g. params.use_cuda is 0, that means CUDA will not be used at all.
  // But if it's N, then it means run on CUDA device N-1.
  while (argIdx < argc) {
    // How many flags to delete from argc/argv
    // 0: not a common option, so leave it
    // 1: a bool parameter like '-h'
    // 2: a parameter followed by a value, like "--cuda 0"
    int remove_flags = 0;
    if (check_arg_int(argIdx, argc, argv, "--threads", params.use_threads)) {
      remove_flags = 2;
    } else if (check_arg_int(argIdx, argc, argv, "--openmp", params.use_openmp)) {
      remove_flags = 2;
    } else if (check_arg_int(argIdx, argc, argv, "--cuda", params.use_cuda)) {
      params.use_cuda++;
      remove_flags = 2;
    } else if (check_arg_int(argIdx, argc, argv, "--hip", params.use_hip)) {
      params.use_hip++;
      remove_flags = 2;
    } else if (check_arg_int(argIdx, argc, argv, "--sycl", params.use_sycl)) {
      params.use_sycl++;
      remove_flags = 2;
    } else if (check_arg_int(argIdx, argc, argv, "--repeat", params.repeat)) {
      remove_flags = 2;
    } else if (check_arg_bool(argIdx, argc, argv, "-h", params.print_help) ||
               check_arg_bool(argIdx, argc, argv, "--help", params.print_help)) {
      remove_flags = 1;
    }

    if (remove_flags) {
      // Shift the remainder of the argv list left by the number of flags
      // removed. Note that argv has (argc + 1) arguments, the last one always
      // being nullptr.  The following loop moves the trailing nullptr element
      // as well
      for (int k = argIdx + remove_flags; k <= argc; ++k) {
        argv[k - remove_flags] = argv[k];
      }
      argc -= remove_flags;
    } else {
      ++argIdx;
    }
  }
}  // parse_common_options()

}  // namespace perf_test

#endif  // KOKKOSKERNELS_PERF_TEST_UTILITIES_HPP
