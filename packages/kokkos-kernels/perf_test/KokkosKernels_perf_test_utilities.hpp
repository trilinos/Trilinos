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

  int repeat = 0;
};

std::string list_common_options() {
  std::ostringstream common_options;
  common_options
      << "\t[Required] BACKEND:\n"
      << "\t\t'--threads [numThreads]' |\n"
      << "\t\t'--openmp [numThreads]' |\n"
      << "\t\t'--cuda [deviceIndex]' |\n"
      << "\t\t'--hip [deviceIndex]' |\n"
      << "\t\t'--sycl [deviceIndex]'\n\n"
      << "\tIf no parallel backend is requested, Serial will be used "
         "(if enabled)\n\n";

  return common_options.str();
}

void process_arg_int(char const* str_val, int& val) {
  errno = 0;
  char* ptr_end;
  val = std::strtol(str_val, &ptr_end, 10);

  if (str_val == ptr_end) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << str_val
       << "' to an integer.\n";
    throw std::invalid_argument(ss.str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for command line argument '" << str_val
       << "' falls out of range.\n";
    throw std::invalid_argument(ss.str());
  }
}

void process_arg_double(char const* str_val, double& val) {
  errno = 0;
  char* ptr_end;
  val = std::strtod(str_val, &ptr_end);

  if (str_val == ptr_end) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << str_val
       << "' to a double.\n";
    throw std::invalid_argument(ss.str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for command line argument '" << str_val
       << "' falls out of range.\n";
    throw std::invalid_argument(ss.str());
  }
}

bool check_arg_int(int const i, int const argc, char** argv, char const* name,
                   int& val) {
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

bool check_arg_double(int const i, int const argc, char** argv,
                      char const* name, double& val) {
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

bool check_arg_bool(int const i, int const /*argc*/, char** argv,
                    char const* name, bool& val) {
  if (0 != Test::string_compare_no_case(argv[i], name)) {
    return false;
  }
  val = true;
  return true;
}

bool check_arg_str(int const i, int const argc, char** argv, char const* name,
                   std::string& val) {
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
    bool remove_flag = false;
    if (check_arg_int(argIdx, argc, argv, "--threads", params.use_threads)) {
      remove_flag = true;
    } else if (check_arg_int(argIdx, argc, argv, "--openmp",
                             params.use_openmp)) {
      remove_flag = true;
    } else if (check_arg_int(argIdx, argc, argv, "--cuda", params.use_cuda)) {
      params.use_cuda++;
      remove_flag = true;
    } else if (check_arg_int(argIdx, argc, argv, "--hip", params.use_hip)) {
      params.use_hip++;
      remove_flag = true;
    } else if (check_arg_int(argIdx, argc, argv, "--sycl", params.use_sycl)) {
      params.use_sycl++;
      remove_flag = true;
    } else if (check_arg_int(argIdx, argc, argv, "--repeat", params.repeat)) {
      remove_flag = true;
    }

    if (remove_flag) {
      // Shift the remainder of the argv list by one.  Note that argv has
      // (argc + 1) arguments, the last one always being nullptr.  The following
      // loop moves the trailing nullptr element as well
      for (int k = argIdx; k < argc - 1; ++k) {
        argv[k]     = argv[k + 2];
        argv[k + 1] = argv[k + 3];
      }
      argc = argc - 2;
    } else {
      ++argIdx;
    }
  }
}  // parse_common_options()

}  // namespace perf_test

#endif  // KOKKOSKERNELS_PERF_TEST_UTILITIES_HPP
