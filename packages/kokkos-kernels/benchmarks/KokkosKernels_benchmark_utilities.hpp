// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
//
// Created by Berger-Vergiat, Luc on 2/6/23.
//

#ifndef KOKKOSKERNELS_BENCHMARK_UTILITIES_HPP
#define KOKKOSKERNELS_BENCHMARK_UTILITIES_HPP

#include "KokkosKernels_TestStringUtils.hpp"  // for string_compare_no_case

// Namepsace that defines common utilities
// for benchmarks
namespace benchmark {

struct CommonInputParams {
  bool use_cuda    = false;
  bool use_hip     = false;
  bool use_sycl    = false;
  bool use_openmp  = false;
  bool use_threads = false;
  bool use_serial  = false;

  int repeat      = 0;
  bool print_help = false;
};

std::string list_common_options() {
  std::ostringstream common_options;
  common_options << "\t[Required]:\n"
                 << "\t[Optional]:\n"
                 << "\t The following backends are available because Kokkos was "
                    "configured with them:\n"
#ifdef KOKKOS_ENABLE_THREADS
                 << "\t\t'--threads'\n"
#endif
#ifdef KOKKOS_ENABLE_OPENMP
                 << "\t\t'--openmp'\n"
#endif
#ifdef KOKKOS_ENABLE_CUDA
                 << "\t\t'--cuda'\n"
#endif
#ifdef KOKKOS_ENABLE_HIP
                 << "\t\t'--hip'\n"
#endif
#ifdef KOKKOS_ENABLE_SYCL
                 << "\t\t'--sycl'\n"
#endif
#ifdef KOKKOS_ENABLE_SERIAL
                 << "\t\t'--serial'\n"
#endif
                 << "\n"
                 << "\t The following backends are not available because Kokkos was not "
                    "configured with them:\n"
#ifndef KOKKOS_ENABLE_THREADS
                 << "\t\t'--threads'\n"
#endif
#ifndef KOKKOS_ENABLE_OPENMP
                 << "\t\t'--openmp'\n"
#endif
#ifndef KOKKOS_ENABLE_CUDA
                 << "\t\t'--cuda'\n"
#endif
#ifndef KOKKOS_ENABLE_HIP
                 << "\t\t'--hip'\n"
#endif
#ifndef KOKKOS_ENABLE_SYCL
                 << "\t\t'--sycl'\n"
#endif
#ifndef KOKKOS_ENABLE_SERIAL
                 << "\t\t'--serial'\n"
#endif
                 << "\n"
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
  while (argIdx < argc) {
    // How many flags to delete from argc/argv
    // 0: not a common option, so leave it
    // 1: a bool parameter like '-h or --cuda'
    // 2: a parameter followed by a value, like "--repeat 3"
    int remove_flags = 0;
    if (check_arg_bool(argIdx, argc, argv, "--threads", params.use_threads)) {
      remove_flags = 1;
    } else if (check_arg_bool(argIdx, argc, argv, "--openmp", params.use_openmp)) {
      remove_flags = 1;
    } else if (check_arg_bool(argIdx, argc, argv, "--cuda", params.use_cuda)) {
      remove_flags = 1;
    } else if (check_arg_bool(argIdx, argc, argv, "--hip", params.use_hip)) {
      remove_flags = 1;
    } else if (check_arg_bool(argIdx, argc, argv, "--sycl", params.use_sycl)) {
      remove_flags = 1;
    } else if (check_arg_bool(argIdx, argc, argv, "--serial", params.use_serial)) {
      remove_flags = 1;
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

}  // namespace benchmark

#endif  // KOKKOSKERNELS_BENCHMARK_UTILITIES_HPP
