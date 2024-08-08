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
#ifndef KOKKOSBLAS3_COMMON_H_
#define KOKKOSBLAS3_COMMON_H_
#include "KokkosKernels_default_types.hpp"

#include <iostream>
#include <fstream>

/************************ perf test default value definitions **********/
#define DEFAULT_TEST BLAS
#define DEFAULT_LOOP SERIAL
#define DEFAULT_MATRIX_START 10
#define DEFAULT_MATRIX_STOP 2430
#define DEFAULT_STEP 3
#define DEFAULT_WARM_UP_N 100
#define DEFAULT_N 100
#define DEFAULT_K 1024
#define DEFAULT_OUT &std::cout
#define DEFAULT_BLAS_ROUTINES "trmm,gemm,"
#define DEFAULT_TEAM_SIZE 1
#define DEFAULT_VECTOR_LEN 1
#define DEFAULT_USE_AUTO 0
#define DEFAULT_BATCH_SIZE_LAST_DIM 0
#define DEFAULT_VERIFY 1
#define DEFAULT_NINTER 4
#define DEFAULT_USE_SIMD 0

/************************ blas routine structure definitions **********/
struct perf_test_trmm_args {
  std::string trmm_args;
  default_scalar alpha;
};
typedef struct perf_test_trmm_args pt_trmm_args_t;

struct perf_test_gemm_args {
  std::string gemm_args;  //[N,T,C][N,T,C] for transA and transB
  default_scalar alpha;
  default_scalar beta;
};
typedef struct perf_test_gemm_args pt_gemm_args_t;
// ADD MORE BLAS3 ROUTINE ARG STRUCTS HERE.

struct blas_args {
  pt_trmm_args_t trmm;
  pt_gemm_args_t gemm;
  // ADD MORE BLAS3 ROUTINES HERE
  int team_size;
  int vector_len;
  bool use_auto, batch_size_last_dim;
  // ADD MORE COMMON BLAS3 OPTIONS HERE
};
typedef struct blas_args blas_args_t;

typedef enum BLAS_ROUTINES {
  TRMM,
  GEMM,
  // ADD MORE BLAS3 ROUTINES HERE
  BLAS_ROUTINES_N
} blas_routines_e;

static std::string blas_routines_e_str[BLAS_ROUTINES_N] = {
    "trmm", "gemm"
    // ADD MORE BLAS3 ROUTINES HERE
};

/************************ perf test type definitions ************************/
/**
 * @var SERIAL:   Run the blas routine iteratively, within a for-loop
 * @var PARALLEL: Run the blas routine iteratively, within a
 * Kokkos::parallel_for-loop
 */
typedef enum LOOP {
  SERIAL,
  PARALLEL,
  // ADD MORE LOOP TYPES HERE
  LOOP_N
} loop_e;

static std::string loop_e_str[LOOP_N] = {"serial", "parallel"};

/**
 * @var BLAS:                          Run the blas routine through the
 *                                     KokkosBlas namespace.
 * @var BATCHED_SERIAL{_BLOCKED}:      Run the serial blas routine through the
 *                                     KokkosBatched namespace.
 * @var BATCHED_SERIAL_SIMD{_BLOCKED}: Run the serial blas routine through the
 *                                     KokkosBatched namespace using SIMD views.
 * @var BATCHED_SERIAL_COMPACT_MKL:    Run the serial blas mkl routine through
 *                                     the KokkosBatched namespace.
 * @var BATCHED_TEAM{_BLOCKED}:        Run the team blas routine through the
 *                                     KokkosBatched namespace.
 * @var BATCHED_TEAM_VECTOR{_BLOCKED}: Run the team vector blas routine through
 *                                     the KokkosBatched namespace.
 * @var BATCHED_TEAM_SIMD{_BLOCKED}:   Run the team vector blas routine through
 * the KokkosBatched namespace using SIMD views.
 * @var EXPERIMENT:                    Run the blas routine as a custom
 * experiment.
 */
typedef enum TEST {
  BLAS,
  BATCHED_HEURISTIC,
  BATCHED_SERIAL,
  BATCHED_SERIAL_BLOCKED,
  BATCHED_SERIAL_SIMD,
  BATCHED_SERIAL_SIMD_BLOCKED,
  BATCHED_SERIAL_COMPACT_MKL,
  BATCHED_TEAM,
  BATCHED_TEAM_BLOCKED,
  BATCHED_TEAM_VECTOR,
  BATCHED_TEAM_VECTOR_BLOCKED,
  BATCHED_TEAM_SIMD,
  BATCHED_TEAM_SIMD_BLOCKED,
  // ADD MORE TEST TYPES HERE
  EXPERIMENT,
  TEST_N
} test_e;

static std::string test_e_str[TEST_N]{
    "blas", "batched_heuristic", "batched_serial", "batched_serial_blocked", "batched_serial_simd",
    "batched_serial_simd_blocked", "batched_serial_compact_mkl", "batched_team", "batched_team_blocked",
    "batched_team_vector", "batched_team_vector_blocked", "batched_team_simd", "batched_team_simd_blocked",
    // ADD MORE TEST TYPES HERE
    "experiment"};

/**
 * @var k: Number of 2D matrices.
 * @var m: Number of rows.
 * @var n: Number of columns.
 */
struct matrix_dim {
  int k, m, n;
};
typedef struct matrix_dim matrix_dim_t;

struct matrix_dims {
  matrix_dim_t a, b, c;
};
typedef struct matrix_dims matrix_dims_t;

/**
 * @var test:          Selects which namespace to test.
 * @var loop:          Selects how to invoke the blas routine.
 * @var start:         Selects which matrix dimensions to start with.
 * @var stop:          Selects which matrix dimensions to end with.
 * @var step:          Selects the multiplier to increase start by.
 * @var warm_up_n:     Selects how many untimed runs to perform for each matrix
 * dimension.
 * @var n:             Selects how many timed runs to perform for each matrix
 * dimension.
 * @var out:           Selects where to write the csv data for each matrix
 * dimension.
 * @var out_file:      The file to write csv data to. Defaults to stdout.
 * @var blas_args:     Arguments for each supported blas routine.
 * @var blas_routines: Selects which supported blas routines to test.
 * @var verify:        Performs verification of the blas routine for each input
 *                     before timing it.
 * @var ninter:        The number of interleaved matrices for armpl.
 */
struct perf_test_options {
  test_e test;
  loop_e loop;
  matrix_dims_t start;
  matrix_dims_t stop;
  uint32_t step;
  uint32_t warm_up_n;
  uint32_t n;
  std::ostream* out;
  std::string out_file;
  blas_args_t blas_args;
  std::string blas_routines;
  bool verify;
  int ninter;
  bool use_simd;
};
typedef struct perf_test_options options_t;

/*************************** Print macros **************************/
// #define PERF_TEST_DEBUG
#ifdef PERF_TEST_DEBUG
#define STATUS printf("STATUS: %s:%d.\n", __func__, __LINE__);
#else
#define STATUS
#endif  // PERF_TEST_DEBUG
#define FATAL_ERROR(msg) printf("FATAL_ERROR: %s:%s:%d %s\n", __FILE__, __func__, __LINE__, (msg));
#endif  // KOKKOSBLAS3_COMMON_H_
