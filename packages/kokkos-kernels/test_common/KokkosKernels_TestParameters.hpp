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

#ifndef KK_TESTPARAMS_H
#define KK_TESTPARAMS_H

namespace KokkosKernels {

namespace Experiment {

struct Parameters {
  int algorithm;
  int accumulator;
  int repeat;
  int chunk_size;
  int multi_color_scale;
  int shmemsize;
  int team_size;
  int use_dynamic_scheduling;
  int verbose;
  int spgemm_step;
  int vector_size;
  int check_output;
  int mkl_sort_option;
  int mkl_keep_output;
  int calculate_read_write_cost;
  char *coloring_input_file;
  char *coloring_output_file;

  int minhashscale;
  int use_threads;
  int use_openmp;
  int use_cuda;
  int use_hip;
  int use_serial;
  int a_mem_space, b_mem_space, c_mem_space, work_mem_space;

  char *a_mtx_bin_file, *b_mtx_bin_file, *c_mtx_bin_file;
  bool compression2step;
  int left_lower_triangle, right_lower_triangle;
  int left_sort, right_sort;

  int triangle_options;
  bool apply_compression;
  int sort_option;
  // 0 - triangle_count
  // 1 - first count then instantiate
  // 2- more options.
  int cache_flush;
  double first_level_hash_cut_off;
  double compression_cut_off;
  size_t MaxColDenseAcc;
  // 0 - no flush
  // 1 - soft flush
  // 2 - hard flush with rand.
  Parameters() {
    algorithm                 = 0;
    accumulator               = 0;
    repeat                    = 6;
    chunk_size                = -1;
    multi_color_scale         = 1;
    shmemsize                 = 16128;
    team_size                 = -1;
    use_dynamic_scheduling    = 0;
    verbose                   = 0;
    spgemm_step               = '0';
    vector_size               = -1;
    check_output              = 0;
    mkl_sort_option           = 7;
    mkl_keep_output           = 1;
    calculate_read_write_cost = 0;
    coloring_input_file       = NULL;
    coloring_output_file      = NULL;
    minhashscale              = 1;
    use_threads               = 0;
    use_openmp                = 0;
    use_cuda                  = 0;
    use_hip                   = 0;
    use_serial                = 0;
    a_mem_space = b_mem_space = c_mem_space = work_mem_space = 1;
    a_mtx_bin_file = b_mtx_bin_file = c_mtx_bin_file = NULL;
    compression2step                                 = true;

    left_lower_triangle  = 0;
    right_lower_triangle = 0;
    left_sort            = 0;
    right_sort           = 2;  // algorithm decides
    triangle_options     = 0;
    apply_compression    = true;
    sort_option          = -1;
    cache_flush          = 1;

    first_level_hash_cut_off = 0.50;
    compression_cut_off      = 0.85;
    MaxColDenseAcc           = 250000;
  }
};
}  // namespace Experiment
}  // namespace KokkosKernels

#endif
