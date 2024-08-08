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
#include "KokkosBlas3_common.hpp"
#include "KokkosBlas3_trmm_perf_test.hpp"
#include "KokkosBlas3_gemm_perf_test.hpp"

#include <cstdlib>
#include <memory>
#include <unistd.h>
#include <getopt.h>

static struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                       {"test", required_argument, 0, 't'},
                                       {"warm_up_loop", required_argument, 0, 'w'},
                                       {"trmm_options", required_argument, 0, 'o'},
                                       {"trmm_alpha", required_argument, 0, 'a'},
                                       {"gemm_options", required_argument, 0, 'g'},
                                       {"gemm_scalars", required_argument, 0, 'p'},
                                       {"team_size", required_argument, 0, 'z'},
                                       {"vector_len", required_argument, 0, 'n'},
                                       {"use_auto", required_argument, 0, 'u'},
                                       {"batch_size", required_argument, 0, 'k'},
                                       {"batch_size_last_dim", required_argument, 0, 'd'},
                                       {"loop_type", required_argument, 0, 'l'},
                                       {"matrix_size_start", required_argument, 0, 'b'},
                                       {"matrix_size_stop", required_argument, 0, 'e'},
                                       {"matrix_size_step", required_argument, 0, 's'},
                                       {"iter", required_argument, 0, 'i'},
                                       {"csv", required_argument, 0, 'c'},
                                       {"routines", required_argument, 0, 'r'},
                                       {"verify", required_argument, 0, 'v'},
                                       {"ninter", required_argument, 0, 'j'},
                                       {"use_simd", required_argument, 0, 'f'},
                                       {0, 0, 0, 0}};

static void __print_help_blas3_perf_test() {
  printf("Options:\n");

  printf("\t-h, --help\n");
  printf("\t\tPrint this help menu.\n");

  printf("\t-t, --test=OPTION\n");
  printf("\t\tAlgorithm selection.\n");
  printf("\t\t\tValid values for OPTION:\n");
  for (int i = 0; i < TEST_N; i++) {
    printf("%c[1m", 27);
    printf("\t\t\t\t%s", test_e_str[i].c_str());
    printf("%c[0m", 27);
    printf("\n");
  }

  printf("\t-o, --trmm_options=OPTION_STRING\n");
  printf("\t\tTRMM side, uplo, trans, and diag options.\n");
  printf(
      "\t\t\tValid format for OPTION_STRING is \"%%c%%c%%c%%c\". (default: "
      "%s)\n",
      DEFAULT_TRMM_ARGS);

  printf("\t-a, --trmm_alpha=SCALAR_VALUE\n");
  printf("\t\tTRMM alpha value.\n");
  printf("\t\t\tThe value of alpha in floating point. (default: %lf)\n", DEFAULT_TRMM_ALPHA);

  printf("\t-g, --gemm_options=OPTION_STRING\n");
  printf("\t\tGEMM transA, and transB options.\n");
  printf(
      "\t\t\tValid format for OPTION_STRING is \"%%c%%c\". (default: "
      "%s)\n",
      DEFAULT_GEMM_ARGS);

  printf("\t-p, --gemm_scalars=ALPHA_SCALAR_VALUE,BETA_SCALAR_VALUE\n");
  printf("\t\tGEMM alpha and beta values.\n");
  printf(
      "\t\t\tThe value of alpha and beta in floating point. (default: "
      "%lf,%lf)\n",
      DEFAULT_GEMM_ALPHA, DEFAULT_GEMM_BETA);

  printf("\t-z, --team_size=SIZE\n");
  printf("\t\tKokkos team size.\n");
  printf("\t\t\tThe value of SIZE as an integer. (default: %d)\n", DEFAULT_TEAM_SIZE);

  printf("\t-n, --vector_len=LEN\n");
  printf("\t\tKokkos vector length (Heirarchical parallelism).\n");
  printf("\t\t\tThe value of LEN as an integer. (default: %d)\n", DEFAULT_VECTOR_LEN);

  printf("\t-u, --use_auto=AUTO\n");
  printf(
      "\t\tWhether to use Kokkos::AUTO for vector_len and team_size "
      "(Heirarchical parallelism).\n");
  printf(
      "\t\t\tValid values for AUTO are 1 to use Kokkos::AUTO and 0 to use "
      "--vector_len and --team_size "
      "instead. (default: %d)\n",
      DEFAULT_USE_AUTO);

  printf("\t-k, --batch_size=LEN\n");
  printf("\t\tBatch size. Adds third dimension to matrices A, B, and C.\n");
  printf("\t\t\tThe value of LEN as an integer. (default: %d)\n", DEFAULT_K);

  printf("\t-d, --batch_size_last_dim=LAST_DIM\n");
  printf("\t\tHow to allocate the batch_size in the matrices.\n");
  printf(
      "\t\t\tValid values for LAST_DIM are 1 make the batch_size the last "
      "dimension and 0 to make the batch_size "
      "the first dimension (default: %d)\n",
      DEFAULT_BATCH_SIZE_LAST_DIM);

  printf("\t-l, --loop_type=OPTION\n");
  printf("\t\tLoop selection.\n");
  printf("\t\t\tValid values for OPTION:\n");
  printf("%c[1m", 27);
  printf("\t\t\t\tserial:");
  printf("%c[0m", 27);
  printf(" invoke blas routine in a serial for-loop. (default)\n");
  printf("%c[1m", 27);
  printf("\t\t\t\tparallel:");
  printf("%c[0m", 27);
  printf(" invoke blas routine in a Kokkos::parallel_for-loop.\n");

  printf("\t-b, --matrix_size_start=MxN,IxJ,PxQ\n");
  printf(
      "\t\tMatrix size selection where A is MxN, B is IxJ, and C is PxQ "
      "(start)\n");
  printf(
      "\t\t\tValid values for M and N are any non-negative 32-bit integers. "
      "(default: %dx%d,%dx%d,%dx%d)\n",
      DEFAULT_MATRIX_START, DEFAULT_MATRIX_START, DEFAULT_MATRIX_START, DEFAULT_MATRIX_START, DEFAULT_MATRIX_START,
      DEFAULT_MATRIX_START);

  printf("\t-e, --matrix_size_stop=SxT,LxK,OxR\n");
  printf(
      "\t\tMatrix size selection where A is SxT, B is LxK, and C is OxR "
      "(stop)\n");
  printf(
      "\t\t\tValid dimension values are any non-negative 32-bit integers. "
      "(default: %dx%d,%dx%d,%dx%d)\n",
      DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP,
      DEFAULT_MATRIX_STOP);

  printf("\t-s, --matrix_size_step=K\n");
  printf("\t\tMatrix step selection.\n");
  printf(
      "\t\t\tValid value for K is any non-negative 32-bit integer. (default: "
      "%d)\n",
      DEFAULT_STEP);

  printf("\t-w, --warm_up_loop=LOOP\n");
  printf("\t\tWarm up loop selection. (untimed)\n");
  printf(
      "\t\t\tValid value for LOOP is any non-negative 32-bit integer that's <= "
      "ITER. (default: %d)\n",
      DEFAULT_WARM_UP_N);

  printf("\t-i, --iter=ITER\n");
  printf("\t\tIteration selection. (timed)\n");
  printf(
      "\t\t\tValid value for ITER is any non-negative 32-bit integer. "
      "(default: %d)\n",
      DEFAULT_N);

  printf("\t-c, --csv=/path/to/file.csv\n");
  printf("\t\tCsv output file selection.\n");
  printf(
      "\t\t\tValid value for /path/to/file.csv is any valid file name. "
      "(default: stdout)\n");

  printf("\t-r, --routines=ROUTINES\n");
  printf("\t\tRoutine selection.\n");
  printf(
      "\t\t\tValid value for ROUTINES is one of more valid blas3 routines "
      "delimited by a comma. (default: %s)\n",
      DEFAULT_BLAS_ROUTINES);

  printf("\t-v, --verify=VERIFY\n");
  printf("\t\tVerification selection. (untimed)\n");
  printf(
      "\t\t\tValid values for VERIFY are either 0 to skip verification or 1 to "
      "verify before timing. "
      "(default: %d)\n",
      DEFAULT_VERIFY);

  printf("\t-j, --ninter=NINTER\n");
  printf("\t\tInterleaving size for armpl. (untimed)\n");
  printf(
      "\t\t\tValid values for NINTER is any positive integer "
      "that evenly divides the batch size. "
      "(default: %d)\n",
      DEFAULT_NINTER);

  printf("\t-u, --use_simd=SIMD\n");
  printf("\t\tWhether to use SIMD views.\n");
  printf(
      "\t\t\tValid values for SIMD are 1 to use SIMD views and 0 to use "
      "non-SIMD"
      "views instead. (default: %d)\n",
      DEFAULT_USE_SIMD);
}

static void __blas3_perf_test_input_error(char ** /*argv*/, char short_opt, char *getopt_optarg) {
  fprintf(stderr, "ERROR: invalid option \"-%c %s\". Try --help.\n", short_opt, getopt_optarg);
  exit(-EINVAL);
}

int main(int argc, char **argv) {
  options_t options;
  int option_idx = 0, ret, i;
  char *n_str = nullptr, *adim = nullptr, *bdim = nullptr, *cdim = nullptr;
  std::filebuf fb;
  std::ostream out(&fb);
  char *out_file                          = nullptr;
  using rt_type                           = decltype(do_trmm_invoke);
  rt_type *routine_table[BLAS_ROUTINES_N] = {
      &do_trmm_invoke, &do_gemm_invoke
      // ADD MORE BLAS3 ROUTINES HERE
  };

  /* set default options */
  options.test                          = DEFAULT_TEST;
  options.loop                          = DEFAULT_LOOP;
  options.start.a.k                     = DEFAULT_K;
  options.start.a.m                     = DEFAULT_MATRIX_START;
  options.start.a.n                     = DEFAULT_MATRIX_START;
  options.stop.a.k                      = DEFAULT_K;
  options.stop.a.m                      = DEFAULT_MATRIX_STOP;
  options.stop.a.n                      = DEFAULT_MATRIX_STOP;
  options.start.b.k                     = DEFAULT_K;
  options.start.b.m                     = DEFAULT_MATRIX_START;
  options.start.b.n                     = DEFAULT_MATRIX_START;
  options.stop.b.k                      = DEFAULT_K;
  options.stop.b.m                      = DEFAULT_MATRIX_STOP;
  options.stop.b.n                      = DEFAULT_MATRIX_STOP;
  options.start.c.k                     = DEFAULT_K;
  options.start.c.m                     = DEFAULT_MATRIX_START;
  options.start.c.n                     = DEFAULT_MATRIX_START;
  options.stop.c.k                      = DEFAULT_K;
  options.stop.c.m                      = DEFAULT_MATRIX_STOP;
  options.stop.c.n                      = DEFAULT_MATRIX_STOP;
  options.step                          = DEFAULT_STEP;
  options.warm_up_n                     = DEFAULT_WARM_UP_N;
  options.n                             = DEFAULT_N;
  options.out                           = DEFAULT_OUT;
  options.blas_routines                 = std::string(DEFAULT_BLAS_ROUTINES);
  options.blas_args.team_size           = DEFAULT_TEAM_SIZE;
  options.blas_args.vector_len          = DEFAULT_VECTOR_LEN;
  options.blas_args.use_auto            = DEFAULT_USE_AUTO;
  options.blas_args.batch_size_last_dim = DEFAULT_BATCH_SIZE_LAST_DIM;
  options.verify                        = DEFAULT_VERIFY;
  options.ninter                        = DEFAULT_NINTER;
  options.use_simd                      = DEFAULT_USE_SIMD;

  options.blas_args.trmm.trmm_args = DEFAULT_TRMM_ARGS;
  options.blas_args.trmm.alpha     = DEFAULT_TRMM_ALPHA;

  options.blas_args.gemm.gemm_args = DEFAULT_GEMM_ARGS;
  options.blas_args.gemm.alpha     = DEFAULT_GEMM_ALPHA;
  options.blas_args.gemm.beta      = DEFAULT_GEMM_BETA;

  while ((ret = getopt_long(argc, argv, "ht:l:b:e:s:w:i:o:a:c:r:g:z:n:k:u:p:d:v:j:f:", long_options, &option_idx)) !=
         -1) {
    switch (ret) {
      case 'h': __print_help_blas3_perf_test(); return 0;
      case 't':
        for (i = 0; i < TEST_N; i++) {
          if (!test_e_str[i].compare(optarg)) {
            options.test = (test_e)i;
            break;
          }
        }
        if (i == TEST_N) {
          __blas3_perf_test_input_error(argv, ret, optarg);
        }
        break;
      case 'o':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        if (strlen(optarg) != 4) {
          __blas3_perf_test_input_error(argv, ret, optarg);
        }
        options.blas_args.trmm.trmm_args = optarg;
        break;
      case 'g':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        if (strlen(optarg) != 2) {
          __blas3_perf_test_input_error(argv, ret, optarg);
        }
        options.blas_args.gemm.gemm_args = optarg;
        break;
      case 'p':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        double alpha, beta;
        if (sscanf(optarg, "%lf,%lf", &alpha, &beta) != 2) __blas3_perf_test_input_error(argv, ret, optarg);

        options.blas_args.gemm.alpha = static_cast<default_scalar>(alpha);
        options.blas_args.gemm.beta  = static_cast<default_scalar>(beta);
        break;
      case 'a':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        options.blas_args.trmm.alpha = (default_scalar)atof(optarg);
        break;
      case 'l':
        for (i = 0; i < LOOP_N; i++) {
          if (!loop_e_str[i].compare(optarg)) {
            options.loop = (loop_e)i;
            break;
          }
        }
        if (i == LOOP_N) {
          __blas3_perf_test_input_error(argv, ret, optarg);
        }
        break;
      case 'b':
        adim    = optarg;
        bdim    = strcasestr(optarg, ",");
        bdim[0] = '\0';
        bdim    = &bdim[1];
        cdim    = strcasestr(bdim, ",");
        cdim[0] = '\0';
        cdim    = &cdim[1];

        n_str = strcasestr(adim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]          = '\0';
        options.start.a.m = atoi(adim);
        options.start.a.n = atoi(&n_str[1]);

        n_str = strcasestr(bdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]          = '\0';
        options.start.b.m = atoi(bdim);
        options.start.b.n = atoi(&n_str[1]);

        n_str = strcasestr(cdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]          = '\0';
        options.start.c.m = atoi(cdim);
        options.start.c.n = atoi(&n_str[1]);
        break;
      case 'e':
        adim    = optarg;
        bdim    = strcasestr(optarg, ",");
        bdim[0] = '\0';
        bdim    = &bdim[1];
        cdim    = strcasestr(bdim, ",");
        cdim[0] = '\0';
        cdim    = &cdim[1];

        n_str = strcasestr(adim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]         = '\0';
        options.stop.a.m = atoi(adim);
        options.stop.a.n = atoi(&n_str[1]);

        n_str = strcasestr(bdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]         = '\0';
        options.stop.b.m = atoi(bdim);
        options.stop.b.n = atoi(&n_str[1]);

        n_str = strcasestr(cdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, ret, optarg);

        n_str[0]         = '\0';
        options.stop.c.m = atoi(cdim);
        options.stop.c.n = atoi(&n_str[1]);
        break;
      case 's': options.step = atoi(optarg); break;
      case 'w': options.warm_up_n = atoi(optarg); break;
      case 'i': options.n = atoi(optarg); break;
      case 'k':
        options.start.a.k = options.start.b.k = options.start.c.k = options.stop.a.k = options.stop.b.k =
            options.stop.c.k                                                         = atoi(optarg);
        break;
      case 'd': options.blas_args.batch_size_last_dim = atoi(optarg); break;
      case 'v': options.verify = atoi(optarg); break;
      case 'j': options.ninter = atoi(optarg); break;
      case 'z': options.blas_args.team_size = atoi(optarg); break;
      case 'n': options.blas_args.vector_len = atoi(optarg); break;
      case 'u': options.blas_args.use_auto = atoi(optarg); break;
      case 'f': options.use_simd = atoi(optarg); break;
      case 'c':
        out_file         = optarg;
        options.out_file = std::string(out_file);
        break;
      case 'r': options.blas_routines = optarg; break;
      case '?':
      default: __blas3_perf_test_input_error(argv, ret, optarg);
    }
  }

  if (out_file != nullptr) {
    fb.open(out_file, std::ios::out);
    options.out = &out;
  }

  if (options.warm_up_n > options.n) {
    fprintf(stderr, "ERROR: warm_up_n=%d > n=%d. Try --help.\n", options.warm_up_n, options.n);
    exit(-EINVAL);
  }

  Kokkos::initialize(argc, argv);
  atexit(Kokkos::finalize);

  int err = 0;
  for (i = 0; i < BLAS_ROUTINES_N; i++) {
    if (options.blas_routines.find(blas_routines_e_str[i]) != std::string::npos) {
      std::cout << "Testing " << blas_routines_e_str[i] << "..." << std::endl;

      auto routine = routine_table[i];

      if (!routine || !routine[0][options.loop][options.test]) {
        std::cerr << "do_" << blas_routines_e_str[i] << "_invoke[";
        err = 1;
        break;
      }
      routine[0][options.loop][options.test](options);
    }
  }

  if (err) {
    std::cerr << loop_e_str[options.loop] << "][" << test_e_str[options.test] << "] not yet implemented!" << std::endl;
    exit(-EINVAL);
  }

  if (out_file != nullptr) fb.close();

  return 0;
}
