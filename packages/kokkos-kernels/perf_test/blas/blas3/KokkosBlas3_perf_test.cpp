/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include "KokkosBlas3_common.hpp"
#include "KokkosBlas3_trmm_perf_test.hpp"

#include <cstdlib>
#include <unistd.h>
#include <getopt.h>

static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"test", required_argument, 0, 't'},
    {"loop_type", required_argument, 0, 'l'},
    {"matrix_size_start", required_argument, 0, 'b'},
    {"matrix_size_stop", required_argument, 0, 'e'},
    {"matrix_size_step", required_argument, 0, 's'},
    {"warm_up_loop", required_argument, 0, 'w'},
    {"iter", required_argument, 0, 'i'},
    {"csv", required_argument, 0, 'c'},
    {"routines", required_argument, 0, 'r'},
    {"trmm_options", required_argument, 0, 'o'},
    {"trmm_alpha", required_argument, 0, 'a'},
    {0, 0, 0, 0}};

static void __print_help_blas3_perf_test() {
  printf("Options:\n");

  printf("\t-h, --help\n");
  printf("\t\tPrint this help menu.\n\n");

  printf("\t-t, --test=OPTION\n");
  printf("\t\tAlgorithm selection.\n");
  printf("\t\t\tValid values for OPTION:\n");
  printf("%c[1m", 27);
  printf("\t\t\t\tblas:");
  printf("%c[0m", 27);
  printf(" invoke Kokkos::trmm the loop-body. (default)\n");
  printf("%c[1m", 27);
  printf("\t\t\t\tbatched:");
  printf("%c[0m", 27);
  printf(" invoke KokkosBatched::SerialTrmm in the loop-body.\n\n");

  printf("\t-o, --trmm_options=OPTION_STRING\n");
  printf("\t\tTRMM side, uplo, trans, and diag options.\n");
  printf(
      "\t\t\tValid format for OPTION_STRING is \"%%c%%c%%c%%c\". (default: "
      "%s)\n",
      DEFAULT_TRMM_ARGS);

  printf("\t-a, --trmm_alpha=SCALAR_VALUE\n");
  printf("\t\tTRMM alpha value.\n");
  printf("\t\t\tThe value of alpha in floating point. (default: %lf)\n",
         DEFAULT_TRMM_ALPHA);

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
  printf(" invoke blas routine in a Kokkos::parallel_for-loop.\n\n");

  printf("\t-b, --matrix_size_start=MxN,IxJ\n");
  printf("\t\tMatrix size selection where A is MxN and B is IxJ (start)\n");
  printf(
      "\t\t\tValid values for M and N are any non-negative 32-bit integers. "
      "(default: %dx%d,%dx%d)\n\n",
      DEFAULT_MATRIX_START, DEFAULT_MATRIX_START, DEFAULT_MATRIX_START,
      DEFAULT_MATRIX_START);

  printf("\t-e, --matrix_size_stop=PxQ,SxT\n");
  printf("\t\tMatrix size selection where A is PxQ and B is SxT (stop)\n");
  printf(
      "\t\t\tValid values for P and Q are any non-negative 32-bit integers. "
      "(default: %dx%d,%dx%d)\n\n",
      DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP, DEFAULT_MATRIX_STOP,
      DEFAULT_MATRIX_STOP);

  printf("\t-s, --matrix_size_step=K\n");
  printf("\t\tMatrix step selection.\n");
  printf(
      "\t\t\tValid value for K is any non-negative 32-bit integer. (default: "
      "%d)\n\n",
      DEFAULT_STEP);

  printf("\t-w, --warm_up_loop=LOOP\n");
  printf("\t\tWarm up loop selection. (untimed)\n");
  printf(
      "\t\t\tValid value for LOOP is any non-negative 32-bit integer that's <= "
      "ITER. (default: %d)\n\n",
      DEFAULT_WARM_UP_N);

  printf("\t-i, --iter=ITER\n");
  printf("\t\tIteration selection. (timed)\n");
  printf(
      "\t\t\tValid value for ITER is any non-negative 32-bit integer. "
      "(default: %d)\n\n",
      DEFAULT_N);

  printf("\t-c, --csv=/path/to/file.csv\n");
  printf("\t\tCsv output file selection.\n");
  printf(
      "\t\t\tValid value for /path/to/file.csv is any valid file name. "
      "(default: stdout)\n\n");

  printf("\t-r, --routines=ROUTINES\n");
  printf("\t\tRoutine selection.\n");
  printf(
      "\t\t\tValid value for ROUTINES is one of more valid blas3 routines "
      "delimited by a comma. (default: %s)\n\n",
      DEFAULT_BLAS_ROUTINES);
}

static void __blas3_perf_test_input_error(char **argv, int option_idx) {
  fprintf(stderr, "ERROR: invalid option \"%s %s\".\n", argv[option_idx],
          argv[option_idx + 1]);
  __print_help_blas3_perf_test();
  exit(-EINVAL);
}

int main(int argc, char **argv) {
  options_t options;
  int option_idx = 0, ret;
  char *n_str = nullptr, *adim = nullptr, *bdim = nullptr;
  std::filebuf fb;
  char *out_file = nullptr;

  /* set default options */
  options.test          = DEFAULT_TEST;
  options.loop          = DEFAULT_LOOP;
  options.start.a.m     = DEFAULT_MATRIX_START;
  options.start.a.n     = DEFAULT_MATRIX_START;
  options.stop.a.m      = DEFAULT_MATRIX_STOP;
  options.stop.a.n      = DEFAULT_MATRIX_STOP;
  options.start.b.m     = DEFAULT_MATRIX_START;
  options.start.b.n     = DEFAULT_MATRIX_START;
  options.stop.b.m      = DEFAULT_MATRIX_STOP;
  options.stop.b.n      = DEFAULT_MATRIX_STOP;
  options.step          = DEFAULT_STEP;
  options.warm_up_n     = DEFAULT_WARM_UP_N;
  options.n             = DEFAULT_N;
  options.out           = DEFAULT_OUT;
  options.blas_routines = std::string(DEFAULT_BLAS_ROUTINES);

  options.blas_args.trmm.trmm_args = DEFAULT_TRMM_ARGS;
  options.blas_args.trmm.alpha     = DEFAULT_TRMM_ALPHA;

  while ((ret = getopt_long(argc, argv, "ht:l:b:e:s:w:i:o:a:c:r:", long_options,
                            &option_idx)) != -1) {
    switch (ret) {
      case 'h': __print_help_blas3_perf_test(); return 0;
      case 't':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        if (!strncasecmp(optarg, "blas", 4)) {
          options.test = BLAS;
        } else if (!strncasecmp(optarg, "batched", 6)) {
          options.test = BATCHED;
        } else {
          __blas3_perf_test_input_error(argv, option_idx);
        }
        break;
      case 'o':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        if (strlen(optarg) != 4) {
          __blas3_perf_test_input_error(argv, option_idx);
        }
        options.blas_args.trmm.trmm_args = optarg;
        break;
      case 'a':
        // printf("optarg=%s. %d\n", optarg, strncasecmp(optarg, "blas", 4));
        options.blas_args.trmm.alpha = (default_scalar)atof(optarg);
        break;
      case 'l':
        if (!strncasecmp(optarg, "serial", 6)) {
          options.loop = SERIAL;
        } else if (!strncasecmp(optarg, "parallel", 8)) {
          options.loop = PARALLEL;
        } else {
          __blas3_perf_test_input_error(argv, option_idx);
        }
        break;
      case 'b':
        adim    = optarg;
        bdim    = strcasestr(optarg, ",");
        bdim[0] = '\0';
        bdim    = &bdim[1];

        n_str = strcasestr(adim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, option_idx);

        n_str[0]          = '\0';
        options.start.a.m = atoi(adim);
        options.start.a.n = atoi(&n_str[1]);

        n_str = strcasestr(bdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, option_idx);

        n_str[0]          = '\0';
        options.start.b.m = atoi(bdim);
        options.start.b.n = atoi(&n_str[1]);
        break;
      case 'e':
        adim    = optarg;
        bdim    = strcasestr(optarg, ",");
        bdim[0] = '\0';
        bdim    = &bdim[1];

        n_str = strcasestr(adim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, option_idx);

        n_str[0]         = '\0';
        options.stop.a.m = atoi(adim);
        options.stop.a.n = atoi(&n_str[1]);

        n_str = strcasestr(bdim, "x");
        if (n_str == NULL) __blas3_perf_test_input_error(argv, option_idx);

        n_str[0]         = '\0';
        options.stop.b.m = atoi(bdim);
        options.stop.b.n = atoi(&n_str[1]);
        break;
      case 's': options.step = atoi(optarg); break;
      case 'w': options.warm_up_n = atoi(optarg); break;
      case 'i': options.n = atoi(optarg); break;
      case 'c':
        out_file         = optarg;
        options.out_file = std::string(out_file);
        break;
      case 'r': options.blas_routines = std::string(optarg); break;
      case '?':
      default: __blas3_perf_test_input_error(argv, option_idx);
    }
  }

  if (out_file != nullptr) {
    fb.open(out_file, std::ios::out);
    std::ostream out(&fb);
    options.out = &out;
  }

  if (options.warm_up_n > options.n)
    __blas3_perf_test_input_error(argv, option_idx);

  Kokkos::initialize(argc, argv);

  for (int i = 0; i < BLAS_ROUTINES_N; i++) {
    if (options.blas_routines.find(blas_routines_e_str[TRMM]) !=
        std::string::npos)
      do_trmm_invoke[options.loop][options.test](options);
    // ADD MORE BLAS3 ROUTINES HERE
  }

  if (out_file != nullptr) fb.close();

  Kokkos::finalize();

  return 0;
}
