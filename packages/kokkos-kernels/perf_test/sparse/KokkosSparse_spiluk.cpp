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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <iomanip>  // std::setprecision

// cuSPARSE ILU and IC factorizations were removed
// completely in cuSPARSE 12.5
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && (CUSPARSE_VERSION < 12500)
#define USE_CUSPARSE_ILU
#endif

#ifdef USE_CUSPARSE_ILU
#include <cusparse.h>
#endif

#include <Kokkos_Core.hpp>

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_spiluk.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

enum { DEFAULT, CUSPARSE, LVLSCHED_RP, LVLSCHED_TP1 /*, LVLSCHED_TP2*/ };

int test_spiluk_perf(std::vector<int> tests, std::string afilename, int kin, int team_size, int /*vector_length*/,
                     /*int idx_offset,*/ int loop) {
  typedef default_scalar scalar_t;
  typedef int lno_t;
  typedef int size_type;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef typename execution_space::memory_space memory_space;

  typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, execution_space, void, size_type> crsmat_t;
  typedef typename crsmat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsmat_t::values_type::non_const_type scalar_view_t;

  typedef Kokkos::View<scalar_t *, memory_space> ValuesType;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, execution_space, memory_space,
                                                           memory_space>
      KernelHandle;
  printf("Execution space: %s, Memory space: %s\n", typeid(execution_space).name(), typeid(memory_space).name());
  scalar_t ZERO             = scalar_t(0);
  scalar_t ONE              = scalar_t(1);
  scalar_t MONE             = scalar_t(-1);
  constexpr int EXPAND_FACT = 6;

  // Read amtx
  // Run all requested algorithms

  std::cout << "\n\n" << std::endl;
  if (!afilename.empty()) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && !defined(USE_CUSPARSE_ILU)
    std::cout << "** Note: cuSPARSE is enabled, but the cusparseXcsrilu*\n";
    std::cout << "   functions were removed in cuSPARSE 12.5.\n";
    std::cout << "   Only KokkosKernels spiluk will be run.\n\n";
#endif
    std::cout << "ILU(K) Begin: Read matrix filename " << afilename << std::endl;
    crsmat_t A            = KokkosSparse::Impl::read_kokkos_crst_matrix<crsmat_t>(afilename.c_str());  // in_matrix
    graph_t graph         = A.graph;                                                                   // in_graph
    const size_type nrows = graph.numRows();
    const int nnz         = A.nnz();
    const typename KernelHandle::const_nnz_lno_t fill_lev = lno_t(kin);

#ifdef USE_CUSPARSE_ILU
    // std::cout << "  cusparse: create handle" << std::endl;
    cusparseStatus_t status;
    cusparseHandle_t handle = 0;
    status                  = cusparseCreate(&handle);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "handle create status error name " << (status) << std::endl;
    cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
    cusparseMatDescr_t descr = 0;
    csrilu02Info_t info      = 0;
    int pBufferSize;
    void *pBuffer = 0;
    int structural_zero;
    int numerical_zero;
    const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;

    // step 1: create a descriptor
    status = cusparseCreateMatDescr(&descr);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "MatDescr create status error name " << (status) << std::endl;
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);

    // step 2: create a empty info structure
    status = cusparseCreateCsrilu02Info(&info);
    if (CUSPARSE_STATUS_SUCCESS != status)
      std::cout << "Csrilu02Info create status error name " << (status) << std::endl;

    // step 3: query how much memory used in csrsv2, and allocate the buffer
    cusparseDcsrilu02_bufferSize(handle, nrows, nnz, descr, A.values.data(), A.graph.row_map.data(),
                                 A.graph.entries.data(), info, &pBufferSize);
    // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
    cudaMalloc((void **)&pBuffer, pBufferSize);
#endif

    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle kh;

      // std::cout << "Create handle" << std::endl;
      switch (test) {
        case LVLSCHED_TP1:
          kh.create_spiluk_handle(SPILUKAlgorithm::SEQLVLSCHD_TP1, nrows, EXPAND_FACT * nnz * (fill_lev + 1),
                                  EXPAND_FACT * nnz * (fill_lev + 1));
          kh.get_spiluk_handle()->print_algorithm();
          kh.get_spiluk_handle()->set_team_size(team_size);
          break;
        // case LVLSCHED_TP2:
        //  kh.create_spiluk_handle(SPILUKAlgorithm::SEQLVLSCHED_TP2, nrows,
        //  EXPAND_FACT*nnz*(fill_lev+1), EXPAND_FACT*nnz*(fill_lev+1));
        //  kh.get_spiluk_handle()->print_algorithm();
        //  break;
        default:
          kh.create_spiluk_handle(SPILUKAlgorithm::SEQLVLSCHD_TP1, nrows, EXPAND_FACT * nnz * (fill_lev + 1),
                                  EXPAND_FACT * nnz * (fill_lev + 1));
          kh.get_spiluk_handle()->print_algorithm();
          kh.get_spiluk_handle()->set_team_size(team_size);
      }

      lno_view_t L_row_map("L_row_map", nrows + 1);
      lno_nnz_view_t L_entries("L_entries", kh.get_spiluk_handle()->get_nnzL());
      scalar_view_t L_values("L_values", kh.get_spiluk_handle()->get_nnzL());
      lno_view_t U_row_map("U_row_map", nrows + 1);
      lno_nnz_view_t U_entries("U_entries", kh.get_spiluk_handle()->get_nnzU());
      scalar_view_t U_values("U_values", kh.get_spiluk_handle()->get_nnzU());

      // Init run to clear cache etc.
      Kokkos::Timer timer;

      timer.reset();
      spiluk_symbolic(&kh, fill_lev, A.graph.row_map, A.graph.entries, L_row_map, L_entries, U_row_map, U_entries);
      std::cout << "ILU(" << fill_lev << ") Symbolic Time: " << timer.seconds() << std::endl;

      Kokkos::resize(L_entries, kh.get_spiluk_handle()->get_nnzL());
      Kokkos::resize(L_values, kh.get_spiluk_handle()->get_nnzL());
      Kokkos::resize(U_entries, kh.get_spiluk_handle()->get_nnzU());
      Kokkos::resize(U_values, kh.get_spiluk_handle()->get_nnzU());

      std::cout << "num levels: " << kh.get_spiluk_handle()->get_num_levels() << std::endl;
      std::cout << "max num rows levels: " << kh.get_spiluk_handle()->get_level_maxrows() << std::endl;
      std::cout << "team size: " << kh.get_spiluk_handle()->get_team_size() << std::endl;
      std::cout << "vector size: " << kh.get_spiluk_handle()->get_vector_size() << std::endl;
      std::cout << "nnzL: " << kh.get_spiluk_handle()->get_nnzL() << std::endl;
      std::cout << "nnzU: " << kh.get_spiluk_handle()->get_nnzU() << std::endl;

      timer.reset();
      spiluk_numeric(&kh, fill_lev, A.graph.row_map, A.graph.entries, A.values, L_row_map, L_entries, L_values,
                     U_row_map, U_entries, U_values);
      Kokkos::fence();
      std::cout << "ILU(" << fill_lev << ") Numeric Time: " << timer.seconds() << std::endl;

      crsmat_t L("L", nrows, nrows, kh.get_spiluk_handle()->get_nnzL(), L_values, L_row_map, L_entries);
      crsmat_t U("U", nrows, nrows, kh.get_spiluk_handle()->get_nnzU(), U_values, U_row_map, U_entries);
      ValuesType e_one("e_one", nrows);
      ValuesType bb("bb", nrows);
      ValuesType bb_tmp("bb_tmp", nrows);

      Kokkos::deep_copy(e_one, scalar_t(1));

      KokkosSparse::spmv("N", ONE, A, e_one, ZERO, bb);
      KokkosSparse::spmv("N", ONE, U, e_one, ZERO, bb_tmp);
      KokkosSparse::spmv("N", ONE, L, bb_tmp, MONE, bb);

      scalar_t bb_nrm = KokkosBlas::nrm2(bb);

      std::cout << "nrm2(A*e-L*U*e) = " << std::setprecision(15) << bb_nrm << std::endl;

#ifdef USE_CUSPARSE_ILU
      if (fill_lev == 0) {
        std::cout << "CUSPARSE: No KK interface added yet" << std::endl;

        lno_view_t A_row_map("A_row_map", nrows + 1);
        lno_nnz_view_t A_entries("A_entries", nnz);
        scalar_view_t A_values("A_values", nnz);

        Kokkos::deep_copy(A_row_map, A.graph.row_map);
        Kokkos::deep_copy(A_entries, A.graph.entries);
        Kokkos::deep_copy(A_values, A.values);

        // step 4: perform analysis
        timer.reset();
        status = cusparseDcsrilu02_analysis(handle, nrows, A.nnz(), descr, A_values.data(), A_row_map.data(),
                                            A_entries.data(), info, policy, pBuffer);
        Kokkos::fence();
        std::cout << "cuSPARSE ILU(0) Symbolic Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;

        status = cusparseXcsrilu02_zeroPivot(handle, info, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
        }

        // step 5: M = L*U
        timer.reset();
        status = cusparseDcsrilu02(handle, nrows, A.nnz(), descr, A_values.data(), A_row_map.data(), A_entries.data(),
                                   info, policy, pBuffer);
        Kokkos::fence();
        std::cout << "cuSPARSE ILU(0) Numeric Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "numeric status error name " << (status) << std::endl;

        status = cusparseXcsrilu02_zeroPivot(handle, info, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("U(%d,%d) is zero\n", numerical_zero, numerical_zero);
        }

        // Error Check
        // KK
        auto h_L_row_map = Kokkos::create_mirror_view(L_row_map);
        auto h_L_entries = Kokkos::create_mirror_view(L_entries);
        auto h_L_values  = Kokkos::create_mirror_view(L_values);
        auto h_U_row_map = Kokkos::create_mirror_view(U_row_map);
        auto h_U_entries = Kokkos::create_mirror_view(U_entries);
        auto h_U_values  = Kokkos::create_mirror_view(U_values);

        Kokkos::deep_copy(h_L_row_map, L_row_map);
        Kokkos::deep_copy(h_L_entries, L_entries);
        Kokkos::deep_copy(h_L_values, L_values);
        Kokkos::deep_copy(h_U_row_map, U_row_map);
        Kokkos::deep_copy(h_U_entries, U_entries);
        Kokkos::deep_copy(h_U_values, U_values);
        // cuSPARSE
        auto h_A_row_map = Kokkos::create_mirror_view(A_row_map);
        auto h_A_entries = Kokkos::create_mirror_view(A_entries);
        auto h_A_values  = Kokkos::create_mirror_view(A_values);

        Kokkos::deep_copy(h_A_row_map, A_row_map);
        Kokkos::deep_copy(h_A_entries, A_entries);
        Kokkos::deep_copy(h_A_values, A_values);

        for (size_type i = 0; i < nrows; ++i) {
          auto a_row_start = h_A_row_map(i);
          auto a_row_end   = h_A_row_map(i + 1);
          auto l_row_start = h_L_row_map(i);
          auto l_row_end   = h_L_row_map(i + 1);
          auto u_row_start = h_U_row_map(i);
          auto u_row_end   = h_U_row_map(i + 1);
          if ((a_row_end - a_row_start) != ((l_row_end - l_row_start) + (u_row_end - u_row_start) - 1)) {
            std::cout << "ILU(0) FAILURE: nnz on row " << i
                      << " do not match -- KK = " << (l_row_end - l_row_start) + (u_row_end - u_row_start) - 1
                      << ", cuSPARSE = " << a_row_end - a_row_start << std::endl;
            return 1;
          } else {
            Kokkos::View<lno_t *, Kokkos::LayoutLeft, Kokkos::HostSpace> h_tmp_entries("h_tmp_entries",
                                                                                       a_row_end - a_row_start);
            Kokkos::View<scalar_t *, Kokkos::LayoutLeft, Kokkos::HostSpace> h_tmp_values("h_tmp_values",
                                                                                         a_row_end - a_row_start);

            Kokkos::deep_copy(subview(h_tmp_entries, Kokkos::make_pair(0, l_row_end - l_row_start)),
                              subview(h_L_entries, Kokkos::make_pair(l_row_start, l_row_end)));  // L part
            Kokkos::deep_copy(
                subview(h_tmp_entries, Kokkos::make_pair(l_row_end - l_row_start - 1, a_row_end - a_row_start)),
                subview(h_U_entries, Kokkos::make_pair(u_row_start, u_row_end + 1)));  // U part

            Kokkos::deep_copy(subview(h_tmp_values, Kokkos::make_pair(0, l_row_end - l_row_start)),
                              subview(h_L_values, Kokkos::make_pair(l_row_start, l_row_end)));  // L part
            Kokkos::deep_copy(
                subview(h_tmp_values, Kokkos::make_pair(l_row_end - l_row_start - 1, a_row_end - a_row_start)),
                subview(h_U_values, Kokkos::make_pair(u_row_start, u_row_end + 1)));  // U part

            for (size_type k = 0; k < (a_row_end - a_row_start); ++k) {
              if (h_tmp_entries(k) != h_A_entries(a_row_start + k)) {
                if (h_A_entries(a_row_start + k) < i)
                  std::cout << "ILU(0) FAILURE: non-zero col idx on row " << i
                            << " do not match -- KK (L part) = " << h_tmp_entries(k)
                            << ", cuSPARSE = " << h_A_entries(a_row_start + k) << std::endl;
                else
                  std::cout << "ILU(0) FAILURE: non-zero col idx on row " << i
                            << " do not match -- KK (U part) = " << h_tmp_entries(k)
                            << ", cuSPARSE = " << h_A_entries(a_row_start + k) << std::endl;
                return 1;
              } else if (abs(h_tmp_values(k) - h_A_values(a_row_start + k)) > 1e-3) {
                if (h_A_entries(a_row_start + k) < i)
                  std::cout << "ILU(0) FAILURE: non-zero entry on row " << i
                            << " do not match -- KK (L part) = " << h_tmp_values(k) << " at col " << h_tmp_entries(k)
                            << ", cuSPARSE = " << h_A_values(a_row_start + k) << " at col "
                            << h_A_entries(a_row_start + k) << std::endl;
                else
                  std::cout << "ILU(0) FAILURE: non-zero entry on row " << i
                            << " do not match -- KK (U part) = " << h_tmp_values(k) << " at col " << h_tmp_entries(k)
                            << ", cuSPARSE = " << h_A_values(a_row_start + k) << " at col "
                            << h_A_entries(a_row_start + k) << std::endl;
                return 1;
              }
            }  // end col
          }
        }  // end row
        std::cout << "ILU(0) SUCCESS!" << std::endl;
      }  // fill_lev=0
#endif

      // Benchmark
      Kokkos::fence();
      double min_time = std::numeric_limits<double>::infinity();
      double max_time = 0.0;
      double ave_time = 0.0;

      for (int i = 0; i < loop; i++) {
        timer.reset();
        spiluk_numeric(&kh, fill_lev, A.graph.row_map, A.graph.entries, A.values, L_row_map, L_entries, L_values,
                       U_row_map, U_entries, U_values);
        Kokkos::fence();
        double time = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }
      std::cout << "LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
      std::cout << "LOOP_MAX_TIME:  " << max_time << std::endl;
      std::cout << "LOOP_MIN_TIME:  " << min_time << std::endl;

#ifdef USE_CUSPARSE_ILU
      if (fill_lev == 0) {
        lno_view_t A_row_map("A_row_map", nrows + 1);
        lno_nnz_view_t A_entries("A_entries", nnz);
        scalar_view_t A_values("A_values", nnz);

        min_time = std::numeric_limits<double>::infinity();
        max_time = 0.0;
        ave_time = 0.0;

        for (int i = 0; i < loop; i++) {
          Kokkos::deep_copy(A_row_map, A.graph.row_map);
          Kokkos::deep_copy(A_entries, A.graph.entries);
          Kokkos::deep_copy(A_values, A.values);

          timer.reset();
          cusparseDcsrilu02(handle, nrows, A.nnz(), descr, A_values.data(), A_row_map.data(), A_entries.data(), info,
                            policy, pBuffer);
          Kokkos::fence();
          double time = timer.seconds();
          ave_time += time;
          if (time > max_time) max_time = time;
          if (time < min_time) min_time = time;
        }
        std::cout << "LOOP_AVG_TIME (cuSPARSE):  " << ave_time / loop << std::endl;
        std::cout << "LOOP_MAX_TIME (cuSPARSE):  " << max_time << std::endl;
        std::cout << "LOOP_MIN_TIME (cuSPARSE):  " << min_time << std::endl;
      }  // fill_lev=0
#endif
    }  // end tests

#ifdef USE_CUSPARSE_ILU
    // step 6: free resources
    cudaFree(pBuffer);
    cusparseDestroyCsrilu02Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);
#endif
  }  // end if (!afilename.empty())

  std::cout << "\n\n" << std::endl;

  return 0;
}

void print_help_spiluk() {
  printf("Options:\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                      lvlrp, lvltp1, lvltp2\n\n");
  printf(
      "  -f [file]       : Read in Matrix Market formatted text file "
      "'file'.\n");
  //  printf("  -s [N]          : generate a semi-random banded (band size
  //  0.01xN) NxN matrix\n"); printf("                    with average of 10
  //  entries per row.\n"); printf("  --schedule [SCH]: Set schedule for kk
  //  variant (static,dynamic,auto [ default ]).\n"); printf("  -afb [file] :
  //  Read in binary Matrix files 'file'.\n"); printf("  --write-binary  : In
  //  combination with -f, generate binary files.\n"); printf("  --offset [O] :
  //  Subtract O from every index.\n"); printf("                    Useful in
  //  case the matrix market file is not 0 based.\n\n");
  printf("  -k [K]          : Fill level (default: 0)\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf(
      "  -vl [V]         : Vector-length (i.e. how many Cuda threads are a "
      "Kokkos 'thread').\n");
  printf(
      "  --loop [LOOP]   : How many spiluk to run to aggregate average time. "
      "\n");
}

int main(int argc, char **argv) {
  std::vector<int> tests;

  std::string afilename;

  int kin           = 0;
  int vector_length = -1;
  int team_size     = -1;
  // int idx_offset = 0;
  int loop = 1;
  // int schedule=AUTO;

  if (argc == 1) {
    print_help_spiluk();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if ((strcmp(argv[i], "lvlrp") == 0)) {
        tests.push_back(LVLSCHED_RP);
      }
      if ((strcmp(argv[i], "lvltp1") == 0)) {
        tests.push_back(LVLSCHED_TP1);
      }
      /*
            if((strcmp(argv[i],"lvltp2")==0)) {
              tests.push_back( LVLSCHED_TP2 );
            }
      */
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0)) {
      afilename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-k") == 0)) {
      kin = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ts") == 0)) {
      team_size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-vl") == 0)) {
      vector_length = atoi(argv[++i]);
      continue;
    }
    // if((strcmp(argv[i],"--offset")==0)) {idx_offset = atoi(argv[++i]);
    // continue;}
    if ((strcmp(argv[i], "--loop") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    /*
        if((strcmp(argv[i],"-afb")==0)) {afilename = argv[++i]; binaryfile =
       true; continue;} if((strcmp(argv[i],"--schedule")==0)) { i++;
          if((strcmp(argv[i],"auto")==0))
            schedule = AUTO;
          if((strcmp(argv[i],"dynamic")==0))
            schedule = DYNAMIC;
          if((strcmp(argv[i],"static")==0))
            schedule = STATIC;
          continue;
        }
    */
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help_spiluk();
      return 0;
    }
  }

  if (tests.size() == 0) {
    tests.push_back(DEFAULT);
  }
  for (size_t i = 0; i < tests.size(); ++i) {
    std::cout << "tests[" << i << "] = " << tests[i] << std::endl;
  }

  Kokkos::initialize(argc, argv);
  {
    int total_errors = test_spiluk_perf(tests, afilename, kin, team_size, vector_length, /*idx_offset,*/ loop);

    if (total_errors == 0)
      printf("Kokkos::SPILUK Test: Passed\n");
    else
      printf("Kokkos::SPILUK Test: Failed\n");
  }
  Kokkos::finalize();
  return 0;
}
