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
#include <limits.h>
#include <cmath>
#include <unordered_map>

#include <sstream>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_Test_Structured_Matrix.hpp>
#include <KokkosSparse_spmv_struct_impl.hpp>
#include <KokkosSparse_spmv_impl.hpp>

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse.h>
#endif

enum { STRUCT, UNSTR };

void print_help() {
  printf("SPMV_struct benchmark code written by Luc Berger-Vergiat.\n");
  printf("Options:\n");
  printf(
      "  --check-errors     : Determine if the result of spmv_struct is "
      "compared to serial unstructured spmv.\n");
  printf(
      "  --compare          : Compare results efficiency of spmv_struct and "
      "spmv.\n");
  printf(
      "  --compare-cusparse : Compare results efficiency of spmv_struct and "
      "cusparseDcsrmv.\n");
  printf(
      "  -l [LOOP]          : How many spmv to run to aggregate average time. "
      "\n");
  printf("  -nx                : How many nodes in x direction. \n");
  printf("  -ny                : How many nodes in y direction. \n");
  printf("  -nz                : How many nodes in z direction. \n");
  printf(
      "  -st                : The stencil type used for discretization: 1 -> "
      "FD, 2 -> FE.\n");
  printf(
      "  -dim               : Number of spacial dimensions used in the "
      "problem: 1, 2 or 3\n");
  printf("  -ws                : Worksets. \n");
  printf("  -ts                : Team Size. \n");
  printf("  -vl                : Vector length. \n");
  printf("  -wsu               : Worksets for unstructured spmv. \n");
  printf("  -tsu               : Team Size for unstructured spmv. \n");
  printf("  -vlu               : Vector length for unstructured spmv. \n");
  printf("  --print-lp         : Print launch parameters to screen.\n");
}

template <typename graph_type>
struct copy_crs_data {
  using execution_space = typename graph_type::device_type::execution_space;
  using cusparse_int_type =
      typename Kokkos::View<int*, typename graph_type::entries_type::array_layout, typename graph_type::device_type>;

  // Dispatch tags
  struct rowPtrTag {};
  struct colIndTag {};

  typename graph_type::row_map_type::const_type Arowptr;
  typename graph_type::entries_type::const_type Acolind;
  cusparse_int_type cusparse_Arowptr, cusparse_Acolind;

  copy_crs_data(typename graph_type::row_map_type::const_type Arowptr_,
                typename graph_type::entries_type::const_type Acolind_, cusparse_int_type cusparse_Arowptr_,
                cusparse_int_type cusparse_Acolind_)
      : Arowptr(Arowptr_),
        Acolind(Acolind_),
        cusparse_Arowptr(cusparse_Arowptr_),
        cusparse_Acolind(cusparse_Acolind_){};

  void doCopy() {
    Kokkos::RangePolicy<execution_space, rowPtrTag> rowPtrPolicy(0, Arowptr.extent(0));
    Kokkos::parallel_for("copy rowPtr to cusparse", rowPtrPolicy, *this);

    Kokkos::RangePolicy<execution_space, colIndTag> colIndPolicy(0, Acolind.extent(0));
    Kokkos::parallel_for("copy colInd to cusparse", colIndPolicy, *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const rowPtrTag&, const size_t idx) const { cusparse_Arowptr(idx) = int(Arowptr(idx)); }

  KOKKOS_INLINE_FUNCTION
  void operator()(const colIndTag&, const size_t idx) const { cusparse_Acolind(idx) = int(Acolind(idx)); }
};

template <class AMatrix, class XVector, class YVector>
void struct_matvec(const int stencil_type,
                   const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
                   typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
                   typename YVector::const_value_type& beta, YVector& y, int team_size_int, int vector_length,
                   int64_t rows_per_thread_int, const bool print_lp) {
  typedef typename AMatrix::ordinal_type ordinal_type;
  typedef typename AMatrix::execution_space execution_space;
  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  int nnzPerRow             = -1;
  int64_t numInteriorPts    = 0;
  int64_t numExteriorPts    = 0;
  int vector_length_default = 0;

  if (structure.extent(0) == 1) {
    numInteriorPts        = structure(0) - 2;
    numExteriorPts        = 2;
    vector_length_default = 1;
  } else if (structure.extent(0) == 2) {
    numInteriorPts = (structure(1) - 2) * (structure(0) - 2);
    numExteriorPts = 2 * (structure(1) + structure(0) - 2);
    if (stencil_type == 1) {
      vector_length_default = 2;
    } else if (stencil_type == 2) {
      vector_length_default = 4;
    }
  } else if (structure.extent(0) == 3) {
    numInteriorPts = (structure(2) - 2) * (structure(1) - 2) * (structure(0) - 2);
    numExteriorPts = structure(2) * structure(1) * structure(0) - numInteriorPts;
    if (stencil_type == 1) {
      vector_length_default = 2;
    } else if (stencil_type == 2) {
      vector_length_default = 8;
    }
  }
  vector_length = (vector_length == -1 ? vector_length_default : vector_length);

  int64_t rows_per_team_int = KokkosSparse::Impl::spmv_struct_launch_parameters<execution_space>(
      numInteriorPts, A.nnz(), nnzPerRow, rows_per_thread_int, team_size_int, vector_length);

  int64_t worksets_int = (numInteriorPts + rows_per_team_int - 1) / rows_per_team_int;

  int rows_per_thread_ext   = -1;
  int team_size_ext         = -1;
  int64_t rows_per_team_ext = KokkosSparse::Impl::spmv_struct_launch_parameters<execution_space>(
      numExteriorPts, A.nnz(), nnzPerRow, rows_per_thread_ext, team_size_ext, vector_length);

  int64_t worksets_ext = (numInteriorPts + rows_per_team_ext - 1) / rows_per_team_ext;

  KokkosSparse::Impl::SPMV_Struct_Functor<execution_space, AMatrix, XVector, YVector, 1, false> spmv_struct(
      structure, stencil_type, alpha, A, x, beta, y, rows_per_team_int, rows_per_team_ext);

  if (print_lp) {
    std::cout << "worksets=" << worksets_ext << ", team_size=" << team_size_ext << ", vector_length=" << vector_length
              << std::endl;
  }

  spmv_struct.compute_interior(execution_space{}, worksets_int, team_size_int, vector_length);
  spmv_struct.compute_exterior(execution_space{}, worksets_ext, team_size_ext, vector_length);

}  // struct_matvec

template <class AMatrix, class XVector, class YVector>
void matvec(typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
            typename YVector::const_value_type& beta, YVector& y, int team_size, int vector_length,
            int64_t rows_per_thread, const bool print_lp) {
  typedef typename AMatrix::ordinal_type ordinal_type;
  typedef typename AMatrix::execution_space execution_space;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  int64_t rows_per_team = KokkosSparse::Impl::spmv_launch_parameters<execution_space>(
      A.numRows(), A.nnz(), rows_per_thread, team_size, vector_length);
  int64_t worksets = (y.extent(0) + rows_per_team - 1) / rows_per_team;

  KokkosSparse::Impl::SPMV_Functor<execution_space, AMatrix, XVector, YVector, 1, false> func(alpha, A, x, beta, y,
                                                                                              rows_per_team);

  if (print_lp) {
    std::cout << "worksets=" << worksets << ", team_size=" << team_size << ", vector_length=" << vector_length
              << std::endl;
  }

  Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> > policy(1, 1);
  if (team_size == -1) {
    policy =
        Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> >(worksets, Kokkos::AUTO, vector_length);
  } else {
    policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> >(worksets, team_size, vector_length);
  }
  Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Static>", policy, func);

}  // matvec

int main(int argc, char** argv) {
  typedef double Scalar;

  int nx                = 100;
  int ny                = 100;
  int nz                = 100;
  int stencil_type      = 1;
  int numDimensions     = 2;
  int numVecs           = 1;
  bool check_errors     = false;
  bool compare          = false;
  bool compare_cusparse = false;
  bool print_lp         = false;
  int loop              = 100;
  int vl                = -1;
  int ts                = -1;
  int ws                = -1;
  int vlu               = -1;
  int tsu               = -1;
  int wsu               = -1;

  if (argc == 1) {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-nx") == 0)) {
      nx = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ny") == 0)) {
      ny = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-nz") == 0)) {
      nz = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-st") == 0)) {
      stencil_type = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-dim") == 0)) {
      numDimensions = atoi(argv[++i]);
      continue;
    }
    // if((strcmp(argv[i],"-mv" )==0)) {numVecs=atoi(argv[++i]); continue;}
    if ((strcmp(argv[i], "-l") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-vl") == 0)) {
      vl = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ts") == 0)) {
      ts = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ws") == 0)) {
      ws = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-vlu") == 0)) {
      vlu = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-tsu") == 0)) {
      tsu = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-wsu") == 0)) {
      wsu = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--check-errors") == 0)) {
      check_errors = true;
      continue;
    }
    if ((strcmp(argv[i], "--compare") == 0)) {
      compare = true;
      continue;
    }
    if ((strcmp(argv[i], "--compare-cusparse") == 0)) {
      compare_cusparse = true;
      continue;
    }
    if ((strcmp(argv[i], "--print-lp") == 0)) {
      print_lp = true;
      continue;
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help();
      return 0;
    }
  }

  Kokkos::initialize(argc, argv);

  {
    using matrix_type = KokkosSparse::CrsMatrix<Scalar, int, Kokkos::DefaultExecutionSpace, void, int>;
    using mv_type     = typename Kokkos::View<Scalar*, Kokkos::LayoutLeft>;
    using h_mv_type   = typename mv_type::HostMirror;

    int leftBC = 1, rightBC = 1, frontBC = 1, backBC = 1, bottomBC = 1, topBC = 1;

    Kokkos::View<int*, Kokkos::HostSpace> structure("Spmv Structure", numDimensions);
    Kokkos::View<int* [3], Kokkos::HostSpace> mat_structure("Matrix Structure", numDimensions);
    if (numDimensions == 1) {
      structure(0)        = nx;
      mat_structure(0, 0) = nx;
    } else if (numDimensions == 2) {
      structure(0)        = nx;
      structure(1)        = ny;
      mat_structure(0, 0) = nx;
      mat_structure(1, 0) = ny;
      if (leftBC == 1) {
        mat_structure(0, 1) = 1;
      }
      if (rightBC == 1) {
        mat_structure(0, 2) = 1;
      }
      if (bottomBC == 1) {
        mat_structure(1, 1) = 1;
      }
      if (topBC == 1) {
        mat_structure(1, 2) = 1;
      }
    } else if (numDimensions == 3) {
      structure(0)        = nx;
      structure(1)        = ny;
      structure(2)        = nz;
      mat_structure(0, 0) = nx;
      mat_structure(1, 0) = ny;
      mat_structure(2, 0) = nz;
      if (leftBC == 1) {
        mat_structure(0, 1) = 1;
      }
      if (rightBC == 1) {
        mat_structure(0, 2) = 1;
      }
      if (frontBC == 1) {
        mat_structure(1, 1) = 1;
      }
      if (backBC == 1) {
        mat_structure(1, 2) = 1;
      }
      if (bottomBC == 1) {
        mat_structure(2, 1) = 1;
      }
      if (topBC == 1) {
        mat_structure(2, 2) = 1;
      }
    }

    std::string discrectization_stencil;
    if (stencil_type == 1) {
      discrectization_stencil = "FD";
    } else if (stencil_type == 2) {
      discrectization_stencil = "FE";
    }

    matrix_type A;
    if (numDimensions == 1) {
      A = Test::generate_structured_matrix1D<matrix_type>(mat_structure);
    } else if (numDimensions == 2) {
      A = Test::generate_structured_matrix2D<matrix_type>(discrectization_stencil, mat_structure);
    } else if (numDimensions == 3) {
      A = Test::generate_structured_matrix3D<matrix_type>(discrectization_stencil, mat_structure);
    }

    mv_type x("X", A.numCols());
    // mv_random_read_type t_x(x);
    mv_type y("Y", A.numRows());
    h_mv_type h_x = Kokkos::create_mirror_view(x);
    h_mv_type h_y = Kokkos::create_mirror_view(y);
    h_mv_type h_y_compare;

    for (int rowIdx = 0; rowIdx < A.numCols(); ++rowIdx) {
      for (int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
        h_x(rowIdx) = static_cast<Scalar>(1.0 * (rand() % 40) - 20.0);
        h_y(rowIdx) = static_cast<Scalar>(1.0 * (rand() % 40) - 20.0);
      }
    }

    if (check_errors) {
      h_y_compare                                                  = Kokkos::create_mirror(y);
      typename matrix_type::StaticCrsGraphType::HostMirror h_graph = Kokkos::create_mirror(A.graph);
      typename matrix_type::values_type::HostMirror h_values       = Kokkos::create_mirror_view(A.values);

      // Error Check Gold Values
      for (int rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        int start = h_graph.row_map(rowIdx);
        int end   = h_graph.row_map(rowIdx + 1);

        for (int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
          h_y_compare(rowIdx) = 0;
        }

        for (int entryIdx = start; entryIdx < end; ++entryIdx) {
          // Scalar tmp_val = h_graph.entries(entryIdx) + i;
          int colIdx = h_graph.entries(entryIdx);
          for (int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
            h_y_compare(rowIdx) += h_values(entryIdx) * h_x(colIdx);
          }
        }
      }
    }

    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);
    Kokkos::View<Scalar*, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> x1("X1", A.numCols());
    Kokkos::View<Scalar*, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> y1("Y1", A.numRows());
    Kokkos::deep_copy(x1, h_x);

    printf(
        "Type NNZ NumRows NumCols ProblemSize(MB) AveBandwidth(GB/s) "
        "MinBandwidth(GB/s) MaxBandwidth(GB/s) AveGFlop MinGFlop MaxGFlop "
        "aveTime(ms) maxTime(ms) minTime(ms)\n");
    {
      Kokkos::Profiling::pushRegion("Structured spmv test");
      // Benchmark
      bool print_lp_struct = print_lp;
      double min_time      = 1.0e32;
      double max_time      = 0.0;
      double ave_time      = 0.0;
      for (int i = 0; i < loop; i++) {
        Kokkos::Timer timer;
        struct_matvec(stencil_type, structure, 1.0, A, x1, 1.0, y1, ts, vl, ws, print_lp_struct);
        Kokkos::fence();
        print_lp_struct = false;
        double time     = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      // Performance Output
      double matrix_size = 1.0 * ((A.nnz() * (sizeof(Scalar) + sizeof(int)) + A.numRows() * sizeof(int))) / 1024 / 1024;
      double struct_matrix_size = 1.0 * ((A.nnz() * sizeof(Scalar) + A.numRows() * sizeof(int))) / 1024 / 1024;
      double vector_size        = 2.0 * A.numRows() * sizeof(Scalar) / 1024 / 1024;
      double vector_readwrite   = (A.nnz() + A.numCols()) * sizeof(Scalar) / 1024 / 1024;

      double problem_size = matrix_size + vector_size;
      printf(
          "Struct %i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf "
          "%6.3lf ) ( %8.5lf %8.5lf %8.5lf )\n",
          A.nnz(), A.numRows(), A.numCols(), problem_size,
          (struct_matrix_size + vector_readwrite) / ave_time * loop / 1024,
          (struct_matrix_size + vector_readwrite) / max_time / 1024,
          (struct_matrix_size + vector_readwrite) / min_time / 1024, 2.0 * A.nnz() * loop / ave_time / 1e9,
          2.0 * A.nnz() / max_time / 1e9, 2.0 * A.nnz() / min_time / 1e9, ave_time / loop * 1000, max_time * 1000,
          min_time * 1000);
      Kokkos::Profiling::popRegion();
    }

    if (compare) {
      Kokkos::Profiling::pushRegion("Unstructured spmv test");
      // Benchmark
      bool print_lp_unstruct = print_lp;
      double min_time        = 1.0e32;
      double max_time        = 0.0;
      double ave_time        = 0.0;
      for (int i = 0; i < loop; i++) {
        Kokkos::Timer timer;
        matvec(1.0, A, x1, 1.0, y1, tsu, vlu, wsu, print_lp_unstruct);
        Kokkos::fence();
        print_lp_unstruct = false;
        double time       = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      // Performance Output
      double matrix_size = 1.0 * ((A.nnz() * (sizeof(Scalar) + sizeof(int)) + A.numRows() * sizeof(int))) / 1024 / 1024;
      double vector_size = 2.0 * A.numRows() * sizeof(Scalar) / 1024 / 1024;
      double vector_readwrite = (A.nnz() + A.numCols()) * sizeof(Scalar) / 1024 / 1024;

      double problem_size = matrix_size + vector_size;
      printf(
          "Unstr  %i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf "
          "%6.3lf ) ( %8.5lf %8.5lf %8.5lf )\n",
          A.nnz(), A.numRows(), A.numCols(), problem_size, (matrix_size + vector_readwrite) / ave_time * loop / 1024,
          (matrix_size + vector_readwrite) / max_time / 1024, (matrix_size + vector_readwrite) / min_time / 1024,
          2.0 * A.nnz() * loop / ave_time / 1e9, 2.0 * A.nnz() / max_time / 1e9, 2.0 * A.nnz() / min_time / 1e9,
          ave_time / loop * 1000, max_time * 1000, min_time * 1000);
      Kokkos::Profiling::popRegion();
    }

    if (compare_cusparse) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
#ifdef CUSPARSE_VERSION
      KokkosKernels::Experimental::Controls controls;

      cusparseIndexType_t myCusparseOffsetType = CUSPARSE_INDEX_32I;
      cusparseIndexType_t myCusparseEntryType  = CUSPARSE_INDEX_32I;
      cudaDataType myCudaDataType              = CUDA_R_64F;

      /* create matrix */
      cusparseSpMatDescr_t A_cusparse;
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
          &A_cusparse, A.numRows(), A.numCols(), A.nnz(), (void*)A.graph.row_map.data(), (void*)A.graph.entries.data(),
          (void*)A.values.data(), myCusparseOffsetType, myCusparseEntryType, CUSPARSE_INDEX_BASE_ZERO, myCudaDataType));

      /* create lhs and rhs */
      cusparseDnVecDescr_t vecX, vecY;
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecX, x1.extent_int(0), (void*)x1.data(), myCudaDataType));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecY, y1.extent_int(0), (void*)y1.data(), myCudaDataType));

      const double alpha = 1.0, beta = 1.0;
      size_t bufferSize = 0;
      void* dBuffer     = NULL;

// CUSPARSE_MM_ALG_DEFAULT was deprecated in CUDA 11.2.1 a.k.a cuSPARSE 11.4.0
#if CUSPARSE_VERSION >= 11400
      cusparseSpMVAlg_t alg = CUSPARSE_SPMV_ALG_DEFAULT;
#else
      cusparseSpMVAlg_t alg = CUSPARSE_MV_ALG_DEFAULT;
#endif
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV_bufferSize(controls.getCusparseHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                        &alpha, A_cusparse, vecX, &beta, vecY, myCudaDataType, alg,
                                                        &bufferSize));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&dBuffer, bufferSize));

      /* perform SpMV */
      Kokkos::Profiling::pushRegion("cuSparse spmv test");
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;
      for (int i = 0; i < loop; i++) {
        Kokkos::Timer timer;
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV(controls.getCusparseHandle(), CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                                               A_cusparse, vecX, &beta, vecY, myCudaDataType, alg, dBuffer));
        Kokkos::fence();
        double time = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      // Performance Output
      double matrix_size = 1.0 * ((A.nnz() * (sizeof(Scalar) + sizeof(int)) + A.numRows() * sizeof(int))) / 1024 / 1024;
      double vector_size = 2.0 * A.numRows() * sizeof(Scalar) / 1024 / 1024;
      double vector_readwrite = (A.nnz() + A.numCols()) * sizeof(Scalar) / 1024 / 1024;

      double problem_size = matrix_size + vector_size;
      printf(
          "cusp   %i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf "
          "%6.3lf ) ( %8.5lf %8.5lf %8.5lf )\n",
          A.nnz(), A.numRows(), A.numCols(), problem_size, (matrix_size + vector_readwrite) / ave_time * loop / 1024,
          (matrix_size + vector_readwrite) / max_time / 1024, (matrix_size + vector_readwrite) / min_time / 1024,
          2.0 * A.nnz() * loop / ave_time / 1e9, 2.0 * A.nnz() / max_time / 1e9, 2.0 * A.nnz() / min_time / 1e9,
          ave_time / loop * 1000, max_time * 1000, min_time * 1000);
      Kokkos::Profiling::popRegion();

      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dBuffer));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(A_cusparse));
#else
      // The data needs to be reformatted for cusparse before launching the
      // kernel. Step one, extract raw data
      using graph_type        = typename matrix_type::StaticCrsGraphType;
      using cusparse_int_type = typename Kokkos::View<int*, typename graph_type::entries_type::array_layout,
                                                      typename graph_type::device_type>;

      typename graph_type::row_map_type::const_type Arowptr   = A.graph.row_map;
      typename graph_type::entries_type::const_type Acolind   = A.graph.entries;
      typename matrix_type::values_type::non_const_type Avals = A.values;
      cusparse_int_type Arowptr_cusparse("Arowptr", Arowptr.extent(0));
      cusparse_int_type Acolind_cusparse("Acolind", Acolind.extent(0));
      copy_crs_data<graph_type> myCopyFunctor(Arowptr, Acolind, Arowptr_cusparse, Acolind_cusparse);
      myCopyFunctor.doCopy();

      int* rows          = reinterpret_cast<int*>(Arowptr_cusparse.data());
      int* cols          = reinterpret_cast<int*>(Acolind_cusparse.data());
      double* vals       = reinterpret_cast<double*>(Avals.data());
      double* x_cusparse = reinterpret_cast<double*>(x1.data());
      double* y_cusparse = reinterpret_cast<double*>(y1.data());

      /* Get handle to the CUSPARSE context */
      cusparseHandle_t cusparseHandle;
      cusparseStatus_t cusparseStatus;
      cusparseStatus = cusparseCreate(&cusparseHandle);
      if (cusparseStatus != CUSPARSE_STATUS_SUCCESS) {
        printf("Error while initialize the cusparse handle!\n");
      }

      cusparseMatDescr_t descrA = 0;
      cusparseStatus            = cusparseCreateMatDescr(&descrA);
      if (cusparseStatus != CUSPARSE_STATUS_SUCCESS) {
        printf("Error while creating the matrix descriptor!\n");
      }

      cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

      const double alpha = 1.0;
      const double beta  = 1.0;
      Kokkos::Profiling::pushRegion("cuSparse spmv test");
      // Benchmark
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;
      for (int i = 0; i < loop; i++) {
        Kokkos::Timer timer;
        cusparseStatus = cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, static_cast<int>(A.numRows()),
                                        static_cast<int>(A.numCols()), static_cast<int>(A.nnz()), &alpha, descrA, vals,
                                        rows, cols, x_cusparse, &beta, y_cusparse);
        Kokkos::fence();
        double time = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }
      if (cusparseStatus != CUSPARSE_STATUS_SUCCESS) {
        printf("Error during the cusparse SpMV!\n");
      }

      // Performance Output
      double matrix_size = 1.0 * ((A.nnz() * (sizeof(Scalar) + sizeof(int)) + A.numRows() * sizeof(int))) / 1024 / 1024;
      double vector_size = 2.0 * A.numRows() * sizeof(Scalar) / 1024 / 1024;
      double vector_readwrite = (A.nnz() + A.numCols()) * sizeof(Scalar) / 1024 / 1024;

      double problem_size = matrix_size + vector_size;
      printf(
          "cusp   %i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf "
          "%6.3lf ) ( %8.5lf %8.5lf %8.5lf )\n",
          A.nnz(), A.numRows(), A.numCols(), problem_size, (matrix_size + vector_readwrite) / ave_time * loop / 1024,
          (matrix_size + vector_readwrite) / max_time / 1024, (matrix_size + vector_readwrite) / min_time / 1024,
          2.0 * A.nnz() * loop / ave_time / 1e9, 2.0 * A.nnz() / max_time / 1e9, 2.0 * A.nnz() / min_time / 1e9,
          ave_time / loop * 1000, max_time * 1000, min_time * 1000);
      Kokkos::Profiling::popRegion();

      // Clean-up cusparse and cublas contexts
      cusparseDestroy(cusparseHandle);
      // cublasDestroy(cublasHandle);
#endif
#else
      printf(
          "Kokkos was not configure with cusparse, the comparison with "
          "cusparse_matvec is not perfromed!\n");
#endif
    }

    if (check_errors) {
      // Error Check
      Kokkos::View<Scalar*, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> x_check("X_check", A.numCols());
      Kokkos::View<Scalar*, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> y_check("Y_check", A.numRows());
      Kokkos::deep_copy(x_check, h_x);
      KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, 1.0, A, x_check, 1.0, y_check);
      Kokkos::fence();

      Kokkos::deep_copy(h_y, y_check);
      Scalar error = 0;
      for (int rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        for (int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
          error += (h_y_compare(rowIdx) - h_y(rowIdx)) * (h_y_compare(rowIdx) - h_y(rowIdx));
        }
      }

      double total_error = 0;
      total_error += error;

      if (total_error == 0) {
        printf("Kokkos::MultiVector Test: Passed\n");
      } else {
        printf("Kokkos::MultiVector Test: Failed\n");
      }
    }
  }
  Kokkos::finalize();
}
