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

#include <fstream>

#define KOKKOSKERNELS_DEBUG_LEVEL 0

/// Kokkos headers
#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "Kokkos_Sort.hpp"

//
// #define KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#include <cusparse.h>
#include "cusolverSp.h"
#include "cusolverDn.h"

#include "KokkosKernels_config.h"
#include "KokkosKernels_SparseUtils_cusparse.hpp"

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "KokkosBatched_Test_Sparse_Helper.hpp"

#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_CrsMatrix.hpp"
#include "KokkosBatched_Krylov_Handle.hpp"
#include "KokkosBatched_Gesv.hpp"
#include "KokkosBatched_JacobiPrec.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Dot_Internal.hpp"
#include "KokkosBatched_Spmv_Serial_Impl.hpp"
#include "KokkosBatched_Copy_Decl.hpp"

typedef Kokkos::DefaultExecutionSpace exec_space;
typedef typename exec_space::memory_space memory_space;
typedef Kokkos::DefaultHostExecutionSpace host_space;
typedef typename Kokkos::Device<exec_space, memory_space> device;

template <typename DeviceType, typename MatrixViewType, typename IntView, typename VectorViewType>
struct Functor_Test_SparseCuSolveQR {
  const MatrixViewType _A;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_Test_SparseCuSolveQR(const MatrixViewType &A, const IntView &r, const IntView &c, const VectorViewType &X,
                               const VectorViewType &B)
      : _A(A), _r(r), _c(c), _X(X), _B(B) {}

  inline double run() {
    std::string name("KokkosBatched::Test::TeamGESV");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());
    cusolverSpHandle_t handle = NULL;
    cusolverSpCreate(&handle);

    const size_t N   = _A.extent(0);
    const size_t nnz = _c.extent(0);
    const size_t m   = _r.extent(0) - 1;

    cusparseMatDescr_t descrA = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

    double tol  = 1e-18;
    int reorder = 0;
    int singularity[1];
    exec_space().fence();
    timer.reset();

    for (size_t i = 0; i < N; ++i) {
      auto csrValA = Kokkos::subview(_A, i, Kokkos::ALL).data();
      auto b       = Kokkos::subview(_B, i, Kokkos::ALL).data();
      auto x       = Kokkos::subview(_X, i, Kokkos::ALL).data();

      cusolverSpDcsrlsvqr(handle, m, nnz, descrA, csrValA, _r.data(), _c.data(), b, tol, reorder, x, singularity);
      if (singularity[0] != -1) std::cout << " Error ! " << singularity[0] << " " << m << std::endl;
    }

    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};

template <typename DeviceType, typename MatrixViewType, typename IntView, typename VectorViewType>
struct Functor_Test_Block_SparseCuSolveQR {
  const MatrixViewType _A;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_Test_Block_SparseCuSolveQR(const MatrixViewType &A, const IntView &r, const IntView &c,
                                     const VectorViewType &X, const VectorViewType &B)
      : _A(A), _r(r), _c(c), _X(X), _B(B) {}

  inline double run() {
    std::string name("KokkosBatched::Test::TeamGESV");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());
    cusolverSpHandle_t handle = NULL;
    cusolverSpCreate(&handle);

    const size_t N   = _A.extent(0);
    const size_t nnz = _c.extent(0);
    const size_t m   = _r.extent(0) - 1;

    const size_t block_nnz = N * nnz;
    const size_t block_m   = N * m;

    cusparseMatDescr_t descrA = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

    double tol  = 1e-18;
    int reorder = 0;
    int singularity[1];

    IntView rowOffsets("values", block_m + 1);
    IntView colIndices("values", block_nnz);

    auto rowOffsets_host = Kokkos::create_mirror_view(rowOffsets);
    auto colIndices_host = Kokkos::create_mirror_view(colIndices);
    auto _c_host         = Kokkos::create_mirror_view(_c);
    auto _r_host         = Kokkos::create_mirror_view(_r);

    Kokkos::deep_copy(_c_host, _c);
    Kokkos::deep_copy(_r_host, _r);

    exec_space().fence();

    rowOffsets_host(0) = 0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t row = 0; row < m; ++row) {
        const size_t current_row_index         = i * m + row;
        const size_t row_length                = _r_host(row + 1) - _r_host(row);
        rowOffsets_host(current_row_index + 1) = rowOffsets_host(current_row_index) + row_length;
        for (size_t nnz_row = 0; nnz_row < row_length; ++nnz_row) {
          const size_t current_block_nnz_index     = rowOffsets_host(current_row_index) + nnz_row;
          const size_t current_block_col_index     = _c_host(_r_host(row) + nnz_row) + i * m;
          colIndices_host(current_block_nnz_index) = current_block_col_index;
        }
      }
    }

    Kokkos::deep_copy(rowOffsets, rowOffsets_host);
    Kokkos::deep_copy(colIndices, colIndices_host);

    exec_space().fence();
    timer.reset();

    auto csrValA = _A.data();
    auto b       = _B.data();
    auto x       = _X.data();

    cusolverSpDcsrlsvqr(handle, block_m, block_nnz, descrA, csrValA, rowOffsets.data(), colIndices.data(), b, tol,
                        reorder, x, singularity);

    if (singularity[0] != -1) std::cout << " Error ! " << singularity[0] << " " << m << std::endl;

    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};

template <typename DeviceType, typename MatrixViewType, typename IntView, typename VectorViewType>
struct Functor_Test_SparseCuSolveChol {
  const MatrixViewType _A;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_Test_SparseCuSolveChol(const MatrixViewType &A, const IntView &r, const IntView &c, const VectorViewType &X,
                                 const VectorViewType &B)
      : _A(A), _r(r), _c(c), _X(X), _B(B) {}

  inline double run() {
    std::string name("KokkosBatched::Test::TeamGESV");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());
    cusolverSpHandle_t handle = NULL;
    cusolverSpCreate(&handle);

    const size_t N   = _A.extent(0);
    const size_t nnz = _c.extent(0);
    const size_t m   = _r.extent(0) - 1;

    cusparseMatDescr_t descrA = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

    double tol  = 1e-18;
    int reorder = 0;
    int singularity[1];
    exec_space().fence();
    timer.reset();

    for (size_t i = 0; i < N; ++i) {
      auto csrValA = Kokkos::subview(_A, i, Kokkos::ALL).data();
      auto b       = Kokkos::subview(_B, i, Kokkos::ALL).data();
      auto x       = Kokkos::subview(_X, i, Kokkos::ALL).data();

      cusolverSpDcsrlsvchol(handle, m, nnz, descrA, csrValA, _r.data(), _c.data(), b, tol, reorder, x, singularity);
      if (singularity[0] != -1) std::cout << " Error ! " << singularity[0] << " " << m << std::endl;
    }

    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};

template <typename DeviceType, typename MatrixViewType, typename IntView, typename VectorViewType>
struct Functor_Test_Block_SparseCuSolveChol {
  const MatrixViewType _A;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_Test_Block_SparseCuSolveChol(const MatrixViewType &A, const IntView &r, const IntView &c,
                                       const VectorViewType &X, const VectorViewType &B)
      : _A(A), _r(r), _c(c), _X(X), _B(B) {}

  inline double run() {
    std::string name("KokkosBatched::Test::TeamGESV");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());
    cusolverSpHandle_t handle = NULL;
    cusolverSpCreate(&handle);

    const size_t N   = _A.extent(0);
    const size_t nnz = _c.extent(0);
    const size_t m   = _r.extent(0) - 1;

    const size_t block_nnz = N * nnz;
    const size_t block_m   = N * m;

    cusparseMatDescr_t descrA = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

    double tol  = 1e-18;
    int reorder = 0;
    int singularity[1];

    IntView rowOffsets("values", block_m + 1);
    IntView colIndices("values", block_nnz);

    auto rowOffsets_host = Kokkos::create_mirror_view(rowOffsets);
    auto colIndices_host = Kokkos::create_mirror_view(colIndices);
    auto _c_host         = Kokkos::create_mirror_view(_r);
    auto _r_host         = Kokkos::create_mirror_view(_c);

    Kokkos::deep_copy(_c_host, _c);
    Kokkos::deep_copy(_r_host, _r);

    exec_space().fence();

    rowOffsets_host(0) = 0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t row = 0; row < m; ++row) {
        const size_t current_row_index         = i * m + row;
        const size_t row_length                = _r_host(row + 1) - _r_host(row);
        rowOffsets_host(current_row_index + 1) = rowOffsets_host(current_row_index) + row_length;
        for (size_t nnz_row = 0; nnz_row < row_length; ++nnz_row) {
          const size_t current_block_nnz_index     = rowOffsets_host(current_row_index) + nnz_row;
          const size_t current_block_col_index     = _c_host(_r_host(row) + nnz_row) + i * m;
          colIndices_host(current_block_nnz_index) = current_block_col_index;
        }
      }
    }

    Kokkos::deep_copy(rowOffsets, rowOffsets_host);
    Kokkos::deep_copy(colIndices, colIndices_host);

    exec_space().fence();
    timer.reset();

    auto csrValA = _A.data();
    auto b       = _B.data();
    auto x       = _X.data();

    cusolverSpDcsrlsvchol(handle, block_m, block_nnz, descrA, csrValA, rowOffsets.data(), colIndices.data(), b, tol,
                          reorder, x, singularity);
    if (singularity[0] != -1) std::cout << " Error ! " << singularity[0] << " " << m << std::endl;

    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};

int main(int argc, char *argv[]) {
  std::cout << " start " << std::endl;
  Kokkos::initialize(argc, argv);
  {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
    cudaProfilerStop();
#endif
    Kokkos::print_configuration(std::cout);

    // typedef Kokkos::ArithTraits<value_type> ats;

    ///
    /// input arguments parsing
    ///
    int n_rep_1       = 10;  // # of repetitions
    int team_size     = 8;
    int n_impl        = 1;
    bool layout_left  = true;
    bool layout_right = false;
    int vector_length = 8;

    std::string name_A = "A.mm";
    std::string name_B = "B.mm";

    std::string name_timer = "timers";
    std::string name_X     = "X";

    std::vector<int> impls;
    for (int i = 1; i < argc; ++i) {
      const std::string &token = argv[i];
      if (token == std::string("-A")) name_A = argv[++i];
      if (token == std::string("-B")) name_B = argv[++i];
      if (token == std::string("-X")) name_X = argv[++i];
      if (token == std::string("-timers")) name_timer = argv[++i];
      if (token == std::string("-team_size")) team_size = std::atoi(argv[++i]);
      if (token == std::string("-vector_length")) vector_length = std::atoi(argv[++i]);
      if (token == std::string("-n_implementations")) n_impl = std::atoi(argv[++i]);
      if (token == std::string("-implementation")) impls.push_back(std::atoi(argv[++i]));
      if (token == std::string("-l")) {
        layout_left  = true;
        layout_right = false;
      }
      if (token == std::string("-r")) {
        layout_left  = false;
        layout_right = true;
      }
    }

    int N, Blk, nnz, ncols;

    readSizesFromMM(name_A, Blk, ncols, nnz, N);

    std::cout << "n = " << Blk << ", N = " << N << ", team_size = " << team_size
              << ", vector_length = " << vector_length << std::endl;

    if (impls.size() == 0)
      for (int i = 0; i < n_impl; ++i) impls.push_back(i);

    // V100 L2 cache 6MB per core
    constexpr size_t LLC_CAPACITY = 80 * 6 * 1024 * 1024;
    KokkosBatched::Flush<LLC_CAPACITY, exec_space> flush;

    printf(" :::: CusolverSp Testing (N = %d, Blk = %d, vl = %d, n = %d)\n", N, Blk, vector_length, n_rep_1);

    typedef Kokkos::LayoutRight LR;
    typedef Kokkos::LayoutLeft LL;

    using IntView            = Kokkos::View<int *, LR>;
    using AMatrixValueViewLR = Kokkos::View<double **, LR>;
    using AMatrixValueViewLL = Kokkos::View<double **, LL>;
    using XYTypeLR           = Kokkos::View<double **, LR>;
    using XYTypeLL           = Kokkos::View<double **, LL>;

    IntView rowOffsets("values", Blk + 1);
    IntView colIndices("values", nnz);
    AMatrixValueViewLR valuesLR("values", N, nnz);
    AMatrixValueViewLL valuesLL("values", N, nnz);

    XYTypeLR xLR("values", N, Blk);
    XYTypeLR yLR("values", N, Blk);

    XYTypeLL xLL("values", N, Blk);
    XYTypeLL yLL("values", N, Blk);

    if (layout_left) printf(" :::: Testing left layout (team_size = %d)\n", team_size);
    if (layout_right) printf(" :::: Testing right layout (team_size = %d)\n", team_size);

    if (layout_left) {
      readCRSFromMM(name_A, valuesLL, rowOffsets, colIndices);
      readArrayFromMM(name_B, yLL);
    }
    if (layout_right) {
      readCRSFromMM(name_A, valuesLR, rowOffsets, colIndices);
      readArrayFromMM(name_B, yLR);
    }

    printf(" :::: GO\n");
    for (auto i_impl : impls) {
      std::vector<double> timers;

      double t_spmv = 0;
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      exec_space().fence();
      Kokkos::deep_copy(xLL, 0.0);
      Kokkos::deep_copy(xLR, 0.0);

      exec_space().fence();

      if (i_impl == 0) {
        if (layout_right) {
          t_spmv = Functor_Test_SparseCuSolveQR<exec_space, AMatrixValueViewLR, IntView, XYTypeLR>(valuesLR, rowOffsets,
                                                                                                   colIndices, xLR, yLR)
                       .run();
        }
      }
      if (i_impl == 1) {
        if (layout_right) {
          t_spmv = Functor_Test_SparseCuSolveChol<exec_space, AMatrixValueViewLR, IntView, XYTypeLR>(
                       valuesLR, rowOffsets, colIndices, xLR, yLR)
                       .run();
        }
      }
      if (i_impl == 2) {
        if (layout_right) {
          t_spmv = Functor_Test_Block_SparseCuSolveQR<exec_space, AMatrixValueViewLR, IntView, XYTypeLR>(
                       valuesLR, rowOffsets, colIndices, xLR, yLR)
                       .run();
        }
      }
      if (i_impl == 3) {
        if (layout_right) {
          t_spmv = Functor_Test_Block_SparseCuSolveChol<exec_space, AMatrixValueViewLR, IntView, XYTypeLR>(
                       valuesLR, rowOffsets, colIndices, xLR, yLR)
                       .run();
        }
      }
      exec_space().fence();

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif

      timers.push_back(t_spmv);

      {
        std::ofstream myfile;
        std::string name;
        if (layout_left) name = name_timer + "_" + std::to_string(i_impl) + "_left.txt";
        if (layout_right) name = name_timer + "_" + std::to_string(i_impl) + "_right.txt";

        myfile.open(name);

        for (size_t i = 0; i < timers.size(); ++i) myfile << timers[i] << " ";

        myfile << std::endl;

        myfile.close();
      }

      double average_time = 0.;

      for (size_t i = 0; i < timers.size(); ++i) average_time += timers[i] / timers.size();

      if (layout_left) printf("Left layout: Implementation %d: solve time = %f\n", i_impl, average_time);
      if (layout_right) printf("Right layout: Implementation %d: solve time = %f\n", i_impl, average_time);

      if (layout_left) {
        writeArrayToMM(name_X + std::to_string(i_impl) + "_l.mm", xLL);
      }
      if (layout_right) {
        writeArrayToMM(name_X + std::to_string(i_impl) + "_r.mm", xLR);
      }
    }
  }
  Kokkos::finalize();

  return 0;
}

#else
int main() { return 0; }
#endif
