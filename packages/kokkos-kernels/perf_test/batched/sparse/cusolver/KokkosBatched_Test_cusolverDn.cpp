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

template <typename DeviceType, typename MatrixViewType, typename VectorViewType>
struct Functor_Test_BatchedDenseCuSolve {
  const MatrixViewType _A;
  const VectorViewType _X;
  const VectorViewType _B;

  KOKKOS_INLINE_FUNCTION
  Functor_Test_BatchedDenseCuSolve(const MatrixViewType &A, const VectorViewType &X, const VectorViewType &B)
      : _A(A), _X(X), _B(B) {}

  inline double run() {
    std::string name("KokkosBatched::Test::TeamGESV");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());
    cusolverDnHandle_t handle = NULL;
    cusolverDnCreate(&handle);

    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

    const int batchSize = _A.extent(0);
    const int m         = _A.extent(1);
    const int lda       = m;
    const int ldb       = m;

    MatrixViewType A(" ", batchSize, m, m);
    Kokkos::deep_copy(A, _A);
    Kokkos::deep_copy(_X, _B);

    int *d_infoArray = NULL;
    int *info        = NULL;

    double **d_Aarray = nullptr;
    double **d_Barray = nullptr;

    cudaMalloc(reinterpret_cast<void **>(&d_Aarray), sizeof(double *) * batchSize);
    cudaMalloc(reinterpret_cast<void **>(&d_Barray), sizeof(double *) * batchSize);

    std::vector<double *> Aarray(batchSize, nullptr);
    std::vector<double *> Barray(batchSize, nullptr);
    for (int i = 0; i < batchSize; ++i) {
      Aarray[i] = Kokkos::subview(A, i, Kokkos::ALL, Kokkos::ALL).data();
      Barray[i] = Kokkos::subview(_X, i, Kokkos::ALL).data();
    }

    cudaMemcpyAsync(d_Aarray, Aarray.data(), sizeof(double *) * batchSize, cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_Barray, Barray.data(), sizeof(double *) * batchSize, cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    exec_space().fence();
    timer.reset();
    auto status1 = cusolverDnDpotrfBatched(handle, uplo, m, d_Aarray, lda, d_infoArray, batchSize);
    if (status1 != CUSOLVER_STATUS_SUCCESS)
      std::cout << "Error in cusolverDnDpotrfBatched with batchSize = " << batchSize << " and m = " << m << std::endl;
    cudaDeviceSynchronize();
    auto status2 = cusolverDnDpotrsBatched(handle, uplo, m, 1, d_Aarray, lda, d_Barray, ldb, info, batchSize);
    if (status2 != CUSOLVER_STATUS_SUCCESS) {
      if (status2 == CUSOLVER_STATUS_NOT_INITIALIZED)
        std::cout << "Error in cusolverDnDpotrsBatched with batchSize = " << batchSize << " and m = " << m
                  << " CUSOLVER_STATUS_NOT_INITIALIZED " << std::endl;
      if (status2 == CUSOLVER_STATUS_INVALID_VALUE)
        std::cout << "Error in cusolverDnDpotrsBatched with batchSize = " << batchSize << " and m = " << m
                  << " CUSOLVER_STATUS_INVALID_VALUE " << std::endl;
      if (status2 == CUSOLVER_STATUS_INTERNAL_ERROR)
        std::cout << "Error in cusolverDnDpotrsBatched with batchSize = " << batchSize << " and m = " << m
                  << " CUSOLVER_STATUS_INTERNAL_ERROR " << std::endl;
      cudaDeviceSynchronize();
      exec_space().fence();
      Kokkos::Profiling::popRegion();
      return 1e8;
    }
    cudaDeviceSynchronize();
    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};

int main(int argc, char *argv[]) {
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

    printf(" :::: CusolverDn Testing (N = %d, Blk = %d, vl = %d, n = %d)\n", N, Blk, vector_length, n_rep_1);

    typedef Kokkos::LayoutRight LR;
    typedef Kokkos::LayoutLeft LL;

    using IntView       = Kokkos::View<int *, LR>;
    using AMatrixViewLR = Kokkos::View<double ***, LR>;
    using AMatrixViewLL = Kokkos::View<double ***, LL>;
    using XYTypeLR      = Kokkos::View<double **, LR>;
    using XYTypeLL      = Kokkos::View<double **, LL>;

    AMatrixViewLR aLR("values", N, Blk, Blk);
    AMatrixViewLL aLL("values", N, Blk, Blk);

    XYTypeLR xLR("values", N, Blk);
    XYTypeLR yLR("values", N, Blk);

    XYTypeLL xLL("values", N, Blk);
    XYTypeLL yLL("values", N, Blk);

    if (layout_left) printf(" :::: Testing left layout (team_size = %d)\n", team_size);
    if (layout_right) printf(" :::: Testing right layout (team_size = %d)\n", team_size);

    if (layout_left) {
      readDenseFromMM(name_A, aLL);
      readArrayFromMM(name_B, yLL);
    }
    if (layout_right) {
      readDenseFromMM(name_A, aLR);
      readArrayFromMM(name_B, yLR);
    }

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
          t_spmv = Functor_Test_BatchedDenseCuSolve<exec_space, AMatrixViewLR, XYTypeLR>(aLR, xLR, yLR).run();
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
