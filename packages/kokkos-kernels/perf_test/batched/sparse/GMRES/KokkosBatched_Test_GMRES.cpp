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

#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Test_Sparse_Helper.hpp"

#include "KokkosBatched_Spmv.hpp"
#include "KokkosBatched_CrsMatrix.hpp"
#include "KokkosBatched_Krylov_Handle.hpp"
#include "KokkosBatched_GMRES.hpp"
#include "KokkosBatched_JacobiPrec.hpp"

typedef Kokkos::DefaultExecutionSpace exec_space;
typedef typename exec_space::memory_space memory_space;
typedef Kokkos::DefaultHostExecutionSpace host_space;
typedef typename Kokkos::Device<exec_space, memory_space> device;

#include "Functor_TestBatchedTeamVectorGMRES_1.hpp"
#include "Functor_TestBatchedTeamVectorGMRES_2.hpp"
#include "Functor_TestBatchedTeamVectorGMRES_3.hpp"

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
    int n_rep_1              = 10;    // # of repetitions
    int n_rep_2              = 1000;  // # of repetitions
    int team_size            = 8;
    int n_impl               = 1;
    int n_iterations         = 10;
    double tol               = 1e-8;
    bool layout_left         = true;
    bool layout_right        = false;
    bool use_preconditioner  = false;
    bool monitor_convergence = false;
    int vector_length        = 8;
    int N_team_potential     = 1;
    int ortho_strategy       = 0;
    int other_level          = 0;
    int arnoldi_level        = 0;

    std::string name_A = "A.mm";
    std::string name_B = "B.mm";

    std::string name_timer = "timers";
    std::string name_X     = "X";
    std::string name_conv  = "res";

    std::vector<int> impls;
    for (int i = 1; i < argc; ++i) {
      const std::string &token = argv[i];
      if (token == std::string("--help") || token == std::string("-h")) {
        std::cout << "Kokkos Batched GMRES performance test options:" << std::endl
                  << "-A                :  Filename of the input batched matrix." << std::endl
                  << "-B                :  Filename of the input batched right-hand "
                     "side."
                  << std::endl
                  << "-X                :  Filename of the output batched solution." << std::endl
                  << "-res              :  Filename of the output residual history." << std::endl
                  << "-timers           :  Filename of the output timers." << std::endl
                  << "-ortho_strategy   :  Select the orthogonalization strategy." << std::endl
                  << "-arnoldi_level    :  Select the scratch pad level (if used) "
                     "where Arnoldi vectors are stored."
                  << std::endl
                  << "-other_level      :  Select the scratch pad level (if used) "
                     "where everything except the Arnoldi vectors are stored."
                  << std::endl
                  << "-n1               :  Number of repetitions of the experience." << std::endl
                  << "-n2               :  Number of the kernel calls inside one "
                     "experience."
                  << std::endl
                  << "-team_size        :  Used team size." << std::endl
                  << "-n_implementations:  Number of implementations to use: test "
                     "all "
                     "implementations [0, specified number -1]."
                  << std::endl
                  << "-implementation   :  Specify only one implementation at a time." << std::endl
                  << "                     Note: implementation 0 : does not use "
                     "scratch pad."
                  << std::endl
                  << "                     Note: implementation 1 : use scratch pad "
                     "for the graph and for the diagonal entries of the matrices."
                  << std::endl
                  << "                     Note: implementation 2 : use scratch pad "
                     "for the graph and for the diagonal entries of the matrices and "
                     "for the temporary variable but not for the Arnoldi vectors."
                  << std::endl
                  << "-l                :  Specify left layout." << std::endl
                  << "-r                :  Specify right layout." << std::endl
                  << "-P                :  Specify if a Jacobi preconditioner is "
                     "used."
                  << std::endl
                  << "-C                :  Specify if the convergence is monitored." << std::endl
                  << "-N_team           :  Specify the number of systems per team." << std::endl
                  << "-vector_length    :  Specify the vector length." << std::endl
                  << std::endl;
        return 0;
      }
      if (token == std::string("-A")) name_A = argv[++i];
      if (token == std::string("-B")) name_B = argv[++i];
      if (token == std::string("-X")) name_X = argv[++i];
      if (token == std::string("-res")) name_conv = argv[++i];
      if (token == std::string("-timers")) name_timer = argv[++i];
      if (token == std::string("-ortho_strategy")) ortho_strategy = std::atoi(argv[++i]);
      if (token == std::string("-arnoldi_level")) arnoldi_level = std::atoi(argv[++i]);
      if (token == std::string("-other_level")) other_level = std::atoi(argv[++i]);
      if (token == std::string("-n1")) n_rep_1 = std::atoi(argv[++i]);
      if (token == std::string("-n2")) n_rep_2 = std::atoi(argv[++i]);
      if (token == std::string("-n_iterations")) n_iterations = std::atoi(argv[++i]);
      if (token == std::string("-tol")) tol = std::stod(argv[++i]);
      if (token == std::string("-team_size")) team_size = std::atoi(argv[++i]);
      if (token == std::string("-N_team")) N_team_potential = std::atoi(argv[++i]);
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
      if (token == std::string("-P")) use_preconditioner = true;
      if (token == std::string("-C")) monitor_convergence = true;
    }

    int N, Blk, nnz, ncols;

    readSizesFromMM(name_A, Blk, ncols, nnz, N);

    std::cout << "N_team_potential = " << N_team_potential << ", n = " << Blk << ", N = " << N
              << ", team_size = " << team_size << ", vector_length = " << vector_length << std::endl;

    if (impls.size() == 0)
      for (int i = 0; i < n_impl; ++i) impls.push_back(i);

    // V100 L2 cache 6MB per core
    constexpr size_t LLC_CAPACITY = 80 * 6 * 1024 * 1024;
    KokkosBatched::Flush<LLC_CAPACITY, exec_space> flush;

    printf(" :::: GMRES Testing (N = %d, Blk = %d, nnz = %d, vl = %d, n = %d)\n", N, Blk, nnz, vector_length, n_rep_1);

    typedef Kokkos::LayoutRight LR;
    typedef Kokkos::LayoutLeft LL;

    using IntView            = Kokkos::View<int *, LR>;
    using AMatrixValueViewLR = Kokkos::View<double **, LR>;
    using AMatrixValueViewLL = Kokkos::View<double **, LL>;
    using XYTypeLR           = Kokkos::View<double **, LR>;
    using XYTypeLL           = Kokkos::View<double **, LL>;

    using alphaViewType = Kokkos::View<double *>;
    alphaViewType alphaV("alpha", N);
    alphaViewType betaV("alpha", N);

    IntView rowOffsets("values", Blk + 1);
    IntView colIndices("values", nnz);
    AMatrixValueViewLR valuesLR("values", N, nnz);
    AMatrixValueViewLL valuesLL("values", N, nnz);

    AMatrixValueViewLR diagLR("values", N, Blk);
    AMatrixValueViewLL diagLL("values", N, Blk);

    XYTypeLR xLR("values", N, Blk);
    XYTypeLR yLR("values", N, Blk);

    XYTypeLL xLL("values", N, Blk);
    XYTypeLL yLL("values", N, Blk);

    if (layout_left) printf(" :::: Testing left layout (team_size = %d)\n", team_size);
    if (layout_right) printf(" :::: Testing right layout (team_size = %d)\n", team_size);

    if (layout_left) {
      readCRSFromMM(name_A, valuesLL, rowOffsets, colIndices);
      readArrayFromMM(name_B, yLL);
      if (use_preconditioner) getInvDiagFromCRS(valuesLL, rowOffsets, colIndices, diagLL);
    }
    if (layout_right) {
      readCRSFromMM(name_A, valuesLR, rowOffsets, colIndices);
      readArrayFromMM(name_B, yLR);
      if (use_preconditioner) getInvDiagFromCRS(valuesLR, rowOffsets, colIndices, diagLR);
    }

    for (auto i_impl : impls) {
      std::vector<double> timers;

      int n_skip = 2;

      int N_team = N_team_potential;

      using ScalarType = typename AMatrixValueViewLL::non_const_value_type;
      using Layout     = typename AMatrixValueViewLL::array_layout;
      using EXSP       = typename AMatrixValueViewLL::execution_space;

      using MagnitudeType = typename Kokkos::ArithTraits<ScalarType>::mag_type;

      using Norm2DViewType   = Kokkos::View<MagnitudeType **, Layout, EXSP>;
      using Scalar3DViewType = Kokkos::View<ScalarType ***, Layout, EXSP>;
      using IntViewType      = Kokkos::View<int *, Layout, EXSP>;

      using KrylovHandleType = KokkosBatched::KrylovHandle<Norm2DViewType, IntViewType, Scalar3DViewType>;
      KrylovHandleType handle(N, N_team, n_iterations, true);
      handle.Arnoldi_view = Scalar3DViewType("", N, n_iterations, Blk + n_iterations + 3);
      // handle.tmp_view = typename KrylovHandleType::TemporaryViewType(
      //  "", N, Blk + n_iterations + 3);

      handle.set_max_iteration(n_iterations);
      handle.set_tolerance(tol);
      handle.set_ortho_strategy(ortho_strategy);
      handle.set_scratch_pad_level(0);
      handle.set_compute_last_residual(false);

      for (int i_rep = 0; i_rep < n_rep_1 + n_skip; ++i_rep) {
        double t_spmv = 0;
        for (int j_rep = 0; j_rep < n_rep_2; ++j_rep) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
          cudaProfilerStart();
#endif
          exec_space().fence();
          Kokkos::deep_copy(xLL, 0.0);
          Kokkos::deep_copy(xLR, 0.0);
          flush.run();
          exec_space().fence();

          if (i_impl == 0 && layout_left) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_1<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, true>(
                            valuesLL, diagLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_1<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, false>(
                            valuesLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          if (i_impl == 1 && layout_left) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_2<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, true>(
                            valuesLL, diagLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_2<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, false>(
                            valuesLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          if (i_impl == 2 && layout_left) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_3<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, true>(
                            valuesLL, diagLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_3<exec_space, AMatrixValueViewLL, IntView, XYTypeLL,
                                                             KrylovHandleType, false>(
                            valuesLL, rowOffsets, colIndices, xLL, yLL, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          if (i_impl == 0 && layout_right) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_1<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, true>(
                            valuesLR, diagLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_1<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, false>(
                            valuesLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          if (i_impl == 1 && layout_right) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_2<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, true>(
                            valuesLR, diagLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_2<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, false>(
                            valuesLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          if (i_impl == 2 && layout_right) {
            if (use_preconditioner)
              t_spmv += Functor_TestBatchedTeamVectorGMRES_3<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, true>(
                            valuesLR, diagLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length,
                            n_iterations, tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
            else
              t_spmv += Functor_TestBatchedTeamVectorGMRES_3<exec_space, AMatrixValueViewLR, IntView, XYTypeLR,
                                                             KrylovHandleType, false>(
                            valuesLR, rowOffsets, colIndices, xLR, yLR, N_team, team_size, vector_length, n_iterations,
                            tol, ortho_strategy, arnoldi_level, other_level, handle)
                            .run();
          }
          exec_space().fence();

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
          cudaProfilerStop();
#endif
        }
        if (i_rep > n_skip) timers.push_back(t_spmv / n_rep_2);
      }

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
      if (monitor_convergence) {
        writeArrayToMM(name_conv + std::to_string(i_impl) + ".mm", handle.residual_norms);
      }
    }
  }
  Kokkos::finalize();

  return 0;
}
