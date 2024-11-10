// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_EXAMPLE_EXTERNALINTERFACE_HPP__
#define __TACHO_EXAMPLE_EXTERNALINTERFACE_HPP__
#include "Tacho_config.h"

#if defined(TACHO_USE_INT_INT)

#include "Tacho.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_CommandLineParser.hpp"
#include "Tacho_Solver.hpp"

namespace tacho {

enum TACHO_PARAM_INDICES {
  USEDEFAULTSOLVERPARAMETERS,
  VERBOSITY,
  SMALLPROBLEMTHRESHOLDSIZE,
  SOLUTION_METHOD,
  TASKING_OPTION_BLOCKSIZE,
  TASKING_OPTION_PANELSIZE,
  TASKING_OPTION_MAXNUMSUPERBLOCKS,
  LEVELSET_OPTION_SCHEDULING,
  LEVELSET_OPTION_DEVICE_LEVEL_CUT,
  LEVELSET_OPTION_DEVICE_FACTOR_THRES,
  LEVELSET_OPTION_DEVICE_SOLVE_THRES,
  LEVELSET_OPTION_NSTREAMS,
  LEVELSET_OPTION_VARIANT,
  INDEX_LENGTH
};

using device_type = typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type;
using exec_space = typename device_type::execution_space;

using host_device_type = typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
using host_space = typename host_device_type::execution_space;

using ViewVectorType = Kokkos::View<double *, device_type>;
using ViewVectorTypeInt = Kokkos::View<int *, device_type>;

template <class SX> class tachoSolver {
public:
  typedef Tacho::Solver<SX, device_type> solver_type;

  typedef Tacho::ordinal_type ordinal_type;
  typedef Tacho::size_type size_type;

  typedef typename solver_type::value_type value_type;
  typedef typename solver_type::ordinal_type_array ordinal_type_array;
  typedef typename solver_type::size_type_array size_type_array;
  typedef typename solver_type::value_type_array value_type_array;
  typedef typename solver_type::value_type_matrix value_type_matrix;

  // host typedefs
  typedef typename solver_type::ordinal_type_array_host ordinal_type_array_host;
  typedef typename solver_type::size_type_array_host size_type_array_host;
  typedef typename solver_type::value_type_array_host value_type_array_host;
  typedef typename solver_type::value_type_matrix_host value_type_matrix_host;

  tachoSolver(const int *solverParams) : m_numRows(0), m_Solver() { setSolverParameters(solverParams); }

  ~tachoSolver() {
    for (auto &solver : m_SolverArray) {
      solver.release();
    }
    m_Solver.release();
  }

  int Initialize(int numRows,
                 /// with TACHO_ENABLE_INT_INT, size_type is "int"
                 int *rowBegin, int *columns, SX *values, int numGraphRows = 0, int *graphRowBegin = nullptr,
                 int *graphColumns = nullptr, int *graphWeights = nullptr) {
    m_numRows = numRows;
    if (m_numRows == 0)
      return 0;

    const int numTerms = rowBegin[numRows];
    size_type_array_host ap_host((size_type *)rowBegin, numRows + 1);
    ordinal_type_array_host aj_host((ordinal_type *)columns, numTerms);

    size_type_array_host graph_ap_host;
    ordinal_type_array_host graph_aj_host;
    ordinal_type_array_host graph_aw_host;

    if (numGraphRows > 0) {
      if (graphRowBegin == nullptr || graphColumns == nullptr || graphWeights == nullptr) {
        std::cout << "ExternalInterface::Error, with non-zero numGraphRows, graph pointers should not be nullptr\n";
        std::logic_error("Error: one of graph pointers is nullptr");
      } else {
        const size_type nnz_graph = graph_ap_host(numGraphRows);
        graph_ap_host = size_type_array_host((size_type *)graphRowBegin, numGraphRows + 1);
        graph_aj_host = ordinal_type_array_host((ordinal_type *)graphColumns, nnz_graph);
        graph_aw_host = ordinal_type_array_host((ordinal_type *)graphWeights, nnz_graph);
      }
    }

#if defined(KOKKOS_ENABLE_CUDA)
    /// transfer A into device
    value_type_array ax(Kokkos::ViewAllocateWithoutInitializing("ax"), numTerms);
    value_type_array_host ax_host(values, numTerms);
    Kokkos::deep_copy(ax, ax_host);
#else
    /// wrap pointer on host
    value_type_array ax(values, numTerms);
#endif

    Kokkos::Timer timer;
    {
      timer.reset();
      if (numGraphRows > 0) {
        m_Solver.analyze(numRows, ap_host, aj_host, numGraphRows, graph_ap_host, graph_aj_host, graph_aw_host);
      } else {
        m_Solver.analyze(numRows, ap_host, aj_host);
      }
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: analyze time " << t << std::endl;
    }

    {
      timer.reset();
      m_Solver.initialize();
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: initialize time " << t << std::endl;
    }

    /// I recommend to separate factorization from the solver initialization
    {
      timer.reset();
      m_Solver.factorize(ax);
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: factorize time " << t << std::endl;
    }

    return 0;
  }

  int Initialize(int numSolvers, int numRows,
                 /// with TACHO_ENABLE_INT_INT, size_type is "int"
                 int *rowBegin, int *columns, SX *values, int numGraphRows = 0, int *graphRowBegin = nullptr,
                 int *graphColumns = nullptr, int *graphWeights = nullptr) {
    m_numRows = numRows;
    if (m_numRows == 0)
      return 0;
    if (numSolvers > 0) {
      /// this is okay
    } else {
      std::cout << "ExternalInterface::Error, nSolver is not a positive number\n";
      std::logic_error("Error: nSolver must be a positive number");
    }

    const int numTerms = rowBegin[numRows];
    size_type_array_host ap_host((size_type *)rowBegin, numRows + 1);
    ordinal_type_array_host aj_host((ordinal_type *)columns, numTerms);

    size_type_array_host graph_ap_host;
    ordinal_type_array_host graph_aj_host;
    ordinal_type_array_host graph_aw_host;

    if (numGraphRows > 0) {
      if (graphRowBegin == nullptr || graphColumns == nullptr || graphWeights == nullptr) {
        std::cout << "ExternalInterface::Error, with non-zero numGraphRows, graph pointers should not be nullptr\n";
        std::logic_error("Error: one of graph pointers is nullptr");
      } else {
        const size_type nnz_graph = graph_ap_host(numGraphRows);
        graph_ap_host = size_type_array_host((size_type *)graphRowBegin, numGraphRows + 1);
        graph_aj_host = ordinal_type_array_host((ordinal_type *)graphColumns, nnz_graph);
        graph_aw_host = ordinal_type_array_host((ordinal_type *)graphWeights, nnz_graph);
      }
    }

#if defined(KOKKOS_ENABLE_CUDA)
    /// transfer A into device
    value_type_array ax(Kokkos::ViewAllocateWithoutInitializing("ax"), numSolvers * numTerms);
    value_type_array_host ax_host(values, numSolvers * numTerms);
    Kokkos::deep_copy(ax, ax_host);
#else
    /// wrap pointer on host
    value_type_array ax(values, numSolvers * numTerms);
#endif

    Kokkos::Timer timer;
    {
      timer.reset();
      /// m_Solver holds symbolic factorization to be shared with array solver
      if (numGraphRows > 0) {
        m_Solver.analyze(numRows, ap_host, aj_host, numGraphRows, graph_ap_host, graph_aj_host, graph_aw_host);
      } else {
        m_Solver.analyze(numRows, ap_host, aj_host);
      }
      /// array solver soft copy
      m_SolverArray.resize(numSolvers);
      for (auto &solver : m_SolverArray) {
        /// duplicate perform soft copy of symbolic factors and
        /// nullify numeric tools object
        solver = m_Solver.duplicate();
      }
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: analyze time " << t << std::endl;
    }

    /// the solver objects in the array still need to be initalized
    /// to allocate required work space
    {
      timer.reset();
      for (auto &solver : m_SolverArray) {
        solver.initialize();
      }
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: initialize time (" << numSolvers << ") " << t << std::endl;
    }

    /// I recommend to separate factorization from the solver initialization
    {
      timer.reset();
      for (ordinal_type i = 0, iend = m_SolverArray.size(); i < iend; ++i) {
        m_SolverArray[i].factorize(value_type_array(ax.data() + numTerms * i, numTerms));
      }
      const double t = timer.seconds();
      std::cout << "ExternalInterface:: factorize time (" << numSolvers << ") " << t << std::endl;
    }

    return 0;
  }

  void exportSupernodes(std::vector<int> &supernodes) {
    const auto supernodes_device = m_Solver.getSupernodes();
    auto supernodes_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), supernodes_device);
    {
      const int n = supernodes_host.extent(0);
      supernodes.resize(n);
      std::copy(supernodes_host.data(), supernodes_host.data() + n, supernodes.data());
    }
  }

  void exportUpperTriangularFactorsToCrsMatrix(const int iSolver, std::vector<int> &rowBeginU,
                                               std::vector<int> &columnsU, std::vector<SX> &valuesU,
                                               std::vector<int> &perm) {
    solver_type *solver = nullptr;
    if (iSolver < 0) {
      solver = &m_Solver;
    } else {
      const int solver_array_size(m_SolverArray.size());
      if (iSolver < solver_array_size) {
        solver = &m_SolverArray[iSolver];
      } else {
        std::cout << "ExternalInterface::Error, non-zero iSolver (" << iSolver
                  << ") is selected where m_SolverArray is sized by (" << solver_array_size << ")\n";
        std::logic_error("Error: iSolver is out of range in m_SolverArray");
      }
    }

    {
      const auto perm_device = solver->getPermutationVector();
      auto perm_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), perm_device);
      {
        const int n = perm_host.extent(0);
        perm.resize(n);
        std::copy(perm_host.data(), perm_host.data() + n, perm.data());
      }
    }

    {
      typename solver_type::crs_matrix_type A;
      solver->exportFactorsToCrsMatrix(A);

      typename solver_type::crs_matrix_type_host A_host;
      A_host.createMirror(A);
      A_host.copy(A);
      {
        const auto ap = A_host.RowPtr();
        const int n = ap.extent(0);
        rowBeginU.resize(n);
        std::copy(ap.data(), ap.data() + n, rowBeginU.data());
      }
      {
        const auto aj = A_host.Cols();
        const int n = aj.extent(0);
        columnsU.resize(n);
        std::copy(aj.data(), aj.data() + n, columnsU.data());
      }
      {
        const auto ax = A_host.Values();
        const int n = ax.extent(0);
        valuesU.resize(n);
        std::copy(ax.data(), ax.data() + n, valuesU.data());
      }
    }
  }

  void exportUpperTriangularFactorsToCrsMatrix(std::vector<int> &rowBeginU, std::vector<int> &columnsU,
                                               std::vector<SX> &valuesU, std::vector<int> &perm) {
    exportUpperTriangularFactorsToCrsMatrix(-1, rowBeginU, columnsU, valuesU, perm);
  }

  void MySolve(int NRHS, value_type_matrix &b, value_type_matrix &x) {
    if (m_numRows == 0)
      return;

    const int m = m_numRows;
    if (static_cast<int>(m_TempRhs.extent(0)) < m || static_cast<int>(m_TempRhs.extent(1)) < NRHS)
      m_TempRhs = value_type_matrix("temp rhs", m, NRHS);

    if (m_SolverArray.empty()) {
      /// solve for a single instance m_Solver
      m_Solver.solve(x, b, m_TempRhs);
    } else {
      /// solve multiple with m_SolverArray
      using range_type = Kokkos::pair<int, int>;
      const auto range_b = range_type(0, NRHS);
      for (int i = 0, iend = m_SolverArray.size(); i < iend; ++i) {
        const auto range_a = range_type(m * i, m * i + m);
        m_SolverArray[i].solve(Kokkos::subview(x, range_a, range_b), Kokkos::subview(b, range_a, range_b), m_TempRhs);
      }
    }
  }

private:
  int m_numRows;
  solver_type m_Solver;
  std::vector<solver_type> m_SolverArray;
  value_type_matrix m_TempRhs;

  void setSolverParameters(const int *solverParams) {
    if (solverParams[USEDEFAULTSOLVERPARAMETERS])
      return;
    // common options
    m_Solver.setVerbose(solverParams[VERBOSITY]);
    m_Solver.setSmallProblemThresholdsize(solverParams[SMALLPROBLEMTHRESHOLDSIZE]);

    // solution method
    m_Solver.setSolutionMethod(solverParams[SOLUTION_METHOD]);

    // tasking options
    m_Solver.setBlocksize(solverParams[TASKING_OPTION_BLOCKSIZE]);
    m_Solver.setPanelsize(solverParams[TASKING_OPTION_PANELSIZE]);
    m_Solver.setMaxNumberOfSuperblocks(solverParams[TASKING_OPTION_MAXNUMSUPERBLOCKS]);

    // levelset options
    m_Solver.setLevelSetScheduling(solverParams[LEVELSET_OPTION_SCHEDULING]);
    m_Solver.setLevelSetOptionDeviceLevelCut(solverParams[LEVELSET_OPTION_DEVICE_LEVEL_CUT]);
    m_Solver.setLevelSetOptionDeviceFunctionThreshold(solverParams[LEVELSET_OPTION_DEVICE_FACTOR_THRES],
                                                      solverParams[LEVELSET_OPTION_DEVICE_SOLVE_THRES]);
    m_Solver.setLevelSetOptionNumStreams(solverParams[LEVELSET_OPTION_NSTREAMS]);
    m_Solver.setLevelSetOptionAlgorithmVariant(solverParams[LEVELSET_OPTION_VARIANT]);
  }
};

} // namespace tacho

#endif
#endif // INTERFACE_TACHO_SOLVER_H
