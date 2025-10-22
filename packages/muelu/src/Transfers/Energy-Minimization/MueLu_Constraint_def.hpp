// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CONSTRAINT_DEF_HPP
#define MUELU_CONSTRAINT_DEF_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include "KokkosBatched_Copy_Internal.hpp"
#include "Teuchos_Assert.hpp"
#include "Xpetra_Access.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_MatrixMatrix.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_ProductOperator.hpp"
#include "MueLu_Constraint_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"
#include "KokkosBatched_ApplyQ_Decl.hpp"
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"
#include "Kokkos_DualView.hpp"
#include "Kokkos_Pair.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "KokkosBatched_QR_Decl.hpp"
#include "KokkosBatched_QR_Serial_Impl.hpp"
#include "KokkosBatched_QR_TeamVector_Impl.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Team_Impl.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Trsv_TeamVector_Impl.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"
#include "KokkosBlas1_set.hpp"

namespace MueLu {

template <class LocalGraph, class LocalVector>
class MinSpmMV {
 private:
  using local_ordinal_type = typename LocalVector::value_type;

  LocalGraph lclGraph;
  LocalVector lhs;
  LocalVector rhs;

  const local_ordinal_type MAX_VAL = Kokkos::ArithTraits<local_ordinal_type>::max();

 public:
  MinSpmMV(LocalGraph lclGraph_, LocalVector lhs_, LocalVector rhs_)
    : lclGraph(lclGraph_)
    , lhs(lhs_)
    , rhs(rhs_) {}

  KOKKOS_INLINE_FUNCTION
  void init(bool& dst) {
    dst = false;
  }

  KOKKOS_INLINE_FUNCTION
  void join(bool& dst, const bool& src) {
    dst = dst || src;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type i, bool& changed) const {
    local_ordinal_type val = MAX_VAL;
    for (local_ordinal_type jj = lclGraph.row_map(i); jj < (local_ordinal_type)lclGraph.row_map(i + 1); ++jj) {
      auto j = lclGraph.entries(jj);
      val    = Kokkos::min(val, lhs(j));
    }
    auto prev = rhs(i);
    rhs(i)    = val;
    changed   = changed || (prev != val);
  }
};

template <class LocalGraph, class LocalVector>
class MinSpmMVT {
 private:
  using local_ordinal_type = typename LocalVector::value_type;

  LocalGraph lclGraph;
  LocalVector lhs;
  LocalVector rhs;

  const local_ordinal_type MAX_VAL = Kokkos::ArithTraits<local_ordinal_type>::max();

 public:
  MinSpmMVT(LocalGraph lclGraph_, LocalVector lhs_, LocalVector rhs_)
    : lclGraph(lclGraph_)
    , lhs(lhs_)
    , rhs(rhs_) {
    Kokkos::deep_copy(rhs, MAX_VAL);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const local_ordinal_type i) const {
    local_ordinal_type val = lhs(i);
    for (local_ordinal_type jj = lclGraph.row_map(i); jj < lclGraph.row_map(i + 1); ++jj) {
      auto j = lclGraph.entries(jj);
      Kokkos::atomic_min(&rhs(j), val);
    }
  }
};

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type
FindBlocks(RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& X) {
  using execution_space = typename Node::execution_space;
  using range_type      = Kokkos::RangePolicy<LocalOrdinal, execution_space>;

  auto numRows = X->getRowMap()->getLocalNumElements();

  Kokkos::View<LocalOrdinal*> blockIds("blockIds", numRows);
  Kokkos::View<LocalOrdinal*> blockIds2("blockIds2", numRows);

  auto lclGraph = X->getLocalGraphDevice();

  auto range = range_type(0, numRows);

  Kokkos::parallel_for(
      "init", range,
      KOKKOS_LAMBDA(const LocalOrdinal i) {
        blockIds(i) = i;
      });

  bool changed    = true;
  bool resultsIn2 = false;
  while (changed) {
    changed = false;
    if (!resultsIn2) {
      MinSpmMV functor(lclGraph, blockIds, blockIds2);
      Kokkos::parallel_reduce("backward_pass", range, functor, changed);
      resultsIn2 = true;
    } else {
      MinSpmMV functor(lclGraph, blockIds2, blockIds);
      Kokkos::parallel_reduce("backward_pass", range, functor, changed);
      resultsIn2 = false;
    }
  }
  if (resultsIn2)
    Kokkos::deep_copy(blockIds, blockIds2);

  Kokkos::View<bool*> blockStatus("", numRows);
  Kokkos::View<LocalOrdinal*> newBlockIds("", numRows);
  Kokkos::parallel_for(
      "init", range,
      KOKKOS_LAMBDA(LocalOrdinal i) {
        blockStatus(blockIds(i)) = true;
      });
  Kokkos::fence();

  LocalOrdinal numBlocks = 0;
  Kokkos::parallel_scan(
      "init", range,
      KOKKOS_LAMBDA(LocalOrdinal i, LocalOrdinal & count, const bool final) {
        if (final)
          newBlockIds(i) = count;
        if (blockStatus(i))
          ++count;
      },
      numBlocks);

  // Build graph with block info
  using graph_type = typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type;
  typename graph_type::row_map_type::non_const_type rowptr("", numBlocks + 1);
  typename graph_type::entries_type::non_const_type indices("", numRows);

  Kokkos::parallel_for(
      "", range, KOKKOS_LAMBDA(const LocalOrdinal i) {
    auto blockId = newBlockIds(blockIds(i));
    if (blockId+2 < numBlocks+1)
      ++rowptr(blockId + 2); });
  Kokkos::fence();

  auto block_range = range_type(0, numBlocks);

  Kokkos::parallel_scan(
      "", block_range, KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& sum, const bool final) {
    sum += rowptr(i+1);
    if (final) {
      rowptr(i+1) = sum;
    } });

  Kokkos::parallel_for(
      "", range, KOKKOS_LAMBDA(const LocalOrdinal i) {
    auto blockId = newBlockIds(blockIds(i));
    indices(rowptr(blockId + 1)) = i;
    ++rowptr(blockId + 1); });

  graph_type blocks(indices, rowptr);

  return blocks;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class BlockInverseFunctor {
 public:
  using CrsGraph          = typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using CrsMatrix         = typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_graph_type  = typename CrsGraph::local_graph_type;
  using local_matrix_type = typename CrsMatrix::local_matrix_type;
  using scalar_type       = typename local_matrix_type::value_type;
  using ATS               = Kokkos::ArithTraits<scalar_type>;

  using shared_matrix = Kokkos::View<scalar_type**, typename Node::execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_vector = Kokkos::View<scalar_type*, typename Node::execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;

  BlockInverseFunctor(local_matrix_type A_, local_graph_type blocks_, LocalOrdinal maxBlocksize_, local_matrix_type invA_, const Kokkos::View<bool*>& singular_)
    : A(A_)
    , blocks(blocks_)
    , maxBlocksize(maxBlocksize_)
    , invA(invA_)
    , singular(singular_) {}

 private:
  const scalar_type zero = ATS::zero();
  const scalar_type one  = ATS::one();

  local_matrix_type A;
  local_graph_type blocks;
  LocalOrdinal maxBlocksize;
  local_matrix_type invA;
  const Kokkos::View<bool*>& singular;

 public:
  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<typename Node::execution_space>::member_type& thread) const {
    using member_type = typename Kokkos::TeamPolicy<typename Node::execution_space>::member_type;

    auto blockId   = thread.league_rank();
    auto blockRow  = blocks.rowConst(blockId);
    auto blockSize = blockRow.length;

    shared_matrix lclA(thread.team_shmem(), blockSize, blockSize);
    shared_matrix lclInvA(thread.team_shmem(), blockSize, blockSize);

    const bool PseudoInverse = (!(singular.extent(0) == 0)) && singular(blockId);

    // Initialize lclA
    // If PseudoInverse, we shift the constant mode.
    KokkosBlas::TeamSet<member_type>::invoke(thread, PseudoInverse ? one : zero, lclA);
    thread.team_barrier();

    // extract block from A
    for (LocalOrdinal ii = 0; ii < blockSize; ++ii) {
      auto i   = blockRow.colidx(ii);
      auto row = A.rowConst(i);
      for (LocalOrdinal jj = 0; jj < row.length; ++jj) {
        auto j = row.colidx(jj);
        auto d = row.value(jj);
        for (LocalOrdinal kk = 0; kk < blockSize; ++kk)
          if (blockRow.colidx(kk) == j) {
            lclA(ii, kk) += d;
            break;
          }
      }
    }

    // LU
    {
      // LU factorization: lclA = L * U
      KokkosBatched::TeamLU<member_type, KokkosBlas::Algo::QR::Unblocked>::invoke(thread, lclA);

      // set lclInvA to identity matrix
      KokkosBatched::TeamSetIdentity<member_type>::invoke(thread, lclInvA);
      thread.team_barrier();

      // // lclInvA = L^{-1}*lclInvA
      for (LocalOrdinal j = 0; j < blockSize; ++j)
        KokkosBatched::TeamVectorTrsv<member_type, KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(thread, one, lclA, Kokkos::subview(lclInvA, Kokkos::ALL(), j));
      thread.team_barrier();

      // // lclInvA = R^{-1}*lclInvA
      for (LocalOrdinal j = 0; j < blockSize; ++j)
        KokkosBatched::TeamVectorTrsv<member_type, KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(thread, one, lclA, Kokkos::subview(lclInvA, Kokkos::ALL(), j));
      thread.team_barrier();
    }

    // The QR in Kokkos Kernels is broken. Once it gets fixed we can use it and remove LU.
    //
    // QR
    // {

    //   shared_vector tau(thread.team_shmem(), blockSize);
    //   shared_vector work(thread.team_shmem(), blockSize);

    //   // QR factorization: lclA = Q * R
    //   KokkosBatched::TeamVectorQR<member_type, KokkosBlas::Algo::QR::Unblocked>::invoke(thread, lclA, tau, work);

    //   // set lclInvA to identity matrix
    //   KokkosBatched::TeamSetIdentity<member_type>::invoke(thread, lclInvA);
    //   thread.team_barrier();

    //   // lclInvA = Q^T*lclInvA
    //   KokkosBatched::TeamVectorApplyQ<member_type, KokkosBatched::Side::Left, KokkosBlas::Trans::Transpose, KokkosBlas::Algo::ApplyQ::Unblocked>::invoke(thread, lclA, tau, lclInvA, work);
    //   thread.team_barrier();

    //   // lclInvA = R^{-1}*lclInvA
    //   for (LocalOrdinal j = 0; j < blockSize; ++j)
    //     KokkosBatched::TeamVectorTrsv<member_type, KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(thread, one, lclA, Kokkos::subview(lclInvA, Kokkos::ALL(), j));
    //   thread.team_barrier();
    // }

    if (PseudoInverse) {
      // Multiply with projection that removes constant vector

      // Set up projection matrix
      for (LocalOrdinal ii = 0; ii < blockSize; ++ii) {
        for (LocalOrdinal kk = 0; kk < blockSize; ++kk) {
          if (ii == kk) {
            lclA(ii, kk) = one - one / (scalar_type)blockSize;
          } else {
            lclA(ii, kk) = -one / (scalar_type)blockSize;
          }
        }
      }
      // Copy lclInvA to temp
      shared_matrix temp(thread.team_shmem(), blockSize, blockSize);
      KokkosBatched::TeamCopy<member_type, KokkosBatched::Trans::NoTranspose>::invoke(thread, lclInvA, temp);
      thread.team_barrier();

      // lclInvA = proj * lclInvA
      KokkosBatched::TeamGemm<member_type, KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Gemm::Unblocked>::invoke(thread, one, lclA, temp, zero, lclInvA);
      thread.team_barrier();
    }

    // write inverse of block to invA
    for (LocalOrdinal ii = 0; ii < blockSize; ++ii) {
      auto i   = blockRow.colidx(ii);
      auto row = invA.row(i);
      for (LocalOrdinal jj = 0; jj < row.length; ++jj) {
        auto j = row.colidx(jj);
        for (LocalOrdinal kk = 0; kk < blockSize; ++kk)
          if (blockRow.colidx(kk) == j) {
            row.value(jj) = lclInvA(ii, kk);
            break;
          }
      }
    }
  }

  // amount of shared memory
  size_t team_shmem_size(int /* team_size */) const {
    return 3 * shared_matrix::shmem_size(maxBlocksize, maxBlocksize) + 2 * shared_vector::shmem_size(maxBlocksize);
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrepareLeastSquaresSolveBelos(const bool detect_singular_blocks) {
  Monitor m(*this, "PrepareLeastSquaresSolveBelos");

  TEUCHOS_ASSERT(!detect_singular_blocks);

  problem_ = rcp(new Belos::LinearProblem<Scalar, MV, OP>());

  std::vector<RCP<Operator>> ops      = {X_, X_};
  std::vector<Teuchos::ETransp> modes = {Teuchos::NO_TRANS, Teuchos::TRANS};
  RCP<Operator> XXt                   = rcp(new ProductOperator(ops, modes));
  auto belosXXt                       = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(XXt));

  problem_->setOperator(belosXXt);
  problem_->setLabel("LeastSquares");

  auto belosList = rcp(new Teuchos::ParameterList());
  belosList->set("Implicit Residual Scaling", "None");
  belosList->set("Convergence Tolerance", 1e-16);
  auto out = GetMueLuOStream();
  belosList->set("Output Stream", out->getOStream());
  // belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  // belosList->set("Output Frequency", 1);
  // belosList->set("Output Style", Belos::Brief);

  Belos::SolverFactory<Scalar, MV, OP> solverFactory;
  solver_ = solverFactory.create("Pseudo Block CG", belosList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
allocateBlockDiagonalMatrix(RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map,
                            const typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type blocks) {
  using graph_type  = typename Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type;
  using matrix_type = typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;

  auto numRows   = map->getLocalNumElements();
  auto numBlocks = blocks.numRows();
  typename graph_type::row_map_type::non_const_type rowptr("rowptr", numRows + 1);

  LocalOrdinal nnz = 0;
  Kokkos::parallel_reduce(
      "", numBlocks, KOKKOS_LAMBDA(const LocalOrdinal blockId, LocalOrdinal& count) {
    auto blockRow = blocks.rowConst(blockId);
    auto blockSize = blockRow.length;

    for (LocalOrdinal k = 0; k < blockSize; ++k) {
      auto rowId = blockRow.colidx(k);
      if ((decltype(numRows))rowId+2<numRows+1)
        rowptr(rowId+2) += blockSize;
      count += blockSize;
    } }, nnz);

  LocalOrdinal temp;
  Kokkos::parallel_scan(
      "", numRows, KOKKOS_LAMBDA(const LocalOrdinal rowId, LocalOrdinal& sum, const bool is_final) {
    sum += rowptr(rowId+1);
    if (is_final) {
      rowptr(rowId+1) = sum;
    } }, temp);

  typename graph_type::entries_type::non_const_type indices("lclInvXXt_indices", nnz);

  Kokkos::parallel_for(
      "", numBlocks, KOKKOS_LAMBDA(const LocalOrdinal blockId) {
    auto blockRow = blocks.rowConst(blockId);
    auto blockSize = blockRow.length;

    for (LocalOrdinal k = 0; k < blockSize; ++k) {
      auto rowId = blockRow.colidx(k);
      for (LocalOrdinal jj = 0; jj < blockSize; ++jj) {
        auto j = blockRow.colidx(jj);
        indices(rowptr(rowId + 1)) = j;
        ++rowptr(rowId + 1);
      }
    } });

  auto lclInvXXtGraph = graph_type(indices, rowptr);
  typename matrix_type::values_type::non_const_type values("lclInvXXt_values", nnz);
  auto lclInvXXt = matrix_type("lclInvXXt", numRows, values, lclInvXXtGraph);
  return Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclInvXXt, map, map, map, map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrepareLeastSquaresSolveDirect(const bool detect_singular_blocks) {
  Monitor m(*this, "PrepareLeastSquaresSolveDirect");

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> XXt;
  {
    SubMonitor m2(*this, "XXt");
    XXt = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*X_, false, *X_, true, XXt, GetOStream(Runtime0), true, true);
  }

  auto XXtgraph = XXt->getCrsGraph();
  auto blocks   = FindBlocks(XXtgraph);
  invXXt_       = allocateBlockDiagonalMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(XXt->getRowMap(), blocks);

  LocalOrdinal numBlocks    = blocks.numRows();
  LocalOrdinal maxBlocksize = invXXt_->getLocalMaxNumRowEntries();

  // If we pass a view of size 0 to the functor, all blocks are assumed to be non-singular.
  Kokkos::View<bool*> block_is_singular;
  if (detect_singular_blocks) {
    block_is_singular = Kokkos::View<bool*>("block_is_singular", numBlocks);
    Kokkos::deep_copy(block_is_singular, true);
  }

  {
    SubMonitor m2(*this, "inversion");

    BlockInverseFunctor<Scalar, LocalOrdinal, GlobalOrdinal, Node> functor(XXt->getLocalMatrixDevice(), blocks, maxBlocksize, invXXt_->getLocalMatrixDevice(), block_is_singular);
    Kokkos::parallel_for("", Kokkos::TeamPolicy<typename Node::execution_space>(numBlocks, 1), functor);
  }

  if (IsPrint(Statistics0)) {
    // print some stats

    auto comm = invXXt_->getRowMap()->getComm();
    GlobalOrdinal globalNumBlocks;
    MueLu_sumAll(comm, (GlobalOrdinal)numBlocks, globalNumBlocks);

    GetOStream(Statistics0) << "Least-squares problem:\n maximum block size: " << invXXt_->getGlobalMaxNumRowEntries() << "\n Number of blocks: " << globalNumBlocks << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrepareLeastSquaresSolveDirect(const Kokkos::View<bool*> block_is_singular) {
  Monitor m(*this, "PrepareLeastSquaresSolveDirect");

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> XXt;
  {
    SubMonitor m2(*this, "XXt");
    XXt = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*X_, false, *X_, true, XXt, GetOStream(Runtime0), true, true);
  }

  auto XXtgraph = XXt->getCrsGraph();
  auto blocks   = FindBlocks(XXtgraph);
  invXXt_       = allocateBlockDiagonalMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(XXt->getRowMap(), blocks);

  LocalOrdinal numBlocks    = blocks.numRows();
  LocalOrdinal maxBlocksize = invXXt_->getLocalMaxNumRowEntries();

  {
    SubMonitor m2(*this, "inversion");

    BlockInverseFunctor<Scalar, LocalOrdinal, GlobalOrdinal, Node> functor(XXt->getLocalMatrixDevice(), blocks, maxBlocksize, invXXt_->getLocalMatrixDevice(), block_is_singular);
    Kokkos::parallel_for("", Kokkos::TeamPolicy<typename Node::execution_space>(numBlocks, 1), functor);
  }

  if (IsPrint(Statistics0)) {
    // print some stats

    auto comm = invXXt_->getRowMap()->getComm();
    GlobalOrdinal globalNumBlocks;
    MueLu_sumAll(comm, (GlobalOrdinal)numBlocks, globalNumBlocks);

    GetOStream(Statistics0) << "Least-squares problem:\n maximum block size: " << invXXt_->getGlobalMaxNumRowEntries() << "\n Number of blocks: " << globalNumBlocks << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrepareLeastSquaresSolve(const std::string& solverType, const bool detect_singular_blocks) {
  if (solverType == "Belos")
    PrepareLeastSquaresSolveBelos(detect_singular_blocks);
  else if (solverType == "direct")
    PrepareLeastSquaresSolveDirect(detect_singular_blocks);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "solverType must be one of (Belos|direct), not \"" << solverType << "\".");
  solverType_ = solverType;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrepareLeastSquaresSolve(const std::string& solverType, const Kokkos::View<bool*> block_is_singular) {
  if (solverType == "Belos") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "solverType must be direct, not \"" << solverType << "\".");
  }
  else if (solverType == "direct")
    PrepareLeastSquaresSolveDirect(block_is_singular);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "solverType must be one of (Belos|direct), not \"" << solverType << "\".");
  solverType_ = solverType;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::LeastSquaresSolveBelos(const MultiVector& B, MultiVector& C) const {
  // Solve (X * X^T) * C = B

  problem_->setLHS(rcpFromRef(C));
  problem_->setRHS(rcpFromRef(B));
  TEUCHOS_ASSERT(problem_->setProblem());

  solver_->setProblem(problem_);
  solver_->solve();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::LeastSquaresSolveDirect(const MultiVector& B, MultiVector& C) const {
  // Solve (X * X^T) * C = B
  invXXt_->apply(B, C);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::LeastSquaresSolve(const MultiVector& B, MultiVector& C) const {
  if (solverType_ == "Belos")
    LeastSquaresSolveBelos(B, C);
  else if (solverType_ == "direct")
    LeastSquaresSolveDirect(B, C);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "solverType must be one of (Belos|direct), not \"" << solverType_ << "\".");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& P,
                                                                  MultiVector& Projected,
                                                                  Teuchos::ETransp mode,
                                                                  Scalar alpha,
                                                                  Scalar beta) const {
  const auto one  = Teuchos::ScalarTraits<Scalar>::one();
  const auto zero = Teuchos::ScalarTraits<Scalar>::zero();

  TEUCHOS_ASSERT(mode == Teuchos::NO_TRANS);
  TEUCHOS_ASSERT(alpha == one);
  TEUCHOS_ASSERT(beta == zero);

  // Projected = P - X^T * (X * X^T)^{-1} * X * P
  Projected = P;
  X_->apply(P, *temp1_, Teuchos::NO_TRANS);
  LeastSquaresSolve(*temp1_, *temp2_);
  X_->apply(*temp2_, *temp3_, Teuchos::TRANS);
  Projected.update(-one, *temp3_, one);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AssignMatrixEntriesToVector(const Matrix& P,
                                                                                        const RCP<const CrsGraph>& pattern,
                                                                                        MultiVector& vecP) const {
  auto lclMat     = P.getLocalMatrixDevice();
  auto lclPattern = pattern->getLocalGraphDevice();
  auto lclVec     = vecP.getLocalViewDevice(Xpetra::Access::OverwriteAll);
  TEUCHOS_ASSERT(lclPattern.entries.extent(0) == lclVec.extent(0));
  Kokkos::deep_copy(lclVec, 0.);

  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
  Kokkos::parallel_for(
      "MueLu::Constraint::filter", range_type(0, lclPattern.numRows()), KOKKOS_LAMBDA(const size_t i) {
        auto row_mat = lclMat.rowConst(i);

        if (row_mat.length == 0)
          return;

        for (size_t jj = lclPattern.row_map(i); jj < lclPattern.row_map(i + 1); ++jj) {
          auto offset     = lclPattern.entries(jj);
          LocalOrdinal kk = 0;
          while ((kk + 1 < row_mat.length) && (row_mat.colidx(kk) != offset))
            ++kk;
          if (row_mat.colidx(kk) == offset) {
            lclVec(jj, 0) = row_mat.value(kk);
          }
        }
      });
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AssignMatrixEntriesToVector(const Matrix& P,
                                                                                        MultiVector& vecP) const {
  AssignMatrixEntriesToVector(P, GetPattern(), vecP);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetMatrixWithEntriesFromVector(MultiVector& vecP) const {
  auto Ppattern = GetPattern();
  RCP<Matrix> P = MatrixFactory::Build(Ppattern);
  {
    auto lclMat = P->getLocalMatrixDevice();
    auto lclVec = vecP.getLocalViewDevice(Xpetra::Access::ReadOnly);
    TEUCHOS_ASSERT(lclMat.values.extent(0) == lclVec.extent(0));
    Kokkos::deep_copy(lclMat.values, Kokkos::subview(lclVec, Kokkos::ALL(), 0));
  }
  P->fillComplete(Ppattern->getDomainMap(), Ppattern->getRowMap());
  return P;
}

}  // namespace MueLu

#endif  // ifndef MUELU_CONSTRAINT_DEF_HPP
