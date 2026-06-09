// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_
#define MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>

#include "Kokkos_Sort.hpp"
#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_QR_Decl.hpp"
#include "KokkosBatched_ApplyQ_Decl.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Util.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_InverseApproximationFactory_decl.hpp"

#if KOKKOSKERNELS_VERSION < 50102
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"
#endif

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  using Magnitude                   = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  validParamList->set<RCP<const FactoryBase>>("A", NoFactory::getRCP(), "Matrix to build the approximate inverse on.\n");

  validParamList->set<std::string>("inverse: approximation type", "diagonal", "Method used to approximate the inverse.");
  validParamList->set<Magnitude>("inverse: drop tolerance", 0.0, "Values below this threshold  are dropped from the matrix (or fixed if diagonal fixing is active).");
  validParamList->set<bool>("inverse: fixing", false, "Keep diagonal and fix small entries with 1.0");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  using STS       = Teuchos::ScalarTraits<SC>;
  const SC one    = STS::one();
  using Magnitude = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  const ParameterList& pL = GetParameterList();
  const bool fixing       = pL.get<bool>("inverse: fixing");

  // check which approximation type to use
  const std::string method = pL.get<std::string>("inverse: approximation type");
  TEUCHOS_TEST_FOR_EXCEPTION(method != "diagonal" && method != "lumping" && method != "sparseapproxinverse", Exceptions::RuntimeError,
                             "MueLu::InverseApproximationFactory::Build: Approximation type can be 'diagonal' or 'lumping' or "
                             "'sparseapproxinverse'.");

  RCP<Matrix> A            = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
  const bool isBlocked     = (bA == Teuchos::null ? false : true);

  // if blocked operator is used, defaults to A(0,0)
  if (isBlocked) A = bA->getMatrix(0, 0);

  const Magnitude tol = pL.get<Magnitude>("inverse: drop tolerance");
  RCP<Matrix> Ainv    = Teuchos::null;

  if (method == "diagonal") {
    const auto diag = VectorFactory::Build(A->getRangeMap(), true);
    A->getLocalDiagCopy(*diag);
    const RCP<const Vector> D = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, tol, one));
    Ainv                      = MatrixFactory::Build(D);
  } else if (method == "lumping") {
    const auto diag           = Utilities::GetLumpedMatrixDiagonal(*A);
    const RCP<const Vector> D = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, tol, one));
    Ainv                      = MatrixFactory::Build(D);
  } else if (method == "sparseapproxinverse") {
    RCP<CrsGraph> sparsityPattern = Utilities::GetThresholdedGraph(A, tol);
    if (IsPrint(Statistics1)) {
      sparsityPattern->computeGlobalConstants();
      GetOStream(Statistics1) << "NNZ Graph(A): " << A->getCrsGraph()->getGlobalNumEntries() << " , NNZ Tresholded Graph(A): " << sparsityPattern->getGlobalNumEntries() << std::endl;
    }
    RCP<Matrix> pAinv = GetSparseInverse(A, sparsityPattern);
    Ainv              = Utilities::GetThresholdedMatrix(pAinv, tol, fixing);
    if (IsPrint(Statistics1)) {
      rcp_const_cast<CrsGraph>(Ainv->getCrsGraph())->computeGlobalConstants();
      GetOStream(Statistics1) << "NNZ Ainv: " << pAinv->getGlobalNumEntries() << ", NNZ Tresholded Ainv (parameter: " << tol << "): " << Ainv->getGlobalNumEntries() << std::endl;
    }
  }

  GetOStream(Statistics1) << "Approximate inverse calculated by: " << method << "." << std::endl;
  GetOStream(Statistics1) << "Ainv has " << Ainv->getGlobalNumRows() << "x" << Ainv->getGlobalNumCols() << " rows and columns." << std::endl;

  Set(currentLevel, "Ainv", Ainv);
}

#if KOKKOSKERNELS_VERSION >= 50102

template <class local_matrix_type>
class LocalSPAIFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using execution_space    = typename local_matrix_type::execution_space;
  using impl_scalar_type   = typename KokkosKernels::ArithTraits<scalar_type>::val_type;
  using impl_ATS           = KokkosKernels::ArithTraits<impl_scalar_type>;

 public:
  using shared_matrix    = Kokkos::View<impl_scalar_type**, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_vector    = Kokkos::View<impl_scalar_type*, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_lo_vector = Kokkos::View<local_ordinal_type*, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;

 private:
  const local_matrix_type lclA;
  local_matrix_type lclAinv;
  const local_ordinal_type maxUniqueColEntries;
  const int scratchLevel;

 public:
  LocalSPAIFunctor(const local_matrix_type& lclA_, local_matrix_type& lclAinv_, local_ordinal_type maxUniqueColEntries_, int scratchLevel_)
    : lclA(lclA_)
    , lclAinv(lclAinv_)
    , maxUniqueColEntries(maxUniqueColEntries_)
    , scratchLevel(scratchLevel_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) const {
    auto rlid    = thread.league_rank();
    auto rowAinv = lclAinv.row(rlid);

    // Loop over entries in row rlid of Ainv and collect all of A's column indices.
    shared_lo_vector column_indices(thread.team_scratch(scratchLevel), maxUniqueColEntries);
    local_ordinal_type numColEntries = 0;
    for (local_ordinal_type ii = 0; ii < rowAinv.length; ++ii) {
      auto i    = rowAinv.colidx(ii);
      auto rowA = lclA.rowConst(i);
      for (local_ordinal_type jj = 0; jj < rowA.length; ++jj) {
        auto j                        = rowA.colidx(jj);
        column_indices(numColEntries) = j;
        ++numColEntries;
      }
    }

    // Get merged list of column indices.
    local_ordinal_type numUniqeColEntries = 0;
    local_ordinal_type diagOffset         = 0;
    {
      // Sort
      Kokkos::Experimental::sort_team(thread, Kokkos::subview(column_indices, Kokkos::make_pair(0, numColEntries)));
      // Merge
      if (numColEntries > 0)
        ++numUniqeColEntries;
      local_ordinal_type pos = 0;
      for (local_ordinal_type m = 1; m < numColEntries; ++m) {
        if (column_indices(pos) != column_indices(m)) {
          column_indices(pos + 1) = column_indices(m);
          ++pos;
          ++numUniqeColEntries;
          if (column_indices(pos) == rlid)
            diagOffset = pos;
        }
      }
    }

    // Extract local part of A into a dense view.
    shared_matrix localA(thread.team_scratch(scratchLevel), numUniqeColEntries, rowAinv.length);
    KokkosBlas::SerialSet::invoke(impl_ATS::zero(), localA);

    // Now fill localA.
    for (local_ordinal_type ii = 0; ii < rowAinv.length; ++ii) {
      auto i    = rowAinv.colidx(ii);
      auto rowA = lclA.rowConst(i);
      for (local_ordinal_type jj = 0; jj < rowA.length; ++jj) {
        auto j = rowA.colidx(jj);
        auto v = rowA.value(jj);
        // Determine local index.
        // Sequential search might not be a great idea.
        for (local_ordinal_type m = 0; m < numUniqeColEntries; ++m) {
          if (column_indices(m) == j) {
            localA(m, ii) = v;
            break;
          }
        }
      }
    }

    shared_matrix ek(thread.team_scratch(scratchLevel), numUniqeColEntries, 1);
    // set to zero, set diagonal entry to one
    for (local_ordinal_type i = 0; i < numUniqeColEntries; ++i) {
      ek(i, 0) = (i == diagOffset) ? impl_ATS::one() : impl_ATS::zero();
    }

    // QR solve
    shared_vector tau(thread.team_scratch(scratchLevel), rowAinv.length);
    shared_vector work(thread.team_scratch(scratchLevel), numUniqeColEntries);
    // factorize localA = Q*R in-place
    KokkosBatched::SerialQR<KokkosBatched::Algo::QR::Unblocked>::invoke(localA, tau, work);
    // ek := Q^T ek
    KokkosBatched::SerialApplyQ<KokkosBatched::Side::Left, KokkosBatched::Trans::Transpose, KokkosBatched::Algo::ApplyQ::Unblocked>::invoke(localA, tau, ek, work);
    // ek[:rowLength] := R^{-1} ek[:rowLength]
    auto sub_A  = Kokkos::subview(localA, Kokkos::make_pair(0, rowAinv.length), Kokkos::ALL());
    auto sub_ek = Kokkos::subview(ek, Kokkos::make_pair(0, rowAinv.length), 0);
    KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(impl_ATS::one(), sub_A, sub_ek);

    // Set entries of Ainv.
    for (local_ordinal_type i = 0; i < rowAinv.length; ++i) {
      rowAinv.value(i) = sub_ek(i);
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSparseInverse(const RCP<Matrix>& Aorg, const RCP<const CrsGraph>& sparsityPattern) const {
  using execution_space = typename Node::execution_space;

  // construct the inverse matrix with the given sparsity pattern
  RCP<Matrix> Ainv = MatrixFactory::Build(sparsityPattern);
  Ainv->resumeFill();

  // gather missing rows from other procs to generate an overlapping map
  RCP<Import> rowImport = ImportFactory::Build(sparsityPattern->getRowMap(), sparsityPattern->getColMap());
  RCP<Matrix> A         = MatrixFactory::Build(Aorg, *rowImport);

  auto maxRowEntriesAinv   = Ainv->getLocalMaxNumRowEntries();
  auto maxRowEntriesA      = A->getLocalMaxNumRowEntries();
  auto maxUniqueColEntries = maxRowEntriesAinv * maxRowEntriesA;
  {
    auto lclA    = A->getLocalMatrixDevice();
    auto lclAinv = Ainv->getLocalMatrixDevice();

    Kokkos::TeamPolicy<execution_space> policy(lclAinv.numRows(), 1);

    using spai_functor_type = LocalSPAIFunctor<decltype(lclAinv)>;
    using shared_matrix     = typename spai_functor_type::shared_matrix;
    using shared_vector     = typename spai_functor_type::shared_vector;
    using shared_lo_vector  = typename spai_functor_type::shared_lo_vector;

    int size = shared_matrix::shmem_size(maxUniqueColEntries, maxRowEntriesAinv) + shared_matrix::shmem_size(maxUniqueColEntries, 1) + shared_vector::shmem_size(3 * maxUniqueColEntries) + shared_vector::shmem_size(maxRowEntriesAinv) + shared_lo_vector::shmem_size(maxUniqueColEntries);

    int scratchLevel = -1;
    if (size < policy.scratch_size_max(/*level=*/(int)0)) {
      policy.set_scratch_size(/*level=*/(int)0, Kokkos::PerTeam(size));
      scratchLevel = 0;
    } else if (size < policy.scratch_size_max(/*level=*/(int)1)) {
      policy.set_scratch_size(/*level=*/(int)1, Kokkos::PerTeam(size));
      scratchLevel = 1;
    } else
      throw Exceptions::RuntimeError("Neither L0 scratch memory (max size " + std::to_string(policy.scratch_size_max((int)0)) +
                                     "), nor L1 scratch memory (max size " + std::to_string(policy.scratch_size_max((int)1)) +
                                     ") is large enough for requested allocation of size " + std::to_string(size));

    LocalSPAIFunctor spaiFunctor(lclA, lclAinv, maxUniqueColEntries, scratchLevel);

    Kokkos::parallel_for("MueLu::InverseFactory::LocalSpai", policy, spaiFunctor);
  }

  Ainv->fillComplete();

  return Ainv;
}

#else

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSparseInverse(const RCP<Matrix>& Aorg, const RCP<const CrsGraph>& sparsityPattern) const {
  // construct the inverse matrix with the given sparsity pattern
  RCP<Matrix> Ainv = MatrixFactory::Build(sparsityPattern);
  Ainv->resumeFill();

  // gather missing rows from other procs to generate an overlapping map
  RCP<Import> rowImport = ImportFactory::Build(sparsityPattern->getRowMap(), sparsityPattern->getColMap());
  RCP<Matrix> A         = MatrixFactory::Build(Aorg, *rowImport);

  // loop over all rows of the inverse sparsity pattern (this can be done in parallel)
  for (size_t k = 0; k < sparsityPattern->getLocalNumRows(); k++) {
    // 1. get column indices Ik of local row k
    ArrayView<const LO> Ik;
    sparsityPattern->getLocalRowView(k, Ik);

    // 2. get all local A(Ik,:) rows
    Array<ArrayView<const LO>> J(Ik.size());
    Array<ArrayView<const SC>> Ak(Ik.size());
    Array<LO> Jk;
    for (LO i = 0; i < Ik.size(); i++) {
      A->getLocalRowView(Ik[i], J[i], Ak[i]);
      for (LO j = 0; j < J[i].size(); j++)
        Jk.append(J[i][j]);
    }
    // set of unique column indices Jk
    std::sort(Jk.begin(), Jk.end());
    Jk.erase(std::unique(Jk.begin(), Jk.end()), Jk.end());
    // create map
    std::map<LO, LO> G;
    for (LO i = 0; i < Jk.size(); i++) G.insert(std::pair<LO, LO>(Jk[i], i));

    // 3. merge rows together
    Teuchos::SerialDenseMatrix<LO, SC> localA(Jk.size(), Ik.size(), true);
    for (LO i = 0; i < Ik.size(); i++) {
      for (LO j = 0; j < J[i].size(); j++) {
        localA(G.at(J[i][j]), i) = Ak[i][j];
      }
    }

    // 4. get direction-vector
    // diagonal needs an entry!
    Teuchos::SerialDenseVector<LO, SC> ek(Jk.size(), true);
    ek[std::find(Jk.begin(), Jk.end(), k) - Jk.begin()] = Teuchos::ScalarTraits<Scalar>::one();
    ;

    // 5. solve linear system for x
    Teuchos::SerialDenseVector<LO, SC> localX(Ik.size());
    Teuchos::SerialQRDenseSolver<LO, SC> qrSolver;
    qrSolver.setMatrix(Teuchos::rcp(&localA, false));
    qrSolver.setVectors(Teuchos::rcp(&localX, false), Teuchos::rcp(&ek, false));
    const int err = qrSolver.solve();
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::GetSparseInverse: Error in serial QR solve.");

    // 6. set calculated row into Ainv
    ArrayView<const SC> Mk(localX.values(), localX.length());
    Ainv->replaceLocalValues(k, Ik, Mk);
  }
  Ainv->fillComplete();

  return Ainv;
}

#endif

}  // namespace MueLu

#endif /* MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_ */
