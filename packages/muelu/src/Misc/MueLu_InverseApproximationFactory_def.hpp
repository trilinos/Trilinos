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
  TEUCHOS_TEST_FOR_EXCEPTION(method != "diagonal" && method != "lumping" && method != "sparseapproxinverse" && method != "factoredsparseapproxinverse", Exceptions::RuntimeError,
                             "MueLu::InverseApproximationFactory::Build: Approximation type can be 'diagonal' or 'lumping' or "
                             "'sparseapproxinverse' or 'factoredsparseapproxinverse'.");

  RCP<Matrix> A            = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
  const bool isBlocked     = (bA == Teuchos::null ? false : true);

  // if blocked operator is used, defaults to A(0,0)
  if (isBlocked) A = bA->getMatrix(0, 0);

  Magnitude tol    = pL.get<Magnitude>("inverse: drop tolerance");
  RCP<Matrix> Ainv = Teuchos::null;

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
    auto maxRowEntries            = sparsityPattern->getLocalMaxNumRowEntries();
    // too many nonzeros per row will cause us to exceed L1 memory space size
    // lets increase threshold to reduce nnzs per row
    while (maxRowEntries > 200) {
      if (tol == 0.0)
        tol = .005;
      else
        tol *= 1.5;
      if (IsPrint(Statistics1)) {
        GetOStream(Statistics1) << "MueLu_InverseApproximationFactory: Too many NNZs per row (" << maxRowEntries << "). Reducing NNZs per row by increasing drop threshold to " << tol << std::endl;
      }
      sparsityPattern = Utilities::GetThresholdedGraph(A, tol);
      maxRowEntries   = sparsityPattern->getLocalMaxNumRowEntries();
    }
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
  } else if (method == "factoredsparseapproxinverse") {
    RCP<CrsGraph> sparsityPattern = Utilities::GetThresholdedLowerTriangularGraph(A, tol);
    auto maxRowEntries            = sparsityPattern->getLocalMaxNumRowEntries();
    // too many nonzeros per row will cause us to exceed L1 memory space size
    // lets increase threshold to reduce nnzs per row
    while (maxRowEntries > 200) {
      if (tol == 0.0)
        tol = .005;
      else
        tol *= 1.5;
      if (IsPrint(Statistics1)) {
        GetOStream(Statistics1) << "MueLu_InverseApproximationFactory: Too many NNZs per row (" << maxRowEntries << "). Reducing NNZs per row by increasing drop threshold to " << tol << std::endl;
      }
      sparsityPattern = Utilities::GetThresholdedLowerTriangularGraph(A, tol);
      maxRowEntries   = sparsityPattern->getLocalMaxNumRowEntries();
    }

    if (IsPrint(Statistics1)) {
      sparsityPattern->computeGlobalConstants();
      GetOStream(Statistics1) << "NNZ Graph(A): " << A->getCrsGraph()->getGlobalNumEntries() << " , NNZ Tresholded Graph(triLower(A)): " << sparsityPattern->getGlobalNumEntries() << std::endl;
    }
    RCP<Matrix> pLinvFactor = GetFactoredSparseInverse(A, sparsityPattern);
    RCP<Matrix> LinvFactor  = Utilities::GetThresholdedMatrix(pLinvFactor, tol, fixing);
    // To create the inverse from the inverse factor, we need to multiply Linv' * LinvFactor. Of course, we could
    // save a fair amount of storage by delaying this to when we actually need it as Linv' * LinvFactor has
    // many more nonzeros than just Linv. One other thing, I'm explicitly using the Transpose computation as opposed to
    // setting one of the booleans to true in the MatrixMatrix:Multiply. There are some comments about true/false combinations
    // that don't work, so I decided not to push my luck here.

    RCP<Matrix> LinvTrans = Utilities::Transpose(*LinvFactor, true);
    Ainv                  = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*LinvTrans, false, *LinvFactor, false, GetOStream(Statistics2), true, true, std::string("Ainv"));

    if (IsPrint(Statistics1)) {
      rcp_const_cast<CrsGraph>(Ainv->getCrsGraph())->computeGlobalConstants();
      GetOStream(Statistics1) << "NNZ Linv: " << LinvFactor->getGlobalNumEntries() << ", NNZ Tresholded Linv (parameter: " << tol << "): " << pLinvFactor->getGlobalNumEntries() << std::endl;
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
    // create a unique version of the column indicies that has the correct length (as opposed
    // to column_indices).  Can we instead resize column_indices with MemoryUnmanaged?
    // This is so that we can use binary search later on sorted list
    shared_lo_vector uniqueColIndicies(thread.team_scratch(scratchLevel), numUniqeColEntries);
    for (local_ordinal_type m = 0; m < numUniqeColEntries; ++m) {
      uniqueColIndicies(m) = column_indices(m);
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

        // do binary search to find column in uniqueColIndices
        auto it              = std::lower_bound(Kokkos::Experimental::begin(uniqueColIndicies), Kokkos::Experimental::end(uniqueColIndicies), j);
        ptrdiff_t newIndex   = it - Kokkos::Experimental::begin(uniqueColIndicies);
        localA(newIndex, ii) = v;
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

  // Transpose needed to match published paper algorithms as the inverse is not symmetric
  // However, non-transposed version seems to work better in row-oriented MinvA algorithm
  // RCP<Matrix> actualSpai = Utilities::Transpose(*Ainv, true); // , label, Tparams);

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

  // Transpose needed to match published paper algorithms as the inverse is not symmetric
  // However, non-transposed version seems to work better in row-oriented MinvA algorithm
  // RCP<Matrix> actualSpai = Utilities::Transpose(*Ainv, true); // , label, Tparams);

  return Ainv;
}

#endif

template <class local_matrix_type, typename global_ordinal_type>
class LocalFSAIFunctor {
 private:
  using scalar_type        = typename local_matrix_type::value_type;
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using execution_space    = typename local_matrix_type::execution_space;
  using impl_scalar_type   = typename KokkosKernels::ArithTraits<scalar_type>::val_type;
  using impl_ATS           = KokkosKernels::ArithTraits<impl_scalar_type>;
  using device_type        = typename local_matrix_type::device_type;
  using local_map_type     = Tpetra::Details::LocalMap<local_ordinal_type, global_ordinal_type, device_type>;

 public:
  using shared_matrix    = Kokkos::View<impl_scalar_type**, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_vector    = Kokkos::View<impl_scalar_type*, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_lo_vector = Kokkos::View<local_ordinal_type*, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;

 private:
  const local_matrix_type lclA;
  local_matrix_type lclAinv;
  local_map_type lclARowMap;
  local_map_type lclAColMap;
  local_map_type lclAinvRowMap;
  local_map_type lclAinvColMap;
  const int scratchLevel;

 public:
  LocalFSAIFunctor(const local_matrix_type& lclA_, local_matrix_type& lclAinv_, const local_map_type& lclARowMap_, const local_map_type& lclAColMap_,
                   const local_map_type& lclAinvRowMap_, const local_map_type& lclAinvColMap_, int scratchLevel_)
    : lclA(lclA_)
    , lclAinv(lclAinv_)
    , lclARowMap(lclARowMap_)
    , lclAColMap(lclAColMap_)
    , lclAinvRowMap(lclAinvRowMap_)
    , lclAinvColMap(lclAinvColMap_)
    , scratchLevel(scratchLevel_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread) const {
    auto rlid    = thread.league_rank();
    auto rowAinv = lclAinv.row(rlid);

    // matlab version of algorithm
    //  n = size(A,1);    nzs = nnz(S);
    //  newrows = zeros(nzs,1); newcols = zeros(nzs,1); newvals= zeros(nzs,1);
    //  count = 0;
    //  for i=1:n,
    //    [~,subCols,~] = find(S(i,:));
    //    diagLocation = find(subCols == i);
    //    submat = A(subCols,subCols);
    //    subn = size(submat,1);
    //    if diagLocation ~= subn, fprintf('pattern not lower triangular?\n'); keyboard; end;
    //    identCol = zeros(subn,1); identCol(diagLocation) = 1;
    //    AinvFactorRow = submat\identCol;
    //    normalizedFactor = AinvFactorRow/sqrt(AinvFactorRow(diagLocation));
    //    newrows(count+1:count+subn) = i;
    //    newcols(count+1:count+subn) = subCols;
    //    newvals(count+1:count+subn) = normalizedFactor;
    //    count = count + subn;
    // end;
    // Lfactor  = sparse(newrows,newcols,newvals,n,n);
    //

    auto A_rowGid = lclARowMap.getGlobalElement(rlid);

    auto numRowEntries                    = rowAinv.length;
    local_ordinal_type diagOffset         = -1;
    scalar_type diagValue                 = 0.0;
    local_ordinal_type A_lclRowIndForDiag = -1;

    // Loop over entries in row rlid of Ainv and collect all of A's column indices.
    shared_lo_vector column_indices(thread.team_scratch(scratchLevel), numRowEntries);
    local_ordinal_type numColEntries = rowAinv.length;
    for (local_ordinal_type ii = 0; ii < rowAinv.length; ++ii) {
      auto i                         = rowAinv.colidx(ii);
      auto Ainv_colGid               = lclAinvColMap.getGlobalElement(i);  // for debugging
      local_ordinal_type A_lclRowInd = lclARowMap.getLocalElement(Ainv_colGid);
      local_ordinal_type A_lclColInd = lclAColMap.getLocalElement(Ainv_colGid);
      column_indices(ii)             = A_lclColInd;
      if (A_rowGid == Ainv_colGid) {
        A_lclRowIndForDiag = A_lclRowInd;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(A_lclRowIndForDiag == -1, Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::GetSparseInverse: no diagonal entry found in A.");

    Kokkos::Experimental::sort_team(thread, column_indices);  // in order to apply binary search later
    for (int kkk = 0; kkk < numColEntries; kkk++) {
      if (column_indices(kkk) == A_lclRowIndForDiag) diagOffset = kkk;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(diagOffset == -1, Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::GetSparseInverse: no diagonal entry found in A.");

    // Extract local part of A into a dense view.
    shared_matrix localA(thread.team_scratch(scratchLevel), numRowEntries, rowAinv.length);
    KokkosBlas::SerialSet::invoke(impl_ATS::zero(), localA);

    // Now fill localA.
    for (local_ordinal_type ii = 0; ii < rowAinv.length; ++ii) {
      auto i                         = rowAinv.colidx(ii);
      auto Ainv_colGid               = lclAinvColMap.getGlobalElement(i);  // for debugging
      local_ordinal_type A_lclRowInd = lclARowMap.getLocalElement(Ainv_colGid);
      TEUCHOS_TEST_FOR_EXCEPTION(A_lclRowInd == -1, Exceptions::RuntimeError, "MueLu::InverseApproximationFactory: Column global ID in Ainv not found in A rowmap\n");
      auto rowA = lclA.rowConst(A_lclRowInd);

      for (local_ordinal_type jj = 0; jj < rowA.length; ++jj) {
        auto j = rowA.colidx(jj);
        // do binary search to find column in column_indices, but first check that it is
        // in lower triangular portion of matrix (because this might be faster?)
        auto A_colGid = lclAColMap.getGlobalElement(j);
        if (A_colGid <= A_rowGid) {
          auto it = std::lower_bound(Kokkos::Experimental::begin(column_indices), Kokkos::Experimental::end(column_indices), j);
          if (it != Kokkos::Experimental::end(column_indices) && *it == j) {
            ptrdiff_t newIndex   = it - Kokkos::Experimental::begin(column_indices);
            auto v               = rowA.value(jj);
            localA(newIndex, ii) = v;
          }
        }
      }
    }
    shared_matrix ek(thread.team_scratch(scratchLevel), numRowEntries, 1);
    // set to zero, set diagonal entry to one
    for (local_ordinal_type i = 0; i < numRowEntries; ++i) {
      ek(i, 0) = (i == diagOffset) ? impl_ATS::one() : impl_ATS::zero();
    }

    // QR solve
    shared_vector tau(thread.team_scratch(scratchLevel), rowAinv.length);
    shared_vector work(thread.team_scratch(scratchLevel), numRowEntries);
    // factorize localA = Q*R in-place
#define QRway
#ifdef QRway
    KokkosBatched::SerialQR<KokkosBatched::Algo::QR::Unblocked>::invoke(localA, tau, work);
#else
    KokkosBatched::SerialCholesky<KokkosBatched::Uplo::Lower, KokkosBatched::Algo::Cholesky::Unblocked>::invoke(localA);                                                                                   // use Cholesky
#endif
    // ek := Q^T ek
#ifdef QRway
    KokkosBatched::SerialApplyQ<KokkosBatched::Side::Left, KokkosBatched::Trans::Transpose, KokkosBatched::Algo::ApplyQ::Unblocked>::invoke(localA, tau, ek, work);
    // ek[:rowLength] := R^{-1} ek[:rowLength]
    auto sub_A = Kokkos::subview(localA, Kokkos::make_pair(0, rowAinv.length), Kokkos::ALL());
#else
    auto sub_A = Kokkos::subview(localA, Kokkos::make_pair(0, rowAinv.length), Kokkos::make_pair(0, rowAinv.length));                                                                                      // use Cholesky
#endif
    auto sub_ek = Kokkos::subview(ek, Kokkos::make_pair(0, rowAinv.length), 0);
#ifdef QRway
    KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(impl_ATS::one(), sub_A, sub_ek);
#else
    KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(impl_ATS::one(), sub_A, sub_ek);  // use Cholesky
    KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Unblocked>::invoke(impl_ATS::one(), sub_A, sub_ek);    // use Cholesky
#endif

    // Set entries of Ainv.

    diagValue = sub_ek(diagOffset);
    TEUCHOS_TEST_FOR_EXCEPTION(impl_ATS::real(diagValue) <= 0.0, Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::GetSparseInverse: non positive diagonal entry.");
    auto scale_factor = impl_ATS::one() / impl_ATS::sqrt(diagValue);
    for (local_ordinal_type i = 0; i < rowAinv.length; ++i) {
//      auto thevalue = sub_ek(i) * scale_factor;
      typename KokkosKernels::ArithTraits<decltype(diagValue)>::val_type thevalue = sub_ek(i) * scale_factor;

      if (thevalue == impl_ATS::zero()) thevalue = impl_ATS::eps();
      rowAinv.value(i) = thevalue;
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFactoredSparseInverse(const RCP<Matrix>& Aorg, const RCP<const CrsGraph>& sparsityPattern) const {
  using execution_space = typename Node::execution_space;

  // construct the inverse matrix with factor the given sparsity pattern
  RCP<Matrix> Ainv = MatrixFactory::Build(sparsityPattern);
  Ainv->resumeFill();

  // gather missing rows from other procs to generate an overlapping map
  RCP<Import> rowImport = ImportFactory::Build(sparsityPattern->getRowMap(), sparsityPattern->getColMap());
  RCP<Matrix> A         = MatrixFactory::Build(Aorg, *rowImport);

  auto maxRowEntriesAinv = Ainv->getLocalMaxNumRowEntries();
  {
    auto lclA          = A->getLocalMatrixDevice();
    auto lclAinv       = Ainv->getLocalMatrixDevice();
    auto lclARowmap    = A->getRowMap()->getLocalMap();
    auto lclAColmap    = A->getColMap()->getLocalMap();
    auto lclAinvRowmap = Ainv->getRowMap()->getLocalMap();
    auto lclAinvColmap = Ainv->getColMap()->getLocalMap();
    auto lclAorgRowmap = Aorg->getRowMap()->getLocalMap();

    Kokkos::TeamPolicy<execution_space> policy(lclAinv.numRows(), 1);

    using fsai_functor_type = LocalFSAIFunctor<decltype(lclAinv), GlobalOrdinal>;
    using shared_matrix     = typename fsai_functor_type::shared_matrix;
    using shared_vector     = typename fsai_functor_type::shared_vector;
    using shared_lo_vector  = typename fsai_functor_type::shared_lo_vector;

    int size = shared_matrix::shmem_size(maxRowEntriesAinv, maxRowEntriesAinv) + shared_matrix::shmem_size(maxRowEntriesAinv, 1) + shared_vector::shmem_size(3 * maxRowEntriesAinv) + shared_vector::shmem_size(maxRowEntriesAinv) + shared_lo_vector::shmem_size(maxRowEntriesAinv);

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

    LocalFSAIFunctor fsaiFunctor(lclA, lclAinv, lclARowmap, lclAColmap, lclAinvRowmap, lclAinvColmap, scratchLevel);

    Kokkos::parallel_for("MueLu::InverseFactory::LocalSpai", policy, fsaiFunctor);
  }

  Ainv->fillComplete();

  return Ainv;
}

}  // namespace MueLu

#endif /* MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_ */
