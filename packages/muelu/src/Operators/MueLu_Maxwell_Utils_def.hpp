// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MAXWELL_UTILS_DEF_HPP
#define MUELU_MAXWELL_UTILS_DEF_HPP

#include <sstream>

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_Maxwell_Utils_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_ThresholdAFilterFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::detectBoundaryConditionsSM(RCP<Matrix>& SM_Matrix_,
                                                                                          RCP<Matrix>& D0_Matrix_,
                                                                                          magnitudeType rowSumTol,
                                                                                          bool useKokkos_,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCrowsKokkos_,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCcolsKokkos_,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCdomainKokkos_,
                                                                                          int& BCedges_,
                                                                                          int& BCnodes_,
                                                                                          Teuchos::ArrayRCP<bool>& BCrows_,
                                                                                          Teuchos::ArrayRCP<bool>& BCcols_,
                                                                                          Teuchos::ArrayRCP<bool>& BCdomain_,
                                                                                          bool& allEdgesBoundary_,
                                                                                          bool& allNodesBoundary_) {
  // clean rows associated with boundary conditions
  // Find rows with only 1 or 2 nonzero entries, record them in BCrows_.
  // BCrows_[i] is true, iff i is a boundary row
  // BCcols_[i] is true, iff i is a boundary column
  int BCedgesLocal = 0;
  int BCnodesLocal = 0;
  if (useKokkos_) {
    BCrowsKokkos_ = Utilities::DetectDirichletRows_kokkos(*SM_Matrix_, Teuchos::ScalarTraits<magnitudeType>::eps(), /*count_twos_as_dirichlet=*/true);

    if (rowSumTol > 0.)
      Utilities::ApplyRowSumCriterion(*SM_Matrix_, rowSumTol, BCrowsKokkos_);

    BCcolsKokkos_   = Kokkos::View<bool*, typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), D0_Matrix_->getColMap()->getLocalNumElements());
    BCdomainKokkos_ = Kokkos::View<bool*, typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletDomains"), D0_Matrix_->getDomainMap()->getLocalNumElements());
    Utilities::DetectDirichletColsAndDomains(*D0_Matrix_, BCrowsKokkos_, BCcolsKokkos_, BCdomainKokkos_);

    auto BCrowsKokkos = BCrowsKokkos_;
    Kokkos::parallel_reduce(
        BCrowsKokkos_.size(), KOKKOS_LAMBDA(int i, int& sum) {
          if (BCrowsKokkos(i))
            ++sum;
        },
        BCedgesLocal);

    auto BCdomainKokkos = BCdomainKokkos_;
    Kokkos::parallel_reduce(
        BCdomainKokkos_.size(), KOKKOS_LAMBDA(int i, int& sum) {
          if (BCdomainKokkos(i))
            ++sum;
        },
        BCnodesLocal);
  } else {
    BCrows_ = Teuchos::arcp_const_cast<bool>(Utilities::DetectDirichletRows(*SM_Matrix_, Teuchos::ScalarTraits<magnitudeType>::eps(), /*count_twos_as_dirichlet=*/true));

    if (rowSumTol > 0.)
      Utilities::ApplyRowSumCriterion(*SM_Matrix_, rowSumTol, BCrows_);

    BCcols_.resize(D0_Matrix_->getColMap()->getLocalNumElements());
    BCdomain_.resize(D0_Matrix_->getDomainMap()->getLocalNumElements());
    Utilities::DetectDirichletColsAndDomains(*D0_Matrix_, BCrows_, BCcols_, BCdomain_);

    for (auto it = BCrows_.begin(); it != BCrows_.end(); ++it)
      if (*it)
        BCedgesLocal += 1;
    for (auto it = BCdomain_.begin(); it != BCdomain_.end(); ++it)
      if (*it)
        BCnodesLocal += 1;
  }

  MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCedgesLocal, BCedges_);
  MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCnodesLocal, BCnodes_);

  allEdgesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCedges_) >= D0_Matrix_->getRangeMap()->getGlobalNumElements();
  allNodesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCnodes_) >= D0_Matrix_->getDomainMap()->getGlobalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::detectBoundaryConditionsSM(RCP<Matrix>& SM_Matrix_,
                                                                                          RCP<Matrix>& D0_Matrix_,
                                                                                          magnitudeType rowSumTol,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCrowsKokkos_,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCcolsKokkos_,
                                                                                          Kokkos::View<bool*, typename Node::device_type::memory_space>& BCdomainKokkos_,
                                                                                          int& BCedges_,
                                                                                          int& BCnodes_,
                                                                                          bool& allEdgesBoundary_,
                                                                                          bool& allNodesBoundary_) {
  // clean rows associated with boundary conditions
  // Find rows with only 1 or 2 nonzero entries, record them in BCrows_.
  // BCrows_[i] is true, iff i is a boundary row
  // BCcols_[i] is true, iff i is a boundary column
  int BCedgesLocal = 0;
  int BCnodesLocal = 0;
  {
    BCrowsKokkos_ = Utilities::DetectDirichletRows_kokkos(*SM_Matrix_, Teuchos::ScalarTraits<magnitudeType>::eps(), /*count_twos_as_dirichlet=*/true);

    if (rowSumTol > 0.)
      Utilities::ApplyRowSumCriterion(*SM_Matrix_, rowSumTol, BCrowsKokkos_);

    BCcolsKokkos_   = Kokkos::View<bool*, typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), D0_Matrix_->getColMap()->getLocalNumElements());
    BCdomainKokkos_ = Kokkos::View<bool*, typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletDomains"), D0_Matrix_->getDomainMap()->getLocalNumElements());
    Utilities::DetectDirichletColsAndDomains(*D0_Matrix_, BCrowsKokkos_, BCcolsKokkos_, BCdomainKokkos_);

    auto BCrowsKokkos = BCrowsKokkos_;
    Kokkos::parallel_reduce(
        BCrowsKokkos_.size(), KOKKOS_LAMBDA(int i, int& sum) {
          if (BCrowsKokkos(i))
            ++sum;
        },
        BCedgesLocal);

    auto BCdomainKokkos = BCdomainKokkos_;
    Kokkos::parallel_reduce(
        BCdomainKokkos_.size(), KOKKOS_LAMBDA(int i, int& sum) {
          if (BCdomainKokkos(i))
            ++sum;
        },
        BCnodesLocal);
  }

  MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCedgesLocal, BCedges_);
  MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCnodesLocal, BCnodes_);

  allEdgesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCedges_) >= D0_Matrix_->getRangeMap()->getGlobalNumElements();
  allNodesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCnodes_) >= D0_Matrix_->getDomainMap()->getGlobalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    removeExplicitZeros(const RCP<Matrix>& A,
                        const magnitudeType tolerance,
                        const bool keepDiagonal,
                        const size_t expectedNNZperRow) {
  Level fineLevel;
  fineLevel.SetFactoryManager(null);
  fineLevel.SetLevelID(0);
  fineLevel.Set("A", A);
  fineLevel.setlib(A->getDomainMap()->lib());
  RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A", tolerance, keepDiagonal, expectedNNZperRow));
  fineLevel.Request("A", ThreshFact.get());
  ThreshFact->Build(fineLevel);
  return fineLevel.Get<RCP<Matrix> >("A", ThreshFact.get());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::removeExplicitZeros(Teuchos::ParameterList& parameterList_,
                                                                                   RCP<Matrix>& D0_Matrix_,
                                                                                   RCP<Matrix>& SM_Matrix_,
                                                                                   RCP<Matrix>& M1_Matrix_,
                                                                                   RCP<Matrix>& Ms_Matrix_) {
  bool defaultFilter = false;

  // Remove zero entries from D0 if necessary.
  // In the construction of the prolongator we use the graph of the
  // matrix, so zero entries mess it up.
  if (parameterList_.get<bool>("refmaxwell: filter D0", true) && D0_Matrix_->getLocalMaxNumRowEntries() > 2) {
    D0_Matrix_ = removeExplicitZeros(D0_Matrix_, 1e-8, false, 2);

    // If D0 has too many zeros, maybe SM and M1 do as well.
    defaultFilter = true;
  }

  if (!M1_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter M1", defaultFilter)) {
    RCP<Vector> diag = VectorFactory::Build(M1_Matrix_->getRowMap());
    // find a reasonable absolute value threshold
    M1_Matrix_->getLocalDiagCopy(*diag);
    magnitudeType threshold = 1.0e-8 * diag->normInf();

    M1_Matrix_ = removeExplicitZeros(M1_Matrix_, threshold, true);
  }

  if (!Ms_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter Ms", defaultFilter)) {
    RCP<Vector> diag = VectorFactory::Build(Ms_Matrix_->getRowMap());
    // find a reasonable absolute value threshold
    Ms_Matrix_->getLocalDiagCopy(*diag);
    magnitudeType threshold = 1.0e-8 * diag->normInf();

    Ms_Matrix_ = removeExplicitZeros(Ms_Matrix_, threshold, true);
  }

  if (!SM_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter SM", defaultFilter)) {
    RCP<Vector> diag = VectorFactory::Build(SM_Matrix_->getRowMap());
    // find a reasonable absolute value threshold
    SM_Matrix_->getLocalDiagCopy(*diag);
    magnitudeType threshold = 1.0e-8 * diag->normInf();

    SM_Matrix_ = removeExplicitZeros(SM_Matrix_, threshold, true);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    thresholdedAbs(const RCP<Matrix>& A,
                   const magnitudeType threshold) {
  using ATS         = Kokkos::ArithTraits<Scalar>;
  using impl_Scalar = typename ATS::val_type;
  using impl_ATS    = Kokkos::ArithTraits<impl_Scalar>;
  using range_type  = Kokkos::RangePolicy<LO, typename NO::execution_space>;

  const impl_Scalar impl_SC_ONE  = impl_ATS::one();
  const impl_Scalar impl_SC_ZERO = impl_ATS::zero();

  {
    auto lclMat = A->getLocalMatrixDevice();
    Kokkos::parallel_for(
        "thresholdedAbs",
        range_type(0, lclMat.nnz()),
        KOKKOS_LAMBDA(const size_t i) {
          if (impl_ATS::magnitude(lclMat.values(i)) > threshold)
            lclMat.values(i) = impl_SC_ONE;
          else
            lclMat.values(i) = impl_SC_ZERO;
        });
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMatvecParams(Matrix& A, RCP<ParameterList> matvecParams) {
  RCP<const Import> xpImporter = A.getCrsGraph()->getImporter();
  if (!xpImporter.is_null())
    xpImporter->setDistributorParameters(matvecParams);
  RCP<const Export> xpExporter = A.getCrsGraph()->getExporter();
  if (!xpExporter.is_null())
    xpExporter->setDistributorParameters(matvecParams);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
Maxwell_Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    PtAPWrapper(const RCP<Matrix>& A, const RCP<Matrix>& P, ParameterList& params, const std::string& label) {
  Level fineLevel, coarseLevel;
  fineLevel.SetFactoryManager(null);
  coarseLevel.SetFactoryManager(null);
  coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
  fineLevel.SetLevelID(0);
  coarseLevel.SetLevelID(1);
  fineLevel.Set("A", A);
  coarseLevel.Set("P", P);
  coarseLevel.setlib(A->getDomainMap()->lib());
  fineLevel.setlib(A->getDomainMap()->lib());
  coarseLevel.setObjectLabel(label);
  fineLevel.setObjectLabel(label);

  RCP<RAPFactory> rapFact = rcp(new RAPFactory());
  ParameterList rapList   = *(rapFact->GetValidParameterList());
  rapList.set("transpose: use implicit", true);
  rapList.set("rap: fix zero diagonals", params.get<bool>("rap: fix zero diagonals", true));
  rapList.set("rap: fix zero diagonals threshold", params.get<double>("rap: fix zero diagonals threshold", Teuchos::ScalarTraits<double>::eps()));
  rapList.set("rap: triple product", params.get<bool>("rap: triple product", false));
  rapFact->SetParameterList(rapList);

  coarseLevel.Request("A", rapFact.get());

  return coarseLevel.Get<RCP<Matrix> >("A", rapFact.get());
}

}  // namespace MueLu

#define MUELU_MAXWELL_UTILS_SHORT
#endif  // ifdef MUELU_MAXWELL_UTILS_DEF_HPP
