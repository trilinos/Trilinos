#ifndef MUELU_FLATOPERATOR_DEF_HPP
#define MUELU_FLATOPERATOR_DEF_HPP

#include "MueLu_FlatOperator_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FlatOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FlatOperator(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mat,
                 const Teuchos::RCP<MueLu::Constraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>> constraint)
  : mat_(mat)
  , constraint_(constraint) {
  CheckMaps();
  AllocateTemporaryMatrix();

  auto pattern = constraint_->GetPattern();
  map_         = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(pattern->getRowMap()->lib(),
                                                                              pattern->getGlobalNumEntries(),
                                                                              pattern->getLocalNumEntries(),
                                                                              pattern->getRowMap()->getIndexBase(),
                                                                              pattern->getComm());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FlatOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
          Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  AllocateTemporaryMatrix();

  TEUCHOS_ASSERT(mode == Teuchos::NO_TRANS);
  TEUCHOS_ASSERT(alpha == Teuchos::ScalarTraits<Scalar>::one());
  TEUCHOS_ASSERT(beta == Teuchos::ScalarTraits<Scalar>::zero());

  {
    auto lclMat = tempMat_->getLocalMatrixDevice();
    auto lclVec = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    TEUCHOS_ASSERT(lclMat.values.extent(0) == lclVec.extent(0));
    Kokkos::deep_copy(lclMat.values, Kokkos::subview(lclVec, Kokkos::ALL(), 0));
  }

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AP;
  AP = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*mat_, false, *tempMat_, false, AP, GetOStream(Runtime0), true, true);
  constraint_->AssignMatrixEntriesToVector(*AP, Y);
}
}  // namespace MueLu
#endif
