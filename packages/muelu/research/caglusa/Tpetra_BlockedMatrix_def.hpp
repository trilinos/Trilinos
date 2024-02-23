#ifndef TPETRA_BLOCKEDMATRIX_DEF_HPP
#define TPETRA_BLOCKEDMATRIX_DEF_HPP

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMatrix(const Teuchos::RCP<matrix_type>& pointA,
                  const Teuchos::RCP<matrix_type>& blockA,
                  const Teuchos::RCP<blocked_map_type>& blockMap,
                  const Teuchos::RCP<blocked_map_type>& ghosted_blockMap)
  : pointA_(pointA)
  , blockA_(blockA)
  , blockMap_(blockMap)
  , ghosted_blockMap_(ghosted_blockMap) {
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRangeMap()));
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRowMap()));
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockMap_->blockMap_));

  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRangeMap()));
  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRowMap()));
  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*blockMap_->pointMap_));
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  pointA_->apply(X, Y, mode, alpha, beta);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    localApply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
               Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
               Teuchos::ETransp mode,
               Scalar alpha,
               Scalar beta) const {
  pointA_->localApply(X, Y, mode, alpha, beta);
}
}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMATRIX_DEF_HPP
