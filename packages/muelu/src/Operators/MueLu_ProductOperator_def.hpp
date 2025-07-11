#ifndef MUELU_PRODUCTOPERATOR_DEF_HPP
#define MUELU_PRODUCTOPERATOR_DEF_HPP

#include "MueLu_ProductOperator.hpp"
namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ProductOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &X,
          Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  AllocateTemporaryMultiVectors(X.getNumVectors());

  const auto one  = Teuchos::ScalarTraits<Scalar>::one();
  const auto zero = Teuchos::ScalarTraits<Scalar>::zero();

#define toggle(opmode) ((opmode) == Teuchos::NO_TRANS) ? Teuchos::TRANS : Teuchos::NO_TRANS

  if (ops_.size() >= 2) {
    if (mode == Teuchos::NO_TRANS) {
      ops_[ops_.size() - 1]->apply(X, *tempVecs_[ops_.size() - 2], modes_[ops_.size() - 1], one, zero);
      for (size_t i = ops_.size() - 2; i > 0; --i) {
        ops_[i]->apply(*tempVecs_[i], *tempVecs_[i - 1], modes_[i], one, zero);
      }
      ops_[0]->apply(*tempVecs_[0], Y, modes_[0], alpha, beta);
    } else {
      ops_[0]->apply(X, *tempVecs_[0], toggle(modes_[0]), one, zero);
      for (size_t i = 1; i < ops_.size() - 2; ++i) {
        ops_[i]->apply(*tempVecs_[i], *tempVecs_[i - 1], toggle(modes_[i]), one, zero);
      }
      ops_[ops_.size() - 1]->apply(*tempVecs_[ops_.size() - 2], Y, toggle(modes_[ops_.size() - 1]), alpha, beta);
    }
  } else {
    if (mode == Teuchos::NO_TRANS)
      ops_[0]->apply(X, Y, modes_[0], alpha, beta);
    else
      ops_[0]->apply(X, Y, toggle(modes_[0]), alpha, beta);
  }

#undef toggle
}

}  // namespace MueLu

#endif
