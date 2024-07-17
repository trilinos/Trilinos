// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_
#define MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_

#include "MueLu_RepartitionBlockDiagonalFactory_decl.hpp"

#include <Teuchos_Utils.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");
}  // DeclareInput()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockCrs;

  RCP<Matrix> originalA = Get<RCP<Matrix> >(currentLevel, "A");
  RCP<BlockCrs> A       = Teuchos::rcp_dynamic_cast<BlockCrs>(originalA);

  //    RCP<BlockCrs> A = Get< RCP<BlockCrs> >(currentLevel, "A");
  TEUCHOS_TEST_FOR_EXCEPTION(A == Teuchos::null, Exceptions::BadCast, "MueLu::RepartitionBlockDiagonalFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

  // Build the block diagonal
  RCP<BlockCrs> DiagonalMatrix = Teuchos::rcp(new BlockCrs(A->getBlockedRangeMap(), A->getBlockedDomainMap(), 0));
  for (size_t i = 0; i < A->Rows(); i++)
    DiagonalMatrix->setMatrix(i, i, A->getMatrix(i, i));

  Set(currentLevel, "A", Teuchos::rcp_dynamic_cast<Matrix>(DiagonalMatrix));

}  // Build()

}  // end namespace MueLu

#endif /* MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_ */
