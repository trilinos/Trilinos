// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP
#define MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_NullspacePresmoothFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory for A");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory for the nonsmoothed nullspace");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "Nullspace");

  if (currentLevel.GetLevelID() == 0)
    Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Nullspace presmooth factory", currentLevel);

  RCP<MultiVector> newB;
  if (currentLevel.GetLevelID() == 0) {
    RCP<Matrix> A      = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<MultiVector> B = Get<RCP<MultiVector> >(currentLevel, "Nullspace");
    newB               = MultiVectorFactory::Build(B->getMap(), B->getNumVectors());

    Teuchos::ArrayRCP<SC> D = Utilities::GetMatrixDiagonal_arcp(*A);

    SC damping = 4. / 3;
    damping /= Utilities::PowerMethod(*A, true, (LO)10, 1e-4);

    A->apply(*B, *newB, Teuchos::NO_TRANS);

    size_t numVec  = newB->getNumVectors();
    LO numElements = newB->getLocalLength();
    for (size_t j = 0; j < numVec; j++) {
      Teuchos::ArrayRCP<const SC> Bj = B->getData(j);
      Teuchos::ArrayRCP<SC> newBj    = newB->getDataNonConst(j);

      for (LO i = 0; i < numElements; i++)
        newBj[i] = Bj[i] - damping * newBj[i] / D[i];
    }
  } else {
    newB = Get<RCP<MultiVector> >(currentLevel, "Nullspace");
  }

  // provide "Nullspace" variable on current level
  Set(currentLevel, "Nullspace", newB);

}  // Build

}  // namespace MueLu

#endif  // MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP
