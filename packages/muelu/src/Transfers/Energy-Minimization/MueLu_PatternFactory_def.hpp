// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PATTERNFACTORY_DEF_HPP
#define MUELU_PATTERNFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_PatternFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
//#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("emin: pattern order");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory for the matrix");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory for the matrix providing nonzero graph");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(coarseLevel, "P");

  const ParameterList& pL = GetParameterList();
  if (pL.get<int>("emin: pattern order") > 0)
    Input(fineLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Ppattern", coarseLevel);

  RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P");

  const ParameterList& pL = GetParameterList();
  int k                   = pL.get<int>("emin: pattern order");

  if (k > 0) {
    RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> AP;

    bool doFillComplete  = true;
    bool optimizeStorage = true;

    for (int i = 0; i < k; i++) {
      AP = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *P, false, GetOStream(Statistics2), doFillComplete, optimizeStorage);
      P.swap(AP);
    }
  }

  Set(coarseLevel, "Ppattern", P->getCrsGraph());
}

}  // namespace MueLu

#endif  // MUELU_PATTERNFACTORY_DEF_HPP
