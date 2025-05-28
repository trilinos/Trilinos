// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
#define MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_EdgeProlongatorPatternFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
//#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("FineD0", Teuchos::null, "Generating factory for the fine discrete gradient");
  validParamList->set<RCP<const FactoryBase> >("CoarseD0", Teuchos::null, "Generating factory for the coarse discrete gradient");
  validParamList->set<RCP<const FactoryBase> >("Pnodal", Teuchos::null, "Generating factory for the nodal prolongator");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "D0", "FineD0");
  Input(coarseLevel, "D0", "CoarseD0");
  Input(coarseLevel, "Pnodal");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "EdgeProlongatorPattern", coarseLevel);

  RCP<Matrix> D  = Get<RCP<Matrix> >(fineLevel, "D0", "FineD0");
  RCP<Matrix> Dc = Get<RCP<Matrix> >(coarseLevel, "D0", "CoarseD0");
  RCP<Matrix> Pn = Get<RCP<Matrix> >(coarseLevel, "Pnodal");

  const auto one = Teuchos::ScalarTraits<Scalar>::one();

  auto absD = MatrixFactory::BuildCopy(D);
  absD->setAllToScalar(one);

  auto absPn = MatrixFactory::BuildCopy(Pn);
  absPn->setAllToScalar(one);

  auto absDc = MatrixFactory::BuildCopy(Dc);
  absDc->setAllToScalar(one);

  RCP<Matrix> temp1 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absD, false, *absPn, false, GetOStream(Statistics2), true, true);
  RCP<Matrix> temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*temp1, false, *absDc, true, GetOStream(Statistics2), true, true);

  Set(coarseLevel, "Ppattern", temp2->getCrsGraph());
}

}  // namespace MueLu

#endif  // MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
