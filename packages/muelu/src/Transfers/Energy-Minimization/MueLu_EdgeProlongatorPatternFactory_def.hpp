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
//#include <Xpetra_IO.hpp>

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
  validParamList->set<RCP<const FactoryBase> >("PnodalEmin", Teuchos::null, "Generating factory for the nodal prolongator");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "D0", "FineD0");
  Input(coarseLevel, "D0", "CoarseD0");
  Input(coarseLevel, "PnodalEmin");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "EdgeProlongatorPattern", coarseLevel);

  RCP<Matrix> D  = Get<RCP<Matrix> >(fineLevel, "D0", "FineD0");
  RCP<Matrix> Dc = Get<RCP<Matrix> >(coarseLevel, "D0", "CoarseD0");
  RCP<Matrix> Pn = Get<RCP<Matrix> >(coarseLevel, "PnodalEmin");

  const auto one = Teuchos::ScalarTraits<Scalar>::one();

  // |FineD| * |Pnodal| * |CoarseD^T|

  auto absD = MatrixFactory::BuildCopy(D);
  absD->setAllToScalar(one);

  auto absPn = MatrixFactory::BuildCopy(Pn);
  absPn->setAllToScalar(one);

  auto absDc = MatrixFactory::BuildCopy(Dc);
  absDc->setAllToScalar(one);

  RCP<Matrix> temp1 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absD, false, *absPn, false, GetOStream(Statistics2), true, true);
  temp1->setAllToScalar(one);
  RCP<Matrix> temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*temp1, false, *absDc, true, GetOStream(Statistics2), true, true);

  RCP<Matrix> filtered;
  filtered = MatrixFactory::Build(temp2->getRowMap(), temp2->getColMap(), temp2->getLocalMaxNumRowEntries());
  ArrayView<const LO> inds2;
  ArrayView<const SC> vals2;

  RCP<MultiVector> oneVec = MultiVectorFactory::Build(absDc->getDomainMap(),1);
  oneVec->putScalar(one );
  RCP<MultiVector> singleParent = MultiVectorFactory::Build(absDc->getRowMap(),1);
  absDc->apply(*oneVec,*singleParent,Teuchos::NO_TRANS);  // I think we need ghosted version
  Teuchos::ArrayRCP<SC> singleParentData = singleParent->getDataNonConst(0);

  Array<LO> inds;
  Array<SC> vals;

  size_t numRows  = temp2->getRowMap()->getLocalNumElements();

  for (LO row = 0; row < (LO) numRows; row++) {
     temp2->getLocalRowView(row, inds2, vals2);
      size_t nnz = inds2.size();
      if (nnz == 0)
        continue;

      inds.resize(inds2.size());
      vals.resize(vals2.size());

      size_t numInds = 0;
      for (size_t j = 0; j < nnz; j++) {
        if ( (vals2[j] == 2.0) || (   (singleParentData[inds2[j]] == 1.0) && (vals2[j] == 1.0)  ) ) {
          inds[numInds] = inds2[j];
          vals[numInds] = vals2[j];
          numInds++;
        }
      }
      if (numInds == 0) continue;
      inds.resize(numInds);
      vals.resize(numInds);

      // Because we used a column map in the construction of the matrix
      // we can just use insertLocalValues here instead of insertGlobalValues
      filtered->insertLocalValues(row, inds, vals);
  }
  RCP<ParameterList> fillCompleteParams(new ParameterList);
  fillCompleteParams->set("No Nonlocal Changes", true);
  filtered->fillComplete(temp2->getDomainMap(), temp2->getRangeMap(), fillCompleteParams);
//  Xpetra::IO<SC, LO, GO, NO>::Write("pattern.code", *filtered);

  Set(coarseLevel, "Ppattern", filtered->getCrsGraph());
}

}  // namespace MueLu

#endif  // MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
