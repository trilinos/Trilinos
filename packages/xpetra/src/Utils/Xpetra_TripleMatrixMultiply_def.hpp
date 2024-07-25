// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP_

#include "Xpetra_TripleMatrixMultiply_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TripleMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiplyRAP(const Matrix& R, bool transposeR,
                                                                                  const Matrix& A, bool transposeA,
                                                                                  const Matrix& P, bool transposeP,
                                                                                  Matrix& Ac,
                                                                                  bool call_FillComplete_on_result,
                                                                                  bool doOptimizeStorage,
                                                                                  const std::string& label,
                                                                                  const RCP<ParameterList>& params) {
  TEUCHOS_TEST_FOR_EXCEPTION(transposeR == false && Ac.getRowMap()->isSameAs(*R.getRowMap()) == false,
                             Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as row map of R");
  TEUCHOS_TEST_FOR_EXCEPTION(transposeR == true && Ac.getRowMap()->isSameAs(*R.getDomainMap()) == false,
                             Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as domain map of R");

  TEUCHOS_TEST_FOR_EXCEPTION(!R.isFillComplete(), Exceptions::RuntimeError, "R is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!P.isFillComplete(), Exceptions::RuntimeError, "P is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (Ac.getRowMap()->lib() == Xpetra::UseEpetra) {
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::TripleMatrixMultiply::MultiplyRAP is only implemented for Tpetra"));
  } else if (Ac.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    using helpers = Xpetra::Helpers<SC, LO, GO, NO>;
    if (helpers::isTpetraCrs(R) && helpers::isTpetraCrs(A) && helpers::isTpetraCrs(P)) {
      // All matrices are Crs
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpR = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(R);
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpP = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(P);
      Tpetra::CrsMatrix<SC, LO, GO, NO>& tpAc      = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(Ac);

      // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
      // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
      Tpetra::TripleMatrixMultiply::MultiplyRAP(tpR, transposeR, tpA, transposeA, tpP, transposeP, tpAc, haveMultiplyDoFillComplete, label, params);
    } else if (helpers::isTpetraBlockCrs(R) && helpers::isTpetraBlockCrs(A) && helpers::isTpetraBlockCrs(P)) {
      // All matrices are BlockCrs (except maybe Ac)
      // FIXME: For the moment we're just going to clobber the innards of Ac, so no reuse. Once we have a reuse kernel,
      // we'll need to think about refactoring BlockCrs so we can do something smarter here.

      if (!A.getRowMap()->getComm()->getRank())
        std::cout << "WARNING: Using inefficient BlockCrs Multiply Placeholder" << std::endl;

      const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpR = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(R);
      const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(A);
      const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpP = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(P);
      //          Tpetra::BlockCrsMatrix<SC,LO,GO,NO> &       tpAc = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraBlockCrs(Ac);

      using CRS           = Tpetra::CrsMatrix<SC, LO, GO, NO>;
      RCP<const CRS> Rcrs = Tpetra::convertToCrsMatrix(tpR);
      RCP<const CRS> Acrs = Tpetra::convertToCrsMatrix(tpA);
      RCP<const CRS> Pcrs = Tpetra::convertToCrsMatrix(tpP);
      //          RCP<CRS> Accrs = Tpetra::convertToCrsMatrix(tpAc);

      // FIXME: The lines below only works because we're assuming Ac is Point
      RCP<CRS> Accrs              = Teuchos::rcp(new CRS(Rcrs->getRowMap(), 0));
      const bool do_fill_complete = true;
      Tpetra::TripleMatrixMultiply::MultiplyRAP(*Rcrs, transposeR, *Acrs, transposeA, *Pcrs, transposeP, *Accrs, do_fill_complete, label, params);

      // Temporary output matrix
      RCP<Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > Ac_t       = Tpetra::convertToBlockCrsMatrix(*Accrs, A.GetStorageBlockSize());
      RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO> > Ac_x = Teuchos::rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(Ac_t));
      RCP<Xpetra::CrsMatrix<SC, LO, GO, NO> > Ac_p            = Ac_x;

      // We can now cheat and replace the innards of Ac
      RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> > Ac_w = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> >(Teuchos::rcpFromRef(Ac));
      Ac_w->replaceCrsMatrix(Ac_p);
    } else {
      // Mix and match
      TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "Mix-and-match Crs/BlockCrs Multiply not currently supported");
    }
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
    fillParams->set("Optimize Storage", doOptimizeStorage);
    Ac.fillComplete((transposeP) ? P.getRangeMap() : P.getDomainMap(),
                    (transposeR) ? R.getDomainMap() : R.getRangeMap(),
                    fillParams);
  }

  // transfer striding information
  RCP<const Map> domainMap = Teuchos::null;
  RCP<const Map> rangeMap  = Teuchos::null;

  const std::string stridedViewLabel("stridedMaps");
  const size_t blkSize = 1;
  std::vector<size_t> stridingInfo(1, blkSize);
  LocalOrdinal stridedBlockId = -1;

  if (R.IsView(stridedViewLabel)) {
    rangeMap = transposeR ? R.getColMap(stridedViewLabel) : R.getRowMap(stridedViewLabel);
  } else {
    rangeMap = transposeR ? R.getDomainMap() : R.getRangeMap();
    rangeMap = StridedMapFactory::Build(rangeMap, stridingInfo, stridedBlockId);
  }

  if (P.IsView(stridedViewLabel)) {
    domainMap = transposeP ? P.getRowMap(stridedViewLabel) : P.getColMap(stridedViewLabel);
  } else {
    domainMap = transposeP ? P.getRangeMap() : P.getDomainMap();
    domainMap = StridedMapFactory::Build(domainMap, stridingInfo, stridedBlockId);
  }
  Ac.CreateView(stridedViewLabel, rangeMap, domainMap);

}  // end Multiply

}  // end namespace Xpetra

#define XPETRA_TRIPLEMATRIXMULTIPLY_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP_ */
