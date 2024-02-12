// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_BLOCKEDRAPFACTORY_DEF_HPP
#define MUELU_BLOCKEDRAPFACTORY_DEF_HPP

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_BlockedRAPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_RAPFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedRAPFactory()
  : checkAc_(false)
  , repairZeroDiagonals_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase> >("R", null, "Restrictor factory");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  const Teuchos::ParameterList &pL = GetParameterList();
  if (pL.get<bool>("transpose: use implicit") == false)
    Input(coarseLevel, "R");

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {  // FIXME make fineLevel const!!
  FactoryMonitor m(*this, "Computing Ac (block)", coarseLevel);

  const ParameterList &pL = GetParameterList();

  RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P");

  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A);
  RCP<BlockedCrsMatrix> bP = rcp_dynamic_cast<BlockedCrsMatrix>(P);
  TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null() || bP.is_null(), Exceptions::BadCast, "Matrices A and P must be of type BlockedCrsMatrix.");

  RCP<BlockedCrsMatrix> bAP;
  RCP<BlockedCrsMatrix> bAc;
  {
    SubFactoryMonitor subM(*this, "MxM: A x P", coarseLevel);

    // Triple matrix product for BlockedCrsMatrixClass
    TEUCHOS_TEST_FOR_EXCEPTION((bA->Cols() != bP->Rows()), Exceptions::BadCast,
                               "Block matrix dimensions do not match: "
                               "A is "
                                   << bA->Rows() << "x" << bA->Cols() << "P is " << bP->Rows() << "x" << bP->Cols());

    bAP = MatrixMatrix::TwoMatrixMultiplyBlock(*bA, false, *bP, false, GetOStream(Statistics2), true, true);
  }

  // If we do not modify matrix later, allow optimization of storage.
  // This is necessary for new faster Epetra MM kernels.
  bool doOptimizeStorage = !checkAc_;

  const bool doTranspose    = true;
  const bool doFillComplete = true;
  if (pL.get<bool>("transpose: use implicit") == true) {
    SubFactoryMonitor m2(*this, "MxM: P' x (AP) (implicit)", coarseLevel);
    bAc = MatrixMatrix::TwoMatrixMultiplyBlock(*bP, doTranspose, *bAP, !doTranspose, GetOStream(Statistics2), doFillComplete, doOptimizeStorage);

  } else {
    RCP<Matrix> R            = Get<RCP<Matrix> >(coarseLevel, "R");
    RCP<BlockedCrsMatrix> bR = rcp_dynamic_cast<BlockedCrsMatrix>(R);
    TEUCHOS_TEST_FOR_EXCEPTION(bR.is_null(), Exceptions::BadCast, "Matrix R must be of type BlockedCrsMatrix.");

    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != bR->Cols(), Exceptions::BadCast,
                               "Block matrix dimensions do not match: "
                               "R is "
                                   << bR->Rows() << "x" << bR->Cols() << "A is " << bA->Rows() << "x" << bA->Cols());

    SubFactoryMonitor m2(*this, "MxM: R x (AP) (explicit)", coarseLevel);
    bAc = MatrixMatrix::TwoMatrixMultiplyBlock(*bR, !doTranspose, *bAP, !doTranspose, GetOStream(Statistics2), doFillComplete, doOptimizeStorage);
  }

  if (checkAc_)
    CheckMainDiagonal(bAc);

  GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*bAc, "Ac (blocked)");

  Set<RCP<Matrix> >(coarseLevel, "A", bAc);

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m1(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;

      GetOStream(Runtime0) << "BlockRAPFactory: call transfer factory: " << fac->description() << std::endl;

      fac->CallBuild(coarseLevel);

      // AP (11/11/13): I am not sure exactly why we need to call Release, but we do need it to get rid
      // of dangling data for CoordinatesTransferFactory
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckMainDiagonal(RCP<BlockedCrsMatrix> &bAc, bool repairZeroDiagonals) {
  RCP<Matrix> c00  = bAc->getMatrix(0, 0);
  RCP<Matrix> Aout = MatrixFactory::Build(c00->getRowMap(), c00->getGlobalMaxNumRowEntries());

  RCP<Vector> diagVec = VectorFactory::Build(c00->getRowMap());
  c00->getLocalDiagCopy(*diagVec);
  ArrayRCP<SC> diagVal = diagVec->getDataNonConst(0);

  // loop over local rows
  for (size_t row = 0; row < c00->getLocalNumRows(); row++) {
    // get global row id
    GO grid = c00->getRowMap()->getGlobalElement(row);  // global row id

    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    c00->getLocalRowView(row, indices, vals);

    // just copy all values in output
    ArrayRCP<GO> indout(indices.size(), Teuchos::OrdinalTraits<GO>::zero());
    ArrayRCP<SC> valout(indices.size(), Teuchos::ScalarTraits<SC>::zero());

    // just copy values
    for (size_t i = 0; i < as<size_t>(indices.size()); i++) {
      GO gcid   = c00->getColMap()->getGlobalElement(indices[i]);  // LID -> GID (column)
      indout[i] = gcid;
      valout[i] = vals[i];
    }

    Aout->insertGlobalValues(grid, indout.view(0, indout.size()), valout.view(0, valout.size()));
    if (diagVal[row] == Teuchos::ScalarTraits<SC>::zero() && repairZeroDiagonals) {
      // always overwrite diagonal entry
      Aout->insertGlobalValues(grid, Teuchos::tuple<GO>(grid), Teuchos::tuple<SC>(1.0));
    }
  }

  Aout->fillComplete(c00->getDomainMap(), c00->getRangeMap());

  bAc->setMatrix(0, 0, Aout);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase> &factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. "
                             "(Note: you can remove this exception if there's a good reason for)");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_BLOCKEDRAPFACTORY_SHORT
#endif  // MUELU_BLOCKEDRAPFACTORY_DEF_HPP

// TODO add plausibility check
// TODO add CheckMainDiagonal for Blocked operator
// Avoid copying block matrix!
// create new empty Operator
