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
#ifndef THYRA_MUELU_TPETRA_Q2Q1PRECONDITIONER_FACTORY_DEF_HPP
#define THYRA_MUELU_TPETRA_Q2Q1PRECONDITIONER_FACTORY_DEF_HPP

#ifdef HAVE_MUELU_EXPERIMENTAL

#include "Thyra_MueLuTpetraQ2Q1PreconditionerFactory_decl.hpp"

#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

#include <Teuchos_Ptr.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerbosityLevel.hpp>

#include <Teko_Utilities.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu.hpp"

#include "../../research/q2q1/MueLu_Q2Q1PFactory.hpp"
#include "../../research/q2q1/MueLu_Q2Q1uPFactory.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_BlockedDirectSolver.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_TrilinosSmoother.hpp"

#include <string>

namespace Thyra {

#define MUELU_GPD(name, type, defaultValue) \
  (paramList.isParameter(name) ? paramList.get<type>(name) : defaultValue)

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// Constructors/initializers/accessors
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MueLuTpetraQ2Q1PreconditionerFactory() {}

// Overridden from PreconditionerFactoryBase
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
  typedef Thyra ::TpetraLinearOp<SC, LO, GO, NO> ThyraTpetraLinOp;
  typedef Tpetra::Operator<SC, LO, GO, NO> TpetraLinOp;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TpetraCrsMat;

  const RCP<const LinearOpBase<SC> > fwdOp           = fwdOpSrc.getOp();
  const RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);
  const RCP<const TpetraLinOp> tpetraFwdOp           = Teuchos::nonnull(thyraTpetraFwdOp) ? thyraTpetraFwdOp->getConstTpetraOperator() : Teuchos::null;
  const RCP<const TpetraCrsMat> tpetraFwdCrsMat      = rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);

  return Teuchos::nonnull(tpetraFwdCrsMat);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<PreconditionerBase<Scalar> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createPrec() const {
  return rcp(new DefaultPreconditioner<SC>);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    initializePrec(const RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc, PreconditionerBase<Scalar>* prec, const ESupportSolveUse supportSolveUse) const {
  // Check precondition
  TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
  TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
  TEUCHOS_ASSERT(prec);

  // Retrieve wrapped concrete Tpetra matrix from FwdOp
  const RCP<const LinearOpBase<SC> > fwdOp = fwdOpSrc->getOp();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

  typedef Thyra::TpetraLinearOp<SC, LO, GO, NO> ThyraTpetraLinOp;
  const RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraTpetraFwdOp));

  typedef Tpetra::Operator<SC, LO, GO, NO> TpetraLinOp;
  const RCP<const TpetraLinOp> tpetraFwdOp = thyraTpetraFwdOp->getConstTpetraOperator();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdOp));

  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TpetraCrsMat;
  const RCP<const TpetraCrsMat> tpetraFwdCrsMat = rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdCrsMat));

  // Retrieve concrete preconditioner object
  const Teuchos::Ptr<DefaultPreconditioner<SC> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<SC>*>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // Workaround since MueLu interface does not accept const matrix as input
  const RCP<TpetraCrsMat> tpetraFwdCrsMatNonConst = rcp_const_cast<TpetraCrsMat>(tpetraFwdCrsMat);

  // Create and compute the initial preconditioner

  // Create a copy, as we may remove some things from the list
  ParameterList paramList = *paramList_;

  typedef Tpetra::MultiVector<SC, LO, GO, NO> MultiVector;
  RCP<MultiVector> coords, nullspace, velCoords, presCoords;
  ArrayRCP<LO> p2vMap;
  Teko::LinearOp thA11, thA12, thA21, thA11_9Pt;
  if (paramList.isType<RCP<MultiVector> >("Coordinates")) {
    coords = paramList.get<RCP<MultiVector> >("Coordinates");
    paramList.remove("Coordinates");
  }
  if (paramList.isType<RCP<MultiVector> >("Nullspace")) {
    nullspace = paramList.get<RCP<MultiVector> >("Nullspace");
    paramList.remove("Nullspace");
  }
  if (paramList.isType<RCP<MultiVector> >("Velcoords")) {
    velCoords = paramList.get<RCP<MultiVector> >("Velcoords");
    paramList.remove("Velcoords");
  }
  if (paramList.isType<RCP<MultiVector> >("Prescoords")) {
    presCoords = paramList.get<RCP<MultiVector> >("Prescoords");
    paramList.remove("Prescoords");
  }
  if (paramList.isType<ArrayRCP<LO> >("p2vMap")) {
    p2vMap = paramList.get<ArrayRCP<LO> >("p2vMap");
    paramList.remove("p2vMap");
  }
  if (paramList.isType<Teko::LinearOp>("A11")) {
    thA11 = paramList.get<Teko::LinearOp>("A11");
    paramList.remove("A11");
  }
  if (paramList.isType<Teko::LinearOp>("A12")) {
    thA12 = paramList.get<Teko::LinearOp>("A12");
    paramList.remove("A12");
  }
  if (paramList.isType<Teko::LinearOp>("A21")) {
    thA21 = paramList.get<Teko::LinearOp>("A21");
    paramList.remove("A21");
  }
  if (paramList.isType<Teko::LinearOp>("A11_9Pt")) {
    thA11_9Pt = paramList.get<Teko::LinearOp>("A11_9Pt");
    paramList.remove("A11_9Pt");
  }

  typedef MueLu::TpetraOperator<SC, LO, GO, NO> MueLuOperator;
  const RCP<MueLuOperator> mueluPrecOp = Q2Q1MkPrecond(paramList, velCoords, presCoords, p2vMap, thA11, thA12, thA21, thA11_9Pt);

  const RCP<LinearOpBase<SC> > thyraPrecOp = Thyra::createLinearOp(RCP<TpetraLinOp>(mueluPrecOp));
  defaultPrec->initializeUnspecified(thyraPrecOp);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    uninitializePrec(PreconditionerBase<Scalar>* prec, RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
  // Check precondition
  TEUCHOS_ASSERT(prec);

  // Retrieve concrete preconditioner object
  const Teuchos::Ptr<DefaultPreconditioner<SC> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<SC>*>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  if (fwdOp) {
    // TODO: Implement properly instead of returning default value
    *fwdOp = Teuchos::null;
  }

  if (supportSolveUse) {
    // TODO: Implement properly instead of returning default value
    *supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED;
  }

  defaultPrec->uninitialize();
}

// Overridden from ParameterListAcceptor
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameterList(const RCP<ParameterList>& paramList) {
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
  paramList_ = paramList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNonconstParameterList() {
  return paramList_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::unsetParameterList() {
  RCP<ParameterList> savedParamList = paramList_;
  paramList_                        = Teuchos::null;
  return savedParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getParameterList() const {
  return paramList_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getValidParameters() const {
  static RCP<const ParameterList> validPL;

  if (validPL.is_null())
    validPL = rcp(new ParameterList());

  return validPL;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Q2Q1MkPrecond(const ParameterList& paramList,
                  const RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& velCoords,
                  const RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& presCoords,
                  const ArrayRCP<LocalOrdinal>& p2vMap,
                  const Teko::LinearOp& thA11, const Teko::LinearOp& thA12, const Teko::LinearOp& thA21, const Teko::LinearOp& thA11_9Pt) const {
  using Teuchos::null;

  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TP_Crs;
  typedef Tpetra::Operator<SC, LO, GO, NO> TP_Op;

  typedef Xpetra::BlockedCrsMatrix<SC, LO, GO, NO> BlockedCrsMatrix;
  typedef Xpetra::CrsMatrix<SC, LO, GO, NO> CrsMatrix;
  typedef Xpetra::CrsMatrixWrap<SC, LO, GO, NO> CrsMatrixWrap;
  typedef Xpetra::MapExtractorFactory<SC, LO, GO, NO> MapExtractorFactory;
  typedef Xpetra::MapExtractor<SC, LO, GO, NO> MapExtractor;
  typedef Xpetra::Map<LO, GO, NO> Map;
  typedef Xpetra::MapFactory<LO, GO, NO> MapFactory;
  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;
  typedef Xpetra::MatrixFactory<SC, LO, GO, NO> MatrixFactory;
  typedef Xpetra::StridedMapFactory<LO, GO, NO> StridedMapFactory;

  typedef MueLu::Hierarchy<SC, LO, GO, NO> Hierarchy;

  const RCP<const Teuchos::Comm<int> > comm = velCoords->getMap()->getComm();

  // Pull out Tpetra matrices
  RCP<Thyra::LinearOpBase<SC> > ThNonConstA11     = rcp_const_cast<Thyra::LinearOpBase<double> >(thA11);
  RCP<Thyra::LinearOpBase<SC> > ThNonConstA21     = rcp_const_cast<Thyra::LinearOpBase<double> >(thA21);
  RCP<Thyra::LinearOpBase<SC> > ThNonConstA12     = rcp_const_cast<Thyra::LinearOpBase<double> >(thA12);
  RCP<Thyra::LinearOpBase<SC> > ThNonConstA11_9Pt = rcp_const_cast<Thyra::LinearOpBase<double> >(thA11_9Pt);

  RCP<TP_Op> TpetA11     = Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, NO>::getTpetraOperator(ThNonConstA11);
  RCP<TP_Op> TpetA21     = Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, NO>::getTpetraOperator(ThNonConstA21);
  RCP<TP_Op> TpetA12     = Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, NO>::getTpetraOperator(ThNonConstA12);
  RCP<TP_Op> TpetA11_9Pt = Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, NO>::getTpetraOperator(ThNonConstA11_9Pt);

  RCP<TP_Crs> TpetCrsA11     = rcp_dynamic_cast<TP_Crs>(TpetA11);
  RCP<TP_Crs> TpetCrsA21     = rcp_dynamic_cast<TP_Crs>(TpetA21);
  RCP<TP_Crs> TpetCrsA12     = rcp_dynamic_cast<TP_Crs>(TpetA12);
  RCP<TP_Crs> TpetCrsA11_9Pt = rcp_dynamic_cast<TP_Crs>(TpetA11_9Pt);

  RCP<Matrix> A_11     = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA11);
  RCP<Matrix> tmp_A_21 = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA21);  // needs map modification
  RCP<Matrix> tmp_A_12 = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA12);  // needs map modification
  RCP<Matrix> A_11_9Pt = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA11_9Pt);

  Xpetra::global_size_t numVel  = A_11->getRowMap()->getLocalNumElements();
  Xpetra::global_size_t numPres = tmp_A_21->getRowMap()->getLocalNumElements();

  // Create new A21 with map so that the global indices of the row map starts
  // from numVel+1 (where numVel is the number of rows in the A11 block)
  RCP<const Map> domainMap2       = tmp_A_12->getDomainMap();
  RCP<const Map> rangeMap2        = tmp_A_21->getRangeMap();
  Xpetra::global_size_t numRows2  = rangeMap2->getLocalNumElements();
  Xpetra::global_size_t numCols2  = domainMap2->getLocalNumElements();
  ArrayView<const GO> rangeElem2  = rangeMap2->getLocalElementList();
  ArrayView<const GO> domainElem2 = domainMap2->getLocalElementList();
  ArrayView<const GO> rowElem1    = tmp_A_12->getRowMap()->getLocalElementList();
  ArrayView<const GO> colElem1    = tmp_A_21->getColMap()->getLocalElementList();

  Xpetra::UnderlyingLib lib = domainMap2->lib();
  GO indexBase              = domainMap2->getIndexBase();

  Array<GO> newRowElem2(numRows2, 0);
  for (Xpetra::global_size_t i = 0; i < numRows2; i++)
    newRowElem2[i] = numVel + rangeElem2[i];

  RCP<const Map> newRangeMap2 = MapFactory::Build(lib, numRows2, newRowElem2, indexBase, comm);

  // maybe should be column map???
  Array<GO> newColElem2(numCols2, 0);
  for (Xpetra::global_size_t i = 0; i < numCols2; i++)
    newColElem2[i] = numVel + domainElem2[i];

  RCP<const Map> newDomainMap2 = MapFactory::Build(lib, numCols2, newColElem2, indexBase, comm);

  RCP<Matrix> A_12 = MatrixFactory::Build(tmp_A_12->getRangeMap(), newDomainMap2, tmp_A_12->getLocalMaxNumRowEntries());
  RCP<Matrix> A_21 = MatrixFactory::Build(newRangeMap2, tmp_A_21->getDomainMap(), tmp_A_21->getLocalMaxNumRowEntries());

  RCP<CrsMatrix> A_11_crs     = rcp_dynamic_cast<CrsMatrixWrap>(A_11)->getCrsMatrix();
  RCP<CrsMatrix> A_12_crs     = rcp_dynamic_cast<CrsMatrixWrap>(A_12)->getCrsMatrix();
  RCP<CrsMatrix> A_21_crs     = rcp_dynamic_cast<CrsMatrixWrap>(A_21)->getCrsMatrix();
  RCP<CrsMatrix> A_11_crs_9Pt = rcp_dynamic_cast<CrsMatrixWrap>(A_11_9Pt)->getCrsMatrix();

#if 0
    RCP<Matrix>    A_22         = MatrixFactory::Build(newRangeMap2,            newDomainMap2,            1);
    RCP<CrsMatrix> A_22_crs     = rcp_dynamic_cast<CrsMatrixWrap>(A_22)    ->getCrsMatrix();

    // FIXME: why do we need to perturb A_22?
    Array<SC> smallVal(1, 1.0e-10);

    // FIXME: could this be sped up using expertStaticFillComplete?
    // There was an attempt on doing it, but it did not do the proper thing
    // with empty columns. See git history
    ArrayView<const LO> inds;
    ArrayView<const SC> vals;
    for (LO row = 0; row < as<LO>(numRows2); ++row) {
      tmp_A_21->getLocalRowView(row, inds, vals);

      size_t nnz = inds.size();
      Array<GO> newInds(nnz, 0);
      for (LO colID = 0; colID < as<LO>(nnz); colID++)
        newInds[colID] = colElem1[inds[colID]];

      A_21_crs->insertGlobalValues(newRowElem2[row], newInds,                        vals);
      A_22_crs->insertGlobalValues(newRowElem2[row], Array<LO>(1, newRowElem2[row]), smallVal);
    }
    A_21_crs->fillComplete(tmp_A_21->getDomainMap(), newRangeMap2);
    A_22_crs->fillComplete(newDomainMap2,            newRangeMap2);
#else
  RCP<Matrix> A_22        = Teuchos::null;
  RCP<CrsMatrix> A_22_crs = Teuchos::null;

  ArrayView<const LO> inds;
  ArrayView<const SC> vals;
  for (LO row = 0; row < as<LO>(numRows2); ++row) {
    tmp_A_21->getLocalRowView(row, inds, vals);

    size_t nnz = inds.size();
    Array<GO> newInds(nnz, 0);
    for (LO colID = 0; colID < as<LO>(nnz); colID++)
      newInds[colID] = colElem1[inds[colID]];

    A_21_crs->insertGlobalValues(newRowElem2[row], newInds, vals);
  }
  A_21_crs->fillComplete(tmp_A_21->getDomainMap(), newRangeMap2);
#endif

  // Create new A12 with map so that the global indices of the ColMap starts
  // from numVel+1 (where numVel is the number of rows in the A11 block)
  for (LO row = 0; row < as<LO>(tmp_A_12->getRowMap()->getLocalNumElements()); ++row) {
    tmp_A_12->getLocalRowView(row, inds, vals);

    size_t nnz = inds.size();
    Array<GO> newInds(nnz, 0);
    for (LO colID = 0; colID < as<LO>(nnz); colID++)
      newInds[colID] = newColElem2[inds[colID]];

    A_12_crs->insertGlobalValues(rowElem1[row], newInds, vals);
  }
  A_12_crs->fillComplete(newDomainMap2, tmp_A_12->getRangeMap());

  RCP<Matrix> A_12_abs = Absolute(*A_12);
  RCP<Matrix> A_21_abs = Absolute(*A_21);

  // =========================================================================
  // Preconditioner construction - I (block)
  // =========================================================================
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out       = *fancy;
  out.setOutputToRootOnly(0);
  RCP<Matrix> BBt     = Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*A_21, false, *A_12, false, out);
  RCP<Matrix> BBt_abs = Xpetra::MatrixMatrix<SC, LO, GO, NO>::Multiply(*A_21_abs, false, *A_12_abs, false, out);

  SC dropTol            = (paramList.get<int>("useFilters") ? paramList.get<double>("tau_1") : 0.00);
  RCP<Matrix> filteredA = FilterMatrix(*A_11, *A_11, dropTol);
  RCP<Matrix> filteredB = FilterMatrix(*BBt, *BBt_abs, dropTol);

  RCP<Matrix> fA_11_crs = rcp_dynamic_cast<CrsMatrixWrap>(filteredA);
  RCP<Matrix> fA_12_crs = Teuchos::null;
  RCP<Matrix> fA_21_crs = Teuchos::null;
  RCP<Matrix> fA_22_crs = rcp_dynamic_cast<CrsMatrixWrap>(filteredB);

  // Build the large filtered matrix which requires strided maps
  std::vector<size_t> stridingInfo(1, 1);
  int stridedBlockId = -1;

  Array<GO> elementList(numVel + numPres);  // Not RCP ...  does this get cleared ?
  Array<GO> velElem  = A_12_crs->getRangeMap()->getLocalElementList();
  Array<GO> presElem = A_21_crs->getRangeMap()->getLocalElementList();

  for (Xpetra::global_size_t i = 0; i < numVel; i++) elementList[i] = velElem[i];
  for (Xpetra::global_size_t i = numVel; i < numVel + numPres; i++) elementList[i] = presElem[i - numVel];
  RCP<const Map> fullMap = StridedMapFactory::Build(Xpetra::UseTpetra, numVel + numPres, elementList(), indexBase, stridingInfo, comm);

  std::vector<RCP<const Map> > partMaps(2);
  partMaps[0] = StridedMapFactory::Build(Xpetra::UseTpetra, numVel, velElem, indexBase, stridingInfo, comm);
  partMaps[1] = StridedMapFactory::Build(Xpetra::UseTpetra, numPres, presElem, indexBase, stridingInfo, comm, stridedBlockId, numVel);

  // Map extractors are necessary for Xpetra's block operators
  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(fullMap, partMaps);
  RCP<BlockedCrsMatrix> fA             = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
  fA->setMatrix(0, 0, fA_11_crs);
  fA->setMatrix(0, 1, fA_12_crs);
  fA->setMatrix(1, 0, fA_21_crs);
  fA->setMatrix(1, 1, fA_22_crs);
  fA->fillComplete();

  // -------------------------------------------------------------------------
  // Preconditioner construction - I.a (filtered hierarchy)
  // -------------------------------------------------------------------------
  MueLu::FactoryManager<SC, LO, GO, NO> M;
  SetDependencyTree(M, paramList);

  RCP<Hierarchy> H              = rcp(new Hierarchy);
  RCP<MueLu::Level> finestLevel = H->GetLevel(0);
  finestLevel->Set("A", rcp_dynamic_cast<Matrix>(fA));
  finestLevel->Set("p2vMap", p2vMap);
  finestLevel->Set("CoordinatesVelocity", Xpetra::toXpetra(velCoords));
  finestLevel->Set("CoordinatesPressure", Xpetra::toXpetra(presCoords));
  finestLevel->Set("AForPat", A_11_9Pt);
  H->SetMaxCoarseSize(MUELU_GPD("coarse: max size", int, 1));

  // The first invocation of Setup() builds the hierarchy using the filtered
  // matrix. This build includes the grid transfers but not the creation of the
  // smoothers.
  // NOTE: we need to indicate what should be kept from the first invocation
  // for the second invocation, which then focuses on building the smoothers
  // for the unfiltered matrix.
  H->Keep("P", M.GetFactory("P").get());
  H->Keep("R", M.GetFactory("R").get());
  H->Keep("Ptent", M.GetFactory("Ptent").get());
  H->Setup(M, 0, MUELU_GPD("max levels", int, 3));

#if 0
    for (int i = 1; i < H->GetNumLevels(); i++) {
      RCP<Matrix>           P     = H->GetLevel(i)->template Get<RCP<Matrix> >("P");
      RCP<BlockedCrsMatrix> Pcrs  = rcp_dynamic_cast<BlockedCrsMatrix>(P);
      RCP<Matrix>           Pp    = Pcrs->getMatrix(1,1);
      RCP<Matrix>           Pv    = Pcrs->getMatrix(0,0);

      Xpetra::IO<SC,LO,GO,NO>::Write("Pp_l" + MueLu::toString(i) + ".mm", *Pp);
      Xpetra::IO<SC,LO,GO,NO>::Write("Pv_l" + MueLu::toString(i) + ".mm", *Pv);
    }
#endif

  // -------------------------------------------------------------------------
  // Preconditioner construction - I.b (smoothers for unfiltered matrix)
  // -------------------------------------------------------------------------
  std::string smootherType = MUELU_GPD("smoother: type", std::string, "vanka");
  ParameterList smootherParams;
  if (paramList.isSublist("smoother: params"))
    smootherParams = paramList.sublist("smoother: params");
  M.SetFactory("Smoother", GetSmoother(smootherType, smootherParams, false /*coarseSolver?*/));

  std::string coarseType = MUELU_GPD("coarse: type", std::string, "direct");
  ParameterList coarseParams;
  if (paramList.isSublist("coarse: params"))
    coarseParams = paramList.sublist("coarse: params");
  M.SetFactory("CoarseSolver", GetSmoother(coarseType, coarseParams, true /*coarseSolver?*/));

#ifdef HAVE_MUELU_DEBUG
  M.ResetDebugData();
#endif

  RCP<BlockedCrsMatrix> A = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
  A->setMatrix(0, 0, A_11);
  A->setMatrix(0, 1, A_12);
  A->setMatrix(1, 0, A_21);
  A->setMatrix(1, 1, A_22);
  A->fillComplete();

  H->GetLevel(0)->Set("A", rcp_dynamic_cast<Matrix>(A));

  H->Setup(M, 0, H->GetNumLevels());

  return rcp(new MueLu::TpetraOperator<SC, LO, GO, NO>(H));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    FilterMatrix(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pattern, Scalar dropTol) const {
  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;
  typedef MueLu::AmalgamationFactory<SC, LO, GO, NO> AmalgamationFactory;
  typedef MueLu::CoalesceDropFactory<SC, LO, GO, NO> CoalesceDropFactory;
  typedef MueLu::FilteredAFactory<SC, LO, GO, NO> FilteredAFactory;
  typedef MueLu::GraphBase<LO, GO, NO> GraphBase;

  RCP<GraphBase> filteredGraph;
  {
    // Get graph pattern for the pattern matrix
    MueLu::Level level;
    level.SetLevelID(1);

    level.Set<RCP<Matrix> >("A", rcpFromRef(Pattern));

    RCP<AmalgamationFactory> amalgFactory = rcp(new AmalgamationFactory());

    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams             = *(dropFactory->GetValidParameterList());
    dropParams.set("lightweight wrap", true);
    dropParams.set("aggregation: drop scheme", "classical");
    dropParams.set("aggregation: drop tol", dropTol);
    // dropParams.set("Dirichlet detection threshold", <>);
    dropFactory->SetParameterList(dropParams);
    dropFactory->SetFactory("UnAmalgamationInfo", amalgFactory);

    // Build
    level.Request("Graph", dropFactory.get());
    dropFactory->Build(level);

    level.Get("Graph", filteredGraph, dropFactory.get());
  }

  RCP<Matrix> filteredA;
  {
    // Filter the original matrix, not the pattern one
    MueLu::Level level;
    level.SetLevelID(1);

    level.Set("A", rcpFromRef(A));
    level.Set("Graph", filteredGraph);
    level.Set("Filtering", true);

    RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
    ParameterList filterParams          = *(filterFactory->GetValidParameterList());
    // We need a graph that has proper structure in it. Therefore, we need to
    // drop older pattern, i.e. not to reuse it
    filterParams.set("filtered matrix: reuse graph", false);
    filterParams.set("filtered matrix: use lumping", false);
    filterFactory->SetParameterList(filterParams);

    // Build
    level.Request("A", filterFactory.get());
    filterFactory->Build(level);

    level.Get("A", filteredA, filterFactory.get());
  }

  // Zero out row sums by fixing the diagonal
  filteredA->resumeFill();
  size_t numRows = filteredA->getRowMap()->getLocalNumElements();
  for (size_t i = 0; i < numRows; i++) {
    ArrayView<const LO> inds;
    ArrayView<const SC> vals;
    filteredA->getLocalRowView(i, inds, vals);

    size_t nnz = inds.size();

    Array<SC> valsNew = vals;

    LO diagIndex = -1;
    SC diag      = Teuchos::ScalarTraits<SC>::zero();
    for (size_t j = 0; j < nnz; j++) {
      diag += vals[j];
      if (inds[j] == Teuchos::as<int>(i))
        diagIndex = j;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(diagIndex == -1, MueLu::Exceptions::RuntimeError,
                               "No diagonal found");
    if (nnz <= 1)
      continue;

    valsNew[diagIndex] -= diag;

    filteredA->replaceLocalValues(i, inds, valsNew);
  }
  filteredA->fillComplete();

  return filteredA;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    SetDependencyTree(MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M, const ParameterList& paramList) const {
  typedef MueLu::BlockedPFactory<SC, LO, GO, NO> BlockedPFactory;
  typedef MueLu::GenericRFactory<SC, LO, GO, NO> GenericRFactory;
  typedef MueLu::BlockedRAPFactory<SC, LO, GO, NO> BlockedRAPFactory;
  typedef MueLu::SmootherFactory<SC, LO, GO, NO> SmootherFactory;
  typedef MueLu::BlockedDirectSolver<SC, LO, GO, NO> BlockedDirectSolver;
  typedef MueLu::FactoryManager<SC, LO, GO, NO> FactoryManager;

  // Pressure and velocity dependency trees are identical. The only
  // difference is that pressure has to go first, so that velocity can use
  // some of pressure data
  RCP<FactoryManager> M11 = rcp(new FactoryManager()), M22 = rcp(new FactoryManager());
  M11->SetKokkosRefactor(paramList.get<bool>("use kokkos refactor"));
  M22->SetKokkosRefactor(paramList.get<bool>("use kokkos refactor"));
  SetBlockDependencyTree(*M11, 0, 0, "velocity", paramList);
  SetBlockDependencyTree(*M22, 1, 1, "pressure", paramList);

  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
  ParameterList pParamList   = *(PFact->GetValidParameterList());
  pParamList.set("backwards", true);  // do pressure first
  PFact->SetParameterList(pParamList);
  PFact->AddFactoryManager(M11);
  PFact->AddFactoryManager(M22);
  M.SetFactory("P", PFact);

  RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
  RFact->SetFactory("P", PFact);
  M.SetFactory("R", RFact);

  RCP<MueLu::Factory> AcFact = rcp(new BlockedRAPFactory());
  AcFact->SetFactory("R", RFact);
  AcFact->SetFactory("P", PFact);
  M.SetFactory("A", AcFact);

  // Smoothers will be set later
  M.SetFactory("Smoother", Teuchos::null);

  RCP<MueLu::Factory> coarseFact = rcp(new SmootherFactory(rcp(new BlockedDirectSolver()), Teuchos::null));
  // M.SetFactory("CoarseSolver", coarseFact);
  M.SetFactory("CoarseSolver", Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    SetBlockDependencyTree(MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M, LocalOrdinal row, LocalOrdinal col, const std::string& mode, const ParameterList& paramList) const {
  typedef MueLu::ConstraintFactory<SC, LO, GO, NO> ConstraintFactory;
  typedef MueLu::EminPFactory<SC, LO, GO, NO> EminPFactory;
  typedef MueLu::GenericRFactory<SC, LO, GO, NO> GenericRFactory;
  typedef MueLu::PatternFactory<SC, LO, GO, NO> PatternFactory;
  typedef MueLu::Q2Q1PFactory<SC, LO, GO, NO> Q2Q1PFactory;
  typedef MueLu::Q2Q1uPFactory<SC, LO, GO, NO> Q2Q1uPFactory;
  typedef MueLu::SubBlockAFactory<SC, LO, GO, NO> SubBlockAFactory;

  RCP<SubBlockAFactory> AFact = rcp(new SubBlockAFactory());
  AFact->SetFactory("A", MueLu::NoFactory::getRCP());
  AFact->SetParameter("block row", Teuchos::ParameterEntry(row));
  AFact->SetParameter("block col", Teuchos::ParameterEntry(col));
  M.SetFactory("A", AFact);

  RCP<MueLu::Factory> Q2Q1Fact;

  const bool isStructured = false;

  if (isStructured) {
    Q2Q1Fact = rcp(new Q2Q1PFactory);

  } else {
    Q2Q1Fact                    = rcp(new Q2Q1uPFactory);
    ParameterList q2q1ParamList = *(Q2Q1Fact->GetValidParameterList());
    q2q1ParamList.set("mode", mode);
    if (paramList.isParameter("dump status"))
      q2q1ParamList.set("dump status", paramList.get<bool>("dump status"));
    if (paramList.isParameter("phase2"))
      q2q1ParamList.set("phase2", paramList.get<bool>("phase2"));
    if (paramList.isParameter("tau_2"))
      q2q1ParamList.set("tau_2", paramList.get<double>("tau_2"));
    Q2Q1Fact->SetParameterList(q2q1ParamList);
  }
  Q2Q1Fact->SetFactory("A", AFact);
  M.SetFactory("Ptent", Q2Q1Fact);

  RCP<PatternFactory> patternFact = rcp(new PatternFactory);
  ParameterList patternParams     = *(patternFact->GetValidParameterList());
  // Our prolongator constructs the exact pattern we are going to use,
  // therefore we do not expand it
  patternParams.set("emin: pattern order", 0);
  patternFact->SetParameterList(patternParams);
  patternFact->SetFactory("A", AFact);
  patternFact->SetFactory("P", Q2Q1Fact);
  M.SetFactory("Ppattern", patternFact);

  RCP<ConstraintFactory> CFact = rcp(new ConstraintFactory);
  CFact->SetFactory("Ppattern", patternFact);
  M.SetFactory("Constraint", CFact);

  RCP<EminPFactory> EminPFact = rcp(new EminPFactory());
  ParameterList eminParams    = *(EminPFact->GetValidParameterList());
  if (paramList.isParameter("emin: num iterations"))
    eminParams.set("emin: num iterations", paramList.get<int>("emin: num iterations"));
  if (mode == "pressure") {
    eminParams.set("emin: iterative method", "cg");
  } else {
    eminParams.set("emin: iterative method", "gmres");
    if (paramList.isParameter("emin: iterative method"))
      eminParams.set("emin: iterative method", paramList.get<std::string>("emin: iterative method"));
  }
  EminPFact->SetParameterList(eminParams);
  EminPFact->SetFactory("A", AFact);
  EminPFact->SetFactory("Constraint", CFact);
  EminPFact->SetFactory("P", Q2Q1Fact);
  M.SetFactory("P", EminPFact);

  if (mode == "velocity" && (!paramList.isParameter("velocity: use transpose") || paramList.get<bool>("velocity: use transpose") == false)) {
    // Pressure system is symmetric, so it does not matter
    // Velocity system may benefit from running emin in restriction mode (with A^T)
    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
    RFact->SetFactory("P", EminPFact);
    M.SetFactory("R", RFact);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::FactoryBase>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetSmoother(const std::string& type, const ParameterList& paramList, bool coarseSolver) const {
  typedef Teuchos::ParameterEntry ParameterEntry;

  typedef MueLu::BlockedDirectSolver<SC, LO, GO, NO> BlockedDirectSolver;
  typedef MueLu::BraessSarazinSmoother<SC, LO, GO, NO> BraessSarazinSmoother;
  typedef MueLu::DirectSolver<SC, LO, GO, NO> DirectSolver;
  typedef MueLu::FactoryManager<SC, LO, GO, NO> FactoryManager;
  typedef MueLu::SchurComplementFactory<SC, LO, GO, NO> SchurComplementFactory;
  typedef MueLu::SmootherFactory<SC, LO, GO, NO> SmootherFactory;
  typedef MueLu::SmootherPrototype<SC, LO, GO, NO> SmootherPrototype;
  typedef MueLu::TrilinosSmoother<SC, LO, GO, NO> TrilinosSmoother;

  RCP<SmootherPrototype> smootherPrototype;
  if (type == "none") {
    return Teuchos::null;

  } else if (type == "vanka") {
    // Set up Vanka smoothing via a combination of Schwarz and block relaxation.
    ParameterList schwarzList;
    schwarzList.set("schwarz: overlap level", as<int>(0));
    schwarzList.set("schwarz: zero starting solution", false);
    schwarzList.set("subdomain solver name", "Block_Relaxation");

    ParameterList& innerSolverList = schwarzList.sublist("subdomain solver parameters");
    innerSolverList.set("partitioner: type", "user");
    innerSolverList.set("partitioner: overlap", MUELU_GPD("partitioner: overlap", int, 1));
    innerSolverList.set("relaxation: type", MUELU_GPD("relaxation: type", std::string, "Gauss-Seidel"));
    innerSolverList.set("relaxation: sweeps", MUELU_GPD("relaxation: sweeps", int, 1));
    innerSolverList.set("relaxation: damping factor", MUELU_GPD("relaxation: damping factor", double, 0.5));
    innerSolverList.set("relaxation: zero starting solution", false);
    // innerSolverList.set("relaxation: backward mode",          MUELU_GPD("relaxation: backward mode",  bool,           true);  NOT SUPPORTED YET

    std::string ifpackType = "SCHWARZ";

    smootherPrototype = rcp(new TrilinosSmoother(ifpackType, schwarzList));

  } else if (type == "schwarz") {
    std::string ifpackType = "SCHWARZ";

    smootherPrototype = rcp(new TrilinosSmoother(ifpackType, paramList));

  } else if (type == "braess-sarazin") {
    // Define smoother/solver for BraessSarazin
    SC omega     = MUELU_GPD("bs: omega", double, 1.0);
    bool lumping = MUELU_GPD("bs: lumping", bool, false);

    RCP<SchurComplementFactory> schurFact = rcp(new SchurComplementFactory());
    schurFact->SetParameter("omega", ParameterEntry(omega));
    schurFact->SetParameter("lumping", ParameterEntry(lumping));
    schurFact->SetFactory("A", MueLu::NoFactory::getRCP());

    // Schur complement solver
    RCP<SmootherPrototype> schurSmootherPrototype;
    std::string schurSmootherType = (paramList.isParameter("schur smoother: type") ? paramList.get<std::string>("schur smoother: type") : "RELAXATION");
    if (schurSmootherType == "RELAXATION") {
      ParameterList schurSmootherParams = paramList.sublist("schur smoother: params");
      // schurSmootherParams.set("relaxation: damping factor", omega);
      schurSmootherPrototype = rcp(new TrilinosSmoother(schurSmootherType, schurSmootherParams));
    } else {
      schurSmootherPrototype = rcp(new DirectSolver());
    }
    schurSmootherPrototype->SetFactory("A", schurFact);

    RCP<SmootherFactory> schurSmootherFact = rcp(new SmootherFactory(schurSmootherPrototype));

    // Define temporary FactoryManager that is used as input for BraessSarazin smoother
    RCP<FactoryManager> braessManager = rcp(new FactoryManager());
    braessManager->SetFactory("A", schurFact);                 // SchurComplement operator for correction step (defined as "A")
    braessManager->SetFactory("Smoother", schurSmootherFact);  // solver/smoother for correction step
    braessManager->SetFactory("PreSmoother", schurSmootherFact);
    braessManager->SetFactory("PostSmoother", schurSmootherFact);
    braessManager->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

    smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetParameter("Sweeps", ParameterEntry(MUELU_GPD("bs: sweeps", int, 1)));
    smootherPrototype->SetParameter("lumping", ParameterEntry(lumping));
    smootherPrototype->SetParameter("Damping factor", ParameterEntry(omega));
    smootherPrototype->SetParameter("q2q1 mode", ParameterEntry(true));
    rcp_dynamic_cast<BraessSarazinSmoother>(smootherPrototype)->AddFactoryManager(braessManager, 0);  // set temporary factory manager in BraessSarazin smoother

  } else if (type == "ilu") {
    std::string ifpackType = "RILUK";

    smootherPrototype = rcp(new TrilinosSmoother(ifpackType, paramList));

  } else if (type == "direct") {
    smootherPrototype = rcp(new BlockedDirectSolver());

  } else {
    throw MueLu::Exceptions::RuntimeError("Unknown smoother type: \"" + type + "\"");
  }

  return coarseSolver ? rcp(new SmootherFactory(smootherPrototype, Teuchos::null)) : rcp(new SmootherFactory(smootherPrototype));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Absolute(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) const {
  typedef Xpetra::CrsMatrix<SC, LO, GO, NO> CrsMatrix;
  typedef Xpetra::CrsMatrixWrap<SC, LO, GO, NO> CrsMatrixWrap;
  typedef Xpetra::Matrix<SC, LO, GO, NO> Matrix;

  const CrsMatrixWrap& Awrap = dynamic_cast<const CrsMatrixWrap&>(A);

  ArrayRCP<const size_t> iaA;
  ArrayRCP<const LO> jaA;
  ArrayRCP<const SC> valA;
  Awrap.getCrsMatrix()->getAllValues(iaA, jaA, valA);

  ArrayRCP<size_t> iaB(iaA.size());
  ArrayRCP<LO> jaB(jaA.size());
  ArrayRCP<SC> valB(valA.size());
  for (int i = 0; i < iaA.size(); i++) iaB[i] = iaA[i];
  for (int i = 0; i < jaA.size(); i++) jaB[i] = jaA[i];
  for (int i = 0; i < valA.size(); i++) valB[i] = Teuchos::ScalarTraits<SC>::magnitude(valA[i]);

  RCP<Matrix> B       = rcp(new CrsMatrixWrap(A.getRowMap(), A.getColMap(), 0));
  RCP<CrsMatrix> Bcrs = rcp_dynamic_cast<CrsMatrixWrap>(B)->getCrsMatrix();
  Bcrs->setAllValues(iaB, jaB, valB);
  Bcrs->expertStaticFillComplete(A.getDomainMap(), A.getRangeMap());

  return B;
}

// Public functions overridden from Teuchos::Describable
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  return "Thyra::MueLuTpetraQ2Q1PreconditionerFactory";
}

}  // namespace Thyra

#endif
#endif  // ifdef THYRA_MUELU_TPETRA_Q2Q1PRECONDITIONER_FACTORY_DEF_HPP
