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

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include <MueLu_FilteredAFactory.hpp>
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include "MueLu_CoalesceDropFactory_decl.hpp"
#include <MueLu_BlockedPFactory.hpp>
#include "MueLu_GenericRFactory.hpp"
#include <MueLu_BlockedRAPFactory.hpp>
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_BlockedDirectSolver.hpp"
#include <MueLu_SubBlockAFactory.hpp>
#include "../../research/q2q1/MueLu_Q2Q1PFactory.hpp"
#include "../../research/q2q1/MueLu_Q2Q1uPFactory.hpp"
#include <MueLu_PatternFactory.hpp>
#include <MueLu_EminPFactory.hpp>
#include <MueLu_ConstraintFactory.hpp>
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include <MueLu_SmootherFactory.hpp>

#include "Xpetra_Matrix.hpp"
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teko_Utilities.hpp"

#include <string>

namespace Thyra {


using Teuchos::RCP;
using Teuchos::ParameterList;


// Constructors/initializers/accessors


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MueLuTpetraQ2Q1PreconditionerFactory()
{}


// Overridden from PreconditionerFactoryBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

  typedef Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyraTpetraLinOp;
  const RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);

  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinOp;
  const RCP<const TpetraLinOp> tpetraFwdOp = Teuchos::nonnull(thyraTpetraFwdOp) ? thyraTpetraFwdOp->getConstTpetraOperator() : Teuchos::null;

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMat;
  const RCP<const TpetraCrsMat> tpetraFwdCrsMat = Teuchos::rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);

  return Teuchos::nonnull(tpetraFwdCrsMat);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<PreconditionerBase<Scalar> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializePrec(
  const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  PreconditionerBase<Scalar> *prec,
  const ESupportSolveUse supportSolveUse
  ) const
{
  // Check precondition

  TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
  TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
  TEUCHOS_ASSERT(prec);

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nEntering Thyra::MueLuTpetraQ2Q1PreconditionerFactory::initializePrec(...) ...\n";
  }

  // Retrieve wrapped concrete Tpetra matrix from FwdOp

  const Teuchos::RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc->getOp();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

  typedef Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyraTpetraLinOp;
  const Teuchos::RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraTpetraFwdOp));

  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinOp;
  const Teuchos::RCP<const TpetraLinOp> tpetraFwdOp = thyraTpetraFwdOp->getConstTpetraOperator();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdOp));

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMat;
  const Teuchos::RCP<const TpetraCrsMat> tpetraFwdCrsMat = Teuchos::rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdCrsMat));

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nCreating a new MueLu::TpetraOperator object...\n";
  }
  timer.start(true);

  // Workaround since MueLu interface does not accept const matrix as input
  const Teuchos::RCP<TpetraCrsMat> tpetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<TpetraCrsMat>(tpetraFwdCrsMat);

  // Create and compute the initial preconditioner

  // Create a copy, as we may remove some things from the list
  Teuchos::ParameterList paramList = *paramList_;

  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
  Teuchos::RCP<MV> coords;
  if (paramList.isType<Teuchos::RCP<MV> >("Coordinates")) {
    coords = paramList.get<Teuchos::RCP<MV> >("Coordinates");
    paramList.remove("Coordinates");
  }

  Teuchos::RCP<MV> null_space;
  if (paramList.isType<Teuchos::RCP<MV> >("Nullspace")) {
    null_space = paramList.get<Teuchos::RCP<MV> >("Nullspace");
    paramList.remove("Nullspace");
  }

  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > coordsVel;
  if (paramList.isType<  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >  >("Velcoords")) {
    coordsVel = paramList.get<Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > ("Velcoords");
    paramList.remove("Velcoords");
  }
  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > coordsPres;
  if (paramList.isType<  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >  >("Prescoords")) {
    coordsPres = paramList.get<Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > ("Prescoords");
    paramList.remove("Prescoords");
  }

  Teuchos::ArrayRCP<LocalOrdinal> p2vMap;
  if (paramList.isType< Teuchos::ArrayRCP<LocalOrdinal> >("p2vMap")) {
    p2vMap = paramList.get< Teuchos::ArrayRCP<LocalOrdinal> > ("p2vMap");
    paramList.remove("p2vMap");
  }
  Teko::LinearOp thA11;
  if (paramList.isType< Teko::LinearOp >("A11")) {
    thA11 = paramList.get<  Teko::LinearOp >("A11");
    paramList.remove("A11");
  }
  Teko::LinearOp thA12;
  if (paramList.isType< Teko::LinearOp >("A12")) {
    thA12 = paramList.get<  Teko::LinearOp >("A12");
    paramList.remove("A12");
  }
  Teko::LinearOp thA21;
  if (paramList.isType< Teko::LinearOp >("A21")) {
    thA21 = paramList.get<  Teko::LinearOp >("A21");
    paramList.remove("A21");
  }

  Teuchos::RCP<MueLu::Hierarchy <Scalar,LocalOrdinal,GlobalOrdinal,Node> > HH = Teuchos::rcp(new MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>());

  const RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  typedef MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> MueLuOperator;
  const Teuchos::RCP<MueLuOperator> mueluPrecOp = Q2Q1MkPrecond(HH, comm, 2, paramList , coordsVel, coordsPres, p2vMap, thA11, thA12, thA21);

  timer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    Teuchos::OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
  }

  const Teuchos::RCP<LinearOpBase<Scalar> > thyraPrecOp = Thyra::createLinearOp(Teuchos::RCP<TpetraLinOp>(mueluPrecOp));
  defaultPrec->initializeUnspecified(thyraPrecOp);

  totalTimer.stop();
  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
    *out << "\nTotal time in Thyra::MueLuTpetraQ2Q1PreconditionerFactory::initializePrec(...) = " << totalTimer.totalElapsedTime() << " sec\n";
  }

  if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
    *out << "\nLeaving Thyra::MueLuTpetraQ2Q1PreconditionerFactory::initializePrec(...) ...\n";
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::uninitializePrec(
  PreconditionerBase<Scalar> *prec,
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > *fwdOp,
  ESupportSolveUse *supportSolveUse
  ) const
{
  // Check precondition

  TEUCHOS_ASSERT(prec);

  // Retrieve concrete preconditioner object

  const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec =
    Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
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
void MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameterList(
  Teuchos::RCP<ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
  paramList_ = paramList;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstParameterList()
{
  return paramList_;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unsetParameterList()
{
  Teuchos::RCP<ParameterList> savedParamList = paramList_;
  paramList_ = Teuchos::null;
  return savedParamList;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getParameterList() const
{
  return paramList_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList>
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;

  if (Teuchos::is_null(validPL)) {
    validPL = Teuchos::rcp(new ParameterList());
  }

  return validPL;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Q2Q1MkPrecond(
           Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & HH,
     const Teuchos::RCP<const Teuchos::Comm<int> > &comm, int maxLevels, const ParameterList paramList,
     const Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & coordsVel,
     const Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & coordsPres,
     const Teuchos::ArrayRCP<LocalOrdinal> & p2vMap, 
     const Teko::LinearOp & thA11, const Teko::LinearOp & thA12, const Teko::LinearOp & thA21) const
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;
  using Teuchos::as;


  typedef Scalar          SC;
  typedef LocalOrdinal    LO;
  typedef GlobalOrdinal   GO;
  typedef Node            NO;

  typedef Tpetra::CrsMatrix<Scalar,int,int> TP_Crs;
  typedef Tpetra::Operator<Scalar,int,int> TP_Op;
  typedef Xpetra::Matrix<SC, LO, GO, NO>   Matrix;
  typedef Xpetra::CrsMatrix<SC,LO,GO,NO>   CrsMatrix;
  typedef MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>        MUtils;

   // Pull out Tpetra matrices
   
   RCP<Thyra::LinearOpBase<double> > ThNonConstA11 = Teuchos::rcp_const_cast< Thyra::LinearOpBase<double> >(thA11);
   RCP<Thyra::LinearOpBase<double> > ThNonConstA21 = Teuchos::rcp_const_cast< Thyra::LinearOpBase<double> >(thA21);
   RCP<Thyra::LinearOpBase<double> > ThNonConstA12 = Teuchos::rcp_const_cast< Thyra::LinearOpBase<double> >(thA12);

   RCP<TP_Op>  TpetA11   = Thyra::TpetraOperatorVectorExtraction<Scalar,int>::getTpetraOperator(ThNonConstA11);
   RCP<TP_Op>  TpetA21   = Thyra::TpetraOperatorVectorExtraction<Scalar,int>::getTpetraOperator(ThNonConstA21);
   RCP<TP_Op>  TpetA12   = Thyra::TpetraOperatorVectorExtraction<Scalar,int>::getTpetraOperator(ThNonConstA12);

   RCP<TP_Crs> TpetCrsA11= rcp_dynamic_cast<TP_Crs>(TpetA11);
   RCP<TP_Crs> TpetCrsA21= rcp_dynamic_cast<TP_Crs>(TpetA21);
   RCP<TP_Crs> TpetCrsA12= rcp_dynamic_cast<TP_Crs>(TpetA12);

   RCP<Matrix> A_11      = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA11);
   RCP<Matrix> XpetA21   = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA21); //Needs map modification
   RCP<Matrix> XpetA12   = MueLu::TpetraCrs_To_XpetraMatrix(TpetCrsA12); //Needs map modification

   //  Create new A21 with map so that the global indices of the RowMap starts
   //  from nv+1 (where nv is the number of rows in the A11 block)

   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > RangeMapBlk2= XpetA21->getRangeMap();
   Xpetra::global_size_t                                         NRowsBlk2= RangeMapBlk2->getNodeNumElements();
   Teuchos::ArrayView<const GO>                           RangeEntriesBlk2= RangeMapBlk2->getNodeElementList();
   Teuchos::ArrayView<const GO>                          MyGlobalCols = XpetA21->getColMap()->getNodeElementList();

   int  nv = A_11->getRangeMap()->getNodeNumElements();
   Teuchos::Array<GO> newGlobalRows(NRowsBlk2, 0.0);
   for (Xpetra::global_size_t i = 0; i < NRowsBlk2; i++) newGlobalRows[i] = RangeEntriesBlk2[i] + nv;

   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal> > newRangeMapBlk2 = 
          rcp( new Xpetra::TpetraMap<LO,GO>(NRowsBlk2,newGlobalRows,RangeMapBlk2->getIndexBase(), comm));

   RCP<Matrix> A_21  = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,Node>(newRangeMapBlk2,
                                          XpetA21->getDomainMap(),XpetA21->getNodeMaxNumRowEntries()));
   RCP<CrsMatrix> A21Crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(A_21)->getCrsMatrix();

   for (LO row = 0; row < Teuchos::as<LO>(XpetA21->getRowMap()->getNodeNumElements()); ++row) {
      size_t nnnz = XpetA21->getNumEntriesInLocalRow(row);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      Teuchos::Array<GO>     newInds(nnnz, 0);
      XpetA21->getLocalRowView(row, indices, vals);
      for (LO colID = 0; colID < Teuchos::as<LO>(nnnz); colID++)
         newInds[colID] = MyGlobalCols[indices[colID]];
      A21Crs->insertGlobalValues(newGlobalRows[row],newInds,vals);
    }
    A21Crs->fillComplete(XpetA21->getDomainMap(),newRangeMapBlk2);
 
#ifdef out
   // This code would be cheaper than the above code, but it did not do the 
   // proper thing with empty columns. I'm leaving it here as it might be
   // useful if fixed up.
   
   Teuchos::ArrayRCP<const size_t>  OrigRowPtr;   Teuchos::ArrayRCP<size_t>  NewRowPtr;
   Teuchos::ArrayRCP<const LO>     OrigColumns;   Teuchos::ArrayRCP<LO>      NewColumns;
   Teuchos::ArrayRCP<const SC>      OrigValues;   Teuchos::ArrayRCP<SC>      NewValues;

   RCP<CrsMatrix> OrigCrs = rcp_dynamic_cast<CrsMatrixWrap>(XpetA21)->getCrsMatrix();
   OrigCrs->getAllValues(OrigRowPtr, OrigColumns, OrigValues);
   RCP<Matrix>    Bad21   = rcp(new CrsMatrixWrap(newRangeMapBlk2, XpetA21->getDomainMap(), 0, Xpetra::StaticProfile));
   RCP<CrsMatrix> Bad21Crs= rcp_dynamic_cast<CrsMatrixWrap>(Bad21)->getCrsMatrix();
   size_t nnz = XpetA21->getNodeNumEntries();

   Bad21Crs->allocateAllValues(nnz,  NewRowPtr, NewColumns, NewValues);
   for (size_t ii = 0; ii <= NRowsBlk2; ii++) NewRowPtr[ii] = OrigRowPtr[ii];
   for (size_t ii = 0; ii <  nnz  ; ii++) NewColumns[ii]= OrigColumns[ii];
   for (size_t ii = 0; ii <  nnz  ; ii++) NewValues[ii] = OrigValues[ii];
   Bad21Crs->setAllValues(NewRowPtr, NewColumns, NewValues);
   Bad21Crs->expertStaticFillComplete(XpetA21->getColMap(),newRangeMapBlk2);
#endif

   //  Create new A12 with map so that the global indices of the ColMap starts
   //  from nv+1 (where nv is the number of rows in the A11 block)

   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >DomainMapBlk2= XpetA12->getDomainMap();
   Xpetra::global_size_t                                         NColsBlk2= DomainMapBlk2->getNodeNumElements();
                                                                      // maybe should be column map???
   Teuchos::ArrayView<const GO>                          DomainEntriesBlk2= DomainMapBlk2->getNodeElementList();
   Teuchos::ArrayView<const GO>                          MyGlobalRows = XpetA12->getRowMap()->getNodeElementList();

   Teuchos::Array<GO> newGlobalCols(NColsBlk2, 0.0);
   for (Xpetra::global_size_t i = 0; i < NColsBlk2; i++) newGlobalCols[i] = DomainEntriesBlk2[i] + nv;

   RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal> > newDomainMapBlk2 = 
          rcp( new Xpetra::TpetraMap<LO,GO>(NColsBlk2,newGlobalCols,DomainMapBlk2->getIndexBase(), comm));

   RCP<Matrix> A_12  = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,Node>(XpetA12->getRangeMap(),newDomainMapBlk2,XpetA12->getNodeMaxNumRowEntries()));
   RCP<CrsMatrix> A12Crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(A_12)->getCrsMatrix();

   for (LO row = 0; row < Teuchos::as<LO>(XpetA12->getRowMap()->getNodeNumElements()); ++row) {
      size_t nnnz = XpetA12->getNumEntriesInLocalRow(row);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      Teuchos::Array<GO>     newInds(nnnz, 0);
      XpetA12->getLocalRowView(row, indices, vals);
      for (LO colID = 0; colID < Teuchos::as<LO>(nnnz); colID++)
         newInds[colID] = newGlobalCols[indices[colID]];
      A12Crs->insertGlobalValues(MyGlobalRows[row],newInds,vals);
    }
    A12Crs->fillComplete(newDomainMapBlk2,XpetA12->getRangeMap());

    RCP<Matrix>    A_22     = Teuchos::null;

    RCP<CrsMatrix> A_11_crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(A_11)->getCrsMatrix();
    RCP<CrsMatrix> A_12_crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(A_12)->getCrsMatrix();
    RCP<CrsMatrix> A_21_crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(A_21)->getCrsMatrix();
    RCP<CrsMatrix> A_22_crs = Teuchos::null;

    // =========================================================================
    // Preconditioner construction - I (block)
    // =========================================================================
    
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;
    out.setOutputToRootOnly(0);
    RCP<Matrix> BBt = MUtils::Multiply(*A_21, false, *A_12, false, out);

    SC dropTol = 0.06;
    RCP<Matrix> filteredA = FilterMatrix(*A_11, dropTol);
    RCP<Matrix> filteredB = FilterMatrix(*BBt,  dropTol);

    RCP<CrsMatrix> fA_11_crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(filteredA)->getCrsMatrix();
    RCP<CrsMatrix> fA_12_crs = Teuchos::null;
    RCP<CrsMatrix> fA_21_crs = Teuchos::null;
    RCP<CrsMatrix> fA_22_crs = rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,Node> >(filteredB)->getCrsMatrix();
    std::vector<Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > > partMaps(2);

    // Build the large filtered matrix which requires strided maps
    
    Xpetra::global_size_t  NumVel = A12Crs->getRangeMap()->getNodeNumElements();
    Xpetra::global_size_t NumPres = A21Crs->getRangeMap()->getNodeNumElements();

    const GO indexBase = 0;
    std::vector<size_t> stridingInfo(1, 1);
    int stridedBlockId = -1;

    Teuchos::Array<GO> elementList(NumVel+NumPres); // Not RCP ...  does this get cleared ?
    Teuchos::Array<GO> VelList = A12Crs->getRangeMap()->getNodeElementList();
    Teuchos::Array<GO> PresList= A21Crs->getRangeMap()->getNodeElementList();

    for (Xpetra::global_size_t i =    0  ; i < NumVel;         i++) elementList[i] = VelList[i];
    for (Xpetra::global_size_t i = NumVel; i < NumVel+NumPres; i++) elementList[i] = PresList[i-NumVel];
    RCP<Xpetra::Map<LO,GO,Node> > fullMap = Xpetra::StridedMapFactory<LO,GO>::Build(Xpetra::UseTpetra, NumVel+NumPres, elementList(), indexBase, stridingInfo, comm);

    partMaps[0] = Xpetra::StridedMapFactory<LO,GO>::Build(Xpetra::UseTpetra,NumVel,VelList,indexBase,stridingInfo,comm);
    partMaps[1] = Xpetra::StridedMapFactory<LO,GO>::Build(Xpetra::UseTpetra, NumPres, PresList, indexBase,
                                           stridingInfo, comm, stridedBlockId, NumVel);
    Teuchos::RCP<const Xpetra::MapExtractor<SC,LO,GO,Node> > mapExtractor = Xpetra::MapExtractorFactory<SC,LO,GO,Node>::Build(fullMap, partMaps);
    RCP<Xpetra::BlockedCrsMatrix<SC,LO,GO,Node> > fA = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<SC,LO,GO,Node>(mapExtractor, mapExtractor, 10));
    fA->setMatrix(0, 0, fA_11_crs);
    fA->setMatrix(0, 1, fA_12_crs);
    fA->setMatrix(1, 0, fA_21_crs);
    fA->setMatrix(1, 1, fA_22_crs);
    fA->fillComplete();

    // -------------------------------------------------------------------------
    // Preconditioner construction - I.a (filtered hierarchy)
    // -------------------------------------------------------------------------
    MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>  M;
    SetDependencyTree(M);

    typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> Hierarchy;
    HH = Teuchos::rcp(new Hierarchy);
    RCP<MueLu::Level> finestLevel = HH->GetLevel(0);
    finestLevel->Set("A",                     rcp_dynamic_cast<Matrix>(fA));
    finestLevel->Set("p2vMap",                p2vMap);
    finestLevel->Set("CoordinatesVelocity",   Xpetra::toXpetra(coordsVel));
    finestLevel->Set("CoordinatesPressure",   Xpetra::toXpetra(coordsPres));
    HH->SetMaxCoarseSize(1);

    // The first invocation of Setup() builds the hierarchy using the filtered
    // matrix. This build includes the grid transfers but not the creation of the
    // smoothers.
    // NOTE: we need to indicate what should be kept from the first invocation
    // for the second invocation, which then focuses on building the smoothers
    // for the unfiltered matrix.
    HH->Keep("P",     M.GetFactory("P")    .get());
    HH->Keep("R",     M.GetFactory("R")    .get());
    HH->Keep("Ptent", M.GetFactory("Ptent").get());
    HH->Setup(M, 0, maxLevels);

    // -------------------------------------------------------------------------
    // Preconditioner construction - I.b (Vanka smoothers for unfiltered matrix)
    // -------------------------------------------------------------------------
    // Set up Vanka smoothing via a combination of Schwarz and block relaxation.
    Teuchos::ParameterList schwarzList;
    schwarzList.set("schwarz: overlap level",                 Teuchos::as<int>(0));
    schwarzList.set("schwarz: zero starting solution",        false);
    schwarzList.set("subdomain solver name",                  "Block_Relaxation");

    Teuchos::ParameterList& innerSolverList = schwarzList.sublist("subdomain solver parameters");
    innerSolverList.set("partitioner: type",                  "user");
    innerSolverList.set("partitioner: overlap",               as<int>(1));
    innerSolverList.set("relaxation: type",                   "Gauss-Seidel");
    innerSolverList.set("relaxation: sweeps",                 as<int>(1));
    innerSolverList.set("relaxation: damping factor",         0.5);
    innerSolverList.set("relaxation: zero starting solution", false);
    // innerSolverList.set("relaxation: backward mode",true);  NOT SUPPORTED YET

    std::string ifpackType = "SCHWARZ";
    typedef MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> SmootherPrototype;
    typedef MueLu::TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node> TrilinosSmoother;
    typedef MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> SmootherFactory;

    Teuchos::RCP<SmootherPrototype> smootherPrototype = Teuchos::rcp(new TrilinosSmoother(ifpackType, schwarzList));
    M.SetFactory("Smoother",     rcp(new SmootherFactory(smootherPrototype)));
    M.SetFactory("CoarseSolver", rcp(new SmootherFactory(smootherPrototype)));


#ifdef HAVE_MUELU_DEBUG
    M.ResetDebugData();
#endif

    RCP<Xpetra::BlockedCrsMatrix<SC,LO,GO,Node> > A = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<SC,LO,GO,Node>(mapExtractor, mapExtractor, 10));
    A->setMatrix(0, 0, A_11_crs);
    A->setMatrix(0, 1, A_12_crs);
    A->setMatrix(1, 0, A_21_crs);
    A->setMatrix(1, 1, A_22_crs);
    A->fillComplete();

    HH->GetLevel(0)->Set("A", rcp_dynamic_cast<Matrix>(A));

    HH->Setup(M, 0, HH->GetNumLevels());

    return rcp(new MueLu::TpetraOperator<SC,LO,GO,Node>(HH));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > 
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
FilterMatrix(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & A, Scalar dropTol) const {

    using Teuchos::RCP;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>   Matrix;

    MueLu::Level level;
    level.SetLevelID(1);
    level.Set<RCP<Matrix> >("A", rcpFromRef(A));

    MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>  M;
    level.SetFactoryManager(rcpFromRef(M));

    RCP<MueLu::CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> > dropFactory = Teuchos::rcp(new MueLu::CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> ());
    ParameterList dropParams = *(dropFactory->GetValidParameterList());
    dropParams.set("lightweight wrap",          true);
    dropParams.set("aggregation: drop scheme",  "classical");
    dropParams.set("aggregation: drop tol",     dropTol);
    // dropParams.set("Dirichlet detection threshold", <>);

    dropFactory->SetParameterList(dropParams);
    M.SetFactory("Graph",     dropFactory);
    M.SetFactory("Filtering", dropFactory);

    RCP<MueLu::FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> > filterFactory = rcp(new MueLu::FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> ());
    ParameterList filterParams = *(filterFactory->GetValidParameterList());
    filterParams.set("filtered matrix: reuse graph", false);
    filterFactory->SetParameterList(filterParams);
    filterFactory->SetFactory("Graph", dropFactory);

    // Build

    level.Request("A", filterFactory.get());
    filterFactory->Build(level);

    RCP<Matrix> filteredA;
    level.Get("A", filteredA, filterFactory.get());

    return filteredA;
  }
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
SetDependencyTree(MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>  & M) const {
    using Teuchos::RCP;
    using Teuchos::rcp;

typedef MueLu::BlockedPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedPFactory;
typedef MueLu::GenericRFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> GenericRFactory;
typedef MueLu::BlockedRAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedRAPFactory;
typedef MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherFactory;
typedef MueLu::BlockedDirectSolver<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedDirectSolver;


    RCP<MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> > M11 = rcp(new MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
    SetBlockDependencyTree(*M11, 0, 0, "velocity");

    RCP<MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> > M22 = rcp(new MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
    SetBlockDependencyTree(*M22, 1, 1, "pressure");

    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
    ParameterList pParamList = *(PFact->GetValidParameterList());
    pParamList.set("backwards", true);      // do pressure first
    PFact->SetParameterList(pParamList);
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);
    M.SetFactory("P", PFact);

    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
    RFact->SetFactory("P", PFact);
    M.SetFactory("R", RFact);

    RCP<MueLu::Factory > AcFact = rcp(new BlockedRAPFactory());
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);
    M.SetFactory("A", AcFact);

    M.SetFactory("Smoother",     Teuchos::null);
    M.SetFactory("CoarseSolver", Teuchos::null);

    RCP<MueLu::Factory> coarseFact =rcp(new SmootherFactory(rcp(new BlockedDirectSolver()), Teuchos::null));

//M.SetFactory("CoarseSolver", coarseFact);
    M.SetFactory("CoarseSolver", Teuchos::null);


}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
SetBlockDependencyTree(MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> & M, LocalOrdinal row, LocalOrdinal col, const std::string& mode)  const {

    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef MueLu::Q2Q1PFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Q2Q1PFactory;
    typedef MueLu::Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Q2Q1uPFactory;
    typedef MueLu::SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node > SubBlockAFactory;
    typedef MueLu::PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node > PatternFactory;
    typedef MueLu::ConstraintFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> ConstraintFactory;
    typedef MueLu::EminPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> EminPFactory;

    RCP<SubBlockAFactory> AFact = rcp(new SubBlockAFactory());
    AFact->SetFactory  ("A",         MueLu::NoFactory::getRCP());
    AFact->SetParameter("block row", Teuchos::ParameterEntry(row));
    AFact->SetParameter("block col", Teuchos::ParameterEntry(col));
    M.SetFactory("A", AFact);

    RCP<MueLu::Factory> Q2Q1Fact;

    bool xSSTRUCTURED = true;

    if (xSSTRUCTURED) {
      Q2Q1Fact = rcp(new Q2Q1PFactory);

    } else {
      Q2Q1Fact = rcp(new Q2Q1uPFactory);
      ParameterList q2q1ParamList = *(Q2Q1Fact->GetValidParameterList());
      q2q1ParamList.set("mode", mode);
      // q2q1ParamList.set("phase2", false);
      Q2Q1Fact->SetParameterList(q2q1ParamList);
    }
    Q2Q1Fact->SetFactory("A", AFact);
    M.SetFactory("Ptent", Q2Q1Fact);

    RCP<PatternFactory> patternFact = rcp(new PatternFactory);
    ParameterList patternParams = *(patternFact->GetValidParameterList());
    patternParams.set("emin: pattern order", 0);
    patternFact->SetParameterList(patternParams);
    patternFact->SetFactory("A", AFact);
    patternFact->SetFactory("P", Q2Q1Fact);
    M.SetFactory("Ppattern", patternFact);

    RCP<ConstraintFactory> CFact = rcp(new ConstraintFactory);
    CFact->SetFactory("Ppattern", patternFact);
    M.SetFactory("Constraint", CFact);

    RCP<EminPFactory> EminPFact = rcp(new EminPFactory());
    EminPFact->SetFactory("A",          AFact);
    EminPFact->SetFactory("Constraint", CFact);
    EminPFact->SetFactory("P",          Q2Q1Fact);
    M.SetFactory("P", EminPFact);

}





// Public functions overridden from Teuchos::Describable

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string MueLuTpetraQ2Q1PreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
{
  return "Thyra::MueLuTpetraQ2Q1PreconditionerFactory";
}

} // namespace Thyra

#endif
#endif // ifdef THYRA_MUELU_TPETRA_Q2Q1PRECONDITIONER_FACTORY_DEF_HPP
