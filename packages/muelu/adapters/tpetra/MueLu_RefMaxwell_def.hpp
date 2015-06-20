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
#ifndef MUELU_REFMAXWELL_DEF_HPP
#define MUELU_REFMAXWELL_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include "MueLu_RefMaxwell_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {

  return Xpetra::toTpetraNonZero(SM_Matrix_->getDomainMap());

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {

  return Xpetra::toTpetraNonZero(SM_Matrix_->getRangeMap());

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(Teuchos::ParameterList& list) {

  disable_addon_  =  list.get("refmaxwell: disable add-on",true);
  MaxCoarseSize_  =  list.get("refmaxwell: max coarse size",1000);
  MaxLevels_      =  list.get("refmaxwell: max levels",5);
  Cycles_         =  list.get("refmaxwell: cycles",1);
  precType11_     =  list.get("refmaxwell: edge smoother","CHEBYSHEV");
  precType22_     =  list.get("refmaxwell: node smoother","CHEBYSHEV");
  mode_           =  list.get("refmaxwell: mode","additive");

  if(list.isSublist("refmaxwell: edge smoother list"))
    precList11_     =  list.sublist("refmaxwell: edge smoother list");

  if(list.isSublist("refmaxwell: node smoother list"))
    precList22_     =  list.sublist("refmaxwell: node smoother list");

  hiptmairPreList_.set("hiptmair: smoother type 1",precType11_);
  hiptmairPreList_.set("hiptmair: smoother type 2",precType22_);
  hiptmairPreList_.set("hiptmair: smoother list 1",precList11_);
  hiptmairPreList_.set("hiptmair: smoother list 2",precList22_);
  hiptmairPreList_.set("hiptmair: pre or post","both");
  hiptmairPreList_.set("hiptmair: zero starting solution",true);

  hiptmairPostList_.set("hiptmair: smoother type 1",precType11_);
  hiptmairPostList_.set("hiptmair: smoother type 2",precType22_);
  hiptmairPostList_.set("hiptmair: smoother list 1",precList11_);
  hiptmairPostList_.set("hiptmair: smoother list 2",precList22_);
  hiptmairPostList_.set("hiptmair: pre or post","both");
  hiptmairPostList_.set("hiptmair: zero starting solution",false);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute() {

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  // clean rows associated with boundary conditions
  Utils::findDirichletRows(SM_Matrix_,BCrows_);
  Utils::findDirichletCols(D0_Matrix_,BCrows_,BCcols_);
  D0_Matrix_->resumeFill();
  Utils::Apply_BCsToMatrixRows(D0_Matrix_,BCrows_);
  Utils::Apply_BCsToMatrixCols(D0_Matrix_,BCcols_);
  D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  //D0_Matrix_->describe(out,Teuchos::VERB_EXTREME);

  // Form TMT_Matrix
  Teuchos::RCP<XMat> C1 = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
  TMT_Matrix_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
  Xpetra::MatrixMatrix::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,*C1,true,true);
  Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*C1,false,*TMT_Matrix_,true,true);
  TMT_Matrix_->resumeFill();
  Utils::Remove_Zeroed_Rows(TMT_Matrix_,1.0e-16);
  TMT_Matrix_->SetFixedBlockSize(1);
  //TMT_Matrix_->describe(out,Teuchos::VERB_EXTREME);

  // build nullspace if necessary
  if(Nullspace_ != Teuchos::null) {
    // no need to do anything - nullspace is built
  }
  else if(Nullspace_ == Teuchos::null && Coords_ != Teuchos::null) {
    Nullspace_ = MultiVectorFactory::Build(SM_Matrix_->getRowMap(),Coords_->getNumVectors());
    D0_Matrix_->apply(*Coords_,*Nullspace_);
  }
  else {
    std::cerr << "MueLu::RefMaxwell::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
  }

  // build special prolongator for (1,1)-block
  if(P11_==Teuchos::null) {
    buildProlongator();
  }

  // build coarse grid operator for (1,1)-block
  formCoarseMatrix();

  // build fine grid operator for (2,2)-block, D0* M1 D0
  Teuchos::RCP<XMat> C = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
  Xpetra::MatrixMatrix::Multiply(*M1_Matrix_,false,*D0_Matrix_,false,*C,true,true);
  A22_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
  Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*C,false,*A22_,true,true);
  A22_->resumeFill();
  Utils::Remove_Zeroed_Rows(A22_,1.0e-16);
  A22_->SetFixedBlockSize(1);

  // build stuff for hierarchies
  Teuchos::RCP<FactoryManager> Manager11 = Teuchos::rcp( new FactoryManager );
  Teuchos::RCP<FactoryManager> Manager22 = Teuchos::rcp( new FactoryManager );
  Teuchos::RCP<SmootherPrototype> SmooProto11
    = Teuchos::rcp( new Ifpack2Smoother(precType11_,precList11_) );
  Teuchos::RCP<SmootherFactory> SmooFact11
    = Teuchos::rcp( new SmootherFactory(SmooProto11) );
  Teuchos::RCP<SmootherPrototype> SmooProto22
    = Teuchos::rcp( new Ifpack2Smoother(precType22_,precList22_) );
  Teuchos::RCP<SmootherFactory> SmooFact22
    = Teuchos::rcp( new SmootherFactory(SmooProto22) );
  Teuchos::RCP<CoalesceDropFactory> Dropfact11
    = Teuchos::rcp( new CoalesceDropFactory() );
  Teuchos::RCP<CoalesceDropFactory> Dropfact22
    = Teuchos::rcp( new CoalesceDropFactory() );
  Teuchos::RCP<UncoupledAggregationFactory> Aggfact11
    = Teuchos::rcp( new UncoupledAggregationFactory() );
  Teuchos::RCP<UncoupledAggregationFactory> Aggfact22
    = Teuchos::rcp( new UncoupledAggregationFactory() );
  Teuchos::ParameterList params;
  params.set("aggregation: drop tol",1.0e-16);
  params.set("aggregation: Dirichlet threshold",1.0e-16);
  Dropfact11->SetParameterList(params);
  Dropfact22->SetParameterList(params);
  Manager11->SetFactory("Aggregates",Aggfact11);
  Manager11->SetFactory("Smoother",SmooFact11);
  Manager11->SetFactory("CoarseSolver",SmooFact11);
  Manager11->SetFactory("Graph",Dropfact11);
  Manager22->SetFactory("Aggregates",Aggfact22);
  Manager22->SetFactory("Smoother",SmooFact22);
  Manager22->SetFactory("CoarseSolver",SmooFact22);
  Manager22->SetFactory("Graph",Dropfact22);
  Hierarchy11_ = Teuchos::rcp( new Hierarchy(A11_) );
  Hierarchy11_ -> SetMaxCoarseSize( MaxCoarseSize_ );
  Hierarchy11_ -> Setup(*Manager11, 0, MaxLevels_ );
  Hierarchy22_ = Teuchos::rcp( new Hierarchy(A22_) );
  Hierarchy22_ -> SetMaxCoarseSize( MaxCoarseSize_ );
  Hierarchy22_ -> Setup(*Manager22, 0, MaxLevels_ );

  // build ifpack2 preconditioners for pre and post smoothing
  Teuchos::RCP<const TCRS> EdgeMatrix = Utils::Op2NonConstTpetraCrs(SM_Matrix_ );
  Teuchos::RCP<const TCRS> NodeMatrix = Utils::Op2NonConstTpetraCrs(TMT_Matrix_);
  Teuchos::RCP<const TCRS> PMatrix    = Utils::Op2NonConstTpetraCrs(D0_Matrix_);
  edgePreSmoother_  = Teuchos::rcp( new Ifpack2::Hiptmair<TCRS>(EdgeMatrix,NodeMatrix,PMatrix) );
  edgePostSmoother_ = Teuchos::rcp( new Ifpack2::Hiptmair<TCRS>(EdgeMatrix,NodeMatrix,PMatrix) );
  edgePreSmoother_ -> setParameters(hiptmairPreList_);
  edgePreSmoother_ -> initialize();
  edgePreSmoother_ -> compute();
  edgePostSmoother_ -> setParameters(hiptmairPostList_);
  edgePostSmoother_ -> initialize();
  edgePostSmoother_ -> compute();

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::buildProlongator() {

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  // build prolongator: algorithm 1 in the reference paper
  // First, aggregate nodal matrix by creating a 2-level hierarchy
  Teuchos::RCP<Hierarchy> auxHierarchy
    = Teuchos::rcp( new Hierarchy(TMT_Matrix_) );
  Teuchos::RCP<FactoryManager> auxManager
    = Teuchos::rcp( new FactoryManager );
  Teuchos::RCP<TentativePFactory> TentPfact
    = Teuchos::rcp( new TentativePFactory );
  Teuchos::RCP<SaPFactory> Pfact
    = Teuchos::rcp( new SaPFactory );
  Teuchos::RCP<UncoupledAggregationFactory> Aggfact
    = Teuchos::rcp( new UncoupledAggregationFactory() );
  Teuchos::ParameterList params;
  params.set("sa: damping factor",0.0);
  Pfact      -> SetParameterList(params);
  auxManager -> SetFactory("P", Pfact);
  auxManager -> SetFactory("Ptent", TentPfact);
  auxManager -> SetFactory("Aggregates", Aggfact);
  auxManager -> SetFactory("Smoother", Teuchos::null);
  auxManager -> SetFactory("CoarseSolver", Teuchos::null);
  auxHierarchy -> Keep("P", Pfact.get());
  auxHierarchy -> SetMaxCoarseSize( MaxCoarseSize_ );
  auxHierarchy -> Setup(*auxManager, 0, 2);

  // pull out tentative P
  Teuchos::RCP<Level> Level1 = auxHierarchy -> GetLevel(1);
  Teuchos::RCP<XMat> P = Level1 -> Get< Teuchos::RCP<XMat> >("P",Pfact.get());

  // make weighting matrix
  Teuchos::RCP<XMat> D0_Matrix_Abs=MatrixFactory2::BuildCopy(D0_Matrix_);
  D0_Matrix_Abs -> resumeFill();
  D0_Matrix_Abs -> setAllToScalar((Scalar)0.5);
  Utils::Apply_BCsToMatrixRows(D0_Matrix_Abs,BCrows_);
  Utils::Apply_BCsToMatrixCols(D0_Matrix_Abs,BCcols_);
  D0_Matrix_Abs -> fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  Teuchos::RCP<XMat> Ptent = MatrixFactory::Build(D0_Matrix_Abs->getRowMap(),0);
  Xpetra::MatrixMatrix::Multiply(*D0_Matrix_Abs,false,*P,false,*Ptent,true,true);

  // put in entries to P11
  size_t dim = Nullspace_->getNumVectors();
  size_t numLocalRows = SM_Matrix_->getNodeNumRows();
  Teuchos::RCP<XMap> BlockColMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getColMap(),dim);
  P11_ = Teuchos::rcp(new XCrsWrap(Ptent->getRowMap(),BlockColMap,0,Xpetra::StaticProfile));

  std::vector< Teuchos::ArrayRCP<const Scalar> > nullspace(dim);
  for(size_t i=0; i<dim; i++) {
    Teuchos::ArrayRCP<const Scalar> datavec = Nullspace_->getData(i);
    nullspace[i]=datavec;
  }

  size_t nnz=0;
  std::vector<size_t>       rowPtrs;
  std::vector<LocalOrdinal> blockCols;
  std::vector<Scalar>       blockVals;
  for(size_t i=0; i<numLocalRows; i++) {
    rowPtrs.push_back(nnz);
    Teuchos::ArrayView<const LocalOrdinal> localCols;
    Teuchos::ArrayView<const Scalar>       localVals;
    Ptent->getLocalRowView(i,localCols,localVals);
    size_t numCols = localCols.size();
    for(size_t j=0; j<numCols; j++) {
      for(size_t k=0; k<dim; k++) {
        blockCols.push_back(localCols[j]*dim+k);
        blockVals.push_back(localVals[j]*nullspace[k][i]);
        nnz++;
      }
    }
  }
  rowPtrs.push_back(nnz);

  ArrayRCP<size_t>       rcpRowPtr;
  ArrayRCP<LocalOrdinal> rcpColumns;
  ArrayRCP<Scalar>       rcpValues;

  RCP<XCRS> TP11 = rcp_dynamic_cast<XCrsWrap>(P11_)->getCrsMatrix();
  TP11->allocateAllValues(nnz, rcpRowPtr, rcpColumns, rcpValues);

  ArrayView<size_t>       rows    = rcpRowPtr();
  ArrayView<LocalOrdinal> columns = rcpColumns();
  ArrayView<Scalar>       values  = rcpValues();

  for (size_t ii = 0; ii < rowPtrs.size();   ii++) rows[ii]    = rowPtrs[ii];
  for (size_t ii = 0; ii < blockCols.size(); ii++) columns[ii] = blockCols[ii];
  for (size_t ii = 0; ii < blockVals.size(); ii++) values[ii]  = blockVals[ii];
  TP11->setAllValues(rcpRowPtr, rcpColumns, rcpValues);
  Teuchos::RCP<XMap> blockCoarseMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getDomainMap(),dim);
  TP11->expertStaticFillComplete(blockCoarseMap,SM_Matrix_->getDomainMap());

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::formCoarseMatrix() {

  // coarse matrix for P11* (M1 + D1* M2 D1) P11
  Teuchos::RCP<XMat> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
  Teuchos::RCP<XMat> Matrix1 = MatrixFactory::Build(P11_->getDomainMap(),0);

  // construct (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix::Multiply(*SM_Matrix_,false,*P11_,false,*C,true,true);

  // construct P11* (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix::Multiply(*P11_,true,*C,false,*Matrix1,true,true);

  if(disable_addon_==true) {
    // if add-on is not chosen
    A11_=Matrix1;
  }
  else {
    // catch a failure
    TEUCHOS_TEST_FOR_EXCEPTION(M0inv_Matrix_==Teuchos::null,std::invalid_argument,
                               "MueLu::RefMaxwell::formCoarseMatrix(): Inverse of "
                               "lumped mass matrix required for add-on (i.e. M0inv_Matrix is null)");

    // coarse matrix for add-on, i.e P11* (M1 D0 M0inv D0* M1) P11
    Teuchos::RCP<XMat> Zaux = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
    Teuchos::RCP<XMat> Z = MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
    Teuchos::RCP<XMat> C2 = MatrixFactory::Build(M0inv_Matrix_->getRowMap(),0);
    // construct M1 P11
    Xpetra::MatrixMatrix::Multiply(*M1_Matrix_,false,*P11_,false,*Zaux,true,true);
    // construct Z = D0* M1 P11
    Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*Zaux,false,*Z,true,true);
    // construct M0inv Z
    Xpetra::MatrixMatrix::Multiply(*M0inv_Matrix_,false,*Z,false,*C2,true,true);
    // construct Z* M0inv Z
    Teuchos::RCP<XMat> Matrix2 = MatrixFactory::Build(Z->getDomainMap(),0);
    Xpetra::MatrixMatrix::Multiply(*Z,true,*C2,false,*Matrix2,true,true);
    // add matrices together
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Utils2::TwoMatrixAdd(*Matrix1,false,(Scalar)1.0,*Matrix2,false,(Scalar)1.0,A11_,*out);
    A11_->fillComplete();
  }

  // set fixed block size for vector nodal matrix
  size_t dim = Nullspace_->getNumVectors();
  A11_->SetFixedBlockSize(dim);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(Teuchos::RCP<TCRS> SM_Matrix_new) {

  // convert Tpetra matrices to Xpetra
  Teuchos::RCP<XCRS> SM_tmp = Teuchos::rcp( new XTCRS(SM_Matrix_new) );
  SM_Matrix_ = Teuchos::rcp( new XCrsWrap(SM_tmp) );

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseAdditive(const XTMV& RHS, XTMV& X) const {

  // compute residuals
  RCP<XMV> residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  RCP<XMV> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);

  // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)
  Hierarchy11_->Iterate(*P11res, *P11x, Cycles_, true);
  Hierarchy22_->Iterate(*D0res,  *D0x,  Cycles_, true);

  // update current solution
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)1.0);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse121(const XTMV& RHS, XTMV& X) const {

  RCP<XMV> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());

  // precondition (1,1)-block
  RCP<XMV> residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, Cycles_, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (2,2)-block
  residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  Cycles_, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (1,1)-block
  residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, Cycles_, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse212(const XTMV& RHS, XTMV& X) const {

  RCP<XMV> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<XMV> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());

  // precondition (2,2)-block
  RCP<XMV> residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  Cycles_, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (1,1)-block
  residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, Cycles_, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (2,2)-block
  residual  = Utils::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  Cycles_, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                                                           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                                                           Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  try {

    TMV& temp_x = const_cast<TMV &>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV       tY(rcpFromRef(Y));
    tY.putScalar(Teuchos::ScalarTraits<Scalar>::zero());

    // apply pre-smoothing
    edgePreSmoother_->apply(X,Y);

    // do solve for the 2x2 block system
    if(mode_=="additive")
      applyInverseAdditive(tX,tY);
    else if(mode_=="121")
      applyInverse121(tX,tY);
    else if(mode_=="212")
      applyInverse212(tX,tY);
    else
      applyInverseAdditive(tX,tY);

    // apply post-smoothing
    edgePostSmoother_->apply(X,Y);

  } catch (std::exception& e) {

    //FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::RefMaxwell::ApplyInverse():" << std::endl
              << e.what() << std::endl;

  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasTransposeApply() const {
  return false;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initialize(const Teuchos::RCP<TCRS> & D0_Matrix,
           const Teuchos::RCP<TCRS> & M0inv_Matrix,
           const Teuchos::RCP<TCRS> & M1_Matrix,
           const Teuchos::RCP<TMV>  & Nullspace,
           const Teuchos::RCP<TMV>  & Coords,
           Teuchos::ParameterList& List)
{
  // some pre-conditions
  TEUCHOS_ASSERT(D0_Matrix!=Teuchos::null);
  TEUCHOS_ASSERT(M1_Matrix!=Teuchos::null);

  Hierarchy11_ = Teuchos::null;
  Hierarchy22_ = Teuchos::null;
  parameterList_ = List;
  disable_addon_ = false;
  MaxCoarseSize_ = 1000;
  MaxLevels_ = 5;
  Cycles_ = 1;
  precType11_ = "CHEBYSHEV";
  precType22_ = "CHEBYSHEV";
  mode_ = "additive";

  // set parameters
  setParameters(List);

  // convert Tpetra matrices to Xpetra
  Teuchos::RCP<XCRS> D0_tmp = Teuchos::rcp( new XTCRS(D0_Matrix) );
  D0_Matrix_ = Teuchos::rcp( new XCrsWrap(D0_tmp) );

  if(M0inv_Matrix != Teuchos::null) {
    Teuchos::RCP<XCRS> M0inv_tmp = Teuchos::rcp( new XTCRS(M0inv_Matrix) );
    M0inv_Matrix_ = Teuchos::rcp( new XCrsWrap(M0inv_tmp) );
  }

  Teuchos::RCP<XCRS> M1_tmp = Teuchos::rcp( new XTCRS(M1_Matrix) );
  M1_Matrix_ = Teuchos::rcp( new XCrsWrap(M1_tmp) );

  // convert Tpetra MultiVector to Xpetra
  if(Coords != Teuchos::null)
    Coords_ = Xpetra::toXpetra(Coords);
  if(Nullspace != Teuchos::null)
    Nullspace_ = Xpetra::toXpetra(Nullspace);
}

} // namespace
#endif //ifdef HAVE_MUELU_TPETRA

#define MUELU_REFMAXWELL_SHORT
#endif //ifdef MUELU_REFMAXWELL_DEF_HPP
