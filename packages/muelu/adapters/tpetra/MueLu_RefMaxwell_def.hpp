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

// TODO move this file to xpetra subfolder

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixMatrix.hpp"

#include "MueLu_RefMaxwell_decl.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_MLParameterListInterpreter.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_HierarchyManager.hpp"

namespace MueLu {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const {
  return SM_Matrix_->getDomainMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const {
  return SM_Matrix_->getRangeMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(Teuchos::ParameterList& list) {

  disable_addon_  =  list.get("refmaxwell: disable add-on",true);
  mode_           =  list.get("refmaxwell: mode","additive");

  if(list.isSublist("refmaxwell: 11list"))
    precList11_     =  list.sublist("refmaxwell: 11list");

  if(list.isSublist("refmaxwell: 22list"))
    precList22_     =  list.sublist("refmaxwell: 22list");

  std::string ref("smoother:");
  std::string replace("coarse:");
  for(Teuchos::ParameterList::ConstIterator i=list.begin(); i !=list.end(); i++) {
    const std::string & pname = list.name(i);
    if(pname.find(ref)!=std::string::npos) {
      smootherList_.setEntry(pname,list.entry(i));
      std::string coarsename(pname);
      coarsename.replace((size_t)0,(size_t)ref.length(),replace);
    }
  }
  if(list.isSublist("smoother: params")) {
    smootherList_.set("coarse: params",list.sublist("smoother: params"));
  }
  smootherList_.set("max levels",1);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute() {

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  // clean rows associated with boundary conditions
  findDirichletRows(SM_Matrix_,BCrows_);
  findDirichletCols(D0_Matrix_,BCrows_,BCcols_);
  D0_Matrix_->resumeFill();
  Apply_BCsToMatrixRows(D0_Matrix_,BCrows_);
  Apply_BCsToMatrixCols(D0_Matrix_,BCcols_);
  D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  //D0_Matrix_->describe(out,Teuchos::VERB_EXTREME);

  // Form TMT_Matrix
  Teuchos::RCP<Matrix> C1 = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
  TMT_Matrix_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,*C1,true,true);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*C1,false,*TMT_Matrix_,true,true);
  TMT_Matrix_->resumeFill();
  Remove_Zeroed_Rows(TMT_Matrix_,1.0e-16);
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
  Teuchos::RCP<Matrix> C = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*D0_Matrix_,false,*C,true,true);
  A22_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*C,false,*A22_,true,true);
  A22_->resumeFill();
  Remove_Zeroed_Rows(A22_,1.0e-16);
  A22_->SetFixedBlockSize(1);

  // Use HierarchyManagers to build 11 & 22 Hierarchies
  typedef MueLu::HierarchyManager<SC,LO,GO,NO>               HierarchyManager;
  RCP<HierarchyManager> Manager11, Manager22, ManagerSmoother;
  std::string syntaxStr = "parameterlist: syntax";
  if (parameterList_.isParameter(syntaxStr) && parameterList_.get<std::string>(syntaxStr) == "ml") {
    parameterList_.remove(syntaxStr);
    Manager11 = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(precList11_));
    Manager22 = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(precList22_));
    ManagerSmoother = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(smootherList_));
  } else {
    Manager11 = rcp(new ParameterListInterpreter  <SC,LO,GO,NO>(precList11_,A11_->getDomainMap()->getComm()));
    Manager22 = rcp(new ParameterListInterpreter  <SC,LO,GO,NO>(precList22_,A22_->getDomainMap()->getComm()));
    ManagerSmoother = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(smootherList_,SM_Matrix_->getDomainMap()->getComm()));
  }

  Hierarchy11_=Manager11->CreateHierarchy();
  Hierarchy11_->setlib(Xpetra::UseTpetra);
  Hierarchy11_->GetLevel(0)->Set("A", A11_);
  Manager11->SetupHierarchy(*Hierarchy11_);

  Hierarchy22_=Manager22->CreateHierarchy();
  Hierarchy22_->setlib(Xpetra::UseTpetra);
  Hierarchy22_->GetLevel(0)->Set("A", A22_);
  Manager22->SetupHierarchy(*Hierarchy22_);

  // build ifpack2 preconditioners for pre and post smoothing
  HierarchySmoother_=ManagerSmoother->CreateHierarchy();
  HierarchySmoother_->setlib(Xpetra::UseTpetra);
  HierarchySmoother_->GetLevel(0)->Set("A", SM_Matrix_);
  ManagerSmoother->SetupHierarchy(*HierarchySmoother_);

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
  auxHierarchy -> SetMaxCoarseSize(1);
  auxHierarchy -> Setup(*auxManager, 0, 2);

  // pull out tentative P
  Teuchos::RCP<Level> Level1 = auxHierarchy -> GetLevel(1);
  Teuchos::RCP<Matrix> P = Level1 -> Get< Teuchos::RCP<Matrix> >("P",Pfact.get());

  // make weighting matrix
  Teuchos::RCP<Matrix> D0_Matrix_Abs=MatrixFactory2::BuildCopy(D0_Matrix_);
  D0_Matrix_Abs -> resumeFill();
  D0_Matrix_Abs -> setAllToScalar((Scalar)0.5);
  Apply_BCsToMatrixRows(D0_Matrix_Abs,BCrows_);
  Apply_BCsToMatrixCols(D0_Matrix_Abs,BCcols_);
  D0_Matrix_Abs -> fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  Teuchos::RCP<Matrix> Ptent = MatrixFactory::Build(D0_Matrix_Abs->getRowMap(),0);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_Abs,false,*P,false,*Ptent,true,true);

  // put in entries to P11
  size_t dim = Nullspace_->getNumVectors();
  size_t numLocalRows = SM_Matrix_->getNodeNumRows();
  Teuchos::RCP<Map> BlockColMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getColMap(),dim);
  P11_ = Teuchos::rcp(new CrsMatrixWrap(Ptent->getRowMap(),BlockColMap,0,Xpetra::StaticProfile));

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

  RCP<CrsMatrix> TP11 = rcp_dynamic_cast<CrsMatrixWrap>(P11_)->getCrsMatrix();
  TP11->allocateAllValues(nnz, rcpRowPtr, rcpColumns, rcpValues);

  ArrayView<size_t>       rows    = rcpRowPtr();
  ArrayView<LocalOrdinal> columns = rcpColumns();
  ArrayView<Scalar>       values  = rcpValues();

  for (size_t ii = 0; ii < rowPtrs.size();   ii++) rows[ii]    = rowPtrs[ii];
  for (size_t ii = 0; ii < blockCols.size(); ii++) columns[ii] = blockCols[ii];
  for (size_t ii = 0; ii < blockVals.size(); ii++) values[ii]  = blockVals[ii];
  TP11->setAllValues(rcpRowPtr, rcpColumns, rcpValues);
  Teuchos::RCP<Map> blockCoarseMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getDomainMap(),dim);
  TP11->expertStaticFillComplete(blockCoarseMap,SM_Matrix_->getDomainMap());

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::formCoarseMatrix() {

  // coarse matrix for P11* (M1 + D1* M2 D1) P11
  Teuchos::RCP<Matrix> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
  Teuchos::RCP<Matrix> Matrix1 = MatrixFactory::Build(P11_->getDomainMap(),0);

  // construct (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*P11_,false,*C,true,true);

  // construct P11* (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*P11_,true,*C,false,*Matrix1,true,true);

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
    Teuchos::RCP<Matrix> Zaux = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
    Teuchos::RCP<Matrix> Z = MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
    Teuchos::RCP<Matrix> C2 = MatrixFactory::Build(M0inv_Matrix_->getRowMap(),0);
    // construct M1 P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*P11_,false,*Zaux,true,true);
    // construct Z = D0* M1 P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*Zaux,false,*Z,true,true);
    // construct M0inv Z
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M0inv_Matrix_,false,*Z,false,*C2,true,true);
    // construct Z* M0inv Z
    Teuchos::RCP<Matrix> Matrix2 = MatrixFactory::Build(Z->getDomainMap(),0);
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Z,true,*C2,false,*Matrix2,true,true);
    // add matrices together
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*Matrix1,false,(Scalar)1.0,*Matrix2,false,(Scalar)1.0,A11_,*out);
    A11_->fillComplete();
  }

  // set fixed block size for vector nodal matrix
  size_t dim = Nullspace_->getNumVectors();
  A11_->SetFixedBlockSize(dim);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(Teuchos::RCP<Matrix> SM_Matrix_new) {
  SM_Matrix_ = SM_Matrix_new;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const {

  // compute residuals
  RCP<MultiVector> residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  RCP<MultiVector> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);

  // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)
  Hierarchy11_->Iterate(*P11res, *P11x, 1, true);
  Hierarchy22_->Iterate(*D0res,  *D0x,  1, true);

  // update current solution
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)1.0);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse121(const MultiVector& RHS, MultiVector& X) const {

  RCP<MultiVector> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());

  // precondition (1,1)-block
  RCP<MultiVector> residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, 1, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (2,2)-block
  residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  1, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (1,1)-block
  residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, 1, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse212(const MultiVector& RHS, MultiVector& X) const {

  RCP<MultiVector> P11res    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> P11x      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0res     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
  RCP<MultiVector> D0x       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());

  // precondition (2,2)-block
  RCP<MultiVector> residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  1, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (1,1)-block
  residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  P11_->apply(*residual,*P11res,Teuchos::TRANS);
  Hierarchy11_->Iterate(*P11res, *P11x, 1, true);
  P11_->apply(*P11x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

  // precondition (2,2)-block
  residual  = Utilities::Residual(*SM_Matrix_, X, RHS);
  D0_Matrix_->apply(*residual,*D0res,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res,  *D0x,  1, true);
  D0_Matrix_->apply(*D0x,*residual,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual, (Scalar) 1.0);

}


/*void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                                                           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                                                           Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {*/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply (const MultiVector& X, MultiVector& Y,
            Teuchos::ETransp mode,
            Scalar alpha,
            Scalar beta) const {
  try {

    /*TMV& temp_x = const_cast<TMV &>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV       tY(rcpFromRef(Y));
    tY.putScalar(Teuchos::ScalarTraits<Scalar>::zero());*/

    Y.putScalar(Teuchos::ScalarTraits<Scalar>::zero());

    // apply pre-smoothing
    HierarchySmoother_->Iterate(X,Y,1);

    // do solve for the 2x2 block system

    if(mode_=="additive")
      applyInverseAdditive(X,Y);
    else if(mode_=="121")
      applyInverse121(X,Y);
    else if(mode_=="212")
      applyInverse212(X,Y);
    else
      applyInverseAdditive(X,Y);

    // apply post-smoothing
    HierarchySmoother_->Iterate(X,Y,1);

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
initialize(const Teuchos::RCP<Matrix> & D0_Matrix,
           const Teuchos::RCP<Matrix> & M0inv_Matrix,
           const Teuchos::RCP<Matrix> & M1_Matrix,
           const Teuchos::RCP<MultiVector>  & Nullspace,
           const Teuchos::RCP<MultiVector>  & Coords,
           Teuchos::ParameterList& List)
{
  // some pre-conditions
  TEUCHOS_ASSERT(D0_Matrix!=Teuchos::null);
  TEUCHOS_ASSERT(M1_Matrix!=Teuchos::null);

  Hierarchy11_ = Teuchos::null;
  Hierarchy22_ = Teuchos::null;
  HierarchySmoother_ = Teuchos::null;
  parameterList_ = List;
  disable_addon_ = false;
  mode_ = "additive";

  // set parameters
  setParameters(List);

  D0_Matrix_ = D0_Matrix;
  M0inv_Matrix_ = M0inv_Matrix;
  M1_Matrix_ = M1_Matrix;
  Coords_ = Coords;
  Nullspace_ = Nullspace;
}

} // namespace

#define MUELU_REFMAXWELL_SHORT
#endif //ifdef MUELU_REFMAXWELL_DEF_HPP
