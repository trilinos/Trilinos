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

#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_TripleMatrixMultiply.hpp"

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

  disable_addon_    = list.get("refmaxwell: disable addon",true);
  mode_             = list.get("refmaxwell: mode","additive");
  dump_matrices_    = list.get("refmaxwell: dump matrices",false);
  read_P_from_file_ = list.get("refmaxwell: read_P_from_file_",false);
  P_filename_       = list.get("refmaxwell: P_filename_",std::string("P.mat"));

  if(list.isSublist("refmaxwell: 11list"))
    precList11_     =  list.sublist("refmaxwell: 11list");

  if(list.isSublist("refmaxwell: 22list"))
    precList22_     =  list.sublist("refmaxwell: 22list");

  if(list.isSublist("smoother: params")) {
    smootherList_ = list.sublist("smoother: params");
  }
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::compute() {

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(false);

  RCP<Teuchos::TimeMonitor> tmCompute = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: compute")));

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("SM.mat"), *SM_Matrix_);

  // clean rows associated with boundary conditions
  findDirichletRows(SM_Matrix_,BCrows_);
  findDirichletCols(D0_Matrix_,BCrows_,BCcols_);

  if (dump_matrices_) {
    std::ofstream outBCrows("BCrows.mat");
    std::copy(BCrows_.begin(), BCrows_.end(), std::ostream_iterator<int>(outBCrows, "\n"));
    std::ofstream outBCcols("BCcols.mat");
    std::copy(BCcols_.begin(), BCcols_.end(), std::ostream_iterator<int>(outBCcols, "\n"));
  }

  // build nullspace if necessary
  if(Nullspace_ != Teuchos::null) {
    // no need to do anything - nullspace is built
  }
  else if(Nullspace_ == Teuchos::null && Coords_ != Teuchos::null) {
    // typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    Teuchos::Array<Scalar> norms(Coords_->getNumVectors());
    Coords_->norm2(norms);
    for (size_t i=0;i<Coords_->getNumVectors();i++) {
      norms[i] = 1./norms[i];
    }
    // Coords_->scale(norms);
    for(size_t j=0;j<Coords_->getNumVectors();j++) {
      Teuchos::ArrayRCP<Scalar> datavec = Coords_->getDataNonConst(j);
      for(size_t i=0;i<static_cast<size_t>(datavec.size());i++)
        datavec[i] *= norms[j];
    }

    if (dump_matrices_)
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("coords.mat"), *Coords_);

    Nullspace_ = MultiVectorFactory::Build(SM_Matrix_->getRowMap(),Coords_->getNumVectors());
    D0_Matrix_->apply(*Coords_,*Nullspace_);
  }
  else {
    std::cerr << "MueLu::RefMaxwell::compute(): either the nullspace or the nodal coordinates must be provided." << std::endl;
  }

  /* Nuke the BC edges */
  size_t dim = Nullspace_->getNumVectors();
  for(size_t j=0;j<dim;j++) {
    Teuchos::ArrayRCP<Scalar> datavec = Nullspace_->getDataNonConst(j);
    for(size_t i=0;i<BCrows_.size();i++)
      datavec[BCrows_[i]]=0;
  }
  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("nullspace.mat"), *Nullspace_);

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("D0_clean.mat"), *D0_Matrix_);
  D0_Matrix_->resumeFill();
  Apply_BCsToMatrixRows(D0_Matrix_,BCrows_);
  Apply_BCsToMatrixCols(D0_Matrix_,BCcols_);
  D0_Matrix_->fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("D0_nuked.mat"), *D0_Matrix_);

  // build special prolongator for (1,1)-block
  if(P11_==Teuchos::null) {
    out << "\n"                                                         \
      "--------------------------------------------------------------------------------\n" \
      "---           build special prolongator for (1,1)-block                      ---\n" \
      "--------------------------------------------------------------------------------\n";

    // Form A_nodal = D0* M1 D0  (aka TMT_agg)
    RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build A_nodal")));
    A_nodal_Matrix_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
    if (parameterList_.get<bool>("rap: triple product", false) == false) {
      Teuchos::RCP<Matrix> C2 = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*D0_Matrix_,false,*C2,true,true);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*C2,false,*A_nodal_Matrix_,true,true);
    } else {
      Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
        MultiplyRAP(*D0_Matrix_, true, *M1_Matrix_, false, *D0_Matrix_, false, *A_nodal_Matrix_, true, true);
    }
    A_nodal_Matrix_->resumeFill();
    Remove_Zeroed_Rows(A_nodal_Matrix_,1.0e-16);
    A_nodal_Matrix_->SetFixedBlockSize(1);
    tm = Teuchos::null;

    if (dump_matrices_)
      Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("A_nodal.mat"), *A_nodal_Matrix_);

    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build special prolongator")));
    buildProlongator();
  }

  // Use HierarchyManagers to build 11 & 22 Hierarchies
  typedef MueLu::HierarchyManager<SC,LO,GO,NO> HierarchyManager;
  std::string syntaxStr = "parameterlist: syntax";

  {
    out << "\n" \
      "--------------------------------------------------------------------------------\n" \
      "---                  build MG for coarse (1,1)-block                         ---\n" \
      "--------------------------------------------------------------------------------\n";
    // build coarse grid operator for (1,1)-block
    RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build AH")));
    formCoarseMatrix();
    tm = Teuchos::null;

    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build coarse (1,1) hierarchy")));
    RCP<HierarchyManager> ManagerH;
    if (parameterList_.isParameter(syntaxStr) && parameterList_.get<std::string>(syntaxStr) == "ml") {
      parameterList_.remove(syntaxStr);
      ManagerH = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(precList11_, AH_->getDomainMap()->getComm()));
    } else {
      ManagerH = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(precList11_,AH_->getDomainMap()->getComm()));
    }
    HierarchyH_=ManagerH->CreateHierarchy();
    HierarchyH_->setlib(Xpetra::UseTpetra);
    HierarchyH_->GetLevel(0)->Set("A", AH_);
    ManagerH->SetupHierarchy(*HierarchyH_);
  }

  {
    out << "\n" \
      "--------------------------------------------------------------------------------\n" \
      "---                  build MG for (2,2)-block                                ---\n" \
      "--------------------------------------------------------------------------------\n";
    // build fine grid operator for (2,2)-block, D0* SM D0  (aka TMT)
    RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build A22")));
    A22_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);
    if (parameterList_.get<bool>("rap: triple product", false) == false) {
      Teuchos::RCP<Matrix> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,*C,true,true);
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*C,false,*A22_,true,true);
    } else {
      Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
        MultiplyRAP(*D0_Matrix_, true, *SM_Matrix_, false, *D0_Matrix_, false, *A22_, true, true);
    }
    A22_->resumeFill();
    Remove_Zeroed_Rows(A22_,1.0e-16);
    A22_->SetFixedBlockSize(1);
    tm = Teuchos::null;

    tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build (2,2) hierarchy")));
    RCP<HierarchyManager> Manager22;
    if (parameterList_.isParameter(syntaxStr) && parameterList_.get<std::string>(syntaxStr) == "ml") {
      parameterList_.remove(syntaxStr);
      Manager22 = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(precList22_, A22_->getDomainMap()->getComm()));
    } else {
      Manager22 = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(precList22_,A22_->getDomainMap()->getComm()));
    }
    Hierarchy22_=Manager22->CreateHierarchy();
    Hierarchy22_->setlib(Xpetra::UseTpetra);
    Hierarchy22_->GetLevel(0)->Set("A", A22_);
    Manager22->SetupHierarchy(*Hierarchy22_);
  }

  {
    out << "\n" \
      "--------------------------------------------------------------------------------\n" \
      "---                  build smoothers for fine (1,1)-block                     ---\n" \
      "--------------------------------------------------------------------------------\n";
    RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MueLu RefMaxwell: Build smoothers for fine (1,1) block")));
    std::string smootherType = parameterList_.get<std::string>("smoother: type", "CHEBYSHEV");
    if (smootherType == "hiptmair" &&
        SM_Matrix_->getDomainMap()->lib() == Xpetra::UseTpetra &&
        A22_->getDomainMap()->lib() == Xpetra::UseTpetra &&
        D0_Matrix_->getDomainMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
      Teuchos::ParameterList hiptmairPreList, hiptmairPostList, smootherPreList, smootherPostList;

      if (smootherList_.isSublist("smoother: pre params"))
        smootherPreList = smootherList_.sublist("smoother: pre params");
      else if (smootherList_.isSublist("smoother: params"))
        smootherPreList = smootherList_.sublist("smoother: params");
      hiptmairPreList.set("hiptmair: smoother type 1",
                           smootherPreList.get<std::string>("hiptmair: smoother type 1", "CHEBYSHEV"));
      hiptmairPreList.set("hiptmair: smoother type 2",
                           smootherPreList.get<std::string>("hiptmair: smoother type 2", "CHEBYSHEV"));
      if(smootherPreList.isSublist("hiptmair: smoother list 1"))
        hiptmairPreList.set("hiptmair: smoother list 1", smootherPreList.sublist("hiptmair: smoother list 1"));
      if(smootherPreList.isSublist("hiptmair: smoother list 2"))
        hiptmairPreList.set("hiptmair: smoother list 2", smootherPreList.sublist("hiptmair: smoother list 2"));
      hiptmairPreList.set("hiptmair: pre or post",
                           smootherPreList.get<std::string>("hiptmair: pre or post", "pre"));
      hiptmairPreList.set("hiptmair: zero starting solution",
                           smootherPreList.get<bool>("hiptmair: zero starting solution", true));

      if (smootherList_.isSublist("smoother: post params"))
        smootherPostList = smootherList_.sublist("smoother: post params");
      else if (smootherList_.isSublist("smoother: params"))
        smootherPostList = smootherList_.sublist("smoother: params");
      hiptmairPostList.set("hiptmair: smoother type 1",
                            smootherPostList.get<std::string>("hiptmair: smoother type 1", "CHEBYSHEV"));
      hiptmairPostList.set("hiptmair: smoother type 2",
                            smootherPostList.get<std::string>("hiptmair: smoother type 2", "CHEBYSHEV"));
      if(smootherPostList.isSublist("hiptmair: smoother list 1"))
        hiptmairPostList.set("hiptmair: smoother list 1", smootherPostList.sublist("hiptmair: smoother list 1"));
      if(smootherPostList.isSublist("hiptmair: smoother list 2"))
        hiptmairPostList.set("hiptmair: smoother list 2", smootherPostList.sublist("hiptmair: smoother list 2"));
      hiptmairPostList.set("hiptmair: pre or post",
                            smootherPostList.get<std::string>("hiptmair: pre or post", "post"));
      hiptmairPostList.set("hiptmair: zero starting solution",
                            smootherPostList.get<bool>("hiptmair: zero starting solution", false));

      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCRS;
      typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TROW;
      Teuchos::RCP<const TCRS> EdgeMatrix = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(SM_Matrix_ );
      Teuchos::RCP<const TCRS> NodeMatrix = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(A22_);
      Teuchos::RCP<const TCRS> PMatrix    = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(D0_Matrix_);
      hiptmairPreSmoother_  = Teuchos::rcp( new Ifpack2::Hiptmair<TROW>(EdgeMatrix,NodeMatrix,PMatrix) );
      hiptmairPreSmoother_ -> setParameters(hiptmairPreList);
      hiptmairPreSmoother_ -> initialize();
      hiptmairPreSmoother_ -> compute();
      hiptmairPostSmoother_ = Teuchos::rcp( new Ifpack2::Hiptmair<TROW>(EdgeMatrix,NodeMatrix,PMatrix) );
      hiptmairPostSmoother_ -> setParameters(hiptmairPostList);
      hiptmairPostSmoother_ -> initialize();
      hiptmairPostSmoother_ -> compute();
      useHiptmairSmoothing_ = true;
#else
      throw(Xpetra::Exceptions::RuntimeError("MueLu must be compiled with Ifpack2 for Hiptmair smoothing."));
#endif  // HAVE_MUELU_IFPACK2
    } else {
      Level level;
      RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
      level.SetFactoryManager(factoryHandler);
      level.SetLevelID(0);
      level.Set("A",SM_Matrix_);
      level.setlib(SM_Matrix_->getDomainMap()->lib());
      RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(smootherType, smootherList_));
      RCP<SmootherFactory> SmootherFact = rcp(new SmootherFactory(smootherPrototype));
      level.Request("PreSmoother",SmootherFact.get());
      SmootherFact->Build(level);
      Smoother_ = level.Get<RCP<SmootherBase> >("PreSmoother",SmootherFact.get());
      useHiptmairSmoothing_ = false;
    }
  }

  // Allocate temporary MultiVectors for solve
  P11res_ = MultiVectorFactory::Build(P11_->getDomainMap(),1);
  P11x_   = MultiVectorFactory::Build(P11_->getDomainMap(),1);
  D0res_  = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),1);
  D0x_    = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),1);
  residual_ = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(),1);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::buildProlongator() {

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out.setShowProcRank(true);

  // build prolongator: algorithm 1 in the reference paper
  // First, aggregate nodal matrix by creating a 2-level hierarchy
  Teuchos::RCP<Matrix> P;
  if (read_P_from_file_) {
    P = Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Read(P_filename_, Xpetra::UseEpetra, P->getDomainMap()->getComm());
  } else {
    Level fineLevel, coarseLevel;
    RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
    fineLevel.SetFactoryManager(factoryHandler);
    coarseLevel.SetFactoryManager(factoryHandler);
    coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
    fineLevel.SetLevelID(0);
    coarseLevel.SetLevelID(1);
    fineLevel.Set("A",A_nodal_Matrix_);

    RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
    RCP<UncoupledAggregationFactory> Aggfact = rcp(new UncoupledAggregationFactory());
    TentativePFact->SetFactory("Aggregates", Aggfact);

    coarseLevel.Request("P",TentativePFact.get());
    TentativePFact->Build(fineLevel,coarseLevel);
    P = coarseLevel.Get<RCP<Matrix> >("P",TentativePFact.get());
  }
  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("P.mat"), *P);

  // make weighting matrix
  Teuchos::RCP<Matrix> D0_Matrix_Abs=MatrixFactory2::BuildCopy(D0_Matrix_);
  D0_Matrix_Abs -> resumeFill();
  D0_Matrix_Abs -> setAllToScalar((Scalar)0.5);
  // Apply_BCsToMatrixRows(D0_Matrix_Abs,BCrows_);
  // Apply_BCsToMatrixCols(D0_Matrix_Abs,BCcols_);
  D0_Matrix_Abs -> fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("AbsD0.mat"), *D0_Matrix_Abs);

  Teuchos::RCP<Matrix> Ptent = MatrixFactory::Build(D0_Matrix_Abs->getRowMap(),0);
  Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_Abs,false,*P,false,*Ptent,true,true);

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("P_intermediate.mat"), *Ptent);

  Ptent->resumeFill();
  Apply_BCsToMatrixRows(Ptent,BCrows_);
  Ptent -> fillComplete(Ptent->getDomainMap(),Ptent->getRangeMap());

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("P_intermediate2.mat"), *Ptent);

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
    size_t nonzeros = 0;
    for(size_t j=0; j<numCols; j++)
      nonzeros += (std::abs(localVals[j])>1.0e-16);

    for(size_t j=0; j<numCols; j++) {
      for(size_t k=0; k<dim; k++) {
        blockCols.push_back(localCols[j]*dim+k);
        if (std::abs(localVals[j]) < 1.0e-16)
          blockVals.push_back(0.);
        else
          blockVals.push_back(nullspace[k][i] / nonzeros);
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

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("P11.mat"), *P11_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::formCoarseMatrix() {

  // coarse matrix for P11* (M1 + D1* M2 D1) P11
  Teuchos::RCP<Matrix> Matrix1 = MatrixFactory::Build(P11_->getDomainMap(),0);
  if (parameterList_.get<bool>("rap: triple product", false) == false) {
    Teuchos::RCP<Matrix> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),0);
    // construct (M1 + D1* M2 D1) P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*SM_Matrix_,false,*P11_,false,*C,true,true);
    // construct P11* (M1 + D1* M2 D1) P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*P11_,true,*C,false,*Matrix1,true,true);
  } else {
    Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
      MultiplyRAP(*P11_, true, *SM_Matrix_, false, *P11_, false, *Matrix1, true, true);
  }

  if(disable_addon_==true) {
    // if add-on is not chosen
    AH_=Matrix1;
  }
  else {
    // catch a failure
    TEUCHOS_TEST_FOR_EXCEPTION(M0inv_Matrix_==Teuchos::null,std::invalid_argument,
                               "MueLu::RefMaxwell::formCoarseMatrix(): Inverse of "
                               "lumped mass matrix required for add-on (i.e. M0inv_Matrix is null)");

    // coarse matrix for add-on, i.e P11* (M1 D0 M0inv D0* M1) P11
    Teuchos::RCP<Matrix> Zaux = MatrixFactory::Build(M1_Matrix_->getRowMap(),0);
    Teuchos::RCP<Matrix> Z = MatrixFactory::Build(D0_Matrix_->getDomainMap(),0);

    // construct M1 P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M1_Matrix_,false,*P11_,false,*Zaux,true,true);
    // construct Z = D0* M1 P11
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*D0_Matrix_,true,*Zaux,false,*Z,true,true);

    // construct Z* M0inv Z
    Teuchos::RCP<Matrix> Matrix2 = MatrixFactory::Build(Z->getDomainMap(),0);
    if (parameterList_.get<bool>("rap: triple product", false) == false) {
      Teuchos::RCP<Matrix> C2 = MatrixFactory::Build(M0inv_Matrix_->getRowMap(),0);
      // construct M0inv Z
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*M0inv_Matrix_,false,*Z,false,*C2,true,true);
      // construct Z* M0inv Z
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Multiply(*Z,true,*C2,false,*Matrix2,true,true);
    }
    else {
      Xpetra::TripleMatrixMultiply<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
        MultiplyRAP(*Z, true, *M0inv_Matrix_, false, *Z, false, *Matrix2, true, true);
    }
    // add matrices together
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*Matrix1,false,(Scalar)1.0,*Matrix2,false,(Scalar)1.0,AH_,*out);
    AH_->fillComplete();
  }

  // set fixed block size for vector nodal matrix
  size_t dim = Nullspace_->getNumVectors();
  AH_->SetFixedBlockSize(dim);

  if (dump_matrices_)
    Xpetra::IO<SC, LO, GlobalOrdinal, Node>::Write(std::string("AH.mat"), *AH_);

}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetMatrix(Teuchos::RCP<Matrix> SM_Matrix_new) {
  SM_Matrix_ = SM_Matrix_new;
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverseAdditive(const MultiVector& RHS, MultiVector& X) const {

  // compute residuals
  Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
  D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);

  // block diagonal preconditioner on 2x2 (V-cycle for diagonal blocks)
  HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
  Hierarchy22_->Iterate(*D0res_,  *D0x_,  1, true);

  // update current solution
  P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
  D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)1.0);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse121(const MultiVector& RHS, MultiVector& X) const {

  // precondition (1,1)-block
  Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
  HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
  P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

  // precondition (2,2)-block
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res_,  *D0x_,  1, true);
  D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

  // precondition (1,1)-block
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
  HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
  P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse212(const MultiVector& RHS, MultiVector& X) const {

  // precondition (2,2)-block
  Scalar one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res_,  *D0x_,  1, true);
  D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

  // precondition (1,1)-block
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  P11_->apply(*residual_,*P11res_,Teuchos::TRANS);
  HierarchyH_->Iterate(*P11res_, *P11x_, 1, true);
  P11_->apply(*P11x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

  // precondition (2,2)-block
  SM_Matrix_->apply(X, *residual_, Teuchos::NO_TRANS, one, zero);
  residual_->update(one, RHS, negone);
  D0_Matrix_->apply(*residual_,*D0res_,Teuchos::TRANS);
  Hierarchy22_->Iterate(*D0res_,  *D0x_,  1, true);
  D0_Matrix_->apply(*D0x_,*residual_,Teuchos::NO_TRANS);
  X.update((Scalar) 1.0, *residual_, (Scalar) 1.0);

}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply (const MultiVector& RHS, MultiVector& X,
                                                                Teuchos::ETransp mode,
                                                                Scalar alpha,
                                                                Scalar beta) const {

  // make sure that we have enough temporary memory
  if (X.getNumVectors() != P11res_->getNumVectors()) {
    P11res_    = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
    P11x_      = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
    D0res_     = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
    D0x_       = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
    residual_  = MultiVectorFactory::Build(SM_Matrix_->getDomainMap(),1);
  }

  // apply pre-smoothing
#ifdef HAVE_MUELU_IFPACK2
  if (useHiptmairSmoothing_) {
    Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX = toTpetra(X);
    Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tRHS = toTpetra(RHS);
    hiptmairPreSmoother_->apply(tRHS, tX);
  }
  else
#endif
    Smoother_->Apply(X, RHS, true);

  // do solve for the 2x2 block system
  if(mode_=="additive")
    applyInverseAdditive(RHS,X);
  else if(mode_=="121")
    applyInverse121(RHS,X);
  else if(mode_=="212")
    applyInverse212(RHS,X);
  else if(mode_=="none") {
    // do nothing
  }
  else
    applyInverseAdditive(RHS,X);

  // apply post-smoothing
#ifdef HAVE_MUELU_IFPACK2
  if (useHiptmairSmoothing_)
    {
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX = toTpetra(X);
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tRHS = toTpetra(RHS);
      hiptmairPostSmoother_->apply(tRHS, tX);
    }
  else
#endif
    Smoother_->Apply(X, RHS, false);

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

  HierarchyH_ = Teuchos::null;
  Hierarchy22_ = Teuchos::null;
  Smoother_ = Teuchos::null;
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
