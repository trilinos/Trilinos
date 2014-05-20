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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const {
  return Xpetra::toTpetraNonZero(SM_Matrix_->getDomainMap());
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const {
  return Xpetra::toTpetraNonZero(SM_Matrix_->getRangeMap());
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setParameters(Teuchos::ParameterList& list) {
  disable_addon_     =  list.get("refmaxwell: disable add-on",false);
  MaxCoarseSize_     =  list.get("refmaxwell: max coarse size",1000);
  MaxLevels_         =  list.get("refmaxwell: max levels",5);
  EdgeSmoother_      =  list.get("refmaxwell: edge smoother","CHEBYSHEV");
  NodeSmoother_      =  list.get("refmaxwell: node smoother","CHEBYSHEV");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::compute() {
  // build special prolongator for (1,1)-block
  if(P11_==Teuchos::null) {
    buildProlongator();
  }
  // build coarse grid operator for (1,1)-block
  if(A11_==Teuchos::null) {
    formCoarseMatrix();
  }
  // build fine grid operator for (2,2)-block, D0* M1 D0
  if(A22_==Teuchos::null) {
    Teuchos::RCP<XMat> C = MatrixFactory::Build(M1_Matrix_->getRowMap(),100);
    Xpetra::MatrixMatrix::Multiply(*M1_Matrix_,false,*D0_Matrix_,false,*C,true,true);
    A22_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),100);
    Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*C,false,*A22_,true,true);
    A22_->SetFixedBlockSize(1);
  }
  // build stuff for hierarchies
  Teuchos::RCP<FactoryManager> Manager11 = Teuchos::rcp( new FactoryManager );
  Teuchos::RCP<FactoryManager> Manager22 = Teuchos::rcp( new FactoryManager );
  Teuchos::RCP<SaPFactory> Pfact = Teuchos::rcp( new SaPFactory );
  std::string precType11 = EdgeSmoother_;
  std::string precType22 = NodeSmoother_;
  Teuchos::ParameterList precList11, precList22;
  precList11.set("chebyshev: degree",3);
  precList22.set("chebyshev: degree",3);
  Teuchos::RCP<SmootherPrototype> SmooProto11
    = Teuchos::rcp( new Ifpack2Smoother(precType11,precList11) );
  Teuchos::RCP<SmootherFactory> SmooFact11
    = Teuchos::rcp( new SmootherFactory(SmooProto11) );
  Teuchos::RCP<SmootherPrototype> SmooProto22
    = Teuchos::rcp( new Ifpack2Smoother(precType22,precList22) );
  Teuchos::RCP<SmootherFactory> SmooFact22
    = Teuchos::rcp( new SmootherFactory(SmooProto22) );
  Teuchos::RCP<UncoupledAggregationFactory> Aggfact
    = Teuchos::rcp( new UncoupledAggregationFactory() );
  Manager11->SetFactory("Aggregates",Aggfact);
  Manager11->SetFactory("Smoother",SmooFact11);
  Manager11->SetFactory("CoarseSolver",SmooFact11);
  Manager22->SetFactory("Aggregates",Aggfact);
  Manager22->SetFactory("Smoother",SmooFact22);
  Manager22->SetFactory("CoarseSolver",SmooFact22);
  Hierarchy11_ = Teuchos::rcp( new Hierarchy(A11_) );
  Hierarchy11_ -> SetMaxCoarseSize( MaxCoarseSize_ );
  Hierarchy11_ -> Setup(*Manager11, 0, MaxLevels_ );
  Hierarchy22_ = Teuchos::rcp( new Hierarchy(A22_) );
  Hierarchy22_ -> SetMaxCoarseSize( MaxCoarseSize_ );
  Hierarchy22_ -> Setup(*Manager22, 0, MaxLevels_ );
  // build ifpack2 preconditioners for Hiptmair
  Teuchos::RCP<const TCRS> EdgeMatrix = Utils::Op2NonConstTpetraCrs(SM_Matrix_ );
  Teuchos::RCP<const TCRS> NodeMatrix = Utils::Op2NonConstTpetraCrs(TMT_Matrix_);
  edgePrec_ = Ifpack2::Factory::create(EdgeSmoother_, EdgeMatrix);
  nodePrec_ = Ifpack2::Factory::create(NodeSmoother_, NodeMatrix);
  edgePrec_->setParameters(precList11); edgePrec_->initialize(); edgePrec_->compute();
  nodePrec_->setParameters(precList22); nodePrec_->initialize(); nodePrec_->compute();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::buildProlongator() {

  // Form TMT_Matrix
  Teuchos::RCP<XMat> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),100);
  TMT_Matrix_=MatrixFactory::Build(D0_Matrix_->getDomainMap(),100);
  Xpetra::MatrixMatrix::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,*C,true,true);
  Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*C,false,*TMT_Matrix_,true,true);
  TMT_Matrix_->SetFixedBlockSize(1);

  // build nullspace
  Nullspace_ = MultiVectorFactory::Build(SM_Matrix_->getRowMap(),Coords_->getNumVectors()); 
  D0_Matrix_->apply(*Coords_,*Nullspace_);

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
  Teuchos::ParameterList params1;
  params1.set("Damping factor",(Scalar)0.0);
  Pfact      -> SetParameterList(params1);
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
  D0_Matrix_Abs -> fillComplete(D0_Matrix_->getDomainMap(),D0_Matrix_->getRangeMap());
  Teuchos::RCP<XMat> Ptent = MatrixFactory::Build(D0_Matrix_Abs->getRowMap(),100);
  Xpetra::MatrixMatrix::Multiply(*D0_Matrix_Abs,false,*P,false,*Ptent,true,true);

  // put in entries to P11
  size_t dim = Nullspace_->getNumVectors();
  size_t numLocalRows = SM_Matrix_->getNodeNumRows();
  Teuchos::RCP<XMap> BlockColMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getColMap(),dim);
  P11_ = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Build(Ptent->getRowMap(),BlockColMap,100);
  
  std::vector< Teuchos::ArrayRCP<const Scalar> > nullspace(dim);
  for(size_t i=0; i<dim; i++) {
    Teuchos::ArrayRCP<const Scalar> datavec = Nullspace_->getData(i);
    nullspace[i]=datavec;
  }
  
  for(size_t i=0; i<numLocalRows; i++) {
    Teuchos::ArrayView<const LocalOrdinal> localCols;
    Teuchos::ArrayView<const Scalar>       localVals;
    Ptent->getLocalRowView(i,localCols,localVals);
    size_t numCols = localCols.size();
    std::vector<LocalOrdinal> blockLocalCols(dim*numCols);
    std::vector<Scalar>       blockLocalVals(dim*numCols);
    for(size_t j=0; j<numCols; j++) {
      for(size_t k=0; k<dim; k++) {
	blockLocalCols[j*dim+k] = localCols[j]*dim+k;
	blockLocalVals[j*dim+k] = localVals[j]*nullspace[k][i];
      }
    }
    P11_ -> insertLocalValues(i,
			      Teuchos::ArrayView<LocalOrdinal>(blockLocalCols),
			      Teuchos::ArrayView<Scalar>(blockLocalVals));
  }
  Teuchos::RCP<XMap> blockCoarseMap
    = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Ptent->getDomainMap(),dim);
  P11_->fillComplete(blockCoarseMap,SM_Matrix_->getDomainMap());

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::formCoarseMatrix() {

  // coarse matrix for P11* (M1 + D1* M2 D1) P11
  Teuchos::RCP<XMat> C = MatrixFactory::Build(SM_Matrix_->getRowMap(),100);
  Teuchos::RCP<XMat> Matrix1 = MatrixFactory::Build(P11_->getDomainMap(),100);
  // construct (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix::Multiply(*SM_Matrix_,false,*P11_,false,*C,true,true);
  // construct P11* (M1 + D1* M2 D1) P11
  Xpetra::MatrixMatrix::Multiply(*P11_,true,*C,false,*Matrix1,true,true);

  if(disable_addon_==true) {
    // if add-on is not chosen
    A11_=Matrix1;
  }
  else {
    // coarse matrix for add-on, i.e P11* (M1 D0 M0inv D0* M1) P11
    Teuchos::RCP<XMat> Zaux = MatrixFactory::Build(M1_Matrix_->getRowMap(),100);
    Teuchos::RCP<XMat> Z = MatrixFactory::Build(D0_Matrix_->getDomainMap(),100);
    Teuchos::RCP<XMat> C2 = MatrixFactory::Build(M0inv_Matrix_->getRowMap(),100);
    // construct M1 P11
    Xpetra::MatrixMatrix::Multiply(*M1_Matrix_,false,*P11_,false,*Zaux,true,true);
    // construct Z = D0* M1 P11
    Xpetra::MatrixMatrix::Multiply(*D0_Matrix_,true,*Zaux,false,*Z,true,true);
    // construct M0inv Z
    Xpetra::MatrixMatrix::Multiply(*M0inv_Matrix_,false,*Z,false,*C2,true,true);
    // construct Z* M0inv Z
    Teuchos::RCP<XMat> Matrix2 = MatrixFactory::Build(Z->getDomainMap(),100);
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyHiptmairSmoother(const XTMV& RHS, XTMV& X) const {

  RCP<XMV> edge_residual, node_residual;
  RCP<TMV> edge_tpetra,   node_tpetra;

  // apply initial relaxation to edge elements
  edge_residual  =  Utils::Residual(*SM_Matrix_, X, RHS);
  edge_tpetra    =  Utils::MV2NonConstTpetraMV(edge_residual);
  edgePrec_->apply(*edge_tpetra,*edge_tpetra);
  edge_residual  =  Xpetra::toXpetra(edge_tpetra);
  X.update((Scalar) 1.0, *edge_residual, (Scalar) 1.0);
  // project to nodal space and smooth
  edge_residual  =  Utils::Residual(*SM_Matrix_, X, RHS);
  node_residual  =  MultiVectorFactory::Build(TMT_Matrix_->getRowMap(),RHS.getNumVectors());
  D0_Matrix_->apply(*edge_residual,*node_residual,Teuchos::TRANS,(Scalar)1.0,(Scalar)0.0);
  node_tpetra    =  Utils::MV2NonConstTpetraMV(node_residual);
  nodePrec_->apply(*node_tpetra,*node_tpetra);
  node_residual  =  Xpetra::toXpetra(node_tpetra);
  D0_Matrix_->apply(*node_residual,*edge_residual,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)0.0);
  X.update((Scalar) 1.0, *edge_residual, (Scalar) 1.0);
  // smooth again on edge elements
  edge_residual  =  Utils::Residual(*SM_Matrix_, X, RHS);
  edge_tpetra    =  Utils::MV2NonConstTpetraMV(edge_residual);
  edgePrec_->apply(*edge_tpetra,*edge_tpetra);
  edge_residual  =  Xpetra::toXpetra(edge_tpetra);
  X.update((Scalar) 1.0, *edge_residual, (Scalar) 1.0);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
									   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
									   Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  try {
    TMV& temp_x = const_cast<TMV &>(X);
    const XTMV tX(rcpFromRef(temp_x));
    XTMV       tY(rcpFromRef(Y));
    tY.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    // apply pre-smoothing
    applyHiptmairSmoother(tX,tY);
    // do solve for the 2x2 block system
    // first, compute residuals
    RCP<XMV> residual = Utils::Residual(*SM_Matrix_, tY, tX);
    RCP<XMV> P11residual = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
    RCP<XMV> P11result = MultiVectorFactory::Build(P11_->getDomainMap(),X.getNumVectors());
    RCP<XMV> D0residual = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
    RCP<XMV> D0result = MultiVectorFactory::Build(D0_Matrix_->getDomainMap(),X.getNumVectors());
    P11_->apply(*residual,*P11residual,Teuchos::TRANS);
    D0_Matrix_->apply(*residual,*D0residual,Teuchos::TRANS);
    // block diagonal preconditioner on 2x2 (V-cycle for each block)
    Hierarchy11_->Iterate(*P11residual, *P11result, 1, true);
    Hierarchy22_->Iterate(*D0residual, *D0result, 1, true);
    // update current solution
    P11_->apply(*P11result,*residual,Teuchos::NO_TRANS);
    D0_Matrix_->apply(*D0result,*residual,Teuchos::NO_TRANS,(Scalar)1.0,(Scalar)1.0);
    tY.update((Scalar) 1.0, *residual, (Scalar) 1.0);
    // apply post-smoothing
    applyHiptmairSmoother(tX,tY);
  } catch (std::exception& e) {
    //FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::TpetraOperator::ApplyInverse():" << std::endl
	      << e.what() << std::endl;
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const {
  return false;
}

} // namespace
#endif //ifdef HAVE_MUELU_TPETRA

#define MUELU_REFMAXWELL_SHORT
#endif //ifdef MUELU_REFMAXWELL_DEF_HPP
