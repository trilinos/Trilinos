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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RAPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RAPFactory(RCP<const FactoryBase> PFact, RCP<const FactoryBase> RFact, RCP<const FactoryBase> AFact)
: PFact_(PFact), RFact_(RFact), AFact_(AFact), implicitTranspose_(false), checkAc_(false), repairZeroDiagonals_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~RAPFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  fineLevel.DeclareInput("A", AFact_.get(), this);  // AFact per default Teuchos::null -> default factory for this
  coarseLevel.DeclareInput("P", PFact_.get(), this); // transfer operators (from PRFactory, not from PFactory and RFactory!)
  coarseLevel.DeclareInput("R", RFact_.get(), this); //TODO: must be request according to (implicitTranspose flag!!!!!
  coarseLevel.DeclareInput("A", PFact_.get(), this); //FIXME hack
  coarseLevel.DeclareInput("Importer", RFact_.get(), this); //FIXME hack, could result in redundant work...

  // call DeclareInput of all user-given transfer factories
  for(std::vector<RCP<const FactoryBase> >::const_iterator it = TransferFacts_.begin(); it!=TransferFacts_.end(); ++it) {
    (*it)->CallDeclareInput(coarseLevel);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PrintMatrixInfo(const Matrix & Ac, const std::string & msgTag) const {
  GetOStream(Statistics0, 0) << msgTag
      << " # global rows = "      << Ac.getGlobalNumRows()
      << ", estim. global nnz = " << Ac.getGlobalNumEntries()
      << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PrintLoadBalancingInfo(const Matrix & Ac, const std::string & msgTag) const {
  // TODO: provide a option to skip this (to avoid global communication)

  //nonzero imbalance
  size_t numMyNnz  = Ac.getNodeNumEntries();
  GO maxNnz, minNnz;
  RCP<const Teuchos::Comm<int> > comm = Ac.getRowMap()->getComm();
  maxAll(comm,(GO)numMyNnz,maxNnz);
  //min nnz over all proc (disallow any processors with 0 nnz)
  minAll(comm, (GO)((numMyNnz > 0) ? numMyNnz : maxNnz), minNnz);
  double imbalance = ((double) maxNnz) / minNnz;

  size_t numMyRows = Ac.getNodeNumRows();
  //Check whether Ac is spread over more than one process.
  GO numActiveProcesses=0;
  sumAll(comm, (GO)((numMyRows > 0) ? 1 : 0), numActiveProcesses);

  //min, max, and avg # rows per proc
  GO minNumRows, maxNumRows;
  double avgNumRows;
  maxAll(comm, (GO)numMyRows, maxNumRows);
  minAll(comm, (GO)((numMyRows > 0) ? numMyRows : maxNumRows), minNumRows);
  assert(numActiveProcesses > 0);
  avgNumRows = Ac.getGlobalNumRows() / numActiveProcesses;

  GetOStream(Statistics1,0) << msgTag << " # processes with rows = " << numActiveProcesses << std::endl;
  GetOStream(Statistics1,0) << msgTag << " min # rows per proc = " << minNumRows << ", max # rows per proc = " << maxNumRows << ", avg # rows per proc = " << avgNumRows << std::endl;
  GetOStream(Statistics1,0) << msgTag << " nonzero imbalance = " << imbalance << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsMatrixClass; // TODO move me

  { //scoping
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);

    int levelID = coarseLevel.GetLevelID();

    //
    // Inputs: A, P
    //

    RCP<Matrix> A = fineLevel.Get< RCP<Matrix> >("A", AFact_.get());
    RCP<Matrix> P = coarseLevel.Get< RCP<Matrix> >("P", PFact_.get());

    //
    // Build Ac = RAP
    //
    RCP<Matrix> Ac;

    const RCP<BlockedCrsMatrixClass> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(A);
    if( bA != Teuchos::null) {
      RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", RFact_.get());
      Ac = BuildRAPBlock(R, A, P, levelID); // Triple matrix product for BlockedCrsMatrixClass
      // TODO add plausibility check
    } else {

      // refactoring needed >>>>>>>>
      // the following code should be moved to his own factory
      //
      if ( coarseLevel.IsAvailable("A",PFact_.get()) && coarseLevel.IsAvailable("Importer",RFact_.get()) ) {
        SubFactoryMonitor m1(*this, "Rebalancing existing Ac", levelID);
        Ac = coarseLevel.Get< RCP<Matrix> >("A", PFact_.get());
        RCP<Matrix> newAc = MatrixFactory::Build(P->getDomainMap(), Ac->getGlobalMaxNumRowEntries());
        RCP<CrsMatrixWrap> crsOp = rcp_dynamic_cast<CrsMatrixWrap>(newAc);
        RCP<CrsMatrix> crsMtx = crsOp->getCrsMatrix();
        RCP<CrsMatrixWrap> origOp = rcp_dynamic_cast<CrsMatrixWrap>(Ac);
        RCP<CrsMatrix> origMtx = origOp->getCrsMatrix();
        RCP<const Import> permImporter = coarseLevel.Get< RCP<const Import> >("Importer",RFact_.get());
        crsMtx->doImport(*origMtx, *permImporter,Xpetra::INSERT);
        crsMtx = Teuchos::null;
        //TODO add plausibility check

        newAc->fillComplete(P->getDomainMap(), P->getDomainMap());
        Ac = newAc;

        PrintMatrixInfo(*Ac, "Ac (rebalanced)");
        PrintLoadBalancingInfo(*Ac, "Ac (rebalanced)");

      } else if (coarseLevel.IsAvailable("A",PFact_.get())) {
        // Ac already built by the load balancing process and no load balancing needed
        SubFactoryMonitor m1(*this, "Ac already computed", levelID);
        Ac = coarseLevel.Get< RCP<Matrix> >("A", PFact_.get());
      }
      // <<<<<<<< refactoring needed

      else {
        // Triple matrix product

        if (implicitTranspose_) {
          // implicit version
          Ac = BuildRAPImplicit(A, P, levelID);
        } else {
          // explicit version
          RCP<Matrix> R = coarseLevel.Get< RCP<Matrix> >("R", RFact_.get());
          Ac = BuildRAPExplicit(R, A, P, levelID);
        }

      }
    } //if( bA != Teuchos::null) ... else ...

    TEUCHOS_TEST_FOR_EXCEPT(Ac == Teuchos::null);

    coarseLevel.Set("A", Ac, this);

  } //scoping

  if (TransferFacts_.begin() != TransferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for(std::vector<RCP<const FactoryBase> >::const_iterator it = TransferFacts_.begin(); it != TransferFacts_.end(); ++it) {
      GetOStream(Runtime0, 0) << "Ac: call transfer factory " << (*it).get() << ": " << (*it)->description() << std::endl;
      (*it)->CallBuild(coarseLevel);
    }
  }
} //Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPExplicit(const RCP<Matrix>& R, const RCP<Matrix>& A, const RCP<Matrix>& P, int levelID) const {
  SubFactoryMonitor m(*this, "Build RAP explicitly", levelID);

  RCP<Matrix> AP;
  {
    SubFactoryMonitor m2(*this, "MxM: A x P", levelID);
    AP = Utils::TwoMatrixMultiply(A, false, P, false);
    //std::string filename="AP.dat";
    //Utils::Write(filename, AP);
  }

  RCP<Matrix> RAP;
  {
    SubFactoryMonitor m2(*this, "MxM: R x (AP)", levelID);
    RAP = Utils::TwoMatrixMultiply(R, false, AP, false, true, false);  // FIXME: no optimization of storage since insertLocalValues cannot be called after fillComplete when values are optimized out (Epetra)
  }

  if(checkAc_) CheckMainDiagonal(RAP);

  PrintMatrixInfo(*RAP, "Ac (explicit)");
  PrintLoadBalancingInfo(*RAP, "Ac (explicit)");

  return RAP;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPImplicit(const RCP<Matrix>& A, const RCP<Matrix>& P, int levelID) const {
  SubFactoryMonitor m(*this, "Build RAP implicitly", levelID);

  GetOStream(Warnings0, 0) << "The implicitTranspose_ flag within RAPFactory for Epetra in parallel produces wrong results" << std::endl;
  RCP<Matrix> AP  = Utils::TwoMatrixMultiply(A, false, P, false);
  RCP<Matrix> RAP = Utils::TwoMatrixMultiply(P, true, AP, false, true, false); // FIXME: no optimization of storage since insertLocalValues cannot be called after fillComplete when values are optimized out (Epetra)

  if(checkAc_) CheckMainDiagonal(RAP);

  PrintMatrixInfo(*RAP, "Ac (implicit)");
  return RAP;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildRAPBlock(const RCP<Matrix>& R, const RCP<Matrix>& A, const RCP<Matrix>& P, int levelID) const {
  SubFactoryMonitor m(*this, "Build RAP block", levelID);

  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsMatrixClass;
  const RCP<BlockedCrsMatrixClass> bR = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(R);
  const RCP<BlockedCrsMatrixClass> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(A);
  const RCP<BlockedCrsMatrixClass> bP = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(P);
  TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix A is not of type BlockedCrsMatrix! error.");
  TEUCHOS_TEST_FOR_EXCEPTION(bP==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix P is not of type BlockedCrsMatrix! error.");
  TEUCHOS_TEST_FOR_EXCEPTION(bR==Teuchos::null, Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: input matrix R is not of type BlockedCrsMatrix! error.");
  if(implicitTranspose_) GetOStream(Warnings0, 0) << "No support for implicitTranspose_ flag within RAPFactory for blocked matrices" << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols()!=bP->Rows(), Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: block matrix dimensions do not match. error.");
  TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows()!=bR->Cols(), Exceptions::BadCast, "MueLu::RAPFactory::BuildRAPblock: block matrix dimensions do not match. error.");

  RCP<BlockedCrsMatrixClass> bAP  = Utils::TwoMatrixMultiplyBlock(bA, false, bP, false, true, true);

  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  //bAP->describe(*fos,Teuchos::VERB_EXTREME);

  RCP<BlockedCrsMatrixClass> bRAP = Utils::TwoMatrixMultiplyBlock(bR, false, bAP, false, true, true);

  // TODO add CheckMainDiagonal for Blocked operator
  // avoid copying block matrix!
  // create new empty Operator
  if(checkAc_) {
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps> > c00 = bRAP->getMatrix(0,0);
    Teuchos::RCP<CrsMatrix> Aout = Xpetra::CrsMatrixFactory<Scalar,LocalOrdinal, GlobalOrdinal, Node,LocalMatOps>::Build(c00->getRowMap(),c00->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile);

    RCP<Vector> diagVec = VectorFactory::Build(c00->getRowMap());
    c00->getLocalDiagCopy(*diagVec);
    Teuchos::ArrayRCP< Scalar > diagVal = diagVec->getDataNonConst(0);

    // loop over local rows
    for(size_t row=0; row<c00->getNodeNumRows(); row++) {
      // get global row id
      GlobalOrdinal grid = c00->getRowMap()->getGlobalElement(row); // global row id

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      c00->getLocalRowView(row, indices, vals);

      // just copy all values in output
      Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
      Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

      //just copy values
      for(size_t i=0; i<(size_t)indices.size(); i++) {
        GlobalOrdinal gcid = c00->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
        indout [i] = gcid;
        valout [i] = vals[i];
      }

      Aout->insertGlobalValues(grid, indout.view(0,indout.size()), valout.view(0,valout.size()));
      if(diagVal[row]==0.0 && repairZeroDiagonals_) {
        // always overwrite diagonal entry
        Aout->insertGlobalValues(grid, Teuchos::tuple<GlobalOrdinal>(grid),Teuchos::tuple<Scalar>(1.0));
      }
    }

    Aout->fillComplete(c00->getDomainMap(), c00->getRangeMap());

    bRAP->setMatrix(0,0,Aout);
  }

  PrintMatrixInfo(*bRAP, "Ac (blocked)");
  return bRAP;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheckMainDiagonal(RCP<Matrix> & Ac) const {
  // plausibility check: no zeros on diagonal
  LO lZeroDiags = 0;
  RCP<Vector> diagVec = VectorFactory::Build(Ac->getRowMap());
  Ac->getLocalDiagCopy(*diagVec);
  Teuchos::ArrayRCP< Scalar > diagVal = diagVec->getDataNonConst(0);
  for (size_t r=0; r<Ac->getRowMap()->getNodeNumElements(); r++) {
    if(diagVal[r]==0.0) {
      lZeroDiags++;
      if(repairZeroDiagonals_) {
        GlobalOrdinal grid = Ac->getRowMap()->getGlobalElement(r);
        LocalOrdinal lcid = Ac->getColMap()->getLocalElement(grid);
        Teuchos::ArrayRCP<LocalOrdinal> indout(1,lcid);
        Teuchos::ArrayRCP<Scalar> valout(1,Teuchos::ScalarTraits<Scalar>::one());
        Ac->insertLocalValues(r, indout.view(0,indout.size()), valout.view(0,valout.size()));
      }
    }
  }

  if(IsPrint(Warnings0)) {
    const RCP<const Teuchos::Comm<int> > & comm = Ac->getRowMap()->getComm();
    GO lZeroDiagsGO = lZeroDiags; /* LO->GO conversion */
    GO gZeroDiags = 0;
    sumAll(comm, lZeroDiagsGO, gZeroDiags);
    if(repairZeroDiagonals_) GetOStream(Warnings0,0) << "RAPFactory (WARNING): repaired " << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
    else                     GetOStream(Warnings0,0) << "RAPFactory (WARNING): found "    << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheckMainDiagonal(RCP<CrsMatrix> & Ac) const {
  // plausibility check: no zeros on diagonal
  LO lZeroDiags = 0;
  RCP<Vector> diagVec = VectorFactory::Build(Ac->getRowMap());
  Ac->getLocalDiagCopy(*diagVec);
  Teuchos::ArrayRCP< Scalar > diagVal = diagVec->getDataNonConst(0);
  for (size_t r=0; r<Ac->getRowMap()->getNodeNumElements(); r++) {
    if(diagVal[r]==0.0) {
      lZeroDiags++;
      //if(repairZeroDiagonals_) {
        GlobalOrdinal grid = Ac->getRowMap()->getGlobalElement(r);
        LocalOrdinal lcid = Ac->getColMap()->getLocalElement(grid);
        Teuchos::ArrayRCP<LocalOrdinal> indout(1,lcid);
        Teuchos::ArrayRCP<Scalar> valout(1,Teuchos::ScalarTraits<Scalar>::one());
        Ac->insertLocalValues(r, indout.view(0,indout.size()), valout.view(0,valout.size()));
      //}
    }
  }

  if(IsPrint(Warnings0)) {
    const RCP<const Teuchos::Comm<int> > & comm = Ac->getRowMap()->getComm();
    GO lZeroDiagsGO = lZeroDiags; /* LO->GO conversion */
    GO gZeroDiags = 0;
    sumAll(comm, lZeroDiagsGO, gZeroDiags);
    if(repairZeroDiagonals_) GetOStream(Warnings0,0) << "RAPFactory (WARNING): repaired " << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
    else                     GetOStream(Warnings0,0) << "RAPFactory (WARNING): found "    << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetImplicitTranspose(bool const &implicit) {
  implicitTranspose_ = implicit;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TransferFacts_.push_back(factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
size_t RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NumTransferFactories() const {
  return TransferFacts_.size();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const FactoryBase> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetPFactory() const {
  return PFact_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<const FactoryBase> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetRFactory() const {
  return RFact_;
}


} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DEF_HPP
