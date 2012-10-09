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
/*
 * MueLu_BraessSarazinSmoother_def.hpp
 *
 *  Created on: Apr 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_BraessSarazinSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"

// include files for default FactoryManager
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BraessSarazinSmoother(LocalOrdinal sweeps, Scalar omega, RCP<const FactoryBase> AFact)
    : type_("Braess Sarazin"), nSweeps_(sweeps), omega_(omega), AFact_(AFact), A_(Teuchos::null)
  {
    RCP<SchurComplementFactory> SchurFact = Teuchos::rcp(new SchurComplementFactory(AFact_,omega));

    // define smoother/solver for BraessSarazin
    Teuchos::ParameterList SCparams;
    std::string SCtype;
    RCP<SmootherPrototype> smoProtoSC     = rcp( new DirectSolver(SCtype,SCparams,SchurFact) );

    RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

    FactManager_ = rcp(new FactoryManager());
    FactManager_->SetFactory("A", SchurFact);
    FactManager_->SetFactory("Smoother", SmooSCFact);
    FactManager_->SetIgnoreUserData(true);

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~BraessSarazinSmoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetFactoryManager(RCP<FactoryManager> FactManager) {
    FactManager_ = FactManager;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
    TEUCHOS_TEST_FOR_EXCEPTION(FactManager_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::DeclareInput: FactManager_ must not be Teuchos::null! error.");
    currentLevel.DeclareInput("PreSmoother",FactManager_->GetFactory("PreSmoother").get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //*********************************************
    // Setup routine can be summarized in 4 steps:
    // - Set the map extractors
    // - Set the blocks
    // - Create and set the inverse of the diagonal of F
    // - Set the smoother for the Schur Complement

    FactoryMonitor m(*this, "Setup blocked Braess-Sarazin Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::BreaessSarazinSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = currentLevel.Get<RCP<Matrix> > ("A", AFact_.get());

    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::BraessSarazinSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

    // store map extractors
    rangeMapExtractor_ = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    // Store the blocks in local member variables
    Teuchos::RCP<CrsMatrix> A00 = bA->getMatrix(0, 0);
    Teuchos::RCP<CrsMatrix> A01 = bA->getMatrix(0, 1);
    Teuchos::RCP<CrsMatrix> A10 = bA->getMatrix(1, 0);
    Teuchos::RCP<CrsMatrix> A11 = bA->getMatrix(1, 1);

    Teuchos::RCP<CrsMatrixWrap> Op00 = Teuchos::rcp(new CrsMatrixWrap(A00));
    Teuchos::RCP<CrsMatrixWrap> Op01 = Teuchos::rcp(new CrsMatrixWrap(A01));
    Teuchos::RCP<CrsMatrixWrap> Op10 = Teuchos::rcp(new CrsMatrixWrap(A10));
    Teuchos::RCP<CrsMatrixWrap> Op11 = Teuchos::rcp(new CrsMatrixWrap(A11));

    F_ = Teuchos::rcp_dynamic_cast<Matrix>(Op00);
    G_ = Teuchos::rcp_dynamic_cast<Matrix>(Op01);
    D_ = Teuchos::rcp_dynamic_cast<Matrix>(Op10);
    Z_ = Teuchos::rcp_dynamic_cast<Matrix>(Op11);

    // Create the inverse of the diagonal of F
    RCP<Vector> diagFVector = VectorFactory::Build(F_->getRowMap());
    F_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F

    ////////// EXPERIMENTAL
    // fix zeros on diagonal
    /*Teuchos::ArrayRCP< Scalar > diagFdata = diagFVector->getDataNonConst(0);
    for(size_t t = 0; t < diagFdata.size(); t++) {
      if(diagFdata[t] == 0.0) {
        std::cout << "fixed zero diagonal entry" << std::endl;
        diagFdata[t] = 1.0;
      }
    }*/
    ////////// EXPERIMENTAL

    diagFVector->reciprocal(*diagFVector);    // build reciprocal
    diagFinv_ = diagFVector;

    // Set the Smoother
    smoo_ = currentLevel.Get<RCP<SmootherBase> > ("PreSmoother", FactManager_->GetFactory("PreSmoother").get());

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::Apply(): Setup() has not been called");

    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    // TODO change notation for v and p
    RCP<MultiVector> vtemp = MultiVectorFactory::Build(F_->getRowMap(),1);
    RCP<MultiVector> ptemp = MultiVectorFactory::Build(Z_->getRowMap(),1);

    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);

    for (LocalOrdinal run = 0; run < nSweeps_; ++run) {
      // Calculate residual
      // note: A_ is the full blocked operator
      // r = B - A*X
      residual->update(1.0,B,0.0); // residual = B
      A_->apply(X, *residual, Teuchos::NO_TRANS, -1.0, 1.0); // residual = B-A*X

      // split the data
      // extract corresponding subvectors from X and residual
      Teuchos::RCP<MultiVector> rvel = rangeMapExtractor_->ExtractVector(residual, 0);
      Teuchos::RCP<MultiVector> rpre = rangeMapExtractor_->ExtractVector(residual, 1);

      Teuchos::RCP<MultiVector> tXvel = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> tXpre = domainMapExtractor_->ExtractVector(rcpX, 1);

      // Calculate qrhs = rpre - D * vtemp (equation 8.14)
      RCP<MultiVector> qrhs = MultiVectorFactory::Build(Z_->getRowMap(),1);
      vtemp->putScalar(0.0);
      vtemp->elementWiseMultiply(1.0/omega_,*diagFinv_,*rvel,0.0); // vtemp = 1/omega*Fhatinv * rvel (equation 8.13)

      D_->apply(*vtemp,*qrhs); // qrhs = D*vtemp (intermediate step)
      qrhs->update(1.0,*rpre,-1.0); // qrhs = rpre - D*vtemp

      //Pressure correction, using the preconditioner
      RCP<MultiVector> q = MultiVectorFactory::Build(Z_->getRowMap(),1); // TODO think about this.
      q->putScalar(0.0);  // just for safety
      smoo_->Apply(*q,*qrhs);

      //Update
      vtemp->putScalar(0.0);
      G_->apply(*q,*vtemp);
      vtemp->update(1.0,*rvel,-1.0); //velres - G*q

      RCP<MultiVector> vx = MultiVectorFactory::Build(F_->getRowMap(),1);
      vx->elementWiseMultiply(1.0/omega_,*diagFinv_,*vtemp,1.0);

      // update with values
      tXvel->update(1.0,*vx,1.0);
      tXpre->update(1.0,*q,1.0);

      domainMapExtractor_->InsertVector(tXvel, 0, rcpX);
      domainMapExtractor_->InsertVector(tXpre, 1, rcpX);
    }

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new BraessSarazinSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#endif /* MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_ */
