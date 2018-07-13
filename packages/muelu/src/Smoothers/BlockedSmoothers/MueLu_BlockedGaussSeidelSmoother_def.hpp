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
/*
 * MueLu_BlockedGaussSeidelSmoother_def.hpp
 *
 *  Created on: 30.01.2012
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_
#define MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedMultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_BlockedGaussSeidelSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedGaussSeidelSmoother()
    : type_("blocked GaussSeidel"), A_(Teuchos::null)
  {
    FactManager_.reserve(10); // TODO fix me!
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BlockedGaussSeidelSmoother() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< Scalar >                ("Damping factor",     1.0, "Damping/Scaling factor in BGS");
    validParamList->set< LocalOrdinal >          ("Sweeps",             1, "Number of BGS sweeps (default = 1)");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
    TEUCHOS_TEST_FOR_EXCEPTION(pos < 0, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::AddFactoryManager: parameter \'pos\' must not be negative! error.");

    size_t myPos = Teuchos::as<size_t>(pos);

    if (myPos < FactManager_.size()) {
      // replace existing entris in FactManager_ vector
      FactManager_.at(myPos) = FactManager;
    } else if( myPos == FactManager_.size()) {
      // add new Factory manager in the end of the vector
      FactManager_.push_back(FactManager);
    } else { // if(myPos > FactManager_.size())
      RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      *out << "Warning: cannot add new FactoryManager at proper position " << pos << ". The FactoryManager is just appended to the end. Check this!" << std::endl;

      // add new Factory manager in the end of the vector
      FactManager_.push_back(FactManager);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //this->Input(currentLevel, "A");
    // TODO: check me: why is this->Input not freeing properly A in release mode?
    currentLevel.DeclareInput("A",this->GetFactory("A").get());

    // loop over all factory managers for the subblocks of blocked operator A
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   *it);

      // request "Smoother" for current subblock row.
      currentLevel.DeclareInput("PreSmoother",(*it)->GetFactory("Smoother").get());

      // request "A" for current subblock row (only needed for Thyra mode)
      currentLevel.DeclareInput("A",(*it)->GetFactory("A").get());
    }

    //RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
    //typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsOMatrix;

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    FactoryMonitor m(*this, "Setup blocked Gauss-Seidel Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true) this->GetOStream(Warnings0) << "MueLu::BlockedGaussSeidelSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A"); // A needed for extracting map extractors
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

    // plausibility check
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != FactManager_.size(), Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Setup: number of block rows of A is " << bA->Rows() << " and does not match number of SubFactoryManagers " << FactManager_.size() << ". error.");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Cols() != FactManager_.size(), Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Setup: number of block cols of A is " << bA->Cols() << " and does not match number of SubFactoryManagers " << FactManager_.size() << ". error.");

    // store map extractors
    rangeMapExtractor_  = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    // loop over all factory managers for the subblocks of blocked operator A
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel), *it);

      // extract Smoother for current block row (BGS ordering)
      RCP<const SmootherBase> Smoo = currentLevel.Get< RCP<SmootherBase> >("PreSmoother",(*it)->GetFactory("Smoother").get());
      Inverse_.push_back(Smoo);

      // store whether subblock matrix is blocked or not!
      RCP<Matrix> Aii = currentLevel.Get< RCP<Matrix> >("A",(*it)->GetFactory("A").get());
      bIsBlockedOperator_.push_back(Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Aii)!=Teuchos::null);
    }

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector& B, bool InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): Setup() has not been called");

#if 0 // def HAVE_MUELU_DEBUG
    // TODO simplify this debug check
    RCP<MultiVector> rcpDebugX = Teuchos::rcpFromRef(X);
    RCP<const MultiVector> rcpDebugB = Teuchos::rcpFromRef(B);
    RCP<BlockedMultiVector> rcpBDebugX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpDebugX);
    RCP<const BlockedMultiVector> rcpBDebugB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpDebugB);
    //RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    if(rcpBDebugB.is_null() == false) {
      //this->GetOStream(Runtime1) << "BlockedGaussSeidel: B is a BlockedMultiVector of size " << B.getMap()->getGlobalNumElements() << " with " << rcpBDebugB->getBlockedMap()->getNumMaps() << " blocks." << std::endl;
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    } else {
      //this->GetOStream(Runtime1) << "BlockedGaussSeidel: B is a MultiVector of size " << B.getMap()->getGlobalNumElements() << std::endl;
      //TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    }
    if(rcpBDebugX.is_null() == false) {
      //this->GetOStream(Runtime1) << "BlockedGaussSeidel: X is a BlockedMultiVector of size " << X.getMap()->getGlobalNumElements() << " with " << rcpBDebugX->getBlockedMap()->getNumMaps() << " blocks." << std::endl;
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    } else {
      //this->GetOStream(Runtime1) << "BlockedGaussSeidel: X is a MultiVector of size " << X.getMap()->getGlobalNumElements() << std::endl;
      //TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    }


#endif
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Input variables used for the rest of the algorithm
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);
    RCP<const MultiVector> rcpB = Teuchos::rcpFromRef(B);

    // make sure that both rcpX and rcpB are BlockedMultiVector objects
    bool bCopyResultX = false;
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    MUELU_TEST_FOR_EXCEPTION(bA.is_null() == true, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): A_ must be a BlockedCrsMatrix");
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
    RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

    if(bX.is_null() == true) {
      RCP<MultiVector> test = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedDomainMap(),rcpX));
      rcpX.swap(test);
      bCopyResultX = true;
    }

    if(bB.is_null() == true) {
      RCP<const MultiVector> test = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedRangeMap(),rcpB));
      rcpB.swap(test);
    }

    // we now can guarantee that X and B are blocked multi vectors
    bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
    bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

    // check the type of operator
    RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbA = Teuchos::rcp_dynamic_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(bA);
    if(rbA.is_null() == false) {
      // A is a ReorderedBlockedCrsMatrix
      Teuchos::RCP<const Xpetra::BlockReorderManager > brm = rbA->getBlockReorderManager();

      // check type of X vector
      if(bX->getBlockedMap()->getNumMaps() != bA->getDomainMapExtractor()->NumMaps()) {
        // X is a blocked multi vector but incompatible to the reordered blocked operator A
        Teuchos::RCP<MultiVector> test =
            buildReorderedBlockedMultiVector(brm, bX);
        rcpX.swap(test);
      }
      if(bB->getBlockedMap()->getNumMaps() != bA->getRangeMapExtractor()->NumMaps()) {
        // B is a blocked multi vector but incompatible to the reordered blocked operator A
        Teuchos::RCP<const MultiVector> test =
            buildReorderedBlockedMultiVector(brm, bB);
        rcpB.swap(test);
      }
    }


    // check the type of operator
    /*RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    MUELU_TEST_FOR_EXCEPTION(bA.is_null() == true, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): A_ must be a BlockedCrsMatrix");
    RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbA = Teuchos::rcp_dynamic_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(bA);
    if(rbA.is_null() == false) {
      // A is a ReorderedBlockedCrsMatrix
      Teuchos::RCP<const Xpetra::BlockReorderManager > brm = rbA->getBlockReorderManager();

      // check type of vectors
      RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
      RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

      // check type of X vector
      if(bX.is_null() == false && bX->getBlockedMap()->getNumMaps() != bA->getDomainMapExtractor()->NumMaps()) {
        RCP<ReorderedBlockedMultiVector> rbX = Teuchos::rcp_dynamic_cast<ReorderedBlockedMultiVector>(bX);
        if(rbX.is_null() == true) {
          // X is a blocked multi vector but not reordered
          // However, A is a reordered blocked operator
          // We have to make sure, that A and X use compatible maps
          Teuchos::RCP<MultiVector> test =
              buildReorderedBlockedMultiVector(brm, bX);
          rcpX.swap(test);
        }
      }
      if(bB.is_null() == false && bB->getBlockedMap()->getNumMaps() != bA->getRangeMapExtractor()->NumMaps()) {
        RCP<const ReorderedBlockedMultiVector> rbB = Teuchos::rcp_dynamic_cast<const ReorderedBlockedMultiVector>(bB);
        if(rbB.is_null() == true) {
          // B is a blocked multi vector but not reordered
          // However, A is a reordered blocked operator
          // We have to make sure, that A and X use compatible maps
          Teuchos::RCP<const MultiVector> test =
              buildReorderedBlockedMultiVector(brm, bB);
          rcpB.swap(test);
        }
      }
    }*/

    // Throughout the rest of the algorithm rcpX and rcpB are used for solution vector and RHS

    RCP<MultiVector> residual = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());
    RCP<MultiVector> tempres = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());

    // extract parameters from internal parameter list
    const ParameterList & pL = Factory::GetParameterList();
    LocalOrdinal nSweeps = pL.get<LocalOrdinal>("Sweeps");
    Scalar omega = pL.get<Scalar>("Damping factor");

    // outer Richardson loop
    for (LocalOrdinal run = 0; run < nSweeps; ++run) {
      // one BGS sweep
      // loop over all block rows
      for(size_t i = 0; i<Inverse_.size(); i++) {

        // calculate block residual r = B-A*X
        // note: A_ is the full blocked operator
        residual->update(1.0,*rcpB,0.0); // r = B
        if(InitialGuessIsZero == false || i > 0 || run > 0)
          bA->bgs_apply(*rcpX, *residual, i, Teuchos::NO_TRANS, -1.0, 1.0);
          //A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -1.0, 1.0);

        // extract corresponding subvectors from X and residual
        bool bRangeThyraMode =  rangeMapExtractor_->getThyraMode();
        bool bDomainThyraMode = domainMapExtractor_->getThyraMode();
        Teuchos::RCP<MultiVector> Xi = domainMapExtractor_->ExtractVector(rcpX, i, bDomainThyraMode);
        Teuchos::RCP<MultiVector> ri = rangeMapExtractor_->ExtractVector(residual, i, bRangeThyraMode);
        Teuchos::RCP<MultiVector> tXi = domainMapExtractor_->getVector(i, X.getNumVectors(), bDomainThyraMode);

        // apply solver/smoother
        Inverse_.at(i)->Apply(*tXi, *ri, false);

        // update vector
        Xi->update(omega,*tXi,1.0);  // X_{i+1} = X_i + omega \Delta X_i

        // update corresponding part of rhs and lhs
        domainMapExtractor_->InsertVector(Xi, i, rcpX, bDomainThyraMode); // TODO wrong! fix me
      }
    }

    if (bCopyResultX == true) {
      RCP<MultiVector> Xmerged = bX->Merge();
      X.update(one, *Xmerged, zero);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    return rcp( new BlockedGaussSeidelSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    // extract parameters from internal parameter list
    const ParameterList & pL = Factory::GetParameterList();
    LocalOrdinal nSweeps = pL.get<LocalOrdinal>("Sweeps");
    Scalar omega = pL.get<Scalar>("Damping factor");

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << " Sweeps: " << nSweeps << " damping: " << omega << std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t BlockedGaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
    // FIXME: This is a placeholder
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }
  


} // namespace MueLu



#endif /* MUELU_BLOCKEDGAUSSSEIDELSMOOTHER_DEF_HPP_ */
