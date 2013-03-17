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
#ifndef MUELU_RAPSHIFTFACTORY_DEF_HPP
#define MUELU_RAPSHIFTFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RAPShiftFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RAPShiftFactory()
    : implicitTranspose_(false), checkAc_(false), repairZeroDiagonals_(false) { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    if (implicitTranspose_ == false) {
      Input(coarseLevel, "R");
    }

    Input(fineLevel,   "K");
    Input(fineLevel,   "M");
    Input(coarseLevel, "P");

    // call DeclareInput of all user-given transfer factories
    for(std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it!=transferFacts_.end(); ++it) {
      (*it)->CallDeclareInput(coarseLevel);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &fineLevel, Level &coarseLevel) const { // FIXME make fineLevel const
    {
      FactoryMonitor m(*this, "Computing Ac", coarseLevel);

      // Inputs: K, M, P
      RCP<Matrix> K = Get< RCP<Matrix> >(fineLevel, "K");
      RCP<Matrix> M = Get< RCP<Matrix> >(fineLevel, "M");
      RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");

      // Build Kc = RKP, Mc = RMP
      RCP<Matrix> KP, MP;

      // Reuse pattern if available (multiple solve)
      if (IsAvailable(coarseLevel, "AP Pattern")) {
        KP = Get< RCP<Matrix> >(coarseLevel, "AP Pattern");
        MP = Get< RCP<Matrix> >(coarseLevel, "AP Pattern");
      }

      {
        SubFactoryMonitor subM(*this, "MxM: K x P", coarseLevel);
        KP = Utils::Multiply(*K, false, *P, false, KP);
        MP = Utils::Multiply(*M, false, *P, false, MP);
        Set(coarseLevel, "AP Pattern", KP);
      }

      // Optimization storage option. If not modifying matrix later (inserting local values),
      // allow optimization of storage.  This is necessary for new faster Epetra MM kernels.
      bool doOptimizedStorage = !checkAc_; 
      
      RCP<Matrix> Ac, Kc, Mc;

      // Reuse pattern if available (multiple solve)
      //     if (IsAvailable(coarseLevel, "RAP Pattern"))
      // Ac = Get< RCP<Matrix> >(coarseLevel, "RAP Pattern");

      if (implicitTranspose_) {
        SubFactoryMonitor m2(*this, "MxM: P' x (KP) (implicit)", coarseLevel);
        Kc = Utils::Multiply(*P, true, *KP, false, Kc, true, doOptimizedStorage);
        Mc = Utils::Multiply(*P, true, *MP, false, Mc, true, doOptimizedStorage);
      }
      else {
        RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");
        SubFactoryMonitor m2(*this, "MxM: R x (KP) (explicit)", coarseLevel);
        Kc = Utils::Multiply(*R, false, *KP, false, Kc, true, doOptimizedStorage);
        Mc = Utils::Multiply(*R, false, *MP, false, Mc, true, doOptimizedStorage);
      }

      // recombine to get K+shift*M
      int level     = coarseLevel.GetLevelID();
      Scalar shift  = shifts_[level];
      Utils2::TwoMatrixAdd(Kc, false, (Scalar) 1.0, Mc, false, shift, Ac);
      Ac->fillComplete();

      if(checkAc_) CheckMainDiagonal(Ac);
      GetOStream(Statistics0, 0) << PrintMatrixInfo(*Ac, "Ac");
      GetOStream(Statistics0, 0) << PrintLoadBalancingInfo(*Ac, "Ac");
      Set(coarseLevel, "A", Ac);
      Set(coarseLevel, "K", Kc);
      Set(coarseLevel, "M", Mc);
      Set(coarseLevel, "RAP Pattern", Ac);
    }

    if (transferFacts_.begin() != transferFacts_.end()) {
      SubFactoryMonitor m(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      for(std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
        GetOStream(Runtime0, 0) << "Ac: call transfer factory " << (*it).get() << ": " << (*it)->description() << std::endl;
        (*it)->CallBuild(coarseLevel);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PrintMatrixInfo(const Matrix & Ac, const std::string & msgTag) {
    std::stringstream ss(std::stringstream::out);
    ss << msgTag
       << " # global rows = "      << Ac.getGlobalNumRows()
       << ", estim. global nnz = " << Ac.getGlobalNumEntries()
       << std::endl;
    return ss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PrintLoadBalancingInfo(const Matrix & Ac, const std::string & msgTag) {
    std::stringstream ss(std::stringstream::out);

    // TODO: provide a option to skip this (to avoid global communication)
      // TODO: skip if nproc == 1

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

    ss << msgTag << " # processes with rows = " << numActiveProcesses << std::endl;
    ss << msgTag << " min # rows per proc = " << minNumRows << ", max # rows per proc = " << maxNumRows << ", avg # rows per proc = " << avgNumRows << std::endl;
    ss << msgTag << " nonzero imbalance = " << imbalance << std::endl;

    return ss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheckMainDiagonal(RCP<Matrix> & Ac) const {
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
      if(repairZeroDiagonals_) GetOStream(Warnings0,0) << "RAPShiftFactory (WARNING): repaired " << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
      else                     GetOStream(Warnings0,0) << "RAPShiftFactory (WARNING): found "    << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "MueLu::RAPShiftFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
    transferFacts_.push_back(factory);
  }

} //namespace MueLu

#define MUELU_RAPSHIFTFACTORY_SHORT
#endif // MUELU_RAPSHIFTFACTORY_DEF_HPP

