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
#ifndef MUELU_RAPSHIFTFACTORY_DEF_HPP
#define MUELU_RAPSHIFTFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RAPShiftFactory_decl.hpp"

#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RAPShiftFactory()
    : implicitTranspose_(false), checkAc_(false), repairZeroDiagonals_(false) { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const { // FIXME make fineLevel const
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
        KP = Utils::Multiply(*K, false, *P, false, KP, GetOStream(Statistics2));
        MP = Utils::Multiply(*M, false, *P, false, MP, GetOStream(Statistics2));
        Set(coarseLevel, "AP Pattern", KP);
      }

      // Optimization storage option. If not modifying matrix later (inserting local values),
      // allow optimization of storage.  This is necessary for new faster Epetra MM kernels.
      bool doOptimizedStorage = !checkAc_;

      RCP<Matrix> Ac, Kc, Mc;

      // Reuse pattern if available (multiple solve)
      //     if (IsAvailable(coarseLevel, "RAP Pattern"))
      // Ac = Get< RCP<Matrix> >(coarseLevel, "RAP Pattern");

      bool doFillComplete=true;
      if (implicitTranspose_) {
        SubFactoryMonitor m2(*this, "MxM: P' x (KP) (implicit)", coarseLevel);
        Kc = Utils::Multiply(*P, true, *KP, false, Kc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
        Mc = Utils::Multiply(*P, true, *MP, false, Mc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
      }
      else {
        RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");
        SubFactoryMonitor m2(*this, "MxM: R x (KP) (explicit)", coarseLevel);
        Kc = Utils::Multiply(*R, false, *KP, false, Kc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
        Mc = Utils::Multiply(*R, false, *MP, false, Mc, GetOStream(Statistics2), doFillComplete, doOptimizedStorage);
      }

      // recombine to get K+shift*M
      int level     = coarseLevel.GetLevelID();
      Scalar shift  = shifts_[level];
      Utils2::TwoMatrixAdd(*Kc, false, (Scalar) 1.0, *Mc, false, shift, Ac, GetOStream(Statistics2));
      Ac->fillComplete();

      if (checkAc_)
        CheckMainDiagonal(Ac);

      RCP<ParameterList> params = rcp(new ParameterList());;
      params->set("printLoadBalancingInfo", true);
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);

      Set(coarseLevel, "A", Ac);
      Set(coarseLevel, "K", Kc);
      Set(coarseLevel, "M", Mc);
      Set(coarseLevel, "RAP Pattern", Ac);
    }

    if (transferFacts_.begin() != transferFacts_.end()) {
      SubFactoryMonitor m(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
        RCP<const FactoryBase> fac = *it;
        GetOStream(Runtime0) << "RAPShiftFactory: call transfer factory: " << fac->description() << std::endl;
        fac->CallBuild(coarseLevel);
        // AP (11/11/13): I am not sure exactly why we need to call Release, but we do need it to get rid
        // of dangling data for CoordinatesTransferFactory
        coarseLevel.Release(*fac);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckMainDiagonal(RCP<Matrix> & Ac) const {
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
      if(repairZeroDiagonals_) GetOStream(Warnings0) << "RAPShiftFactory (WARNING): repaired " << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
      else                     GetOStream(Warnings0) << "RAPShiftFactory (WARNING): found "    << gZeroDiags << " zeros on main diagonal of Ac." << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPShiftFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "MueLu::RAPShiftFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. This is very strange. (Note: you can remove this exception if there's a good reason for)");
    transferFacts_.push_back(factory);
  }

} //namespace MueLu

#define MUELU_RAPSHIFTFACTORY_SHORT
#endif // MUELU_RAPSHIFTFACTORY_DEF_HPP
