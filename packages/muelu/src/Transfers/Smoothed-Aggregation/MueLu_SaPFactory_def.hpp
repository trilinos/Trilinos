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
#ifndef MUELU_SAPFACTORY_DEF_HPP
#define MUELU_SAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("sa: damping factor");
    SET_VALID_ENTRY("sa: calculate eigenvalue estimate");
    SET_VALID_ENTRY("sa: eigenvalue estimate num iterations");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Tentative prolongator factory");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this); // --
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }

    // Level Get
    RCP<Matrix> A     = Get< RCP<Matrix> >(fineLevel, "A");
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    if(restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utils2::Transpose(*A, true); // build transpose of A explicitely
    }

    //Build final prolongator
    RCP<Matrix> finalP; // output

    const ParameterList & pL = GetParameterList();
    Scalar dampingFactor = as<Scalar>(pL.get<double>("sa: damping factor"));
    LO maxEigenIterations = as<LO>(pL.get<int>("sa: eigenvalue estimate num iterations"));
    bool estimateMaxEigen = pL.get<bool>("sa: calculate eigenvalue estimate");
    if (dampingFactor != Teuchos::ScalarTraits<Scalar>::zero()) {

      Scalar lambdaMax;
      {
        SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
        lambdaMax = A->GetMaxEigenvalueEstimate();
        if (lambdaMax == -Teuchos::ScalarTraits<SC>::one() || estimateMaxEigen) {
          GetOStream(Statistics1) << "Calculating max eigenvalue estimate now (max iters = "<< maxEigenIterations << ")" << std::endl;
          Magnitude stopTol = 1e-4;
          lambdaMax = Utils::PowerMethod(*A, true, maxEigenIterations, stopTol);
          A->SetMaxEigenvalueEstimate(lambdaMax);
        } else {
          GetOStream(Statistics1) << "Using cached max eigenvalue estimate" << std::endl;
        }
        GetOStream(Statistics0) << "Prolongator damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;
      }

      {
        SubFactoryMonitor m2(*this, "Fused (I-omega*D^{-1} A)*Ptent", coarseLevel);
        Teuchos::RCP<Vector> invDiag = Utils::GetMatrixDiagonalInverse(*A);

        SC omega = dampingFactor / lambdaMax;

        // finalP = Ptent + (I - \omega D^{-1}A) Ptent
        finalP = Utils::Jacobi(omega, *invDiag, *A, *Ptent, finalP, GetOStream(Statistics2));
      }

    } else {
      finalP = Ptent;
    }

    // Level Set
    if (!restrictionMode_) {
      // prolongation factory is in prolongation mode
      Set(coarseLevel, "P", finalP);

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        finalP->CreateView("stridedMaps", Ptent);

    } else {
      // prolongation factory is in restriction mode
      RCP<Matrix> R = Utils2::Transpose(*finalP, true); // use Utils2 -> specialization for double
      Set(coarseLevel, "R", R);

      // NOTE: EXPERIMENTAL
      if (Ptent->IsView("stridedMaps"))
        R->CreateView("stridedMaps", Ptent, true);
    }

    if (IsPrint(Statistics1)) {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      params->set("printCommInfo",          true);
      GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*finalP, (!restrictionMode_ ? "P" : "R"), params);
    }

  } //Build()

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetDampingFactor(Scalar dampingFactor) {
    SetParameter("sa: damping factor", ParameterEntry(dampingFactor)); // revalidate
  }

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDampingFactor() {
    const ParameterList & pL = GetParameterList();
    return as<Scalar>(pL.get<double>("sa: damping factor"));
  }

} //namespace MueLu

#endif // MUELU_SAPFACTORY_DEF_HPP

//TODO: restrictionMode_ should use the parameter list.
