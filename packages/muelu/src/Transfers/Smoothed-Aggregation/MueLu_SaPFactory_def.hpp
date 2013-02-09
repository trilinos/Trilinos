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
#ifndef MUELU_SAPFACTORY_DEF_HPP
#define MUELU_SAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< Scalar >                ("Damping factor",          4./3, "Smoothed-Aggregation damping factor");
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Tentative prolongator factory");
    // validParamList->set                       ("Diagonal view",      "current", "Diagonal view used during the prolongator smoothing process");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "A");

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    RCP<const FactoryBase> initialPFact = GetFactory("P");
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P", initialPFact.get(), this); // --
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
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
      A = Utils2::Transpose(A, true); // build transpose of A explicitely
    }

    //Build final prolongator
    RCP<Matrix> finalP; // output

    //FIXME Xpetra::Matrix should calculate/stash max eigenvalue
    //FIXME SC lambdaMax = A->GetDinvALambda();

    const ParameterList & pL = GetParameterList();
    Scalar dampingFactor = pL.get<Scalar>("Damping factor");
    if (dampingFactor != Teuchos::ScalarTraits<Scalar>::zero()) {

      //Teuchos::ParameterList matrixList;
      //RCP<Matrix> I = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>("Identity", Get< RCP<Matrix> >(fineLevel, "A")->getRowMap(), matrixList);
      //RCP<Matrix> newPtent = Utils::TwoMatrixMultiply(I, false, Ptent, false);
      //Ptent = newPtent; //I tried a checkout of the original Ptent, and it seems to be gone now (which is good)

      RCP<Matrix> AP;
      {
        SubFactoryMonitor m2(*this, "MxM: A x Ptentative", coarseLevel);
        //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
        //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
        //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
        bool doFillComplete=true;

        bool optimizeStorage=true;

        // FIXME: ADD() need B.getProfileType()==DynamicProfile
        if (A->getRowMap()->lib() == Xpetra::UseTpetra) {
          optimizeStorage=false;
        }

        //FIXME Improved Epetra MM returns error code -1 optimizeStorage==true.
#if !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_EPETRAEXT) || !defined(HAVE_MUELU_ML)
        if (A->getRowMap()->lib() == Xpetra::UseEpetra) {
          optimizeStorage=false;
        }
#endif
        //

        //FIXME but once fixed, reenable the next line.
        //if (A->getRowMap()->lib() == Xpetra::UseTpetra) optimizeStorage=false;
        AP = Utils::Multiply(*A, false, *Ptent, false, doFillComplete, optimizeStorage);
      }

      {
        SubFactoryMonitor m2(*this, "Scaling (A x Ptentative) by D^{-1}", coarseLevel);
        bool doFillComplete=true;
        bool optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(*A);
        Utils::MyOldScaleMatrix(AP, diag, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag
      }

      Scalar lambdaMax;
      {
        SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
        Magnitude stopTol = 1e-4;
        lambdaMax = Utils::PowerMethod(*A, true, (LO) 10, stopTol);
        //Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 50, (Scalar)1e-7, true);
        GetOStream(Statistics1, 0) << "Damping factor = " << dampingFactor/lambdaMax << " (" << dampingFactor << " / " << lambdaMax << ")" << std::endl;
      }

      {
        SubFactoryMonitor m2(*this, "M+M: P = (Ptentative) + (D^{-1} x A x Ptentative)", coarseLevel);

        bool doTranspose=false;
        Utils2::TwoMatrixAdd(Ptent, doTranspose, Teuchos::ScalarTraits<Scalar>::one(), AP, doTranspose, -dampingFactor/lambdaMax, finalP);
      }

      {
        SubFactoryMonitor m2(*this, "FillComplete() of P", coarseLevel);
        finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
      }

    } else {
      finalP = Ptent;
    }

    // Level Set
    if(!restrictionMode_)
      {
        // prolongation factory is in prolongation mode
        Set(coarseLevel, "P", finalP);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) finalP->CreateView("stridedMaps", Ptent);
        ///////////////////////// EXPERIMENTAL
      }
    else
      {
        // prolongation factory is in restriction mode
        RCP<Matrix> R = Utils2::Transpose(finalP, true); // use Utils2 -> specialization for double
        Set(coarseLevel, "R", R);

        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) R->CreateView("stridedMaps", Ptent, true);
        ///////////////////////// EXPERIMENTAL
      }

  } //Build()

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDampingFactor(Scalar dampingFactor) {
    SetParameter("Damping factor", ParameterEntry(dampingFactor)); // revalidate
  }

  // deprecated
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDampingFactor() {
    const ParameterList & pL = GetParameterList();
    return pL.get<Scalar>("Damping factor");
  }

} //namespace MueLu

#endif // MUELU_SAPFACTORY_DEF_HPP

//TODO: restrictionMode_ should use the parameter list.
