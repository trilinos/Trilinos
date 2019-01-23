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
#ifndef MUELU_RAPFACTORY_DEF_HPP
#define MUELU_RAPFACTORY_DEF_HPP


#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RAPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_RAPFactory_decl.hpp"
//#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RAPFactory()
    : hasDeclaredInput_(false) { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("transpose: use implicit");
    SET_VALID_ENTRY("rap: triple product");
    SET_VALID_ENTRY("rap: fix zero diagonals");
#undef  SET_VALID_ENTRY
    validParamList->set< RCP<const FactoryBase> >("A",                   null, "Generating factory of the matrix A used during the prolongator smoothing process");
    validParamList->set< RCP<const FactoryBase> >("P",                   null, "Prolongator factory");
    validParamList->set< RCP<const FactoryBase> >("R",                   null, "Restrictor factory");

    validParamList->set< bool >                  ("CheckMainDiagonal",  false, "Check main diagonal for zeros");
    validParamList->set< bool >                  ("RepairMainDiagonal", false, "Repair zeros on main diagonal");

    // Make sure we don't recursively validate options for the matrixmatrix kernels
    ParameterList norecurse;
    norecurse.disableRecursiveValidation();
    validParamList->set<ParameterList> ("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const Teuchos::ParameterList& pL = GetParameterList();
    if (pL.get<bool>("transpose: use implicit") == false)
      Input(coarseLevel, "R");

    Input(fineLevel,   "A");
    Input(coarseLevel, "P");

    // call DeclareInput of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
      (*it)->CallDeclareInput(coarseLevel);

    hasDeclaredInput_ = true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    const bool doTranspose       = true;
    const bool doFillComplete    = true;
    const bool doOptimizeStorage = true;
    {
      FactoryMonitor m(*this, "Computing Ac", coarseLevel);
      std::ostringstream levelstr;
      levelstr << coarseLevel.GetLevelID();
      std::string labelstr = FormattingHelper::getColonLabel(coarseLevel.getObjectLabel());

      TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
        "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

      const Teuchos::ParameterList& pL = GetParameterList();
      RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel,   "A");
      RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P"), AP, Ac;


      if (pL.get<bool>("rap: triple product") == false) {
        // Reuse pattern if available (multiple solve)
        RCP<ParameterList> APparams;
        if(pL.isSublist("matrixmatrix: kernel params"))
          APparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
        else
          APparams = rcp(new ParameterList);

        // By default, we don't need global constants for A*P
        APparams->set("compute global constants: temporaries",APparams->get("compute global constants: temporaries",false));
        APparams->set("compute global constants",APparams->get("compute global constants",false));

        if (coarseLevel.IsAvailable("AP reuse data", this)) {
          GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous AP data" << std::endl;

          APparams = coarseLevel.Get< RCP<ParameterList> >("AP reuse data", this);

          if (APparams->isParameter("graph"))
            AP = APparams->get< RCP<Matrix> >("graph");
        }

        {
          SubFactoryMonitor subM(*this, "MxM: A x P", coarseLevel);

          AP = MatrixMatrix::Multiply(*A, !doTranspose, *P, !doTranspose, AP, GetOStream(Statistics2),
                                      doFillComplete, doOptimizeStorage, labelstr+std::string("MueLu::A*P-")+levelstr.str(), APparams);
        }

        // Reuse coarse matrix memory if available (multiple solve)
        RCP<ParameterList> RAPparams;
        if(pL.isSublist("matrixmatrix: kernel params"))
          RAPparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
        else
          RAPparams = rcp(new ParameterList);



        if (coarseLevel.IsAvailable("RAP reuse data", this)) {
          GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous RAP data" << std::endl;

          RAPparams = coarseLevel.Get< RCP<ParameterList> >("RAP reuse data", this);

          if (RAPparams->isParameter("graph"))
            Ac = RAPparams->get< RCP<Matrix> >("graph");

          // Some eigenvalue may have been cached with the matrix in the previous run.
          // As the matrix values will be updated, we need to reset the eigenvalue.
          Ac->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
        }

        // We *always* need global constants for the RAP, but not for the temps
        RAPparams->set("compute global constants: temporaries",RAPparams->get("compute global constants: temporaries",false));
        RAPparams->set("compute global constants",true);

        // Allow optimization of storage.
        // This is necessary for new faster Epetra MM kernels.
        // Seems to work with matrix modifications to repair diagonal entries.

        if (pL.get<bool>("transpose: use implicit") == true) {
          SubFactoryMonitor m2(*this, "MxM: P' x (AP) (implicit)", coarseLevel);

          Ac = MatrixMatrix::Multiply(*P,  doTranspose, *AP, !doTranspose, Ac, GetOStream(Statistics2),
                                      doFillComplete, doOptimizeStorage, labelstr+std::string("MueLu::R*(AP)-implicit-")+levelstr.str(), RAPparams);

        } else {
          RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");

          SubFactoryMonitor m2(*this, "MxM: R x (AP) (explicit)", coarseLevel);

          Ac = MatrixMatrix::Multiply(*R, !doTranspose, *AP, !doTranspose, Ac, GetOStream(Statistics2),
                                      doFillComplete, doOptimizeStorage, labelstr+std::string("MueLu::R*(AP)-explicit-")+levelstr.str(), RAPparams);
        }
        bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
        bool checkAc             = pL.get<bool>("CheckMainDiagonal")|| pL.get<bool>("rap: fix zero diagonals"); ;
        if (checkAc || repairZeroDiagonals)
          Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, GetOStream(Warnings1));

        if (IsPrint(Statistics2)) {
          RCP<ParameterList> params = rcp(new ParameterList());;
          params->set("printLoadBalancingInfo", true);
          params->set("printCommInfo",          true);
          GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);
        }

        Set(coarseLevel, "A",         Ac);

        APparams->set("graph", AP);
        Set(coarseLevel, "AP reuse data",  APparams);
        RAPparams->set("graph", Ac);
        Set(coarseLevel, "RAP reuse data", RAPparams);
      } else {
        RCP<ParameterList> RAPparams;
        if(pL.isSublist("matrixmatrix: kernel params"))
          RAPparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
        else
          RAPparams = rcp(new ParameterList);

        // We *always* need global constants for the RAP, but not for the temps
        RAPparams->set("compute global constants: temporaries",RAPparams->get("compute global constants: temporaries",false));
        RAPparams->set("compute global constants",true);

        if (pL.get<bool>("transpose: use implicit") == true) {

          Ac = MatrixFactory::Build(P->getDomainMap(), Teuchos::as<LO>(0));

          SubFactoryMonitor m2(*this, "MxMxM: R x A x P (implicit)", coarseLevel);

          Xpetra::TripleMatrixMultiply<SC,LO,GO,NO>::
            MultiplyRAP(*P, doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                        doOptimizeStorage, labelstr+std::string("MueLu::R*A*P-implicit-")+levelstr.str(),
                        RAPparams);

        } else {
          RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");
          Ac = MatrixFactory::Build(R->getRowMap(), Teuchos::as<LO>(0));

          SubFactoryMonitor m2(*this, "MxMxM: R x A x P (explicit)", coarseLevel);

          Xpetra::TripleMatrixMultiply<SC,LO,GO,NO>::
            MultiplyRAP(*R, !doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                        doOptimizeStorage, labelstr+std::string("MueLu::R*A*P-explicit-")+levelstr.str(),
                        RAPparams);
        }
        bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
        bool checkAc             = pL.get<bool>("CheckMainDiagonal")|| pL.get<bool>("rap: fix zero diagonals"); ;
        if (checkAc || repairZeroDiagonals)
          Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, GetOStream(Warnings1));

        if (IsPrint(Statistics2)) {
          RCP<ParameterList> params = rcp(new ParameterList());;
          params->set("printLoadBalancingInfo", true);
          params->set("printCommInfo",          true);
          GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);
        }

        Set(coarseLevel, "A",         Ac);

        // RAPparams->set("graph", Ac);
        // Set(coarseLevel, "RAP reuse data", RAPparams);
      }


    }

    if (transferFacts_.begin() != transferFacts_.end()) {
      SubFactoryMonitor m(*this, "Projections", coarseLevel);

      // call Build of all user-given transfer factories
      for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
        RCP<const FactoryBase> fac = *it;
        GetOStream(Runtime0) << "RAPFactory: call transfer factory: " << fac->description() << std::endl;
        fac->CallBuild(coarseLevel);
        // Coordinates transfer is marginally different from all other operations
        // because it is *optional*, and not required. For instance, we may need
        // coordinates only on level 4 if we start repartitioning from that level,
        // but we don't need them on level 1,2,3. As our current Hierarchy setup
        // assumes propagation of dependencies only through three levels, this
        // means that we need to rely on other methods to propagate optional data.
        //
        // The method currently used is through RAP transfer factories, which are
        // simply factories which are called at the end of RAP with a single goal:
        // transfer some fine data to coarser level. Because these factories are
        // kind of outside of the mainline factories, they behave different. In
        // particular, we call their Build method explicitly, rather than through
        // Get calls. This difference is significant, as the Get call is smart
        // enough to know when to release all factory dependencies, and Build is
        // dumb. This led to the following CoordinatesTransferFactory sequence:
        // 1. Request level 0
        // 2. Request level 1
        // 3. Request level 0
        // 4. Release level 0
        // 5. Release level 1
        //
        // The problem is missing "6. Release level 0". Because it was missing,
        // we had outstanding request on "Coordinates", "Aggregates" and
        // "CoarseMap" on level 0.
        //
        // This was fixed by explicitly calling Release on transfer factories in
        // RAPFactory. I am still unsure how exactly it works, but now we have
        // clear data requests for all levels.
        coarseLevel.Release(*fac);
      }
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                               "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                               "This is very strange. (Note: you can remove this exception if there's a good reason for)");
    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");
    transferFacts_.push_back(factory);
  }

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif // MUELU_RAPFACTORY_DEF_HPP
