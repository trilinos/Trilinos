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

#ifndef PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_
#define PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_

//! @file
//! @brief Various adapters that will create a MueLu preconditioner that is an Xpetra::Matrix.

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include <MueLu.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

#include <stdlib.h>

namespace MueLu {

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Xpetra.
    @ingroup MueLuAdapters
    Given an Xpetra::Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] inParamList Parameter list
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op,
                             const Teuchos::ParameterList& inParamList) {
    using SC = Scalar;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using NO = Node;

    using HierarchyManager = MueLu::HierarchyManager<SC, LO, GO, NO>;
    using HierarchyUtils = MueLu::HierarchyUtils<SC, LO, GO, NO>;
    using Hierarchy = MueLu::Hierarchy<SC, LO, GO, NO>;
    using MLParameterListInterpreter = MLParameterListInterpreter<SC, LO, GO, NO>;
    using ParameterListInterpreter = ParameterListInterpreter<SC, LO, GO, NO>;

    bool hasParamList = inParamList.numParams();

    RCP<HierarchyManager> mueLuFactory;

    // Rip off non-serializable data before validation
    Teuchos::ParameterList nonSerialList,paramList;
    MueLu::ExtractNonSerializableData(inParamList, paramList, nonSerialList);

    std::string label;
    if (hasParamList && paramList.isParameter("hierarchy label")) {
      label = paramList.get<std::string>("hierarchy label");
      paramList.remove("hierarchy label");
    } else
      label = op->getObjectLabel();

    std::string timerName;
    if (label != "")
      timerName = "MueLu setup time (" + label + ")";
    else
      timerName = "MueLu setup time";
    RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer(timerName);
    tm->start();

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MLParameterListInterpreter(paramList));
    } else {
      mueLuFactory = rcp(new ParameterListInterpreter(paramList,op->getDomainMap()->getComm()));
    }

    // Create Hierarchy
    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy(label);
    H->setlib(op->getDomainMap()->lib());

    // Set fine level operator
    H->GetLevel(0)->Set("A", op);
    H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());

    // Stick the non-serializible data on the hierarchy.
    HierarchyUtils::AddNonSerializableDataToHierarchy(*mueLuFactory,*H, nonSerialList);

    mueLuFactory->SetupHierarchy(*H);

    tm->stop();
    tm->incrementNumCalls();

    if (H->GetVerbLevel() & Statistics0) {
      const bool alwaysWriteLocal = true;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = timerName;
      Teuchos::TimeMonitor::summarize(op->getRowMap()->getComm().ptr(), H->GetOStream(Statistics0), alwaysWriteLocal, writeGlobalStats,
                                      writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    tm->reset();

    return H;
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Xpetra.
    @ingroup MueLuAdapters
    Given an Xpetra::Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] xmlFileName std::string
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op,
                             const std::string& xmlFileName) {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *op->getDomainMap()->getComm());
    return CreateXpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op, paramList);
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Xpetra.
    @ingroup MueLuAdapters
    Given an Xpetra::Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op) {
    Teuchos::ParameterList paramList;
    return CreateXpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op, paramList);
  }

  /*!
    @brief Helper function to reuse an existing MueLu preconditioner.
    @ingroup MueLuAdapters

    @param[in] inA Matrix
    @param[in] Op  Existing MueLu preconditioner.
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReuseXpetraPreconditioner(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                                 Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>>& H) {
    std::string label = H->GetLevel(0)->getObjectLabel();

    std::string timerName;
    if (label != "")
      timerName = "MueLu setup time (" + label + ")";
    else
      timerName = "MueLu setup time";
    RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer(timerName);
    tm->start();

    typedef Scalar          SC;
    typedef LocalOrdinal    LO;
    typedef GlobalOrdinal   GO;
    typedef Node            NO;

    typedef Xpetra::Matrix<SC,LO,GO,NO>     Matrix;
    typedef Xpetra::Operator<SC,LO,GO,NO>   Operator;

    TEUCHOS_TEST_FOR_EXCEPTION(!H->GetNumLevels(), Exceptions::RuntimeError,
                               "MueLu::ReuseXpetraPreconditioner: Hierarchy has no levels in it");
    TEUCHOS_TEST_FOR_EXCEPTION(!H->GetLevel(0)->IsAvailable("A"), Exceptions::RuntimeError,
                               "MueLu::ReuseXpetraPreconditioner: Hierarchy has no fine level operator");
    RCP<Level> level0 = H->GetLevel(0);

    RCP<Operator> O0 = level0->Get<RCP<Operator> >("A");
    RCP<Matrix>   A0 = Teuchos::rcp_dynamic_cast<Matrix>(O0);

    if (!A0.is_null()) {
      // If a user provided a "number of equations" argument in a parameter list
      // during the initial setup, we must honor that settings and reuse it for
      // all consequent setups.
      A->SetFixedBlockSize(A0->GetFixedBlockSize());
    }
    level0->Set("A", A);

    H->SetupRe();

    tm->stop();
    tm->incrementNumCalls();

    if (H->GetVerbLevel() & Statistics0) {
      const bool alwaysWriteLocal = true;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = timerName;
      Teuchos::TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), H->GetOStream(Statistics0), alwaysWriteLocal, writeGlobalStats,
                                      writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    tm->reset();
  }



} //namespace

#define XPETRA_CREATEXPETRAPRECONDITIONER_SHORT

#endif /* PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_ */
