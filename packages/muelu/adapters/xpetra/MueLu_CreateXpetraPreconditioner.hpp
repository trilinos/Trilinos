/*
 * MueLu_CreateXpetraPreconditioner.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: tawiesn
 */

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
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
  */

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op,
                             const Teuchos::ParameterList& inParamList,
                             Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords = Teuchos::null,
                             Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > nullspace = Teuchos::null) {
    typedef MueLu::HierarchyManager<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyManager;
    typedef MueLu::HierarchyUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyUtils;
    typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> Hierarchy;
    typedef MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> MLParameterListInterpreter;
    typedef MueLu::ParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> ParameterListInterpreter;

    std::string timerName = "MueLu setup time";
    RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer(timerName);
    tm->start();

    bool hasParamList = inParamList.numParams();

    RCP<HierarchyManager> mueLuFactory;

    // Rip off non-serializable data before validation
    Teuchos::ParameterList nonSerialList,paramList;
    MueLu::ExtractNonSerializableData(inParamList, paramList, nonSerialList);

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MLParameterListInterpreter(paramList));
    } else {
      mueLuFactory = rcp(new ParameterListInterpreter(paramList,op->getDomainMap()->getComm()));
    }

    // Create Hierarchy
    std::string label;
    if (hasParamList && paramList.isParameter("hierarchy label")) {
      label = paramList.get<std::string>("hierarchy label");
      paramList.remove("hierarchy label");
    } else
      label = op->getObjectLabel();
    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy(label);
    H->setlib(op->getDomainMap()->lib());

    // Stick the non-serializible data on the hierarchy.
    HierarchyUtils::AddNonSerializableDataToHierarchy(*mueLuFactory,*H, nonSerialList);

    // Set fine level operator
    H->GetLevel(0)->Set("A", op);

    // Set coordinates if available
    if (coords != Teuchos::null)
      H->GetLevel(0)->Set("Coordinates", coords);

    // Set nullspace if available
    if (nullspace != Teuchos::null)
      H->GetLevel(0)->Set("Nullspace", nullspace);

    mueLuFactory->SetupHierarchy(*H);

    tm->stop();
    tm->incrementNumCalls();

    if (H->GetVerbLevel() & Statistics0) {
      const bool alwaysWriteLocal = true;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = timerName;
      Teuchos::TimeMonitor::summarize(op->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                                      writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    tm->reset();

    return H;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op,
                             const Teuchos::ParameterList& inParamList,
                             const Teuchos::ParameterList& dummy) {
    typedef MueLu::HierarchyManager<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyManager;
    typedef MueLu::HierarchyUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node> HierarchyUtils;
    typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> Hierarchy;
    typedef MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> MLParameterListInterpreter;
    typedef MueLu::ParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> ParameterListInterpreter;

    std::string timerName = "MueLu setup time";
    RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer(timerName);
    tm->start();

    bool hasParamList = inParamList.numParams();

    RCP<HierarchyManager> mueLuFactory;

    // Rip off non-serializable data before validation
    Teuchos::ParameterList nonSerialList,paramList;
    MueLu::ExtractNonSerializableData(inParamList, paramList, nonSerialList);

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MLParameterListInterpreter(paramList));
    } else {
      mueLuFactory = rcp(new ParameterListInterpreter(paramList,op->getDomainMap()->getComm()));
    }

    // Create Hierarchy
    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();
    H->setlib(op->getDomainMap()->lib());

    // Stick the non-serializible data on the hierarchy.
    HierarchyUtils::AddNonSerializableDataToHierarchy(*mueLuFactory,*H, nonSerialList);

    // Set fine level operator
    H->GetLevel(0)->Set("A", op);

    mueLuFactory->SetupHierarchy(*H);

    tm->stop();
    tm->incrementNumCalls();

    if (H->GetVerbLevel() & Statistics0) {
      const bool alwaysWriteLocal = true;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = timerName;
      Teuchos::TimeMonitor::summarize(op->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                                      writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    tm->reset();

    return H;
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
    std::string timerName = "MueLu setup time";
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
      Teuchos::TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                                      writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    tm->reset();
  }

} //namespace

#define XPETRA_CREATEXPETRAPRECONDITIONER_SHORT

#endif /* PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_ */
