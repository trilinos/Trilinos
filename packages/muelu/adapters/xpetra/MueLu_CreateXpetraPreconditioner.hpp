/*
 * MueLu_CreateXpetraPreconditioner.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: tawiesn
 */

#ifndef PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_
#define PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op, const Teuchos::ParameterList& inParamList, Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > nullspace) {
#include "MueLu_UseShortNames.hpp"
    std::string timerName = "MueLu setup time";
    RCP<Teuchos::Time> tm = Teuchos::TimeMonitor::getNewTimer(timerName);
    tm->start();

    ParameterList paramList = inParamList;

    bool hasParamList = paramList.numParams();

    RCP<HierarchyManager> mueLuFactory;

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MLParameterListInterpreter(paramList));
    } else {
      mueLuFactory = rcp(new ParameterListInterpreter(paramList,op->getDomainMap()->getComm()));
    }

    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();
    H->setlib(op->getDomainMap()->lib());


    // Set fine level operator
    H->GetLevel(0)->Set("A", op);

    // Set coordinates if available
    if (coords != Teuchos::null) {
      H->GetLevel(0)->Set("Coordinates", coords);
    }

    // Wrap nullspace if available, otherwise use constants
    if (nullspace == Teuchos::null) {
      int nPDE = MueLu::MasterList::getDefault<int>("number of equations");
      if (paramList.isSublist("Matrix")) {
        // Factory style parameter list
        const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
        if (operatorList.isParameter("PDE equations"))
          nPDE = operatorList.get<int>("PDE equations");

      } else if (paramList.isParameter("number of equations")) {
        // Easy style parameter list
        nPDE = paramList.get<int>("number of equations");
      }

      nullspace = MultiVectorFactory::Build(op->getDomainMap(), nPDE);
      if (nPDE == 1) {
        nullspace->putScalar(Teuchos::ScalarTraits<Scalar>::one());

      } else {
        for (int i = 0; i < nPDE; i++) {
          Teuchos::ArrayRCP<Scalar> nsData = nullspace->getDataNonConst(i);
          for (int j = 0; j < nsData.size(); j++) {
            GlobalOrdinal GID = op->getDomainMap()->getGlobalElement(j) - op->getDomainMap()->getIndexBase();

            if ((GID-i) % nPDE == 0)
              nsData[j] = Teuchos::ScalarTraits<Scalar>::one();
          }
        }
      }
    }
    H->GetLevel(0)->Set("Nullspace", nullspace);


    Teuchos::ParameterList nonSerialList,dummyList;
    MueLu::ExtractNonSerializableData(paramList, dummyList, nonSerialList);
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

#endif /* PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_ */
