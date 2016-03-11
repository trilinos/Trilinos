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
    typedef MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node>  Hierarchy;
    typedef MueLu::HierarchyManager<Scalar,LocalOrdinal,GlobalOrdinal,Node>  HierarchyManager;

    Teuchos::ParameterList paramList = inParamList;

    bool hasParamList = paramList.numParams();

    RCP<HierarchyManager> mueLuFactory;

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MueLu::MLParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node>(paramList));
    } else {
      mueLuFactory = rcp(new MueLu::ParameterListInterpreter  <Scalar,LocalOrdinal,GlobalOrdinal,Node>(paramList,op->getDomainMap()->getComm()));
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

      nullspace = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(op->getDomainMap(), nPDE);
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
    MueLu::HierarchyUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::AddNonSerializableDataToHierarchy(*mueLuFactory,*H, nonSerialList);

    mueLuFactory->SetupHierarchy(*H);

    return H;
  }

} //namespace

#endif /* PACKAGES_MUELU_ADAPTERS_XPETRA_MUELU_CREATEXPETRAPRECONDITIONER_HPP_ */
