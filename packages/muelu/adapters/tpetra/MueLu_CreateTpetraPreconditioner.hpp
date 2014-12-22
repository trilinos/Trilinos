#ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_TPETRA_PRECONDITIONER_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

//! @file MueLu_CreateTpetraPreconditioner.hpp

namespace MueLu {

  /*! \fn TpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Tpetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  TpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Atpetra) {
    typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTCrsMatrix;
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       XCrsMatrix;
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>   XCrsMatrixWrap;

    RCP<XCrsMatrix> Atmp = rcp(new XTCrsMatrix(Atpetra));
    return rcp(new XCrsMatrixWrap(Atmp));
  }

  /*! \fn TpetraMultiVector_To_XpetraMultiVector
    @brief Helper function to convert a Tpetra::MultiVector to an Xpetra::MultiVector
    TODO move this function to an Xpetra utility file
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  TpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Vtpetra) {
    return rcp(new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Vtpetra));
  }

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] paramList Parameter list
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                             Teuchos::ParameterList& paramList,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inNullspace = Teuchos::null)
  {
    typedef Scalar          SC;
    typedef LocalOrdinal    LO;
    typedef GlobalOrdinal   GO;
    typedef Node            NO;

    typedef Xpetra::MultiVector<SC,LO,GO,NO>            MultiVector;
    typedef Xpetra::Matrix<SC,LO,GO,NO>                 Matrix;
    typedef Hierarchy<SC,LO,GO,NO>                      Hierarchy;
    typedef HierarchyManager<SC,LO,GO,NO>               HierarchyManager;

    bool hasParamList = paramList.numParams();

    RCP<HierarchyManager> mueLuFactory;
    RCP<Hierarchy>        H;

    Teuchos::ParameterList serialList, nonSerialList;
    if (hasParamList) {
      MueLu::ExtractNonSerializableData(paramList,serialList,nonSerialList);

      std::string Plist = paramList.get("parameterlist: syntax","muelu");
      if(!Plist.compare("ml")) {
	serialList.remove("parameterlist: syntax");
	mueLuFactory = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(serialList));
      }
      else mueLuFactory = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(serialList));
      H = mueLuFactory->CreateHierarchy();

    } else {
      H = rcp(new Hierarchy());
    }

    H->setlib(Xpetra::UseTpetra);

    // Wrap A
    RCP<Matrix> A = TpetraCrs_To_XpetraMatrix<SC,LO,GO,NO>(inA);
    H->GetLevel(0)->Set("A", A);

    // Wrap coordinates if available
    if (inCoords != Teuchos::null) {
      RCP<MultiVector> coordinates = TpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
      H->GetLevel(0)->Set("Coordinates", coordinates);
    }

    // Wrap nullspace if available, otherwise use constants
    RCP<MultiVector> nullspace;
    if (inNullspace != Teuchos::null) {
      nullspace = TpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inNullspace);

    } else {
      int nPDE = 1;
      if (paramList.isSublist("Matrix")) {
        // Factory style parameter list
        const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
        if (operatorList.isParameter("PDE equations"))
          nPDE = operatorList.get<int>("PDE equations");

      } else if (paramList.isParameter("number of equations")) {
        // Easy style parameter list
        nPDE = paramList.get<int>("number of equations");
      }

      nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(A->getDomainMap(), nPDE);
      if (nPDE == 1) {
        nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

      } else {
        for (int i = 0; i < nPDE; i++) {
          Teuchos::ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
          for (int j = 0; j < nsData.size(); j++) {
            GO GID = A->getDomainMap()->getGlobalElement(j) - A->getDomainMap()->getIndexBase();

            if ((GID-i) % nPDE == 0)
              nsData[j] = Teuchos::ScalarTraits<SC>::one();
          }
        }
      }
    }
    H->GetLevel(0)->Set("Nullspace", nullspace);
    
    if (hasParamList) {
      MueLu::HierarchyUtils<SC,LO,GO,NO>::AddNonSerializableDataToHierarchy(*mueLuFactory,*H,nonSerialList);
      mueLuFactory->SetupHierarchy(*H);
    }
    else
      H->Setup();

    return rcp(new TpetraOperator<SC,LO,GO,NO>(H));
  }

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inNullspace = Teuchos::null) {
    Teuchos::ParameterList paramList;
    return CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(inA, paramList, inCoords, inNullspace);
  }

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  CreateTpetraPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                             const std::string& xmlFileName,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inNullspace = Teuchos::null)
  {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *inA->getComm());

    return CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(inA, paramList, inCoords, inNullspace);
  }

} //namespace

#endif //ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
