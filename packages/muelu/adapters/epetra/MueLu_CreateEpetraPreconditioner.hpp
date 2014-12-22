#ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_EPETRA_PRECONDITIONER_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

//! @file MueLu_CreateEpetraPreconditioner.hpp

namespace MueLu {

  /*! \fn EpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Epetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  EpetraCrs_To_XpetraMatrix(const Teuchos::RCP<Epetra_CrsMatrix>& A) {
    typedef Xpetra::EpetraCrsMatrix                                            XECrsMatrix;
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       XCrsMatrix;
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>   XCrsMatrixWrap;

    RCP<XCrsMatrix> Atmp = rcp(new XECrsMatrix(A));
    return rcp(new XCrsMatrixWrap(Atmp));
  }

  /*! \fn EpetraMultiVector_To_XpetraMultiVector
    @brief Helper function to convert a Epetra::MultiVector to an Xpetra::MultiVector
    TODO move this function to an Xpetra utility file
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  EpetraMultiVector_To_XpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector>& V) {
    return rcp(new Xpetra::EpetraMultiVector(V));
  }

  /*! \fn CreateEpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.

    Given a Epetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] paramList Parameter list
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>&   inA,
                             Teuchos::ParameterList& paramList,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null)
  {
    typedef double                                                              SC;
    typedef int                                                                 LO;
    typedef int                                                                 GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType                         NO;

    typedef Xpetra::MultiVector<SC, LO, GO, NO>      MultiVector;
    typedef Xpetra::Matrix<SC, LO, GO, NO>           Matrix;
    typedef Hierarchy<SC,LO,GO,NO>                   Hierarchy;
    typedef HierarchyManager<SC,LO,GO,NO>            HierarchyManager;

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
    
    H->setlib(Xpetra::UseEpetra);

    // Wrap A
    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(inA);
    H->GetLevel(0)->Set("A", A);

    // Wrap coordinates if available
    if (inCoords != Teuchos::null) {
      RCP<MultiVector> coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
      H->GetLevel(0)->Set("Coordinates", coordinates);
    }

    // Wrap nullspace if available, otherwise use constants
    RCP<MultiVector> nullspace;
    if (inNullspace != Teuchos::null) {
      nullspace = EpetraMultiVector_To_XpetraMultiVector<SC, LO, GO, NO>(inNullspace);

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

    return rcp(new EpetraOperator(H));
  }

  /*! \fn CreateEpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.

    Given a Epetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & inA,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null) {
    Teuchos::ParameterList paramList;
    return CreateEpetraPreconditioner(inA, paramList, inCoords, inNullspace);
  }

  /*! \fn CreateEpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.

    Given a Epetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & inA,
                             const std::string& xmlFileName,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null)
  {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *Xpetra::toXpetra(inA->Comm()));

    return CreateEpetraPreconditioner(inA, paramList, inCoords, inNullspace);
  }

} //namespace

#endif //ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
