#ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_EPETRA_PRECONDITIONER_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <MueLu.hpp>

#include <MueLu_EpetraOperator.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyHelpers.hpp>

//! @file MueLu_CreateEpetraPreconditioner.hpp

namespace MueLu {

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
                             // FIXME: why is it non-const
                             Teuchos::ParameterList& paramListIn,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null)
  {
    typedef double                                                              SC;
    typedef int                                                                 LO;
    typedef int                                                                 GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType                         NO;

    using   Teuchos::ParameterList;

    typedef Xpetra::MultiVector<SC, LO, GO, NO>     MultiVector;
    typedef Xpetra::Matrix<SC, LO, GO, NO>          Matrix;
    typedef Hierarchy<SC,LO,GO,NO>                  Hierarchy;
    typedef HierarchyManager<SC,LO,GO,NO>           HierarchyManager;

    bool hasParamList = paramListIn.numParams();

    RCP<HierarchyManager> mueLuFactory;
    ParameterList paramList = paramListIn;

    std::string syntaxStr = "parameterlist: syntax";
    if (hasParamList && paramList.isParameter(syntaxStr) && paramList.get<std::string>(syntaxStr) == "ml") {
      paramList.remove(syntaxStr);
      mueLuFactory = rcp(new MLParameterListInterpreter<SC,LO,GO,NO>(paramList));

    } else {
      mueLuFactory = rcp(new ParameterListInterpreter  <SC,LO,GO,NO>(paramList));
    }

    RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();
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
      int nPDE = MasterList::getDefault<int>("number of equations");
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

    Teuchos::ParameterList nonSerialList,dummyList;
    ExtractNonSerializableData(paramList, dummyList, nonSerialList);    
    HierarchyUtils<SC,LO,GO,NO>::AddNonSerializableDataToHierarchy(*mueLuFactory,*H, nonSerialList);

    mueLuFactory->SetupHierarchy(*H);

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
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & A,
                             const std::string& xmlFileName,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null)
  {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *Xpetra::toXpetra(A->Comm()));

    return CreateEpetraPreconditioner(A, paramList, inCoords, inNullspace);
  }

  void ReuseEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& inA, MueLu::EpetraOperator& Op) {
    typedef double                                                              SC;
    typedef int                                                                 LO;
    typedef int                                                                 GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType                         NO;

    typedef Xpetra::Matrix<SC,LO,GO,NO>     Matrix;
    typedef Xpetra::Operator<SC,LO,GO,NO>   Operator;
    typedef MueLu ::Hierarchy<SC,LO,GO,NO>  Hierarchy;

    RCP<Hierarchy> H = Op.GetHierarchy();

    TEUCHOS_TEST_FOR_EXCEPTION(!H->GetNumLevels(), Exceptions::RuntimeError,
                               "ReuseTpetraPreconditioner: Hierarchy has no levels in it");
    TEUCHOS_TEST_FOR_EXCEPTION(!H->GetLevel(0)->IsAvailable("A"), Exceptions::RuntimeError,
                               "ReuseTpetraPreconditioner: Hierarchy has no fine level operator");
    RCP<Level> level0 = H->GetLevel(0);

    RCP<Operator> O0 = level0->Get<RCP<Operator> >("A");
    RCP<Matrix>   A0 = Teuchos::rcp_dynamic_cast<Matrix>(O0);

    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC,LO,GO,NO>(inA);
    if (!A0.is_null()) {
      // If a user provided a "number of equations" argument in a parameter list
      // during the initial setup, we must honor that settings and reuse it for
      // all consequent setups.
      A->SetFixedBlockSize(A0->GetFixedBlockSize());
    }
    level0->Set("A", A);
    H->SetupRe();
  }

} //namespace

#endif //ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
