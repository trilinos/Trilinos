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
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_HierarchyUtils.hpp>

//! @file
//! @brief Various adapters that will create a MueLu preconditioner that is an Epetra_Operator.
#if defined(HAVE_MUELU_EPETRA)
namespace MueLu {

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a EpetraCrs_Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] paramListIn Parameter list
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
    typedef double              SC;
    typedef int                 LO;
    typedef int                 GO;
    typedef Xpetra::EpetraNode  NO;

    using   Teuchos::ParameterList;

    typedef Xpetra::MultiVector<SC, LO, GO, NO>     MultiVector;
    typedef Xpetra::Matrix<SC, LO, GO, NO>          Matrix;
    typedef Hierarchy<SC,LO,GO,NO>                  Hierarchy;
    typedef HierarchyManager<SC,LO,GO,NO>           HierarchyManager;

    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(inA);
    RCP<MultiVector> coordinates = Teuchos::null;
    if (inCoords != Teuchos::null) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    RCP<MultiVector> nullspace = Teuchos::null;
    if (inNullspace != Teuchos::null) {
      nullspace = EpetraMultiVector_To_XpetraMultiVector<SC, LO, GO, NO>(inNullspace);
    }
    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A,paramListIn,coordinates,nullspace);
    return rcp(new EpetraOperator(H));
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
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

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
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
    typedef double             SC;
    typedef int                LO;
    typedef int                GO;
    typedef Xpetra::EpetraNode NO;

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
#endif // HAVE_MUELU_SERIAL and HAVE_MUELU_EPETRA

#endif //ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
