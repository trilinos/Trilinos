#ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_CPP
#define MUELU_CREATE_EPETRA_PRECONDITIONER_CPP

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
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>&   inA,
                             // FIXME: why is it non-const
                             Teuchos::ParameterList& paramListIn)
  {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using Teuchos::ParameterList;

    using MultiVector      = Xpetra::MultiVector<SC, LO, GO, NO>;
    using Matrix           = Xpetra::Matrix<SC, LO, GO, NO>;
    using Hierarchy        = Hierarchy<SC,LO,GO,NO>;
    using HierarchyManager = HierarchyManager<SC,LO,GO,NO>;

    Teuchos::ParameterList& userList = paramListIn.sublist("user data");
    if (userList.isParameter("Coordinates")) {
      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType,LO,GO,NO> > coordinates = Teuchos::null;
      try {
        coordinates = EpetraMultiVector_To_XpetraMultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType,LO,GO,NO>(userList.get<RCP<Epetra_MultiVector> >("Coordinates"));
      } catch(Teuchos::Exceptions::InvalidParameterType&) {
        coordinates = userList.get<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType, LO, GO, NO> > >("Coordinates");
      }
      if(Teuchos::nonnull(coordinates)){
        userList.set<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType,LO,GO,NO> > >("Coordinates", coordinates);
      }
    }
    if (userList.isParameter("Nullspace")) {
      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType,LO,GO,NO> > nullspace = Teuchos::null;
      try {
        nullspace = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(userList.get<RCP<Epetra_MultiVector> >("Nullspace"));
      } catch(Teuchos::Exceptions::InvalidParameterType&) {
        nullspace = userList.get<RCP<Xpetra::MultiVector<SC, LO, GO, NO> > >("Nullspace");
      }
      if(Teuchos::nonnull(nullspace)){
        userList.set<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<SC>::coordinateType,LO,GO,NO> > >("Nullspace", nullspace);
      }
    }

    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(inA);
    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A, paramListIn);
    return rcp(new EpetraOperator(H));
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options.
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & A,
                             const std::string& xmlFileName)
  {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *Xpetra::toXpetra(A->Comm()));

    return CreateEpetraPreconditioner(A, paramList);
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix.
    */
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & A)
  {
    Teuchos::ParameterList paramList;
    return CreateEpetraPreconditioner(A, paramList);
  }

  void ReuseEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& inA, MueLu::EpetraOperator& Op) {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using Teuchos::ParameterList;

    using Matrix           = Xpetra::Matrix<SC, LO, GO, NO>;
    using Hierarchy        = Hierarchy<SC,LO,GO,NO>;

    RCP<Hierarchy> H = Op.GetHierarchy();
    RCP<Matrix>    A = EpetraCrs_To_XpetraMatrix<SC,LO,GO,NO>(inA);

    MueLu::ReuseXpetraPreconditioner<SC,LO,GO,NO>(A, H);
  }

#ifdef HAVE_MUELU_DEPRECATED_CODE

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a EpetraCrs_Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] paramListIn Parameter list
    @param[in] inCoords Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace Near nullspace of the matrix.
    */
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>&   inA,
                             // FIXME: why is it non-const
                             Teuchos::ParameterList& paramListIn,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace)
  {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using Teuchos::ParameterList;

    using MultiVector      = Xpetra::MultiVector<SC, LO, GO, NO>;
    using Matrix           = Xpetra::Matrix<SC, LO, GO, NO>;
    using Hierarchy        = Hierarchy<SC,LO,GO,NO>;
    using HierarchyManager = HierarchyManager<SC,LO,GO,NO>;

    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(inA);

    Teuchos::ParameterList& userParamList = paramListIn.sublist("user data");
    RCP<MultiVector> coordinates = Teuchos::null;
    if (Teuchos::nonnull(inCoords)) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }
    RCP<MultiVector> nullspace = Teuchos::null;
    if (Teuchos::nonnull(inNullspace)) {
      nullspace = EpetraMultiVector_To_XpetraMultiVector<SC, LO, GO, NO>(inNullspace);
    }
    if (Teuchos::nonnull(nullspace)) {
      userParamList.set<RCP<MultiVector> >("Nullspace", nullspace);
    }

    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A, paramListIn);
    return rcp(new EpetraOperator(H));
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a EpetraCrs_Matrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] paramListIn Parameter list
    @param[in] inCoords Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    */
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>&   inA,
                             // FIXME: why is it non-const
                             Teuchos::ParameterList& paramListIn,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords)
  {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using Teuchos::ParameterList;

    using MultiVector      = Xpetra::MultiVector<SC, LO, GO, NO>;
    using Matrix           = Xpetra::Matrix<SC, LO, GO, NO>;
    using Hierarchy        = Hierarchy<SC,LO,GO,NO>;
    using HierarchyManager = HierarchyManager<SC,LO,GO,NO>;

    RCP<Matrix> A = EpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(inA);

    Teuchos::ParameterList& userParamList = paramListIn.sublist("user data");
    RCP<MultiVector> coordinates = Teuchos::null;
    if (Teuchos::nonnull(inCoords)) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }

    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner<SC,LO,GO,NO>(A, paramListIn);
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
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & inA,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace) {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using MultiVector = Xpetra::MultiVector<SC, LO, GO, NO>;

    Teuchos::ParameterList paramList;
    Teuchos::ParameterList& userParamList = paramList.sublist("user data");

    RCP<MultiVector> coordinates = Teuchos::null;
    if (Teuchos::nonnull(inCoords)) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }

    RCP<MultiVector> nullspace = Teuchos::null;
    if (Teuchos::nonnull(inNullspace)) {
      nullspace = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inNullspace);
    }
    if (Teuchos::nonnull(nullspace)) {
      userParamList.set<const Teuchos::RCP<MultiVector> >("Nullspace", nullspace);
    }

    return CreateEpetraPreconditioner(inA, paramList);
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & inA,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords) {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using MultiVector = Xpetra::MultiVector<SC, LO, GO, NO>;

    Teuchos::ParameterList paramList;
    Teuchos::ParameterList& userParamList = paramList.sublist("user data");

    RCP<MultiVector> coordinates = Teuchos::null;
    if (Teuchos::nonnull(inCoords)) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }

    return CreateEpetraPreconditioner(inA, paramList);
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace Near nullspace of the matrix.
    */
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & A,
                             const std::string& xmlFileName,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords,
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace)
  {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using MultiVector = Xpetra::MultiVector<SC, LO, GO, NO>;

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *Xpetra::toXpetra(A->Comm()));
    Teuchos::ParameterList& userParamList = paramList.sublist("user data");

    RCP<MultiVector> coordinates = Teuchos::null;
    if (inCoords != Teuchos::null) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }

    RCP<MultiVector> nullspace = Teuchos::null;
    if (inCoords != Teuchos::null) {
      nullspace = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inNullspace);
    }
    if (Teuchos::nonnull(nullspace)) {
      userParamList.set<const Teuchos::RCP<MultiVector> >("Nullspace", nullspace);
    }

    return CreateEpetraPreconditioner(A, paramList);
  }

  /*!
    @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
    @ingroup MueLuAdapters
    Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    */
  MUELU_DEPRECATED
  Teuchos::RCP<MueLu::EpetraOperator>
  CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>  & A,
                             const std::string& xmlFileName,
                             const Teuchos::RCP<Epetra_MultiVector>& inCoords)
  {
    using SC = double;
    using LO = int;
    using GO = int;
    using NO = Xpetra::EpetraNode;

    using MultiVector = Xpetra::MultiVector<SC, LO, GO, NO>;

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *Xpetra::toXpetra(A->Comm()));
    Teuchos::ParameterList& userParamList = paramList.sublist("user data");

    RCP<MultiVector> coordinates = Teuchos::null;
    if (inCoords != Teuchos::null) {
      coordinates = EpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    }
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set<RCP<MultiVector> >("Coordinates", coordinates);
    }

    return CreateEpetraPreconditioner(A, paramList);
  }

#endif // HAVE_MUELU_DEPRECATED_CODE

} //namespace
#endif // HAVE_MUELU_SERIAL and HAVE_MUELU_EPETRA

#endif //ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_CPP
