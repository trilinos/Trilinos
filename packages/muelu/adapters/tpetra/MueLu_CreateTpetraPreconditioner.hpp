#ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_TPETRA_PRECONDITIONER_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Utilities.hpp>

//! @file MueLu_CreateTpetraPreconditioner.hpp

namespace MueLu {

/*! \fn TpetraCrs_To_XpetraMatrix
    @brief Helper function to convert a Tpetra::CrsMatrix to an Xpetra::Matrix
    TODO move this function to an Xpetra utility file
*/
template <class SC, class LO, class GO, class NO>
RCP<Xpetra::Matrix<SC, LO, GO, NO> >
TpetraCrs_To_XpetraMatrix(Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > const & Atpetra)
{
  RCP<Xpetra::TpetraCrsMatrix<SC, LO, GO, NO> > Atmp = rcp(new Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>(Atpetra));
  RCP<Xpetra::Matrix <SC, LO, GO, NO> > Axpetra  = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(Atmp));
  return Axpetra;
}

/*! \fn TpetraMultiVector_To_XpetraMultiVector
    @brief Helper function to convert a Tpetra::MultiVector to an Xpetra::MultiVector
    TODO move this function to an Xpetra utility file
*/
template <class SC, class LO, class GO, class NO>
RCP<Xpetra::MultiVector<SC, LO, GO, NO> >
TpetraMultiVector_To_XpetraMultiVector(Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const Vtpetra)
{
    RCP<Xpetra::MultiVector<SC, LO, GO, NO> > Vxpetra= rcp(new Xpetra::TpetraMultiVector<SC, LO, GO, NO>(Vtpetra));
    return Vxpetra;
}

/*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] Ain Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
*/
template <class SC, class LO, class GO, class NO>
Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> >
CreateTpetraPreconditioner(Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > const &Ain,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inCoords = Teuchos::null,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inNullspace = Teuchos::null)
{
  RCP<Hierarchy<SC,LO,GO,NO> > H;
  H = rcp(new Hierarchy<SC,LO,GO,NO>());

  RCP<Xpetra::Matrix <SC, LO, GO, NO> > Amuelu  = TpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(Ain);
  H->GetLevel(0)->Set("A", Amuelu);

  //Wrap coordinates if available.
  if (inCoords != Teuchos::null) {
    RCP<Xpetra::MultiVector<SC, LO, GO, NO> > coordinates = TpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    H->GetLevel(0)->Set("Coordinates", coordinates);
  }

  //Wrap nullspace if available, otherwise use constants.
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > nullspace;
  if (inNullspace != Teuchos::null) {
    nullspace = TpetraMultiVector_To_XpetraMultiVector<SC, LO, GO, NO>(inNullspace);
  } else {
    nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Amuelu->getDomainMap(),1);
    nullspace->putScalar( Teuchos::ScalarTraits<SC>::one() );
  }
  H->GetLevel(0)->Set("Nullspace", nullspace);
  H->Setup();

  RCP<TpetraOperator<SC,LO,GO,NO> > tH = rcp(new TpetraOperator<SC,LO,GO,NO>(H));
  return tH;
}


/*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] Ain Matrix
    @param[in] paramList Parameter list
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
*/
template <class SC, class LO, class GO, class NO>
Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> >
CreateTpetraPreconditioner(Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > const &Ain, Teuchos::ParameterList &paramList,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inCoords = Teuchos::null,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inNullspace = Teuchos::null)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using MueLu::Hierarchy;
  using MueLu::ParameterListInterpreter;
  using MueLu::TpetraOperator;

  RCP<Hierarchy<SC,LO,GO,NO> > H;
  //MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(paramList);
  MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(paramList);
  H = mueLuFactory.CreateHierarchy();

  //Wrap A
  RCP<Xpetra::Matrix <SC, LO, GO, NO> > Amuelu  = TpetraCrs_To_XpetraMatrix<SC, LO, GO, NO>(Ain);
  H->GetLevel(0)->Set("A", Amuelu);

  //Wrap coordinates if available.
  if (inCoords != Teuchos::null) {
    RCP<Xpetra::MultiVector<SC, LO, GO, NO> > coordinates = TpetraMultiVector_To_XpetraMultiVector<SC,LO,GO,NO>(inCoords);
    H->GetLevel(0)->Set("Coordinates", coordinates);
  }

  //Wrap nullspace if available, otherwise use constants.
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > nullspace;
  if (inNullspace != Teuchos::null) {
    nullspace = TpetraMultiVector_To_XpetraMultiVector<SC, LO, GO, NO>(inNullspace);
  } else {
    nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Amuelu->getDomainMap(),1);
    nullspace->putScalar( Teuchos::ScalarTraits<SC>::one() );
  }
  H->GetLevel(0)->Set("Nullspace", nullspace);

  mueLuFactory.SetupHierarchy(*H);

  RCP<TpetraOperator<SC,LO,GO,NO> > tH = rcp(new TpetraOperator<SC,LO,GO,NO>(H));

  return tH;
}

/*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] Ain Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
*/
template <class SC, class LO, class GO, class NO>
Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> >
CreateTpetraPreconditioner(Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > const &Ain, std::string const &xmlFileName,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inCoords = Teuchos::null,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inNullspace = Teuchos::null)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Ain->getComm();
  Teuchos::ParameterList paramList;
  //TODO why are we doing this directly and not using the MueLu::ParameterListInterpreter?
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);
  Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > mueluPrecond = CreateTpetraPreconditioner<SC, LO, GO, NO>(Ain, paramList, inCoords, inNullspace);
  return mueluPrecond;
}

} //namespace

#endif //ifndef MUELU_CREATE_TPETRA_PRECONDITIONER_HPP
