#ifndef MUELU_PRECONDITIONERINTERFACE_HPP
#define MUELU_PRECONDITIONERINTERFACE_HPP

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

/*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] Ain Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
*/
template <class SC, class LO, class GO, class NO>
Teuchos::RCP<MueLu::TpetraOperator<SC,LO,GO,NO> >
CreateTpetraPreconditioner(Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > const &Ain, std::string const &xmlFileName = "",
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inCoords = Teuchos::null,
                          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > const &inNullspace = Teuchos::null)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using MueLu::Hierarchy;
  using MueLu::ParameterListInterpreter;
  using MueLu::TpetraOperator;

  RCP<Hierarchy<SC,LO,GO,NO> > H;
  RCP<ParameterListInterpreter<SC,LO,GO,NO> > mueLuFactory;
  if (!(xmlFileName == "")) {
    mueLuFactory = rcp(new ParameterListInterpreter<SC,LO,GO,NO>(xmlFileName));
    H = mueLuFactory->CreateHierarchy();
  } else {
    H = rcp(new Hierarchy<SC,LO,GO,NO>());
  }

  H->SetDefaultVerbLevel(MueLu::High);

  //Wrap A
  RCP<Xpetra::CrsMatrix<SC,LO,GO,NO> > Atmp;
  Atmp = rcp(new Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>(Ain));
  RCP<Xpetra::Matrix <SC, LO, GO, NO> > Amuelu  = 
     rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(Atmp));
  H->GetLevel(0)->Set("A", Amuelu);

  //Wrap nullspace if available, otherwise use constants.
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > nullspace;
  if (inNullspace != Teuchos::null) {
    nullspace = rcp(new Xpetra::TpetraMultiVector<SC, LO, GO, NO>(inNullspace));
  } else {
    nullspace = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(Amuelu->getDomainMap(),1);
    nullspace->putScalar( Teuchos::ScalarTraits<SC>::one() );
  }
  H->GetLevel(0)->Set("Nullspace", nullspace);

  //Wrap coordinates if available.
  if (inCoords != Teuchos::null) {
    RCP<Xpetra::MultiVector<SC, LO, GO, NO> > coordinates = rcp(new Xpetra::TpetraMultiVector<SC, LO, GO, NO>(inCoords));
    H->GetLevel(0)->Set("Coordinates", coordinates);
  }

  if (mueLuFactory != Teuchos::null) mueLuFactory->SetupHierarchy(*H);
  else                               H->Setup();
  RCP<TpetraOperator<SC,LO,GO,NO> > tH = rcp(new TpetraOperator<SC,LO,GO,NO>(H));

  return tH;
}


} //namespace

#endif //ifndef MUELU_PRECONDITIONERINTERFACE_HPP
