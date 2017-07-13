#ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
#define MUELU_CREATE_EPETRA_PRECONDITIONER_HPP

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>

#include <Teuchos_RCP.hpp>

#include <MueLu.hpp>

#include <MueLu_EpetraOperator.hpp>

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
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null);

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
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null);

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
                             const Teuchos::RCP<Epetra_MultiVector>& inNullspace = Teuchos::null);

  void ReuseEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& inA, MueLu::EpetraOperator& Op);

} //namespace
#endif // HAVE_MUELU_SERIAL and HAVE_MUELU_EPETRA

#endif //ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
