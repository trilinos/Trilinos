// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
  */
Teuchos::RCP<MueLu::EpetraOperator>
CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& inA,
                           // FIXME: why is it non-const
                           Teuchos::ParameterList& paramListIn);

/*!
  @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
  @ingroup MueLuAdapters
  Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
  @param[in] inA Matrix
  @param[in] xmlFileName XML file containing MueLu options.
  */
Teuchos::RCP<MueLu::EpetraOperator>
CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& A,
                           const std::string& xmlFileName);

/*!
  @brief Helper function to create a MueLu preconditioner that can be used by Epetra.
  @ingroup MueLuAdapters
  Given a Epetra_CrsMatrix, this function returns a constructed MueLu preconditioner.
  @param[in] inA Matrix
  */
Teuchos::RCP<MueLu::EpetraOperator>
CreateEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& A,
                           const std::string& xmlFileName);

void ReuseEpetraPreconditioner(const Teuchos::RCP<Epetra_CrsMatrix>& inA, MueLu::EpetraOperator& Op);

}  // namespace MueLu
#endif  // HAVE_MUELU_SERIAL and HAVE_MUELU_EPETRA

#endif  // ifndef MUELU_CREATE_EPETRA_PRECONDITIONER_HPP
