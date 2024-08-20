// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_PERFORMSOLVE_HPP
#define PIRO_EPETRA_PERFORMSOLVE_HPP

#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"

class Epetra_Vector;
class Epetra_MultiVector;

namespace Piro {

namespace Epetra {

//! \name Top-level Epetra solve drivers

//@{
//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to non-<tt>const</tt> objects.
//! \sa The \link Piro::PerformSolve Thyra version\endlink of the driver
//! \ingroup Piro_Epetra_solve_driver_grp
void PerformSolve(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Epetra_Vector> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Epetra_MultiVector> > > &sensitivities);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \sa The \link Piro::PerformSolve Thyra version\endlink of the driver
//! \ingroup Piro_Epetra_solve_driver_grp
void PerformSolve(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Epetra_Vector> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Epetra_MultiVector> > > &sensitivities);
//@}

} // namespace Epetra

} // namespace Piro

#endif /* PIRO_EPETRA_PERFORMSOLVE_HPP */
