// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_PERFORMANALYSIS_HPP
#define PIRO_EPETRA_PERFORMANALYSIS_HPP

#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

class Epetra_Vector;

namespace Piro {

namespace Epetra {

//! \name Top-level Epetra analysis driver
//@{
//! \brief Performs analysis of an Epetra-based solved model.
//! \details This function can either call one of the \link Piro_Thyra_analysis_driver_grp package-specific drivers\endlink
//!          or perform a \link Piro_Epetra_solve_driver_grp forward solve\endlink.
//! \sa The \link Piro::PerformAnalysis Thyra version\endlink of the driver
//! \ingroup Piro_Epetra_analysis_driver_grp
int PerformAnalysis(
    EpetraExt::ModelEvaluator &piroModel,
    Teuchos::ParameterList &analysisParams,
    Teuchos::RCP<Epetra_Vector> &p);
//@}

} // namespace Epetra

} // namespace Piro

#endif /* PIRO_EPETRA_PERFORMANALYSIS_HPP */
