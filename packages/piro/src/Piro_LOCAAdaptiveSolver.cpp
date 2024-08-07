// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_LOCAAdaptiveSolver_Def.hpp"

namespace Piro {

// Explicit template instantiations
// LOCA currently only supports Scalar = double

template class LOCAAdaptiveSolver<double>;

template Teuchos::RCP<LOCAAdaptiveSolver<double> > observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> > &adjointModel,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<Piro::ObserverBase<double> > &observer);

} // namespace Piro

