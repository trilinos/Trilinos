// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_SOLVERFACTORY_H
#define PIRO_SOLVERFACTORY_H

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_AdaptiveSolutionManager.hpp"

#include "Piro_ObserverBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Piro {

/*! \brief Factory for creating Thyra-based %Piro solvers
 *
 *  Piro::Epetra::SolverFactory is the counterpart for Epetra-based models.
 */
class SolverFactory {
public:
  /*! \brief Create a solved model
   *
   *  The type of %Piro solver to instantiate is determined by the value of the string entry labeled <tt>"Solver Type"</tt>
   *  and located at the top level of parameter list \c piroParams.
   *
   *  Currently, the following solver types are available (each accompanied by the corresponding token value):
   *  - Piro::NOXSolver (<tt>"NOX"</tt>)
   *  - Piro::LOCASolver (<tt>"LOCA"</tt>)
   *  - Piro::VelocityVerletSolver (<tt>"Velocity Verlet"</tt>)
   *  - Piro::TrapezoidRuleSolver (<tt>"Trapezoid Rule"</tt>)
   *  - Piro::TempusSolver (<tt>"Tempus"</tt>)
   *
   *  For Epetra-based models, additional options are available in Piro::Epetra::SolverFactory.
   */
  template <typename Scalar>
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > createSolverAdaptive(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel = Teuchos::null,
      const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr = Teuchos::null,
      const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null);

  template <typename Scalar>
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > createSolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel = Teuchos::null,
      const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null);

};

} // namespace Piro

#include "Piro_SolverFactory_Def.hpp"

#endif /* PIRO_SOLVERFACTORY_H */
