// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_NOXSOLVER_HPP
#define PIRO_NOXSOLVER_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Piro_ObserverBase.hpp"

#include "NOX.H"
#include "NOX_Thyra.H"

#include "Teuchos_ParameterList.hpp"

namespace Piro {

/** \brief Thyra-based Model Evaluator for NOX solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class NOXSolver
    : public SteadyStateSolver<Scalar>
{
  public:

  /** \name Constructors/initializers */
  //@{
  /** \brief . */
  NOXSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams,
            const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
            const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel = Teuchos::null,
            const Teuchos::RCP<ObserverBase<Scalar> > &observer = Teuchos::null);
  //@}

  void reset(){ if(Teuchos::nonnull(solver)) solver->resetSolver(); }

  Teuchos::RCP<Thyra::NOXNonlinearSolver> getSolver() {return solver;}

  Teuchos::RCP<ObserverBase<Scalar> > getObserver() {return observer;}

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > getSubModel() {return model;}
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > getAdjointSubModel() {return adjointModel;}

  private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  Teuchos::RCP<Teuchos::ParameterList> appParams;
  Teuchos::RCP<ObserverBase<Scalar> > observer;

  Teuchos::RCP<Thyra::NOXNonlinearSolver> solver;

  Teuchos::RCP<Teuchos::FancyOStream> out;

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > adjointModel; 

  /** \brief to write only when the converged solutions or all solutions. */
  bool writeOnlyConvergedSol;

  /** \brief Whether to throw an exception when solve fails. */
  bool exitUponFailedNOXSolve; 

  /** \brief Whether to Solve again the system with zero intial guess after failure. */
  bool reComputeWithZeroInitialGuess;

  mutable bool solveState;

  //Store current iteration of Analysis solver
  mutable int current_iteration;

};

}

#endif /*PIRO_NOXSOLVER_HPP*/
