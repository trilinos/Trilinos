// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
            const Teuchos::RCP<ObserverBase<Scalar> > &observer = Teuchos::null);
  //@}

  void reset(){ if(Teuchos::nonnull(solver)) solver->resetSolver(); }

  Teuchos::RCP<Thyra::NOXNonlinearSolver> getSolver() {return solver;}

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

  /** \brief to write only when the converged solutions or all solutions. */
  bool writeOnlyConvergedSol;

  /** \brief Whether to throw an exception when solve fails. */
  bool exitUponFailedNOXSolve; 

  /** \brief Derivative layouts for Thyra operator"
   * OP:  Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP
   * COL: Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM
   * ROW: Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM
   * */
  enum DerivativeLayout { OP, COL, ROW };

  //Store current iteration of Analysis solver
  mutable int current_iteration;

};

}

#endif /*PIRO_NOXSOLVER_HPP*/
