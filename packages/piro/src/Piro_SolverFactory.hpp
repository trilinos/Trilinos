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
   *  - Piro::RythmosSolver (<tt>"Rythmos"</tt>)
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
