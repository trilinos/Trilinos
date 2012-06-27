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

#ifndef PIRO_RYTHMOSSOLVER_H
#define PIRO_RYTHMOSSOLVER_H

#include <iostream>

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"


/** \brief Thyra-based Model Evaluator subclass for Charon!
 *
 * This class will support a wide number of different types of abstract
 * problem types that will allow NOX, LOCA, Rythmos, Aristos, and MOOCHO to
 * solve different types of problems with Charon.
 * 
 * ToDo: Finish documentation!
 */

namespace Piro {

template <typename Scalar>
class RythmosSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
                Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model,
                Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer = Teuchos::null
                );

  //@}

  ~RythmosSolver();

  Teuchos::RCP<const Rythmos::IntegratorBase<Scalar> > getRythmosIntegrator() const;

  /** \name Overridden from Thyra::ModelEvaluatorBase . */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int i) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int i) const;

  /** \brief . */
  void evalModelImpl( const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidRythmosParameters() const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model;
   Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer;
  
   int num_p;
   int num_g;

   Teuchos::RCP<Rythmos::StepperBase<Scalar> > fwdStateStepper;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;
   Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > fwdStateIntegrator;
   Teuchos::RCP<Thyra::NonlinearSolverBase<double> > fwdTimeStepSolver;
  //Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> > fwdTimeStepSolver;
   Scalar t_final;

};

}

#include "Piro_RythmosSolver_Def.hpp"
#endif
