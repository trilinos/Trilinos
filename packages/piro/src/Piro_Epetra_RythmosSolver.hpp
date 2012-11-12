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

#ifndef PIRO_EPETRA_RYTHMOSSOLVER_H
#define PIRO_EPETRA_RYTHMOSSOLVER_H

#include "EpetraExt_ModelEvaluator.h"

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_StepperBase.hpp"

#include "Piro_RythmosSolver.hpp"

#include "Teuchos_ParameterList.hpp"

/** \brief . */

namespace Piro {
namespace Epetra {

class RythmosSolver
    : public EpetraExt::ModelEvaluator
{
  public:

  typedef double Scalar;
  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  RythmosSolver(Teuchos::RCP<Teuchos::ParameterList> piroParams,
                Teuchos::RCP<EpetraExt::ModelEvaluator> model,
                Teuchos::RCP<Rythmos::IntegrationObserverBase<double> > observer = Teuchos::null
                );

  //@}

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;

  /** \brief . */
  void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;
  //@}

  /** \name Basic information . */
  //@{
  /** \brief Return the number of sets of auxiliary parameters. */
  int Np() const;
  /** \brief Return the number of sets of auxiliary response functions. */
  int Ng() const;
  //@}

  private:
  /** \name Parameter list validation . */
  //@{
  //! Valid list for old "Rythmos" parameter list style
  Teuchos::RCP<const Teuchos::ParameterList> getValidRythmosParameters() const;
  //! Valid list for new "Rythmos Solver" parameter list style
  Teuchos::RCP<const Teuchos::ParameterList> getValidRythmosSolverParameters() const;
  //@}

   Teuchos::RCP<EpetraExt::ModelEvaluator> model;

   int num_p;
   int num_g;

   typedef ::Piro::RythmosSolver<double> ThyraRythmosSolver;
   Teuchos::RCP<ThyraRythmosSolver> thyraImplementation_;

   // To pass to implementation
   Teuchos::RCP<Rythmos::DefaultIntegrator<double> > fwdStateIntegrator;
   Teuchos::RCP<Rythmos::StepperBase<double> > fwdStateStepper;
   Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<double> > fwdTimeStepSolver;
   double t_final;
   Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > fwdStateModel;
   Teuchos::EVerbosityLevel solnVerbLevel;

   Teuchos::RCP<Teuchos::FancyOStream> out;
};

}
}
#endif
