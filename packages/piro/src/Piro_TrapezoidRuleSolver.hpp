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

#ifndef PIRO_TRAPEZOIDRULESOLVER_H
#define PIRO_TRAPEZOIDRULESOLVER_H

#include <iostream>

#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Piro_ObserverBase.hpp"

#include "Piro_NOXSolver.hpp"
#include "Thyra_AdaptiveSolutionManager.hpp"

namespace Piro {

template <typename Scalar>
class TrapezoidDecorator
    : public Thyra::ModelEvaluatorDelegatorBase<Scalar> {

  public:

  /** \name Constructors/initializers */
  //@{

  TrapezoidDecorator(
                const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model
                );

  //@}

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x_dot() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x_dotdot() const;

   void reportFinalPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint, const bool wasSolved);



  //! Method to give info to compute xDotDot(x), so that the
  // NOX solver can treat the time dep problem as steady
  void injectData(const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_a_, Scalar fdt2_,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_v_, Scalar tdt_,
                  Scalar time_ );

  // Resize internal arrays when mesh is adapted
  void resize(const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_);

  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar> > DMEWSF;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_a;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_v;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_save;
   Scalar fdt2;
   Scalar tdt;
   Scalar time;

   Teuchos::RCP<Teuchos::FancyOStream> out;

};

template <typename Scalar>
class TrapezoidRuleSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar> {


  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  TrapezoidRuleSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams,
                      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
                      const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
                      const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null
                      );

  //@}

/** \name Overridden from Thyra::ModelEvaluatorBase. */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  //@}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase. */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  /** \brief . */
  void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidTrapezoidRuleParameters() const;

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<Piro::TrapezoidDecorator<Scalar> > model;
   Teuchos::RCP<Piro::ObserverBase<Scalar> > observer;
   Teuchos::RCP<Thyra::AdaptiveSolutionManager> solMgr;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;

   Teuchos::RCP<Piro::NOXSolver<Scalar> > noxSolver;

   int num_p;
   int num_g;

   int numTimeSteps;
   Scalar t_init, t_final, delta_t;

};

}

#include "Piro_TrapezoidRuleSolver_Def.hpp"

#endif
