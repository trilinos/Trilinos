//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef Rythmos_INTEGRATOR_BUILDER_DECL_H
#define Rythmos_INTEGRATOR_BUILDER_DECL_H

// Rythmos classes:
#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Rythmos_ErrWtVecCalcBase.hpp"
#include "Rythmos_InterpolatorBase.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"

// Thyra classes:
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_ModelEvaluator.hpp"

// Teuchos:
#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"

namespace Rythmos {


/** \brief Concrete integrator builder class.
 *
 * ToDo: Finish documentation!
 *
 */
template<class Scalar>
  class IntegratorBuilder : virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief . */
  IntegratorBuilder();

  /** \brief . */
  virtual ~IntegratorBuilder();

  /** \brief Set a new Integrator factory object. */
  void setIntegratorFactory(
    const RCP<const AbstractFactory<IntegratorBase<Scalar> > > &integratorFactory,
    const std::string &integratorFactoryName
    );

  /** \brief Set a new Integration Control Strategy factory object. */
  void setIntegrationControlFactory(
    const RCP<const AbstractFactory<IntegrationControlStrategyBase<Scalar> > > &integrationControlFactory,
    const std::string &integrationControlName
    );

  /** \brief Set the Stepper Builder object. */
  void setStepperBuilder(
    const RCP<StepperBuilder<Scalar> > &stepperBuilder
    );

  /** \brief Get the Stepper Builder object. */
  RCP<StepperBuilder<Scalar> > getStepperBuilder();

  /** \brief Set the RK Butcher Tableau Builder object. */
  void setRKButcherTableauBuilder(
      const RCP<RKButcherTableauBuilder<Scalar> > & rkbtBuilder
      );

  /** \brief Set a new Step Control Strategy factory object. */
  void setStepControlFactory(
    const RCP<const AbstractFactory<StepControlStrategyBase<Scalar> > > &stepControlStrategyFactory,
    const std::string &stepControlName
    );

  /** \brief Set an InterpolationBuffer factory object. */
  void setInterpolationBufferFactory(
    const RCP<const AbstractFactory<InterpolationBufferBase<Scalar> > > &interpolationBufferFactory,
    const std::string &interpolationBufferName
    );

  /** \brief Set an InterpolationBufferAppender factory object. */
  void setInterpolationBufferAppenderFactory(
    const RCP<const AbstractFactory<InterpolationBufferAppenderBase<Scalar> > > &interpolationBufferAppenderFactory,
    const std::string &interpolationBufferAppenderName
    );

  /** \brief Set an ErrWtVecCalc factory object. */
  void setErrWtVecCalcFactory(
    const RCP<const AbstractFactory<ErrWtVecCalcBase<Scalar> > > &errWtVecCalcFactory,
    const std::string &errWtVecCalcFactoryName
    );
  
  /** \brief Set an Interpolator factory object. */
  void setInterpolatorFactory(
    const RCP<const AbstractFactory<InterpolatorBase<Scalar> > > &interpolatorFactory,
    const std::string &interpolatorFactoryName
    );
  
  /** \brief Set a W factory object. */
  void setWFactoryObject(
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &wFactoryObject
    );
  
  /** \brief Create an fully formed integrator ready to go.
   *
   * \param model [in] The model evaluator object.  What initialization state
   * does this object need to be in before calling this function?
   *
   * \param initialCondition [in] Initial condition for the integration.  If
   * the real initial condition is not known, then a zero IC can be passed in
   * here and then be rest with returnVal->setInitialCondition(...) later.
   *
   * \param nlSolver [in] The nonlinear solver that will be set on an implicit
   * stepper object.  If an explicit stepper will be created, then this can be
   * left null.
   *
   */
  RCP<IntegratorBase<Scalar> > create(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver
    ) const;
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  /** \brief. */
  RCP<ParameterList> getNonconstParameterList();

  /** \brief. */
  RCP<ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getParameterList() const;
 
  //@}

private:

  // //////////////////////////////////////
  // Private data members

  RCP<Teuchos::ObjectBuilder<IntegratorBase<Scalar> > > integratorBuilder_;
  RCP<Teuchos::ObjectBuilder<IntegrationControlStrategyBase<Scalar> > > integrationControlBuilder_;
  RCP<StepperBuilder<Scalar> > stepperBuilder_;
  RCP<RKButcherTableauBuilder<Scalar> > rkbtBuilder_;
  RCP<Teuchos::ObjectBuilder<StepControlStrategyBase<Scalar> > > stepControlBuilder_;
  RCP<Teuchos::ObjectBuilder<InterpolationBufferBase<Scalar> > > interpolationBufferBuilder_;
  RCP<Teuchos::ObjectBuilder<InterpolationBufferAppenderBase<Scalar> > > interpolationBufferAppenderBuilder_;
  RCP<Teuchos::ObjectBuilder<ErrWtVecCalcBase<Scalar> > > errWtVecCalcBuilder_;
  RCP<Teuchos::ObjectBuilder<InterpolatorBase<Scalar> > > interpolatorBuilder_;

  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > wFactoryObject_;

  RCP<ParameterList> paramList_;
  mutable RCP<ParameterList> validPL_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults_();

};


/** \brief Nonmember constructor to uninitialized.
 *
 * \relates IntegratorBuilder
 */
template<class Scalar>
RCP<IntegratorBuilder<Scalar> > integratorBuilder();


/** \brief Nonmember constructor given parameter list.
 *
 * \relates IntegratorBuilder
 */
template<class Scalar>
RCP<IntegratorBuilder<Scalar> >
integratorBuilder(const RCP<ParameterList> &paramList);

/** \brief Creates a Forward Sensitivity Integrator
 *   Initial conditions for the sensitivities are all zero.
 *
 * \relates IntegratorBuilder
 */
template<class Scalar>
RCP<IntegratorBase<Scalar> > createForwardSensitivityIntegrator(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const int& p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& model_ic,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver,
    const RCP<ParameterList>& integratorBuilderPL
    );

} // namespace Rythmos

#endif //Rythmos_INTEGRATOR_BUILDER_DECL_H

