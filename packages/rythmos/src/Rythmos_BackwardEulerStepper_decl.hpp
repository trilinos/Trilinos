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

#ifndef Rythmos_BACKWARD_EULER_STEPPER_DECL_H
#define Rythmos_BACKWARD_EULER_STEPPER_DECL_H

#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_InterpolatorAcceptingObjectBase.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Rythmos {


/** \brief Simple concrete stepper subclass implementing an implicit backward
 * Euler method.
 *
 * This class exists primarily as a simple example of an implicit time stepper
 * and as a vehicle for experimentation.  The <tt>ImplicitBDFStepper</tt> also
 * implements backward Euler and is a more powerful stepper class.  This class
 * does not implement a local truncation error test and therefore also does
 * not handle the automatic step size selection.  Therefore, if you need these
 * features, you should really use the <tt>ImplicitBDFStepper</tt> class.
 */
template<class Scalar>
class BackwardEulerStepper : 
  virtual public SolverAcceptingStepperBase<Scalar>,
  virtual public InterpolatorAcceptingObjectBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \name Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  BackwardEulerStepper();
  
  /** \brief . */
  BackwardEulerStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  //@}
  
  /** \name Overridden from InterpolatorAcceptingObjectBase */
  //@{
  
  /** \brief . */
  void setInterpolator(const RCP<InterpolatorBase<Scalar> >& interpolator);

  /** \brief . */
  RCP<InterpolatorBase<Scalar> >
    getNonconstInterpolator();

  /** \brief . */
  RCP<const InterpolatorBase<Scalar> >
    getInterpolator() const;
  
  /** \brief . */
  RCP<InterpolatorBase<Scalar> > unSetInterpolator();

  //@}

  /** \name Overridden from SolverAcceptingStepperBase */
  //@{

  /** \brief . */
  void setSolver(
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  /** \brief . */
  RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver();

  /** \brief . */
  RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{
 
  /** \brief Returns true. */
  bool supportsCloning() const;

  /** \brief Creates copies of all internal data (including the parameter
   * list) except the model which is assumed to stateless.
   *
   * If a shallow copy of the model is not appropirate for some reasone, then
   * the client can simply reset the model using
   * <tt>returnVal->setModel()</tt>.
   */
  RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

  /** \brief . */
  bool isImplicit() const;

  /** \brief . */
  void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> > &model);
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getModel() const;

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Scalar takeStep(Scalar dt, StepSizeType flag);
  
  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;
  
  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );
  
  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
  
  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;
  
  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;
  
  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;
 
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  void describe(
    Teuchos::FancyOStream  &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  // ///////////////////////
  // Private date members

  bool isInitialized_;
  bool haveInitialCondition_;
  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  RCP<Thyra::VectorBase<Scalar> > scaled_x_old_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_old_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  RCP<Thyra::VectorBase<Scalar> > x_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_;
  Scalar t_;
  Scalar t_old_;

  Scalar dt_;
  int numSteps_;

  RCP<Rythmos::SingleResidualModelEvaluator<Scalar> >  neModel_;

  RCP<Teuchos::ParameterList> parameterList_;

  RCP<InterpolatorBase<Scalar> > interpolator_;


  // //////////////////////////
  // Private member functions

  void defaultInitializeAll_();
  void initialize();

};

/** \brief Nonmember constructor.
 *
 * \relates BackwardEulerStepper
 */
template<class Scalar>
RCP<BackwardEulerStepper<Scalar> >
backwardEulerStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

template<class Scalar>
RCP<BackwardEulerStepper<Scalar> >
backwardEulerStepper();


} // namespace Rythmos

#endif //Rythmos_BACKWARD_EULER_STEPPER_DECL_H
