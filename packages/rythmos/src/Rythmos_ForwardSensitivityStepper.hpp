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

#ifndef RYTHMOS_FORWARD_SENSITIVITY_STEPPER_HPP
#define RYTHMOS_FORWARD_SENSITIVITY_STEPPER_HPP


#include "Rythmos_StepperBase.hpp"
#include "Rythmos_ForwardSensitivityModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_SingleResidualModelEvaluatorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Foward sensitivity stepper concrete subclass.
 *
 * This class provides a very general implemenation of a forward sensitivity
 * stepper.
 *
 * \section Rythmos_ForwardSensitivityStepper_intro_sec Introduction
 *
 * The form of the parameterized state equation is:

 \verbatim

   f(x_dot(t),x(t),p) = 0, over t = [t0,tf]

   x(t0) = x_init(p)

 \endverbatim

 * As shown above, the parameters are assumed to be steady state and can enter
 * through the intial condition and/or through the DAE equation itself.
 *
 * The forward sensitivity equations that are solved along with the state
 * equation, as written in multi-vector form, are:

 \verbatim

   d(f)/d(x_dot)*S_dot + d(f)/d(x)*S + d(f)/d(p) = 0, over t = [t0,tf]

   S(t0) = d(x_init)/d(p)

 \endverbatim

 * where <tt>S</tt> is a multi-vector with <tt>np</tt> columns where each
 * column <tt>S(:,j) = d(x)/d(p_j)</tt> is the sensitivity of <tt>x(t)</tt>
 * with respect to the <tt>p_j</tt> parameter.
 *
 * The forward sensitivity equations are a DAE and must be solved using a time
 * integrator.  Conceptually, the state plus forward sensitivity system can be
 * throught of as a big composite DAE model of the form:

 \verbatim

   f_bar(x_bar_dot(t),x_bar(t)) = 0, over t = [t0,tf]

   x_bar(t0) = x_bar_init

 \endverbatim

 * where

 \verbatim

   x_bar = [ x; s_bar ] 

   s_bar = [ S(:,0); S(:,0); ...; S(:,np-1) ]

 \endverbatim

 * and <tt>f_bar(...)</tt> is the obvious concatenated state and sensitivity
 * system.
 *
 * The <tt>InterpolationBufferBase</tt> interface implemented by this class is
 * defined with respect to the full composite solution vector <tt>x_bar</tt>
 * which is returned as a product vector with two components. The first
 * component is <tt>x</tt>.  The second component is another product vector
 * for the implicit concatenation of the columns of the sensitivities shown as
 * <tt>s_bar</tt> above.  The <tt>s_bar</tt> product vector is really
 * implemented directly as a multi-vector and represents an efficient way to
 * represent the forward sensitivities.  Therefore, through the interpolation
 * buffer interface function <tt>getPoints()</tt> a client can access the
 * state and the sensitivities at any point in the range of the current
 * timestep.
 *
 * Note this class does not implement the generic model supporting functions
 * <tt>setModel()</tt> or <tt>getModel()</tt>.  The general assumption of the
 * <tt>StepperBase</tt> interface is that the space for <tt>x</tt> in the
 * model will be the same as used in the interpolation buffer.  That is not
 * true in this case!
 *
 * \section Rythmos_ForwardSensitivityStepper_details_sec Implementation Details
 *
 * There are a variety of ways that one can go about implementing at state
 * plus forward sensitivity stepper.  Three ways for doing this are described
 * in the report "Design of New DASPK for Sensitivity Analysis" by Shengtai Li
 * and Linda Petzold.  The three ways are the <em>simultaneous corrector</em>,
 * the <tt>staggered direct</em> and the <em>staggered corrector</em> methods.
 *
 * The <em>simultaneous corrector</em> method would be equivalent to forming
 * one big ModelEvaluator for the state and sensitivities where the "state"
 * variables would be <tt>x_bar</tt> describe above and then solving them with
 * a single stepper object and as single nonlinear solver.  The advantage of
 * this approach is that it makes great reuse of all of the timestepping
 * software.  Also, by being able to specialize the nonlinear solver (which
 * you can't do in the Sundials software or in DASPK) you could set up the
 * nonlinear solver to first solve the nonlinear state timestep equation, and
 * then solve the linear sensitivity equations.  The problem with this
 * approach would be that it would be very wasteful if the timestep had to be
 * cut back in order to reduce the local truncation error of the state
 * solution.  This would result in the sensitivity solution being thrown away
 * for each cut-back iteration.  Because of this fundamental problem, we have
 * not implemented the simultaneous corrector method.  Actually, we are not
 * really sure why anyone implements ths method.
 *
 * The <em>staggered direct</em> and <em>staggered corrector</em> methods are
 * similar in several ways.  In each method, the state timestep is first fully
 * solved (including any stepsize reduction iterations that are needed) and
 * then the sensitivities are solved, taking advantage of the pre-computed
 * state timestep Jacobian.  One difference between the two methods is that in
 * the staggered direct method, the sensitivities are solved for directly.
 * This can result in numerical accuracy problems and does not allow the reuse
 * of an inexact solve when a direct factorization is used.  In the staggered
 * corrector method, an implicit corrector solve is used to compute the change
 * in the sensitivity variables using a Newton method.  This results in better
 * numerical stability and allows the reuse of an old Jacobian and
 * factorization in the Newton method.  However, if an exact Jacobian is used
 * for the solve, then the Newton method will converge in one iteration
 * (assuming the linear solver tolerance is made tight enough) with no harm
 * done.
 *
 * Because of the advantages of the staggered corrector method over the other
 * methods, the staggered corrector method is what is implemented in this
 * stepper class.  However, the term "corrector" is not really appropriate to
 * be used in the context of this class since this class does not have to
 * assume anything about how timesteps are computed and does not care if a
 * predictor/corrector method is used or not.
 *
 * This subclass never forms a full ModelEvaluator for the full state plus
 * forward sensitivity DAE <tt>f_bar(x_bar_hat,x_bar)</tt>.  Instead, the step
 * is solved first for the state equation and then a ModelEvaluator for just
 * the forward sensitivity equations is formed and is solved over the same
 * time step as the forward solve.
 *
 * Currently, timestep control is not performed for the forward sensitivity
 * variables.  In the future, however, it would not be too difficult to allow
 * for the timestep to be reduced for the sensitivity variables but this would
 * require a "undoStep()" operation be implimented for the state stepper
 * object and this is not currently supported by the <tt>StepperBase</tt>
 * interface.
 *
 *
 * 2007/15/21: rabart: ToDo: This class only works for implicit models and
 * steppers right now but it would be easy to get this to work for explicit
 * steppers and models later with a little work.  With an examplicit
 * method/model, we don't need to reuse W_tilde!

 * ToDo: Finish documentation!
 */
template<class Scalar> 
class ForwardSensitivityStepper
  : virtual public StepperBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Intializers, Misc. */
  //@{

  /** \brief Constructs to uninitialized. */
  ForwardSensitivityStepper();

  /** \brief Intialize.
   *
   * \param  stateModel
   *           [in,persisting] The ModelEvaluator that defines the
   *           parameterized state model <tt>f(x_dot,x,p)</tt>.
   * \param  p_index
   *           [in] The index of the parameter subvector in <tt>stateModel</tt>
   *           for which sensitivities will be computed for.
   * \param  baseStatePoint
   *           [in] Whatever input arguments are needed to define the state
   *           of the model including the parameters except x, x_dot, and t!
   * \param  stateStepper
   *           [in,persisting] The stepper object that will be used to advance
   *           the state solution <tt>x(t)</tt>.  This stepper need not be
   *           setup with a model or a nonlinear timestep solver.  All this
   *           stepper object needs is to be given its parameters to determine
   *           exactly what timestepping algorithm will be employed.  The
   *           model and the timestep solver objects will be set internally.
   * \param  stateTimeStepSolver
   *           [in,persisting] The nonlinear solver object that is used to 
   *           solve for the state timestep equation.  This is needed to
   *           extract the Jacobian matrix that is used in the sensitivity model.
   *           If the stepper is not an implicit stepper and does not use
   *           an implicit time step solver, then this argument can be left
   *           null.
   * \param  sensStepper
   *           [in,persisting] The stepper object that will be used to advance
   *           the sensitivity solution <tt>S(t)</tt>.  This stepper need not
   *           be setup with a model or a nonlinear timestep solver.  All this
   *           stepper object needs is to be given its parameters to determine
   *           exactly what timestepping algorithm will be employed.  The
   *           model and the timestep solver objects will be set internally.
   *           If this argument is null, then the <tt>stateStepper</tt> object
   *           will be cloned to generate this stepper object.  The most
   *           common use cases should just pass in <tt>Teuchos::null</tt> and
   *           just use the identical stepper as the state stepper.  However,
   *           this argument allows a client to specialize exactly what the
   *           sensitivity stepper does and therefore this hook is allowed.
   * \param  sensTimeStepSolver
   *           [in,persisting] The nonlinear solver object that is used to
   *           solve for the sensitivity timestep equation.  If the stepper is
   *           not an implicit stepper and does not use an implicit time step
   *           solver, then this argument can be left null.  If the stepper is
   *           implicit, and this argument is left null, then this solver
   *           object will be created by cloning the
   *           <tt>stateTimeStepSolver</tt> object.  The most common use cases
   *           should just pass in <tt>Teuchos::null</tt> and just use the
   *           identical nonlinear solver as the state stepper.  However, this
   *           argument allows a client to specialize exactly what the
   *           nonlinear solver in the sensitivity stepper does and therefore
   *           this hook is allowed.
   */
  void initialize(
    const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const Teuchos::RefCountPtr<StepperBase<Scalar> > &stateStepper,
    const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
    const Teuchos::RefCountPtr<StepperBase<Scalar> > &sensStepper = Teuchos::null,
    const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
    );

  /** \brief Return the forward sensitivity model evaluator object that got
   * created internally when <tt>initialize()</tt> was called.
   */
  Teuchos::RefCountPtr<const ForwardSensitivityModelEvaluator<Scalar> >
  getFwdSensModel() const;

  /** \brief Set the initial condition for the forward sensitivities.
   *
   * 2007/05/23: rabartl: ToDo: This should really be handled through the
   * setInitialCondition(...) function.  However, in order for a client to
   * call this function, first an InArgs object must be set up and where would
   * that InArgs object get created from?  Typically it would be created from
   * this->getModel()->createInArgs() but this class currently does not
   * implement such a model.  We could provide a function to just return the
   * sensModel object and we could us it's InArgs object structure to create
   * the InArgs object.  I think the cleanest thing to do is to implement a
   * dummy StateAndForwardSensitivitiesModelEvaluator class that can be used
   * to implement this->getModel() and then this can be used to get the needed
   * InArgs object safely.
   */
  void setFwdSensInitialCondition(
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &s_bar_init
    );

  //@}

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{

  /** \brief Returns false. */
  bool acceptsModel() const;

  /** \brief Throws exception. */
  void setModel(
    const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
    );

  /** \brief Returns the state model ???. */
  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
  getModel() const;

  // RAB: ToDo: 2007/05/15: I need to talk with Todd about potentially
  // removing the setModel() and getModel() functions from the StepperBase
  // interface.  In the case of this forward sensitivity solver, I am not sure
  // that it makes a lot of sense to define a model.  This surely will not be
  // the model that a generic client would expect.  The assumption I am sure
  // would be that this model has the same space for x as the interpolation
  // buffer but that is not true in this case.

  /** \brief Sets the full initial condition for <tt>x_bar</tt> and
   * <tt>x_bar_dot</tt>.
   *
   * The vectors <tt>x_bar</tt> and <tt>x_bar_dot</tt> can be created using
   * the function <tt>this->get_x_space()</tt>.
   */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Scalar takeStep( Scalar dt, StepSizeType stepType );

  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;

  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief Returns the space for <tt>x_bar</tt> and <tt>x_bar_dot</tt>.
   *
   * This space is a nested product vector space as described above.  Dynamic
   * casting is required to get at the <tt>ProductVectorSapceBase</tt> and
   * <tt>ProductVectorBase</tt> intefaces.
   */
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  bool setPoints(
    const std::vector<Scalar>& time_vec,
    const std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >& x_vec,
    const std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >& xdot_vec,
    const std::vector<ScalarMag> & accuracy_vec
    );

  /** \brief . */
  bool setRange(
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar>& IB
    );

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  bool getPoints(
    const std::vector<Scalar>& time_vec,
    std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >* x_vec,
    std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    std::vector<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  bool getNodes(std::vector<Scalar>* time_vec) const;

  /** \brief . */
  bool removeNodes(std::vector<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > stateModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;
  Teuchos::RefCountPtr<StepperBase<Scalar> > stateStepper_;
  Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > stateTimeStepSolver_;
  Teuchos::RefCountPtr<StepperBase<Scalar> > sensStepper_;
  Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > sensTimeStepSolver_;

  bool isSingleResidualStepper_;
  Teuchos::RefCountPtr<ForwardSensitivityModelEvaluator<Scalar> > sensModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_t_;

  Teuchos::RefCountPtr<const Thyra::DefaultProductVectorSpace<Scalar> > x_bar_space_;

  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > s_bar_init_;

  Scalar t_;
  Scalar t_old_;
  Scalar dt_;

  // /////////////////////////
  // Private member functions
  
  // Create a wrapped product vector of the form x_bar = [ x; s_bar ]
  //
  // Note: This does not copy any vector data, it only creates the wrapped
  // product vector.
  Teuchos::RefCountPtr<const Thyra::DefaultProductVector<Scalar> >
  create_x_bar_vec(
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &x_vec,
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &s_bar_vec
    ) const;

};


/** \brief Nonmember constructor.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar> 
inline
Teuchos::RefCountPtr<ForwardSensitivityStepper<Scalar> >
forwardSensitivityStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &stateStepper,
  const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &sensStepper = Teuchos::null,
  const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
  )
{
  Teuchos::RefCountPtr<ForwardSensitivityStepper<Scalar> >
    fwdSensStepper = Teuchos::rcp(new ForwardSensitivityStepper<Scalar>());
  fwdSensStepper->initialize(
    stateModel, p_index, stateBasePoint, stateStepper, stateTimeStepSolver );
  return fwdSensStepper;
}


//
// Implementation
//


// Constructors, Intializers, Misc.


template<class Scalar> 
ForwardSensitivityStepper<Scalar>::ForwardSensitivityStepper()
  :isSingleResidualStepper_(false), t_(0.0), t_old_(0.0), dt_(-1.0)
{}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::initialize(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &stateStepper,
  const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const Teuchos::RefCountPtr<StepperBase<Scalar> > &sensStepper,
  const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver
  )
  
{

  using Teuchos::rcp_dynamic_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::tuple;

  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // Validate input
  //

  TEST_FOR_EXCEPT( is_null(stateModel) );
  TEST_FOR_EXCEPT( is_null(stateStepper) );
  TEST_FOR_EXCEPT( is_null(stateTimeStepSolver) ); // ToDo: allow to be null for explicit methods!


  //
  // Create the sensModel which will do some more validation
  //
  
  Teuchos::RefCountPtr<ForwardSensitivityModelEvaluator<Scalar> >
    sensModel = Teuchos::rcp(new ForwardSensitivityModelEvaluator<Scalar>);
  sensModel->initializeStructure(stateModel,p_index);
  
  //
  // Get the input objects
  //

  stateModel_ = stateModel;

  stateBasePoint_ = stateBasePoint;

  stateStepper_ = stateStepper;
  
  stateTimeStepSolver_ = stateTimeStepSolver;

  sensModel_ = sensModel;

  if (!is_null(sensStepper)) {
    sensStepper_ = sensStepper;
  }
  else {
    sensStepper_ = stateStepper_->cloneStepper();
    TEST_FOR_EXCEPTION(
      is_null(sensStepper_), std::logic_error,
      "Error, if the client does not pass in a stepper for the senitivity\n"
      "equations then the stateStepper object must support cloning to create\n"
      "the sensitivity stepper!"
      );
  }

  if (!is_null(sensTimeStepSolver)) {
    sensTimeStepSolver_ = sensTimeStepSolver;
  }
  else {
    sensTimeStepSolver_ = stateTimeStepSolver_->cloneNonlinearSolver();
    TEST_FOR_EXCEPTION(
      is_null(sensTimeStepSolver_), std::logic_error,
      "Error, if the client does not pass in a time step solver for the senitivity\n"
      "equations then the stateTimeStepSolver object must support cloning to create\n"
      "the needed solver!"
      );
  }

  //
  // Setup the steppers
  //

  isSingleResidualStepper_ = true; // ToDo: Add dynamic cast on
                                   // stateTimeStepSolver to check this!

  stateStepper_->setModel(stateModel_);
  rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(
    stateStepper_)->setSolver(stateTimeStepSolver_);
  sensStepper_->setModel(sensModel_);
  rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(
    sensStepper_)->setSolver(sensTimeStepSolver_);

  stateBasePoint_t_ = stateModel_->createInArgs();

  //
  // Setup the structure of the full vector x_bar = [ x; s_bar ]
  //

  x_bar_space_ = Thyra::productVectorSpace(
    tuple<RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > >(
      stateModel_->get_x_space(), sensModel_->get_x_space()
      )
    );

  //
  // Setup the initial condition
  //

  t_ = stateBasePoint_t_.get_t();

  // 2007/05/18: rabartl: ToDo: Move the above initialization code to give
  // setInitializeCondition(...) a chance to set the initial condition.
  
}

  
template<class Scalar> 
Teuchos::RefCountPtr<const ForwardSensitivityModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getFwdSensModel() const
{
  return sensModel_;
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setFwdSensInitialCondition(
  const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &s_bar_init
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(s_bar_init));
  THYRA_ASSERT_VEC_SPACES(
    "ForwardSensitivityStepper<Scalar>::setFwdSensInitialCondition(s_bar_init)",
    *sensModel_->get_x_space(), *s_bar_init->space() );
#endif
  
  s_bar_init_ = s_bar_init;
  
  MEB::InArgs<Scalar> sensIC = sensModel_->createInArgs();
  sensIC.set_t(stateBasePoint_.get_t());
  sensIC.set_x(s_bar_init_);
  sensStepper_->setInitialCondition(sensIC);
  
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT("Do not accept a parameter list yet!");
}


template<class Scalar> 
Teuchos::RefCountPtr<Teuchos::ParameterList>
ForwardSensitivityStepper<Scalar>::getParameterList()
{
  return Teuchos::null;
}


template<class Scalar> 
Teuchos::RefCountPtr<Teuchos::ParameterList>
ForwardSensitivityStepper<Scalar>::unsetParameterList()
{
  return Teuchos::null;
}


template<class Scalar> 
Teuchos::RefCountPtr<const Teuchos::ParameterList>
ForwardSensitivityStepper<Scalar>::getParameterList() const
{
  return Teuchos::null;
}


template<class Scalar> 
Teuchos::RefCountPtr<const Teuchos::ParameterList>
ForwardSensitivityStepper<Scalar>::getValidParameters() const
{
  return Teuchos::null;
}


// Overridden from StepperBase

template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::acceptsModel() const
{
  return false;
}

template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setModel(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  )
{
  TEST_FOR_EXCEPT("Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
}


template<class Scalar> 
Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getModel() const
{
  TEST_FOR_EXCEPT("Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
  return Teuchos::null;
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  TEST_FOR_EXCEPT(true);
}
 

template<class Scalar> 
Scalar
ForwardSensitivityStepper<Scalar>::takeStep(
  Scalar dt, StepSizeType stepType
  )
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::VerboseObjectTempState<InterpolationBufferBase<Scalar> > VOTSIBB;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nEntering " << Teuchos::TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  //
  // Compute the state timestep
  //

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nTaking state step using stepper : "
      << stateStepper_->description() << "\n";
  }

  VOTSIBB stateStepper_outputTempState(stateStepper_,out,verbLevel);
  
  const Scalar state_dt = stateStepper_->takeStep(dt,stepType);

  if (state_dt < Scalar(-ST::one())) {
    // The state stepper has failed so return a failed timestep!
    return dt_;
  }

  //
  // Setup the sensitivity model for this timestep
  //

  const StepStatus<Scalar> stateStepStatus = stateStepper_->getStepStatus();

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ) {
    *out
      << "\nState step status:\n" << stateStepStatus;
  }
  
  TEST_FOR_EXCEPTION(
    stateStepStatus.stepStatus != STEP_STATUS_CONVERGED, std::logic_error,
    "Error, the status should be converged since a positive step size was returned!"
    );

  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
    x = stateStepStatus.solution;
  
  TEST_FOR_EXCEPTION(
    is_null(x), std::logic_error,
    "Error, the state solution can not be null!"
    );

  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
    x_dot = stateStepStatus.solutionDot;
  
  TEST_FOR_EXCEPTION(
    is_null(x_dot), std::logic_error,
    "Error, the state solution derivative can not be null!"
    );
  
  stateBasePoint_t_ = stateBasePoint_;
  stateBasePoint_t_.set_x_dot( x_dot );
  stateBasePoint_t_.set_x( x );
  stateBasePoint_t_.set_t( t_ + state_dt );

  Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<Scalar> >
    W_tilde = stateTimeStepSolver_->get_W();
  
  TEST_FOR_EXCEPTION(
    !stateTimeStepSolver_->is_W_current() || is_null(W_tilde),
    std::logic_error,
    "Error, the W from the state time step must be current and must be nonnull!"
    );
  
  Teuchos::RefCountPtr<const Rythmos::SingleResidualModelEvaluatorBase<Scalar> >
    singleResidualModel
    = Teuchos::rcp_dynamic_cast<const Rythmos::SingleResidualModelEvaluatorBase<Scalar> >(
      stateTimeStepSolver_->getModel()
      );
  
  const Scalar
    coeff_x_dot = singleResidualModel->get_coeff_x_dot(),
    coeff_x = singleResidualModel->get_coeff_x();

  sensModel_->initializeState(
    stateBasePoint_t_, W_tilde, coeff_x_dot, coeff_x );

  //
  // Compute the sensitivity timestep for the exact same timestep as was used
  // for the state solve.
  //

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nTaking sensitivity step using stepper : "
      << sensStepper_->description() << "\n";
  }

  // 2007/05/18: rabartl: ToDo: Copy the stepper control logic from
  // stateStepper_ to sensStepper_!

  VOTSIBB sensStepper_outputTempState(sensStepper_,out,verbLevel);

  const Scalar sens_dt = sensStepper_->takeStep(state_dt,FIXED_STEP);

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) ) {
    const StepStatus<Scalar> sensStepStatus = sensStepper_->getStepStatus();
    *out
      << "\nSensitivity step status:\n" << sensStepStatus;
  }
  
  TEST_FOR_EXCEPTION(
    sens_dt != state_dt, std::logic_error,
    "Error, the sensitivity step failed for some reason.  We should\n"
    "just return a negative step size and reject the step but currently\n"
    "there is no way to roll back the state timestep it for back to\n"
    "the status before this function was called!"
    );

  // 2007/05/18: rabartl: ToDo: If stepType == VARIABLE_STEP and the state
  // timestep sucessed but the sensitivity timestep failed, then we need to
  // really throw an excpetion because there is nothing that we can really do
  // here!

  //
  // Update the timestep
  //

  t_old_ = t_;
  t_ += state_dt;
  dt_ = state_dt;

  if ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) ) {
    *out
      << "\nLeaving " << Teuchos::TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  return dt_;

  // 2007/05/18: rabartl: ToDo: Replace the above std::logic_error type with
  // a Rythmos::CatastrophicFailure or just use Thyra::CatastrophicFailure!

}


template<class Scalar> 
const StepStatus<Scalar>
ForwardSensitivityStepper<Scalar>::getStepStatus() const
{

  const StepStatus<Scalar>
    stateStepStatus = stateStepper_->getStepStatus();
  const StepStatus<Scalar>
    sensStepStatus = sensStepper_->getStepStatus();

  StepStatus<Scalar> stepStatus;
  
  stepStatus.message = sensStepStatus.message;
  stepStatus.stepStatus = sensStepStatus.stepStatus;
  stepStatus.stepLETStatus = sensStepStatus.stepLETStatus;
  stepStatus.stepSize = sensStepStatus.stepSize;
  stepStatus.order = sensStepStatus.order;
  stepStatus.time = sensStepStatus.time;
  stepStatus.stepLETValue = sensStepStatus.stepLETValue;
  stepStatus.solution = create_x_bar_vec(
    stateStepStatus.solution, sensStepStatus.solution
    );
  stepStatus.solutionDot = create_x_bar_vec(
    stateStepStatus.solutionDot, sensStepStatus.solutionDot
    );
  stepStatus.residual = Teuchos::null;
  stepStatus.extraParameters = sensStepStatus.extraParameters;

  return stepStatus;

}


// Overridden from InterpolationBufferBase


template<class Scalar> 
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityStepper<Scalar>::get_x_space() const
{
  return x_bar_space_;
}


template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::setPoints(
  const std::vector<Scalar>& time_vec,
  const std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >& x_vec,
  const std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >& xdot_vec,
  const std::vector<ScalarMag> & accuracy_vec
  )
{
  TEST_FOR_EXCEPT("Not implemented setPoints(...) yet but we could if we wanted!");
  return false;
}


template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& IB
  )
{
  return false;
}


template<class Scalar> 
TimeRange<Scalar>
ForwardSensitivityStepper<Scalar>::getTimeRange() const
{
  return sensStepper_->getTimeRange();
}


template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::getPoints(
  const std::vector<Scalar>& time_vec,
  std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >* x_bar_vec,
  std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >* x_bar_dot_vec,
  std::vector<ScalarMag>* accuracy_vec
  ) const
{

  using Teuchos::as;
  using Teuchos::RefCountPtr;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( as<int>(time_vec.size()) == 0 );
#endif

  const int numTimePoints = time_vec.size();

  bool result = false;
  
  std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >
    x_vec, x_dot_vec;

  result = stateStepper_->getPoints(
    time_vec,
    x_bar_vec ? &x_vec: 0,
    x_bar_dot_vec ? &x_dot_vec: 0,
    0 // Ignoring accuracy from state for now!
    );
  
  if (!result)
    return false;

  std::vector<Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > >
    s_bar_vec, s_bar_dot_vec;

  result = sensStepper_->getPoints(
    time_vec,
    x_bar_vec ? &s_bar_vec: 0,
    x_bar_dot_vec ? &s_bar_dot_vec: 0,
    accuracy_vec
    );
  
  if (!result)
    return false;

  if ( x_bar_vec ) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      x_bar_vec->push_back(create_x_bar_vec(x_vec[i],s_bar_vec[i]));
    }
  }
  
  if ( x_bar_dot_vec ) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      x_bar_dot_vec->push_back(create_x_bar_vec(x_dot_vec[i],s_bar_dot_vec[i]));
    }
  }

  return true;
  
}


template<class Scalar>
bool ForwardSensitivityStepper<Scalar>::getNodes(
  std::vector<Scalar>* time_vec
  ) const
{
  TEST_FOR_EXCEPT("Not implemented yet but we can!");
  return false;
}


template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::removeNodes(
  std::vector<Scalar>& time_vec
  )
{
  TEST_FOR_EXCEPT("Not implemented yet but we can!");
  return false;
}


template<class Scalar> 
int ForwardSensitivityStepper<Scalar>::getOrder() const
{
  return sensStepper_->getOrder();
  // Note: This assumes that stateStepper will have the same order!
}


// private


template<class Scalar> 
Teuchos::RefCountPtr<const Thyra::DefaultProductVector<Scalar> >
ForwardSensitivityStepper<Scalar>::create_x_bar_vec(
  const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &x_vec,
  const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &s_bar_vec
  ) const
{

  using Teuchos::tuple;
  using Teuchos::RefCountPtr;
  typedef RefCountPtr<const Thyra::VectorBase<Scalar> > RCPCV;

  return Thyra::defaultProductVector<Scalar>(
    x_bar_space_, tuple<RCPCV>(x_vec,s_bar_vec)
    );
  
}


} // namespace Rythmos


#endif //RYTHMOS_FORWARD_SENSITIVITY_STEPPER_HPP
