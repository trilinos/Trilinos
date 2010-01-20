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
#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_ForwardSensitivityModelEvaluatorBase.hpp"
#include "Rythmos_ForwardSensitivityImplicitModelEvaluator.hpp"
#include "Rythmos_ForwardSensitivityExplicitModelEvaluator.hpp"
#include "Rythmos_StateAndForwardSensitivityModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_SingleResidualModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
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
 * The most general form of the parameterized state equation is:

 \verbatim

   f(x_dot(t),x(t),p) = 0, over t = [t0,tf]

   x(t0) = x_init(p)

 \endverbatim

 * As shown above, the parameters <tt>p</tt> are assumed to be steady state
 * parameters and can enter through the intial condition and/or through the
 * DAE equation itself.  This class also supports a form of the problem where
 * parameters <tt>p</tt> are only assumed to bin the initial condition
 * <tt>x_init(p)</tt> and not in the state equation.  In this case, you can
 * just drop out all of the terms <tt>d(f)/d(p)</tt> shown in the equations
 * below because they are zero.
 *
 * The forward sensitivity equations that are solved along with the state
 * equation, written in multi-vector form, are:

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
 * system.  See the class <tt>StateAndForwardSensitivityModelEvaluatorBase</tt>
 * for a description of how to get at the components of <tt>x</tt>,
 * <tt>s_bar</tt>, and <tt>S</tt> given <tt>x_bar</tt>.
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
 * Note that this class does not implement the function <tt>setModel()</tt>
 * since any arbitrary combined state + sensitivity model can not be
 * supported.
 *
 * \section Rythmos_ForwardSensitivityStepper_details_sec Implementation Details
 *
 * There are a variety of ways that one can go about implementing a state
 * plus forward sensitivity stepper.  Three ways for doing this are described
 * in the report "Design of New DASPK for Sensitivity Analysis" by Shengtai Li
 * and Linda Petzold.  The three ways are the <em>simultaneous corrector</em>,
 * the <em>staggered direct</em> and the <em>staggered corrector</em> methods.
 *
 * The <em>simultaneous corrector</em> method would be equivalent to forming
 * one big ModelEvaluator for the state and sensitivities where the "state"
 * variables would be the <tt>x_bar</tt> variables described above and then
 * this big system would be solved with a single stepper object and as single
 * nonlinear solver.  The advantage of this approach is that it makes great
 * reuse of all of the timestepping software.  Also, by being able to
 * specialize the nonlinear solver (which you can't do in most software) you
 * could set up the nonlinear solver to first solve the nonlinear state
 * timestep equation, and then solve the linear sensitivity equations.  The
 * problem with this approach would be that it would be very wasteful if the
 * timestep had to be cut back in order to reduce the local truncation error
 * of the state solution.  This would result in the sensitivity solution being
 * thrown away for each cut-back iteration.  Because of this fundamental
 * problem, we have not implemented the simultaneous corrector method.
 * Actually, we are not really sure why anyone implements ths method.
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
 * While this class does provide a full ModelEvaluator for the full state plus
 * forward sensitivity DAE <tt>f_bar(x_bar_hat,x_bar)</tt> it is not solved
 * all at one time as described above.  Instead, the step is solved first for
 * the state equation and then a ModelEvaluator for just the linear forward
 * sensitivity equations is formed and is solved over the same time step as
 * the forward solve.
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
 * steppers and models later with a little work.  With an explicit method and
 * model, we don't need to reuse W_tilde so this is easier in a way!
 *
 * ToDo: Finish documentation!
 */
template<class Scalar> 
class ForwardSensitivityStepper
  : virtual public StepperBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Intializers, Misc. */
  //@{

  /** \brief Constructs to uninitialized. */
  ForwardSensitivityStepper();

  /** \brief Intialize for synced state and sens steppers.
   *
   * \param stateModel [in,persisting] The ModelEvaluator that defines the
   * parameterized state model <tt>f(x_dot,x,p)</tt>.
   *
   * \param p_index [in] The index of the parameter subvector in
   * <tt>stateModel</tt> for which sensitivities will be computed for.
   *
   * \param baseStatePoint [in] Whatever input arguments are needed to define
   * the state of the model including the parameters except x, x_dot, and t!
   *
   * \param stateStepper [in,persisting] The stepper object that will be used
   * to advance the state solution <tt>x(t)</tt>.  This stepper need not be
   * setup with a model or a nonlinear timestep solver.  All this stepper
   * object needs is to be given its parameters to determine exactly what
   * timestepping algorithm will be employed.  The model and the timestep
   * solver objects will be set internally.
   *
   * \param stateTimeStepSolver [in,persisting] The nonlinear solver object
   * that is used to solve for the state timestep equation.  This is needed to
   * extract the Jacobian matrix that is used in the sensitivity model.  If
   * the stepper is not an implicit stepper and does not use an implicit time
   * step solver, then this argument should be left null.
   *
   * \param sensStepper [in,persisting] The stepper object that will be used
   * to advance the sensitivity solution <tt>S(t)</tt>.  This stepper need not
   * be setup with a model or a nonlinear timestep solver.  All this stepper
   * object needs is to be given its parameters to determine exactly what
   * timestepping algorithm will be employed.  The model and the timestep
   * solver objects will be set internally.  If this argument is null, then
   * the <tt>stateStepper</tt> object will be cloned to generate this stepper
   * object.  The most common use cases should just pass in
   * <tt>Teuchos::null</tt> and just use the identical stepper as the state
   * stepper.  However, this argument allows a client to specialize exactly
   * what the sensitivity stepper does and therefore this hook is allowed.
   *
   * \param sensTimeStepSolver [in,persisting] The nonlinear solver object
   * that is used to solve for the (linear) sensitivity timestep equation.  If
   * the stepper is not an implicit stepper and does not use an implicit
   * timestep solver, then this argument can be left null.  If the stepper is
   * implicit, and this argument is left null, then a
   * <tt>Thyra::LinearNonlinearSolver</tt> object will be created and used.
   * The most common use cases should just pass in <tt>Teuchos::null</tt> and
   * just use the simple linear nonlinear solver to will perform just a single
   * linear solve.  However, this argument allows a client to specialize
   * exactly what the nonlinear solver in the sensitivity stepper does and
   * therefore this hook is exposed to clients.
   *
   * Here <tt>*this</tt> is set up to synchronize the state and sensitivity
   * solvers.  Currently, error control is only done by the state stepper and
   * not the sensitivity stepper but the overall implementation has a high
   * degree of resuse and will therefore compute sensitivities quite fast.
   */
  void initializeSyncedSteppers(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const RCP<StepperBase<Scalar> > &stateStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
    const RCP<StepperBase<Scalar> > &sensStepper = Teuchos::null,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
    );

  /** \brief Intialize for synced state and sens steppers for an
   * initial-condition only parametrized sensitivity problem.
   *
   * \param stateModel [in,persisting] See initializeSyncedSteppers().
   *
   * \param p_space [in] The vector space for the parameterized initial
   * condition parameters.
   *
   * \param baseStatePoint [in]  See initializeSyncedSteppers().
   *
   * \param stateStepper [in,persisting] See initializeSyncedSteppers().
   *
   * \param stateTimeStepSolver [in,persisting] See initializeSyncedSteppers().
   *
   * \param sensStepper [in,persisting] See initializeSyncedSteppers().
   *
   * \param sensTimeStepSolver [in,persisting] See initializeSyncedSteppers().
   *
   * Here <tt>*this</tt> is set up to synchronize the state and sensitivity
   * solvers for an initial-condition only forward sensitivity problem.
   * Currently, error control is only done by the state stepper and not the
   * sensitivity stepper but the overall implementation has a high degree of
   * resuse and will therefore compute sensitivities quite fast.
   */
  void initializeSyncedSteppersInitCondOnly(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
    const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& stateBasePoint,
    const RCP<StepperBase<Scalar> >& stateStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& stateTimeStepSolver,
    const RCP<StepperBase<Scalar> >& sensStepper = Teuchos::null,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& sensTimeStepSolver = Teuchos::null
    );

  /** \brief Intialize for decoupled state and sens steppers.
   *
   * \param stateModel [in,persisting] See <tt>initializeSyncedSteppers()</tt>.
   *
   * \param p_index [in] See <tt>initializeSyncedSteppers()</tt>.
   *
   * \param baseStatePoint [in] See <tt>initializeSyncedSteppers()</tt>.
   *
   * \param stateStepper [in,persisting] See
   * <tt>initializeSyncedSteppers()</tt>.
   *
   * \param stateTimeStepSolver [in,persisting] See
   * <tt>initializeSyncedSteppers()</tt>.
   *
   * \param stateIntegrator [in,persisting] The intergrator that will be used
   * to integrate the state given <tt>stateStepper</tt>.  This integrator must
   * be set up with a trailing interpolation buffer in order to be able to
   * allow for complete flexibility in how the time steps for the sens
   * equations are solved.
   *
   * \param finalTime [in] The final time for the state integrator.  This is
   * needed to initialize <tt>stateIntegrator</tt> with <tt>stateStepper</tt>.
   *
   * \param sensStepper [in,persisting] See
   * <tt>initializeSyncedSteppers()</tt>.
   *
   * \param sensTimeStepSolver [in,persisting] See
   * <tt>initializeSyncedSteppers()</tt>.
   *
   * Here <tt>*this</tt> is set up to run the state and sens steppers
   * completely independently; each with the their own error control
   * strategies.  The state stepper in driven through the state integrator
   * which in turn is driven by the ForwardSensitivityModelEvaluatorBase that is
   * driven by the sens stepper.
   */
  void initializeDecoupledSteppers(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const RCP<StepperBase<Scalar> > &stateStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
    const RCP<IntegratorBase<Scalar> > &stateIntegrator,
    const Scalar &finalTime,
    const RCP<StepperBase<Scalar> > &sensStepper = Teuchos::null,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
    );

  /** \brief Return if the state model is const-only or not. */
  bool stateModelIsConst() const;

  /** \brief Return the state model that was passed into an initialize
   * function.
   */
  RCP<const Thyra::ModelEvaluator<Scalar> >
  getStateModel() const;

  /** \brief Return the state stepper that was passed into an initialize
   * function.
   */
  RCP<StepperBase<Scalar> >
  getNonconstStateStepper();

  /** \brief Return the forward sensitivity model evaluator object that got
   * created internally when the initialize function was called.
   */
  RCP<const ForwardSensitivityModelEvaluatorBase<Scalar> >
  getFwdSensModel() const;

  /** \brief Return the state and forward sensitivity model evaluator object
   * that got created internally when the nitialize function was called.
   *
   * This is also the same model returned by the function <tt>getModel()</tt>,
   * except through it's concrete subclass type.
   */
  RCP<const StateAndForwardSensitivityModelEvaluator<Scalar> >
  getStateAndFwdSensModel() const;

  //@}

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{

  /** \brief Returns false. */
  bool acceptsModel() const;

  /** \brief Throws exception. */
  void setModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model
    );

  /** \brief Throws exception. */
  void setNonconstModel(
    const RCP<Thyra::ModelEvaluator<Scalar> >& model
    );

  /** \brief Returns <tt>getStateAndFwdSensModel()</tt>.
   *
   * Warning, currently the returned model does not implement evalModel(...) 
   * or define a W object.  It is just used for getting the spaces and for
   * creating an InArgs object for setting the initial condition.
   */
  RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> > getNonconstModel();

  // RAB: ToDo: 2007/05/15: I need to talk with Todd about potentially
  // removing the setModel() and getModel() functions from the StepperBase
  // interface.  In the case of this forward sensitivity solver, I am not sure
  // that it makes a lot of sense to define a model.  This surely will not be
  // the model that a generic client would expect.  The assumption I am sure
  // would be that this model has the same space for x as the interpolation
  // buffer but that is not true in this case.

  /** \brief Sets the full initial condition for <tt>x_bar = [ x; s_bar] </tt>
   * and <tt>x_bar_dot = [ x_dot; s_bar_dot ]</tt> as well as the initial time
   * and the parameter values.
   *
   * The InArgs object must be created using
   * <tt>this->getModel()->createInArgs()</tt> and then populated with the
   * initial values.  The product vectors for <tt>x_bar</tt> and
   * <tt>x_bar_dot</tt> can be created using
   * <tt>this->getStateAndFwdSensModel()->create_x_bar_vec(...)</tt>.  All of
   * the input objects in <tt>state_and_sens_ic</tt> will be cloned and
   * therefore no memory of the objects in <tt>state_and_sens_ic</tt> will be
   * retained after calling this function.
   */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &state_and_sens_ic
    );

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getInitialCondition() const;

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
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
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

  /** \name Deprecated. */
  //@{

  /** \brief Deprecated. */
  void initialize(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const RCP<StepperBase<Scalar> > &stateStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
    const RCP<StepperBase<Scalar> > &sensStepper = Teuchos::null,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
    )
    {
      initializeSyncedSteppers(
        stateModel, p_index, stateBasePoint, stateStepper, stateTimeStepSolver,
        sensStepper, sensTimeStepSolver
        );
    }

  //@}

private:
  // ///////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<Thyra::ModelEvaluator<Scalar> > CNCME;

  // /////////////////////////
  // Private data members

  bool forceUpToDateW_;
  CNCME stateModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_;
  RCP<StepperBase<Scalar> > stateStepper_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > stateTimeStepSolver_;
  RCP<IntegratorBase<Scalar> > stateIntegrator_;
  Scalar finalTime_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateAndSensBasePoint_;
  RCP<StepperBase<Scalar> > sensStepper_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > sensTimeStepSolver_;

  bool isSingleResidualStepper_;
  RCP<ForwardSensitivityModelEvaluatorBase<Scalar> > sensModel_;
  RCP<StateAndForwardSensitivityModelEvaluator<Scalar> > stateAndSensModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stateBasePoint_t_;

  static const std::string forceUpToDateW_name_;
  static const bool forceUpToDateW_default_;

  // /////////////////////////
  // Private member functions
  

  // Common initialization helper
  //
  // Preconditions:
  // (*) p_index >=0 or nonnull(p_space) == true
  //
  void initializeCommon(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
    const int p_index,
    const RCP<const Thyra::VectorSpaceBase<Scalar> > &p_space,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
    const RCP<StepperBase<Scalar> > &stateStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
    const RCP<StepperBase<Scalar> > &sensStepper,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver
    );

  Scalar takeSyncedStep( Scalar dt, StepSizeType stepType );

  Scalar takeDecoupledStep( Scalar dt, StepSizeType stepType );

};


// 2009/09/05: rabartl: ToDo: To fix the const and non-const handling of the
// stateModel in this class is going to be a lot of work but here is what
// needs to be done:
//
// (*) Duplicate each function that sets the stateModel, one for const and one
// for non-const.
//
// (*) Create a single a private version for each of these functions that
// accepts a Teuchos::ConstNonconstObjectContainer<> object and will implement
// the guts of the set up same as the existing functions.
//
// (*) Get all of the concrete StepperBase subclasses to implement the
// setModel(const RCP<const ME>&) and modelIsConst() functions and get them to
// use the Teuchos::ConstNonconstObjectContainer<> class as described above.
// This should be pretty easy as the only function that needs to be addressed
// in most cases is just the setModel(...) function.
//


/** \brief Nonmember constructor.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar> 
inline
RCP<ForwardSensitivityStepper<Scalar> >
forwardSensitivityStepper()
{
  return Teuchos::rcp(new ForwardSensitivityStepper<Scalar>());
}


/** \brief Nonmember constructor.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar> 
inline
RCP<ForwardSensitivityStepper<Scalar> >
forwardSensitivityStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const RCP<StepperBase<Scalar> > &sensStepper = Teuchos::null,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver = Teuchos::null
  )
{
  RCP<ForwardSensitivityStepper<Scalar> >
    fwdSensStepper = Teuchos::rcp(new ForwardSensitivityStepper<Scalar>());
  fwdSensStepper->initializeSyncedSteppers(
    stateModel, p_index, stateBasePoint, stateStepper, stateTimeStepSolver );
  return fwdSensStepper;
}


/** \brief Return the index of the parameter subvector in the underlying state
 * model.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar>
int getParameterIndex(
  const ForwardSensitivityStepper<Scalar> &fwdSensStepper
  )
{
  return fwdSensStepper.getFwdSensModel()->get_p_index();
}


/** \brief Set up default initial conditions for the state and sensitivity
 * stepper with default zero initial conditions for the sensitivity
 * quantities.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar> 
inline
Thyra::ModelEvaluatorBase::InArgs<Scalar>
createStateAndSensInitialCondition(
  const ForwardSensitivityStepper<Scalar> &fwdSensStepper,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &state_ic,
  const RCP<const Thyra::MultiVectorBase<Scalar> > S_init = Teuchos::null,
  const RCP<const Thyra::MultiVectorBase<Scalar> > S_dot_init = Teuchos::null
  )
{

  using Teuchos::outArg;
  using Thyra::assign;
  typedef Thyra::ModelEvaluatorBase MEB;
  
  RCP<const Thyra::VectorBase<Scalar> > s_bar_init;
  if (nonnull(S_init)) {
    s_bar_init = create_s_bar_given_S(*fwdSensStepper.getFwdSensModel(), S_init);
  }
  else {
    RCP<Thyra::VectorBase<Scalar> > s_bar_init_loc =
      createMember(fwdSensStepper.getFwdSensModel()->get_x_space());
    assign( outArg(*s_bar_init_loc), 0.0 );
    s_bar_init = s_bar_init_loc;
  }

  RCP<const Thyra::VectorBase<Scalar> > s_bar_dot_init;
  if (nonnull(S_dot_init)) {
    s_bar_dot_init = create_s_bar_given_S(*fwdSensStepper.getFwdSensModel(), S_dot_init);
  }
  else {
    RCP<Thyra::VectorBase<Scalar> > s_bar_dot_init_loc =
      createMember(fwdSensStepper.getFwdSensModel()->get_x_space());
    assign( outArg(*s_bar_dot_init_loc), 0.0 );
    s_bar_dot_init = s_bar_dot_init_loc;
  }
  
  RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
    stateAndSensModel = fwdSensStepper.getStateAndFwdSensModel();
  
  MEB::InArgs<Scalar>
    state_and_sens_ic = fwdSensStepper.getModel()->createInArgs();
  
  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar = [ x; s_bar ]
  state_and_sens_ic.set_x(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x(), s_bar_init)
    );
  // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
  state_and_sens_ic.set_x_dot(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(), s_bar_dot_init)
    );

  return state_and_sens_ic;

}


/** \brief Extract out the initial condition for just the state model given
 * the initial condition for the state and sensitivity model.
 *
 * \relates ForwardSensitivityStepper
 */
template<class Scalar> 
inline
Thyra::ModelEvaluatorBase::InArgs<Scalar>
extractStateInitialCondition(
  const ForwardSensitivityStepper<Scalar> &fwdSensStepper,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &state_and_sens_ic
  )
{

  using Thyra::productVectorBase;
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>
    state_ic = fwdSensStepper.getStateModel()->createInArgs();
  
  // Copy time, parameters etc.
  state_ic.setArgs(state_and_sens_ic);
  state_ic.set_x(
    productVectorBase(state_and_sens_ic.get_x())->getVectorBlock(0));
  state_ic.set_x_dot(
    productVectorBase(state_and_sens_ic.get_x_dot())->getVectorBlock(0));
  
  return state_ic;
  
}


//
// Implementation
//


// Static members


template<class Scalar> 
const std::string ForwardSensitivityStepper<Scalar>::forceUpToDateW_name_
= "Force Up-To-Date Jacobian";

template<class Scalar> 
const bool ForwardSensitivityStepper<Scalar>::forceUpToDateW_default_
= true;


// Constructors, Intializers, Misc.


template<class Scalar> 
ForwardSensitivityStepper<Scalar>::ForwardSensitivityStepper()
  :forceUpToDateW_(false),
   isSingleResidualStepper_(false)
{}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::initializeSyncedSteppers(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const RCP<StepperBase<Scalar> > &sensStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver
  )
  
{
  initializeCommon( stateModel, p_index, Teuchos::null, stateBasePoint, stateStepper,
    stateTimeStepSolver, sensStepper, sensTimeStepSolver );
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::initializeSyncedSteppersInitCondOnly(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const RCP<const Thyra::VectorSpaceBase<Scalar> >& p_space,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& stateBasePoint,
  const RCP<StepperBase<Scalar> >& stateStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& stateTimeStepSolver,
  const RCP<StepperBase<Scalar> >& sensStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& sensTimeStepSolver
  )
{
  initializeCommon(stateModel, -1, p_space, stateBasePoint, stateStepper,
    stateTimeStepSolver, sensStepper, sensTimeStepSolver );
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::initializeDecoupledSteppers(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &stateModel,
  const int p_index,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const RCP<IntegratorBase<Scalar> > &stateIntegrator,
  const Scalar &finalTime,
  const RCP<StepperBase<Scalar> > &sensStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver
  )
{
  TEUCHOS_ASSERT(nonnull(stateIntegrator));
  initializeCommon( stateModel, p_index, Teuchos::null, stateBasePoint, stateStepper,
    stateTimeStepSolver, sensStepper, sensTimeStepSolver );
  stateIntegrator_ = stateIntegrator;
  finalTime_ = finalTime;
}


template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::stateModelIsConst() const
{
  return stateModel_.isConst();
}


template<class Scalar> 
RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getStateModel() const
{
  return stateModel_.getConstObj();
}

  
template<class Scalar> 
RCP<StepperBase<Scalar> >
ForwardSensitivityStepper<Scalar>::getNonconstStateStepper()
{
  return stateStepper_;
}

  
template<class Scalar> 
RCP<const ForwardSensitivityModelEvaluatorBase<Scalar> >
ForwardSensitivityStepper<Scalar>::getFwdSensModel() const
{
  return sensModel_;
}

  
template<class Scalar> 
RCP<const StateAndForwardSensitivityModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getStateAndFwdSensModel() const
{
  return stateAndSensModel_;
}


// Overridden from Teuchos::ParameterListAcceptor


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);
  forceUpToDateW_ = paramList->get(forceUpToDateW_name_,forceUpToDateW_default_);
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const Teuchos::ParameterList>
ForwardSensitivityStepper<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set( forceUpToDateW_name_, forceUpToDateW_default_,
      "If set to true, then the Jacobian matrix W used in the\n"
      "state timestep equation will be forced to be up to date\n"
      "with the final value for x for the nonlinear solve.  If\n"
      "you are willing to live with slightly less accurate sensitivities\n"
      "then set this to false."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from StepperBase

template<class Scalar> 
bool ForwardSensitivityStepper<Scalar>::acceptsModel() const
{
  return false;
}

template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model
  )
{
  TEST_FOR_EXCEPT_MSG( true,
    "Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setNonconstModel(
  const RCP<Thyra::ModelEvaluator<Scalar> >& model
  )
{
  TEST_FOR_EXCEPT_MSG( true,
    "Error, this stepper subclass does not accept a model"
    " as defined by the StepperBase interface!");
}


template<class Scalar> 
RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getModel() const
{
  return stateAndSensModel_;
}


template<class Scalar> 
RCP<Thyra::ModelEvaluator<Scalar> >
ForwardSensitivityStepper<Scalar>::getNonconstModel() 
{
  return stateAndSensModel_;
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &state_and_sens_ic
  )
{
  
  typedef Thyra::ModelEvaluatorBase MEB;

  stateAndSensBasePoint_ = state_and_sens_ic;

  // Get the product vectors for x_bar = [ x; s_bar ] and x_bar_dot

  TEST_FOR_EXCEPTION(
    is_null(state_and_sens_ic.get_x()), std::logic_error,
    "Error, the initial condition for x_bar = [ x; s_bar ] can not be null!" );

  const RCP<const Thyra::ProductVectorBase<Scalar> >
    x_bar_init = Thyra::productVectorBase<Scalar>(
      state_and_sens_ic.get_x()
      );

  RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_dot_init;
  if (state_and_sens_ic.supports(MEB::IN_ARG_x_dot)) {
      x_bar_dot_init = Thyra::productVectorBase<Scalar>(
        state_and_sens_ic.get_x_dot()
        );
  }

  // Remove x and x_dot from state_and_sens_ic_in to avoid cloning x and x dot!
  
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
    state_and_sens_ic_no_x = state_and_sens_ic;
  state_and_sens_ic_no_x.set_x(Teuchos::null);
  if (state_and_sens_ic_no_x.supports(MEB::IN_ARG_x_dot)) {
    state_and_sens_ic_no_x.set_x_dot(Teuchos::null);
  }

  // Set initial condition for the state

  MEB::InArgs<Scalar> state_ic = stateModel_->createInArgs();
  state_ic.setArgs(state_and_sens_ic_no_x,true,true); // Set time, parameters etc.
  state_ic.set_x(x_bar_init->getVectorBlock(0)->clone_v());
  if (state_ic.supports(MEB::IN_ARG_x_dot)) {
    state_ic.set_x_dot(
        !is_null(x_bar_dot_init)
        ? x_bar_dot_init->getVectorBlock(0)->clone_v()
        : Teuchos::null
        );
  }
  stateStepper_->setInitialCondition(state_ic);

  // Set up the integrator if needed
  //if (!is_null(stateIntegrator_)) {
  //  stateIntegrator_->setStepper( stateStepper_, finalTime_ );
  //  sensModel_->setStateIntegrator( stateIntegrator_, state_ic );
  //}

  // Set initial condition for the sensitivities
  
  MEB::InArgs<Scalar> sens_ic = sensModel_->createInArgs();
  sens_ic.setArgs(state_and_sens_ic_no_x,true,true); // Set time etc.
  sens_ic.set_x(x_bar_init->getVectorBlock(1)->clone_v());
  if (sens_ic.supports(MEB::IN_ARG_x_dot)) {
    sens_ic.set_x_dot(
        !is_null(x_bar_dot_init)
        ? x_bar_dot_init->getVectorBlock(1)->clone_v()
        : Teuchos::null
        );
  }
  sensStepper_->setInitialCondition(sens_ic);

}


template<class Scalar> 
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ForwardSensitivityStepper<Scalar>::getInitialCondition() const
{
  return stateAndSensBasePoint_;
}
 

template<class Scalar> 
Scalar
ForwardSensitivityStepper<Scalar>::takeStep(
  Scalar dt, StepSizeType stepType
  )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityStepper::takeStep");
#endif

  if (!is_null(stateIntegrator_)) {
    return takeDecoupledStep(dt,stepType);
  }

  return takeSyncedStep(dt,stepType);

}


template<class Scalar> 
const StepStatus<Scalar>
ForwardSensitivityStepper<Scalar>::getStepStatus() const
{

  const StepStatus<Scalar> sensStepStatus = sensStepper_->getStepStatus();
  StepStatus<Scalar> stepStatus;
  
  stepStatus.message = sensStepStatus.message;
  stepStatus.stepStatus = sensStepStatus.stepStatus;
  stepStatus.stepLETStatus = sensStepStatus.stepLETStatus;
  stepStatus.stepSize = sensStepStatus.stepSize;
  stepStatus.order = sensStepStatus.order;
  stepStatus.time = sensStepStatus.time;
  stepStatus.stepLETValue = sensStepStatus.stepLETValue;
  stepStatus.extraParameters = sensStepStatus.extraParameters;

  if (is_null(stateIntegrator_)) {
    const StepStatus<Scalar> 
      stateStepStatus = stateStepper_->getStepStatus();
    if (!is_null(stateStepStatus.solution) && !is_null(sensStepStatus.solution))
      stepStatus.solution = stateAndSensModel_->create_x_bar_vec(
        stateStepStatus.solution, sensStepStatus.solution
        );
    if (!is_null(stateStepStatus.solutionDot) && !is_null(sensStepStatus.solutionDot))
      stepStatus.solutionDot = stateAndSensModel_->create_x_bar_vec(
        stateStepStatus.solutionDot, sensStepStatus.solutionDot
        );
  }

  return stepStatus;

}


// Overridden from InterpolationBufferBase


template<class Scalar> 
RCP<const Thyra::VectorSpaceBase<Scalar> >
ForwardSensitivityStepper<Scalar>::get_x_space() const
{
  return stateAndSensModel_->get_x_space();
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::addPoints(
  const Array<Scalar>& time_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
  )
{
  TEST_FOR_EXCEPT("Not implemented addPoints(...) yet but we could if we wanted!");
}


template<class Scalar> 
TimeRange<Scalar>
ForwardSensitivityStepper<Scalar>::getTimeRange() const
{
  return sensStepper_->getTimeRange();
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_bar_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_bar_dot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{

  using Teuchos::as;

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPT( as<int>(time_vec.size()) == 0 );
#endif

  const int numTimePoints = time_vec.size();

  if (x_bar_vec)
    x_bar_vec->clear();

  if (x_bar_dot_vec)
    x_bar_dot_vec->clear();
  
  Array<RCP<const Thyra::VectorBase<Scalar> > >
    x_vec, x_dot_vec;

  if (!is_null(stateIntegrator_)) {
    stateIntegrator_->getPoints(
      time_vec,
      x_bar_vec ? &x_vec: 0,
      x_bar_dot_vec ? &x_dot_vec: 0,
      0 // Ignoring accuracy from state for now!
      );
  }
  else {
    stateStepper_->getPoints(
      time_vec,
      x_bar_vec ? &x_vec: 0,
      x_bar_dot_vec ? &x_dot_vec: 0,
      0 // Ignoring accuracy from state for now!
      );
  }
  
  Array<RCP<const Thyra::VectorBase<Scalar> > >
    s_bar_vec, s_bar_dot_vec;

  sensStepper_->getPoints(
    time_vec,
    x_bar_vec ? &s_bar_vec: 0,
    x_bar_dot_vec ? &s_bar_dot_vec: 0,
    accuracy_vec
    );
  
  if ( x_bar_vec ) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      x_bar_vec->push_back(
        stateAndSensModel_->create_x_bar_vec(x_vec[i],s_bar_vec[i])
        );
    }
  }
  
  if ( x_bar_dot_vec ) {
    for ( int i = 0; i < numTimePoints; ++i ) {
      x_bar_dot_vec->push_back(
        stateAndSensModel_->create_x_bar_vec(x_dot_vec[i],s_bar_dot_vec[i])
        );
    }
  }

}


template<class Scalar>
void ForwardSensitivityStepper<Scalar>::getNodes(
  Array<Scalar>* time_vec
  ) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (is_null(stateIntegrator_) && is_null(stateStepper_)) {
    return;
  }
  if (!is_null(stateIntegrator_)) {
    stateIntegrator_->getNodes(time_vec);
  }
  else {
    stateStepper_->getNodes(time_vec);
  }
}


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::removeNodes(
  Array<Scalar>& time_vec
  )
{
  TEST_FOR_EXCEPT("Not implemented yet but we can!");
}


template<class Scalar> 
int ForwardSensitivityStepper<Scalar>::getOrder() const
{
  return sensStepper_->getOrder();
  // Note: This assumes that stateStepper will have the same order!
}


// private


template<class Scalar> 
void ForwardSensitivityStepper<Scalar>::initializeCommon(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& stateModel,
  const int p_index,
  const RCP<const Thyra::VectorSpaceBase<Scalar> > &p_space,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &stateBasePoint,
  const RCP<StepperBase<Scalar> > &stateStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &stateTimeStepSolver,
  const RCP<StepperBase<Scalar> > &sensStepper,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &sensTimeStepSolver
  )
{

  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_dynamic_cast;

  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // Validate input
  //

  TEUCHOS_ASSERT( p_index >= 0 || nonnull(p_space) );
  if (nonnull(p_space)) {
    TEUCHOS_ASSERT_EQUALITY(p_index, -1);
  }
  if (p_index >= 0) {
    TEUCHOS_ASSERT(is_null(p_space));
  }
  TEST_FOR_EXCEPT( is_null(stateModel) );
  TEST_FOR_EXCEPT( is_null(stateStepper) );
  if (stateStepper->isImplicit()) {
    TEST_FOR_EXCEPT( is_null(stateTimeStepSolver) ); // allow to be null for explicit methods
  }

  //
  // Create the sensModel which will do some more validation
  //
  
  RCP<ForwardSensitivityModelEvaluatorBase<Scalar> > sensModel;
  MEB::InArgs<Scalar> stateModelInArgs = stateModel->createInArgs();
  if (stateModelInArgs.supports(MEB::IN_ARG_x_dot)) {
    // Implicit DE formulation
    sensModel = Teuchos::rcp(new ForwardSensitivityImplicitModelEvaluator<Scalar>);
  }
  else {
    // Explicit DE formulation
    sensModel = Teuchos::rcp(new ForwardSensitivityExplicitModelEvaluator<Scalar>);
  }

  if (p_index >= 0) {
    sensModel->initializeStructure(stateModel, p_index);
  }
  else {
    sensModel->initializeStructureInitCondOnly(stateModel, p_space);
  }
  
  //
  // Get the input objects
  //

  stateModel_.initialize(stateModel);

  stateBasePoint_ = stateBasePoint;

  stateStepper_ = stateStepper;
  
  stateTimeStepSolver_ = stateTimeStepSolver;

  sensModel_ = sensModel;

  stateAndSensModel_ = Teuchos::rcp(new StateAndForwardSensitivityModelEvaluator<Scalar>);
  stateAndSensModel_->initializeStructure(sensModel_);

  if (!is_null(sensStepper)) {
    sensStepper_ = sensStepper;
  }
  else {
    sensStepper_ = stateStepper_->cloneStepperAlgorithm();
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
    RCP<Thyra::LinearNonlinearSolver<Scalar> >
      linearNonlinearSolver(new Thyra::LinearNonlinearSolver<Scalar>);
    // ToDo: Set tolerance on the nonlinear solver???
    sensTimeStepSolver_ = linearNonlinearSolver;
  }

  //
  // Setup the steppers
  //

  isSingleResidualStepper_ = true; // ToDo: Add dynamic cast on
                                   // stateTimeStepSolver to check this!

  setStepperModel(Teuchos::inOutArg(*stateStepper_),stateModel_);
  if (stateStepper_->isImplicit()) {
    rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(
        stateStepper_,true)->setSolver(stateTimeStepSolver_);
  }
  sensStepper_->setModel(sensModel_); 
  if (sensStepper_->isImplicit()) {
    rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(
        sensStepper_,true)->setSolver(sensTimeStepSolver_);
  }

  stateBasePoint_t_ = stateModel_->createInArgs();

  // 2007/05/18: rabartl: ToDo: Move the above initialization code to give
  // setInitializeCondition(...) a chance to set the initial condition.

}


template<class Scalar> 
Scalar ForwardSensitivityStepper<Scalar>::takeSyncedStep(
  Scalar dt, StepSizeType stepType
  )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("Rythmos:ForwardSensitivityStepper::takeStep: synced",
    TopLevel);
#endif

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::VerboseObjectTempState<InterpolationBufferBase<Scalar> > VOTSIBB;
  typedef Thyra::ModelEvaluatorBase MEB;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool lowTrace =
    ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) );
  const bool mediumTrace =
    ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) );
  Teuchos::OSTab tab(out);

  if (lowTrace) {
    *out
      << "\nEntering " << TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeSyncedStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  //
  // A) Compute the state timestep
  //

  if (lowTrace) {
    *out
      << "\nTaking state step using stepper : "
      << stateStepper_->description() << "\n";
  }

  Scalar state_dt = -1.0;
  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityStepper::takeStep: stateStep");
#endif
    VOTSIBB stateStepper_outputTempState(stateStepper_,out,verbLevel);
    state_dt = stateStepper_->takeStep(dt,stepType);
  }

  if (state_dt < Scalar(-ST::one())) {
    if (lowTrace)
      *out << "\nThe state stepper has failed so return a failed timestep!\n";
    return state_dt;
  }

  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityStepper::takeStep: updateSensModel");
#endif
    // Set up the sensitivity model for this timestep
    sensModel_->initializePointState(Teuchos::inOutArg(*stateStepper_),forceUpToDateW_);
  } 

  //
  // C) Compute the sensitivity timestep for the exact same timestep as was
  // used for the state solve.
  //

  if (lowTrace) {
    *out
      << "\nTaking sensitivity step using stepper : "
      << sensStepper_->description() << "\n";
  }

  Scalar sens_dt = -1.0;
  {
#ifdef ENABLE_RYTHMOS_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityStepper::takeStep: sensStep");
#endif
    // Copy the step control data to make sure that the sensStepper takes the
    // same type of step that the statStepper took.  This is needed to ensure
    // that the W matrix is the same for one.
    sensStepper_->setStepControlData(*stateStepper_);
    VOTSIBB sensStepper_outputTempState(sensStepper_,out,verbLevel);
    sens_dt = sensStepper_->takeStep(state_dt,STEP_TYPE_FIXED);
  }

  if (mediumTrace) {
    const StepStatus<Scalar> sensStepStatus = sensStepper_->getStepStatus();
    *out << "\nSensitivity step status:\n" << sensStepStatus;
  }
  
  TEST_FOR_EXCEPTION(
    sens_dt != state_dt, std::logic_error,
    "Error, the sensitivity step failed for some reason.  We should\n"
    "just return a negative step size and reject the step but currently\n"
    "there is no way to roll back the state timestep it for back to\n"
    "the status before this function was called!"
    );

  // 2007/05/18: rabartl: ToDo: If stepType == STEP_TYPE_VARIABLE and the state
  // timestep sucessed but the sensitivity timestep failed, then we need to
  // really throw an excpetion because there is nothing that we can really do
  // here!

  // 2007/05/18: rabartl: ToDo: Replace the above std::logic_error type with
  // a Rythmos::CatastrophicFailure or just use Thyra::CatastrophicFailure!

  if (lowTrace) {
    *out
      << "\nLeaving " << TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeSyncedStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  return state_dt;

}


template<class Scalar> 
Scalar ForwardSensitivityStepper<Scalar>::takeDecoupledStep(
  Scalar dt, StepSizeType stepType
  )
{

#ifdef ENABLE_RYTHMOS_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Rythmos:ForwardSensitivityStepper::takeStep: decoupled");
#endif

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::VerboseObjectTempState<InterpolationBufferBase<Scalar> > VOTSIBB;
  typedef Thyra::ModelEvaluatorBase MEB;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool lowTrace =
    ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW) );
  const bool mediumTrace =
    ( !is_null(out) && as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM) );
  Teuchos::OSTab tab(out);

  if (lowTrace) {
    *out
      << "\nEntering " << TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeDecoupledStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }
  
  //
  // A) Take the sens timestep
  //

  if (lowTrace) {
    *out
      << "\nTaking sensitivity step using stepper : "
      << sensStepper_->description() << "\n";
  }

  Scalar sens_dt = -1.0;
  VOTSIBB sensStepper_outputTempState(sensStepper_,out,verbLevel);
  sens_dt = sensStepper_->takeStep(dt,stepType);

  if (mediumTrace) {
    const StepStatus<Scalar> sensStepStatus = sensStepper_->getStepStatus();
    *out << "\nSensitivity step status:\n" << sensStepStatus;
  }

  //
  // B) Wipe out all state interp buffer info before this sens timestep
  //
  
  //TEST_FOR_EXCEPT(true);

  if (lowTrace) {
    *out
      << "\nLeaving " << TypeNameTraits<ForwardSensitivityStepper<Scalar> >::name()
      << "::takeDecoupledStep("<<dt<<","<<toString(stepType)<<") ...\n"; 
  }

  return sens_dt;
  
}


} // namespace Rythmos


#endif //RYTHMOS_FORWARD_SENSITIVITY_STEPPER_HPP
