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

#ifndef Rythmos_IMPLICITBDF_STEPPER_DECL_H
#define Rythmos_IMPLICITBDF_STEPPER_DECL_H

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_StepControlStrategyAcceptingStepperBase.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_SolveSupportTypes.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {

/** \brief . */
template<class Scalar>
class ImplicitBDFStepper 
  : virtual public SolverAcceptingStepperBase<Scalar>
  , virtual public StepControlStrategyAcceptingStepperBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \brief Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  ImplicitBDFStepper();

  /** \brief . */
  ImplicitBDFStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> >  &model
    ,const RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver
    );

  /** \brief . */
  ImplicitBDFStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> >  &model
    ,const RCP<Thyra::NonlinearSolverBase<Scalar> >  &solver
    ,const RCP<Teuchos::ParameterList> &parameterList
    );

  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > get_solution() const;

  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > get_residual() const;

  /** \brief . */
  const Thyra::VectorBase<Scalar>& getxHistory(int index) const;

  /** \brief . */
  void setStepControlData(const StepperBase<Scalar> & stepper);

  //@}

  /** \name Overridden from StepControlStrategyAcceptingStepperBase */
  //@{

  /** \brief . */
  void setStepControlStrategy(
      const RCP<StepControlStrategyBase<Scalar> >& stepControlStrategy
      );

  /** \brief . */
  RCP<StepControlStrategyBase<Scalar> > 
    getNonconstStepControlStrategy();

  /** \brief . */
  RCP<const StepControlStrategyBase<Scalar> > 
    getStepControlStrategy() const;

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

  /** \brief Overridden from StepperBase */
  //@{
 
  /** \brief Returns true. */
  bool isImplicit() const;

  /** \brief Returns true. */
  bool supportsCloning() const;

  /** \brief Creates copies of all internal data (including the parameter
   * list) except the model which is assumed to stateless.
   *
   * If a shallow copy of the model is not appropirate for some reason, then
   * the client can simply reset the model using
   * <tt>returnVal->setModel()</tt>.
   */
  RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

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
    const Array<Scalar>& time_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
    
  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec
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

  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  //
  // Private data members
  //

  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  Rythmos::SingleResidualModelEvaluator<Scalar>   neModel_;

  RCP<Thyra::VectorBase<Scalar> > xn0_;
  RCP<Thyra::VectorBase<Scalar> > xpn0_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_base_;
  Array<RCP<Thyra::VectorBase<Scalar> > > xHistory_;
  RCP<Thyra::VectorBase<Scalar> > ee_;
  RCP<Thyra::VectorBase<Scalar> > residual_;

  RCP<StepControlStrategyBase<Scalar> > stepControl_; 

  Scalar time_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;

  Scalar hh_;        // Current step-size
  int currentOrder_; // Current order of integration
  int maxOrder_;     // maximum order = std::min(5,user option maxord) - see below.
  int usedOrder_;    // order used in current step (used after currentOrder is updated)
  Array<Scalar> alpha_;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  Scalar LETvalue_;   // ck * enorm
  EStepLETStatus stepLETStatus_; // Local Error Test Status
  Array<Scalar> gamma_;    // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  // calculate time derivative of history array for predictor 
  Array<Scalar> beta_;     // coefficients used to evaluate predictor from history array
  Array<Scalar> psi_;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
  // compute $\beta_j(n)$
  Scalar alpha_s_;    // $\alpha_s$ fixed-leading coefficient of this BDF method
  int numberOfSteps_;// number of total time integration steps taken
  int nef_; // number of error failures 
  Scalar usedStep_;
  int nscsco_;
  bool haveInitialCondition_;
  bool isInitialized_;

  int newtonConvergenceStatus_;

  RCP<Teuchos::ParameterList> parameterList_;

  //
  // Private member functions
  //

  void defaultInitializeAll_();
  void getInitialCondition_();
  void obtainPredictor_();
  void interpolateSolution_(
    const Scalar& timepoint,
    Thyra::VectorBase<Scalar>* x_ptr_,
    Thyra::VectorBase<Scalar>* xdot_ptr_,
    ScalarMag* accuracy_ptr_
    ) const;
  void updateHistory_();
  void restoreHistory_();
  void updateCoeffs_();
  void initialize_();
  void completeStep_();

};

template<class Scalar>
RCP<ImplicitBDFStepper<Scalar> > implicitBDFStepper();

template<class Scalar>
RCP<ImplicitBDFStepper<Scalar> > implicitBDFStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver,
  const RCP<Teuchos::ParameterList> &parameterList
  );

} // namespace Rythmos 

#endif // Rythmos_IMPLICITBDF_STEPPER_DECL_H
