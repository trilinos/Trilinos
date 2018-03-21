// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKStartingStepSize_hpp
#define Tempus_StepperRKStartingStepSize_hpp

//#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionHistory.hpp"
//#include "Tempus_StepperState.hpp"
#include "Tempus_StepperExplicitRKObserver.hpp"


namespace Tempus {

/** \brief Starting Step Size for Explicit and Implicit Runge-Kutta Methods
 *
 *   Purpose: Return the optimal starting step-size
 *   Solving Ordinary Differential Equation I - nonstiff¬
 *   Hairer. pg. 169. Starting Step Size Algorithm
 *
 */
template<class Scalar>
class StepperRKStartingStepSize :
   virtual public StepperExplicitRKObserver<Scalar>
{
public:

  /// Constructor
  StepperRKStartingStepSize(){ isStarted_ = false;}

  /// Destructor
  virtual ~StepperRKStartingStepSize(){}

  /** \brief Determine the time step size.*/
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh, 
    StepperExplicitRK<Scalar> & ERK) override
  {
     if (isStarted_) return;

     Teuchos::RCP<SolutionState<Scalar> > workingState=sh->getWorkingState();
     Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
     const Scalar errorAbs = metaData->getTolRel();
     const Scalar errorRel = metaData->getTolAbs();
     int order = metaData->getOrder();
     Scalar time = metaData->getTime();

     // get the model
     Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> inModel = ERK.getModel();

     auto err_func = [] (Teuchos::RCP<Thyra::VectorBase<Scalar> > U,
           const Scalar rtol, const Scalar atol,
           Teuchos::RCP<Thyra::VectorBase<Scalar> > absU)
     {
        // compute err = Norm_{WRMS} with w = Atol + Rtol * | U | 
        Thyra::assign(absU.ptr(), *U);
        Thyra::abs(*U, absU.ptr()); // absU = | X0 |
        Thyra::Vt_S(absU.ptr(), rtol); // absU *= Rtol
        Thyra::Vp_S(absU.ptr(), atol); // absU += Atol
        Thyra::ele_wise_divide(Teuchos::as<Scalar>(1.0), *U, *absU, absU.ptr());
        Scalar err = Thyra::norm_inf(*absU);
        return err;
     };


     Teuchos::RCP<Thyra::VectorBase<Scalar> > stageX_, scratchX;
     stageX_ = Thyra::createMember(inModel->get_f_space());
     scratchX = Thyra::createMember(inModel->get_f_space());
     Thyra::assign(stageX_.ptr(), *(workingState->getX()));

     std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_(2);
     for (int i=0; i<2; ++i) {
        stageXDot_[i] = Thyra::createMember(inModel->get_f_space());
        assign(stageXDot_[i].ptr(), Teuchos::ScalarTraits<Scalar>::zero());
     }

     // A: one functione evaluation at F(t_0, X_0)
     typedef Thyra::ModelEvaluatorBase MEB;
     Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_ = inModel->getNominalValues();
     Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_ = inModel->createOutArgs();
     inArgs_.set_x(stageX_);
     if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);
     if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
     outArgs_.set_f(stageXDot_[0]); // K1
     inModel->evalModel(inArgs_,outArgs_);

     Scalar d0 = err_func(stageX_, errorRel, errorAbs, scratchX);
     Scalar d1 = err_func(stageXDot_[0], errorRel, errorAbs, scratchX);

     // b) first guess for the step size
     Scalar dt = Teuchos::as<Scalar>(0.01)*(d0/d1);

     // c) perform one explicit Euler step (X_1)
     Thyra::Vp_StV(stageX_.ptr(), dt, *(stageXDot_[0]));

     // compute F(t_0 + dt, X_1)
     inArgs_.set_x(stageX_);
     if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time + dt);
     if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
     outArgs_.set_f(stageXDot_[1]); // K2
     inModel->evalModel(inArgs_,outArgs_);

     // d) compute estimate of the second derivative of the solution
     // d2 = || f(t_0 + dt, X_1) - f(t_0, X_0) || / dt
     Teuchos::RCP<Thyra::VectorBase<Scalar> > errX;
     errX = Thyra::createMember(inModel->get_f_space());
     assign(errX.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
     Thyra::V_VmV(errX.ptr(), *(stageXDot_[1]), *(stageXDot_[0]));
     Scalar d2 = err_func(errX, errorRel, errorAbs, scratchX) / dt;

     // e) compute step size h_1 (from m = 0 order Taylor series)
     Scalar max_d1_d2 = std::max(d1, d2);
     Scalar h1 = std::pow((0.01/max_d1_d2),(1.0/(order+1)));

     // f) propse starting step size
     dt = std::min(100*dt, h1);

     // update order and dt
     metaData->setOrder(order);
     metaData->setDt(dt);
     isStarted_ = true;
  } // observeBeginTakeStep

private:
  bool isStarted_ = false;

};
} // namespace Tempus
#endif // Tempus_StepperRKStartingStepSize_hpp
