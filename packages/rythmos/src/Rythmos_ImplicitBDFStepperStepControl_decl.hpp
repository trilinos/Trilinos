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

#ifndef Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DECL_H
#define Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DECL_H

#include "Rythmos_ErrWtVecCalcAcceptingStepControlStrategyBase.hpp"

namespace Rythmos {

/** \brief . */
enum BDFactionFlag { ACTION_UNSET, ACTION_LOWER, ACTION_MAINTAIN, ACTION_RAISE };

// Step Control Strategy object for ImplicitBDFStpper
//
// Order of calls:
// setRequestedStepSize()
// nextStepSize()
// optional:  nextStepOrder()
// setCorrection
// acceptStep
// completeStep or rejectStep
// repeat
//
// 08/16/07 tscoffe:  This order of operations must be enforced through
// preconditions or I need to re-think how to set up the interface for this
// strategy object.
//
template<class Scalar>
class ImplicitBDFStepperStepControl 
  : virtual public ErrWtVecCalcAcceptingStepControlStrategyBase<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** \name Overridden from StepControlStrategyBase */
    //@{
    /** \brief . */
    void setRequestedStepSize(const StepperBase<Scalar>& stepper,
      const Scalar& stepSize, const StepSizeType& stepSizeType);

    /** \brief . */
    void nextStepSize(const StepperBase<Scalar>& stepper, Scalar* stepSize,
      StepSizeType* stepSizeType, int* order);

    /** \brief . */
    void setCorrection(
         const StepperBase<Scalar>& stepper
        ,const RCP<const Thyra::VectorBase<Scalar> >& soln
        ,const RCP<const Thyra::VectorBase<Scalar> >& ee
        ,int solveStatus
        ); 

    /** \brief . */
    bool acceptStep(const StepperBase<Scalar>& stepper, Scalar* LETValue);

    /** \brief . */
    void completeStep(const StepperBase<Scalar>& stepper);
    
    /** \brief . */
    AttemptedStepStatusFlag rejectStep(const StepperBase<Scalar>& stepper);

    /** \brief . */
    StepControlStrategyState getCurrentState();

    /** \brief . */
    int getMinOrder() const;

    /** \brief . */
    int getMaxOrder() const;

    /** \brief . */
    void setStepControlData(const StepperBase<Scalar>& stepper);

    /** \brief . */
    bool supportsCloning() const;

    /** \brief . */
    RCP<StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const;

    //@}

    /** \name Overridden from ErrWtVecCalcAcceptingStepControlStrategyBase */
    //@{

    /** \brief . */
    void setErrWtVecCalc(const RCP<ErrWtVecCalcBase<Scalar> >& errWtVecCalc);

    /** \brief . */
    RCP<const ErrWtVecCalcBase<Scalar> > getErrWtVecCalc() const;

    //@}

    ImplicitBDFStepperStepControl();

    /** \name Overridden from Teuchos::Describable */
    //@{
    /** \brief . */
    void describe(
      Teuchos::FancyOStream &out,
      const Teuchos::EVerbosityLevel verbLevel
      ) const;
    //@}

    /** \name Overridden from ParameterListAcceptor */
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

    /** \brief . */
    void initialize(const StepperBase<Scalar>& stepper);
    

  private:

    // Private data members

    void defaultInitializeAllData_();

    StepControlStrategyState stepControlState_;

    RCP<ErrWtVecCalcBase<Scalar> > errWtVecCalc_;

    Scalar hh_;
    int numberOfSteps_;
    StepSizeType stepSizeType_;
    int minOrder_;
    int maxOrder_;
    int  nef_;
    bool midStep_; // true after setStepSize and nextStepSize and nextStepOrder and setCorrection and acceptStep
                   // false after rejectStep and completeStep
    RCP<const Thyra::VectorBase<Scalar> > ee_; // Newton update
    RCP<Thyra::VectorBase<Scalar> > errWtVec_; // error weight vector
    RCP<Thyra::VectorBase<Scalar> > delta_;

    bool checkReduceOrderCalled_;

    RCP<Teuchos::ParameterList> parameterList_;

    Scalar time_;

    ScalarMag relErrTol_; // relative error tolerance
    ScalarMag absErrTol_; // absolute error tolerance
    Scalar usedStep_;
    int currentOrder_; // Current order of integration
    int usedOrder_;    // order used in current step (used after currentOrder is updated)
    int nscsco_;
    Array<Scalar> alpha_;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
    // note:   $h_n$ = current step size, n = current time step
    Array<Scalar> sigma_;    // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    Array<Scalar> gamma_;    // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
    // calculate time derivative of history array for predictor 
    Array<Scalar> beta_;     // coefficients used to evaluate predictor from history array
    Array<Scalar> psi_;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
    // compute $\beta_j(n)$
    Scalar alpha_s_;    // $\alpha_s$ fixed-leading coefficient of this BDF method
    Scalar alpha_0_;     // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    Scalar cj_ ;        // $-\alpha_s/h_n$ coefficient used in local error test
    Scalar ck_ ;        // local error coefficient gamma[0] = 0; 
    Scalar ck_enorm_;   // ck * enorm

    bool constantStepSize_;
    Scalar Ek_;
    Scalar Ekm1_;
    Scalar Ekm2_;
    Scalar Ekp1_;
    Scalar Est_;
    Scalar Tk_;
    Scalar Tkm1_;
    Scalar Tkm2_;
    Scalar Tkp1_;
    int newOrder_;
    int oldOrder_;
    bool initialPhase_;
    Scalar stopTime_;
    
    // Magic Numbers:
    Scalar h0_safety_;
    Scalar h0_max_factor_;
    Scalar h_phase0_incr_;
    Scalar h_max_inv_;
    Scalar Tkm1_Tk_safety_;
    Scalar Tkp1_Tk_safety_;
    Scalar r_factor_;
    Scalar r_safety_;
    Scalar r_fudge_;
    Scalar r_min_;
    Scalar r_max_;
    Scalar r_hincr_test_;
    Scalar r_hincr_;
    int    max_LET_fail_;
    Scalar minTimeStep_;
    Scalar maxTimeStep_;

    int newtonConvergenceStatus_;
    bool failStepIfNonlinearSolveFails_;

    // Private member functions

    Scalar wRMSNorm_(
        const Thyra::VectorBase<Scalar>& weight, 
        const Thyra::VectorBase<Scalar>& vector
        ) const;

    Scalar checkReduceOrder_(const StepperBase<Scalar>& stepper);
  
    void getFirstTimeStep_(const StepperBase<Scalar>& stepper);

    void setStepControlState_(StepControlStrategyState state);

    void updateCoeffs_();

    void setDefaultMagicNumbers_(Teuchos::ParameterList &magicNumberList);

};

} // namespace Rythmos

#endif // Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DECL_H

