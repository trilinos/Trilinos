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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DECL_H
#define Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DECL_H

#include "Rythmos_StepControlStrategyBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Rythmos {

// Step Control Strategy object for FirstOrderErrorStepControlStrategy
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
/** \brief Step Control Strategy for first-order time integration.
 *
 * This step-control strategy assumes a first-order predictor (Forward
 * Euler) for a first-order time integrator (Backward Euler) and tries
 * to maintain a constant error based on relative and absolute tolerances
 * by adjusting the step size.  Although specifically built from Forward
 * and Backward Euler, one could use this for other first-order schemes.
 * See Grasho and Sani, "Incompressible Flow and the Finite Element Method",
 * Volume One, p. 268.
 *
 * \section Rythmos_FirstOrderErrorStepControlStrategy_Theory_sec Theory
 *
 * To control the error through step-size selection, Forward Euler is
 * used as a predictor
 * \f[ x^P_n = x_{n-1} + \Delta t f(x_{n-1}, t_{n-1}) \f]
 * with the local truncation error (LTE) of
 * \f{eqnarray*}{
 *   d_n &=& x^P_n - x(t_n)
 *       &=& -\frac{\Delta t^2}{2} \ddot{x}_n + O(\Delta t^3)
 * \f}
 * Similarly for Backward Euler
 * \f[ x_n = x_{n-1} + \Delta t f(x_n, t_n) \f]
 * with the local truncation error (LTE) of
 * \f{eqnarray*}{
 *   d_n &=& x_n - x(t_n)
 *       &=& \frac{\Delta t^2}{2} \ddot{x}_n + O(\Delta t^3)
 * \f}
 * Using the Forward Euler as an estimate of the LTE for the Backward Euler,
 * eliminating \f$ \ddot{x}_n \f$, and forming \f$ d_n = x_n - x(t_n) \f$
 * in terms of \f$ (x_n - x^P_n) \f$, we get
 * \f[ d_n = x_n - x(t_n) =  (x_n - x^P_n)/2 \f]
 * which gives us an estimate of the LTE, \f$ d_n \f$ in term of the
 * difference between the Backward Euler solution, \f$ x_n \f$, and the
 * predicted Forward Euler solution, \f$ x^P_n \f$.
 *
 * To select the next step size, we would like the LTE to be meet a
 * user-supplied tolerance,
 * \f[ d_{n+1} \leq \epsilon_r ||x_n|| + \epsilon_a \f]
 * where \f$ \epsilon_r \f$ and \f$ \epsilon_a \f$ are the relative and
 * absolute tolerances, and \f$ ||\cdot|| \f$ is the norm.  Normalizing
 * this with respect to the current step, which allows us the relate the
 * current step size to the next step size, and using the LTE, we find
 * \f[ \frac{d_{n+1}}{d_n} = \frac{\epsilon_r ||x_n|| + \epsilon_a}{d_n}
 *     = \frac{\Delta t^2_{n+1} \ddot{x}_{n+1}}{\Delta t^2_n \ddot{x}_n}
 *     = \left(\frac{\Delta t_{n+1}}{\Delta t_n}\right)^2 + O(\Delta t) \f]
 * leading to the value for the next step size
 * \f[ \Delta t_{n+1} =
 *     \Delta t_n \left(\frac{\epsilon_r ||x_n|| + \epsilon_a}{d_n}\right)^{1/2}
 * \f]
 *
 */
template<class Scalar>
class FirstOrderErrorStepControlStrategy
  : virtual public StepControlStrategyBase<Scalar>
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
    int getMaxOrder() const;

    /** \brief . */
    void setStepControlData(const StepperBase<Scalar>& stepper);

    /** \brief . */
    bool supportsCloning() const;

    /** \brief . */
    RCP<StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const;

    //@}

    FirstOrderErrorStepControlStrategy();

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

    void initialize(const StepperBase<Scalar>& stepper);


  private:

    // Private data members

    void defaultInitializeAllData_();

    StepControlStrategyState stepControlState_;

    RCP<Teuchos::ParameterList> parameterList_;

    Scalar initialStepSize_;
    Scalar requestedStepSize_;
    Scalar currentStepSize_;
    Scalar nextStepSize_;
    Scalar stepSizeFactor_;
    StepSizeType stepSizeType_;

    Scalar minStepSize_;
    Scalar maxStepSize_;
    Scalar maxStepSizeIncreaseFactor_;
    Scalar minStepSizeDecreaseFactor_;
    int numStepFailures_;
    int maxStepFailures_;
    int maxOrder_;
    Scalar errorRelativeTolerance_;
    Scalar errorAbsoluteTolerance_;
    int solveStatus_;

    RCP<const Thyra::VectorBase<Scalar> > x_;
    RCP<const Thyra::VectorBase<Scalar> > dx_;
    RCP<Thyra::VectorBase<Scalar> > errWtVec_;


    static const std::string initialStepSizeName_;
    static const double initialStepSizeDefault_;

    static const std::string minStepSizeName_;
    static const double minStepSizeDefault_;

    static const std::string maxStepSizeName_;
    static const double maxStepSizeDefault_;

    static const std::string maxStepSizeIncreaseFactorName_;
    static const double maxStepSizeIncreaseFactorDefault_;

    static const std::string minStepSizeDecreaseFactorName_;
    static const double minStepSizeDecreaseFactorDefault_;

    static const std::string maxStepFailuresName_;
    static const int maxStepFailuresDefault_;

    static const std::string errorRelativeToleranceName_;
    static const double errorRelativeToleranceDefault_;

    static const std::string errorAbsoluteToleranceName_;
    static const double errorAbsoluteToleranceDefault_;


    // Private member functions

    void setStepControlState_(StepControlStrategyState state);

};

} // namespace Rythmos

#endif // Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DECL_H

