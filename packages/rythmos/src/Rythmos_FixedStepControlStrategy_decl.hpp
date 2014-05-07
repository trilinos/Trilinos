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

#ifndef Rythmos_FIXED_STEP_CONTROL_STRATEGY_DECL_H
#define Rythmos_FIXED_STEP_CONTROL_STRATEGY_DECL_H

#include "Rythmos_StepControlStrategyBase.hpp"

namespace Rythmos {

// Step Control Strategy object for FixedStepControlStrategy
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
template<class Scalar>
class FixedStepControlStrategy
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

    FixedStepControlStrategy();

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

    StepControlStrategyState stepControlState_;
    RCP<Teuchos::ParameterList> parameterList_;

    // Private member functions

    void setStepControlState_(StepControlStrategyState state);

};

} // namespace Rythmos

#endif // Rythmos_FIXED_STEP_CONTROL_STRATEGY_DECL_H

