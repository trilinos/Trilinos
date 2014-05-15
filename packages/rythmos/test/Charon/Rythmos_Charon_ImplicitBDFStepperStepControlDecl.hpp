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

#ifndef RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_DECL_HPP
#define RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_DECL_HPP

#include "Rythmos_StepControlStrategyBase.hpp"

namespace Rythmos {
  template<class Scalar> class ErrWtVecCalcBase;
} // namespace Rythmos

namespace RythmosCharon {

template<class Scalar>
class CharonImplicitBDFStepperStepControl
  : virtual public Rythmos::StepControlStrategyBase<Scalar>
{
  public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  CharonImplicitBDFStepperStepControl();
  virtual ~CharonImplicitBDFStepperStepControl();
  void setErrWtVecCalc(const Teuchos::RCP<Rythmos::ErrWtVecCalcBase<Scalar> >& errWtVecCalc);
  Teuchos::RCP<const Rythmos::ErrWtVecCalcBase<Scalar> > getErrWtVecCalc() const;
  void initialize(const Rythmos::StepperBase<Scalar>& stepper);
  void setRequestedStepSize(
      const Rythmos::StepperBase<Scalar>& stepper
      , const Scalar& stepSize
      , const Rythmos::StepSizeType& stepSizeType
      );
  void nextStepSize(
      const Rythmos::StepperBase<Scalar>& stepper
      , Scalar* stepSize
      , Rythmos::StepSizeType* stepSizeType
      , int* order 
      );
  void setCorrection(
      const Rythmos::StepperBase<Scalar>& stepper
      , const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& soln
      , const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& ee
      , int solveStatus
      );
  bool acceptStep(
      const Rythmos::StepperBase<Scalar>& stepper
      ,Scalar* LETValue
      );
  void completeStep(
      const Rythmos::StepperBase<Scalar>& stepper
      );
  Rythmos::AttemptedStepStatusFlag rejectStep(
      const Rythmos::StepperBase<Scalar>& stepper
      );
  Rythmos::StepControlStrategyState getCurrentState();
  int getMaxOrder() const;
  void setStepControlData(const Rythmos::StepperBase<Scalar>& stepper);
  bool supportsCloning() const;
  Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const;
  // Overridden from Teuchos::ParameterListAcceptor
  void setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList );
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  private:
    Rythmos::StepControlStrategyState stepControlState_;
    Teuchos::RCP<Rythmos::ErrWtVecCalcBase<Scalar> > errWtVecCalc_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > errWtVec_; 
    Teuchos::RCP<Teuchos::ParameterList> parameterList_;
    ScalarMag relErrTol_; 
    ScalarMag absErrTol_; 
};

} // namespace RythmosCharon

#endif // RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_DECL_HPP

