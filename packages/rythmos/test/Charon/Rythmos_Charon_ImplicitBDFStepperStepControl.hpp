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

#include "Rythmos_Types.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
//#include "Rythmos_ErrWtVecCalcBase.hpp"

namespace RythmosCharon {

using Teuchos::RCP;

template<class Scalar>
class CharonImplicitBDFStepperStepControl
  : virtual public Rythmos::StepControlStrategyBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  CharonImplicitBDFStepperStepControl() { }
  virtual ~CharonImplicitBDFStepperStepControl() { }
  void setErrWtVecCalc(const RCP<Rythmos::ErrWtVecCalcBase<Scalar> >& errWtVecCalc)
  { }
  RCP<const Rythmos::ErrWtVecCalcBase<Scalar> > getErrWtVecCalc() const
  {
    return Teuchos::null;
  }
  void initialize(const Rythmos::StepperBase<Scalar>& stepper)
  { }
  void setRequestedStepSize(
      const Rythmos::StepperBase<Scalar>& stepper
      , const Scalar& stepSize
      , const Rythmos::StepSizeType& stepSizeType
      )
  { }
  void nextStepSize(
      const Rythmos::StepperBase<Scalar>& stepper
      , Scalar* stepSize
      , Rythmos::StepSizeType* stepSizeType
      , int* order 
      )
  { }
  void setCorrection(
      const Rythmos::StepperBase<Scalar>& stepper
      , const RCP<const Thyra::VectorBase<Scalar> >& soln
      , const RCP<const Thyra::VectorBase<Scalar> >& ee
      , int solveStatus
      )
  { }
  bool acceptStep(
      const Rythmos::StepperBase<Scalar>& stepper
      ,Scalar* LETValue
      )
  { 
    return false;
  }
  void completeStep(
      const Rythmos::StepperBase<Scalar>& stepper
      )
  { }
  Rythmos::AttemptedStepStatusFlag rejectStep(
      const Rythmos::StepperBase<Scalar>& stepper
      )
  { 
    return Rythmos::REP_ERR_FAIL;
  }
  Rythmos::StepControlStrategyState getCurrentState()
  {
    return Rythmos::UNINITIALIZED;
  }
  int getMaxOrder() const
  { 
    return 0;
  }
  void setStepControlData(const Rythmos::StepperBase<Scalar>& stepper)
  { }
  bool supportsCloning() const
  { 
    return false;
  }
  RCP<Rythmos::StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const
  {
    return Teuchos::null;
  }
  // Overridden from Teuchos::ParameterListAcceptorDefaultBase
  void setParameterList( RCP<Teuchos::ParameterList> const& paramList )
  { }
  RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }
};

} // namespace RythmosCharon


