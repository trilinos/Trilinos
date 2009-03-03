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

#ifndef Rythmos_STEPPER_BUILDER_NEWNEW_H
#define Rythmos_STEPPER_BUILDER_NEWNEW_H

#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"

#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"
#include "Rythmos_ExplicitTaylorPolynomialStepper.hpp"

namespace Rythmos {

template<class Scalar>
  class StepperBuilderNewNew : virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief . */
  StepperBuilderNewNew();

  /** \brief . */
  ~StepperBuilderNewNew();

  /** \brief Set a new Stepper factory object. */
  void setStepperFactory(
    const RCP<const Teuchos::AbstractFactory<StepperBase<Scalar> > > &stepperFactory,
    const std::string &stepperFactoryName
    );
  
  /** \brief Get the name of the Stepper that will be created
   * on the next call to <tt>this->create()</tt>.
   */
  std::string getStepperName() const;

  /** \brief . */
  RCP<StepperBase<Scalar> > create(
    const std::string &stepperName = ""
    ) const;

  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
  
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getParameterList() const;

  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;
 
  //@}
  

private:

  // //////////////////////////////////////
  // Private data members

  Teuchos::ObjectBuilder<StepperBase<Scalar> > builder_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults_();

};

} // namespace Rythmos

namespace Rythmos {


// Nonmember constructor
template<class Scalar>
RCP<StepperBuilderNewNew<Scalar> > stepperBuilderNewNew()
{
  RCP<StepperBuilderNewNew<Scalar> > sb = rcp(new StepperBuilderNewNew<Scalar> );
  return sb;
}

template<class Scalar>
StepperBuilderNewNew<Scalar>::StepperBuilderNewNew()
{
  this->initializeDefaults_();
}

template<class Scalar>
StepperBuilderNewNew<Scalar>::~StepperBuilderNewNew()
{
}

template<class Scalar>
void StepperBuilderNewNew<Scalar>::setStepperFactory(
  const RCP<const Teuchos::AbstractFactory<StepperBase<Scalar> > > &stepperFactory,
  const std::string &stepperName
  )
{
  builder_.setObjectFactory(stepperFactory, stepperName);
}



template<class Scalar>
std::string
StepperBuilderNewNew<Scalar>::getStepperName() const
{
  return builder_.getObjectName();
}

template<class Scalar>
void StepperBuilderNewNew<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  builder_.setParameterList(paramList);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
StepperBuilderNewNew<Scalar>::getNonconstParameterList()
{
  return builder_.getNonconstParameterList();
}


template<class Scalar>
RCP<Teuchos::ParameterList>
StepperBuilderNewNew<Scalar>::unsetParameterList()
{
  return builder_.unsetParameterList();
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
StepperBuilderNewNew<Scalar>::getParameterList() const
{
  return builder_.getParameterList();
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
StepperBuilderNewNew<Scalar>::getValidParameters() const
{
  return builder_.getValidParameters();
}

template<class Scalar>
RCP<StepperBase<Scalar> >
StepperBuilderNewNew<Scalar>::create(
  const std::string &stepperName
  ) const
{
  return builder_.create(stepperName);
}

template<class Scalar>
void StepperBuilderNewNew<Scalar>::initializeDefaults_()
{

  using Teuchos::abstractFactoryStd;

  builder_.setObjectName("Rythmos::Stepper");
  builder_.setObjectTypeName("Stepper Type");

  //
  // Steppers
  //

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, ImplicitBDFStepper<Scalar> >(),
      "Implicit BDF"
      );

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, ForwardEulerStepper<Scalar> >(),
      "Forward Euler"
      );

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, ExplicitRKStepper<Scalar> >(),
      "Explicit RK"
      );

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, ImplicitRKStepper<Scalar> >(),
      "Implicit RK"
      );

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, ExplicitTaylorPolynomialStepper<Scalar> >(),
      "Explicit Taylor Polynomial"
      );

  builder_.setObjectFactory(
      abstractFactoryStd< StepperBase<Scalar>, BackwardEulerStepper<Scalar> >(),
      "Backward Euler"
      );

  // Note: Above, the last Stepper object set will be the default unless we are
  // on multiple processors!

}


} // namespace Rythmos

#endif //Rythmos_STEPPER_BUILDER_NEWNEW_H

