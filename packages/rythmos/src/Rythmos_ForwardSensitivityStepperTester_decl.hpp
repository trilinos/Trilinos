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

#ifndef Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DECL_H
#define Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DECL_H


#include "Rythmos_StepperSupportTypes.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Rythmos {


template<class Scalar> class ForwardSensitivityStepperTester;


namespace ForwardSensitivityStepperTesterUtils {

const std::string FdCalc_name = "FD Calc";

const std::string ErrorTol_name = "Error Tol";
const double ErrorTol_default = 1e-6;

} // namespace ForwardSensitivityStepperTesterUtils


/** \brief Nonmember constructor.
 *
 * \relates ForwardSensitivityStepperTester
 */
template<class Scalar>
RCP<ForwardSensitivityStepperTester<Scalar> >
forwardSensitivityStepperTester();


/** \brief Nonmember constructor.
 *
 * \relates ForwardSensitivityStepperTester
 */
template<class Scalar>
RCP<ForwardSensitivityStepperTester<Scalar> >
forwardSensitivityStepperTester(const RCP<ParameterList> &paramList);


/** \brief Concrete testing class for forward sensitivities.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar> 
class ForwardSensitivityStepperTester
  : virtual public Teuchos::VerboseObject<ForwardSensitivityStepperTester<Scalar> >,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Testing functions */
  //@{

  /** \brief Test a forward sensitivity stepper.
   */
  bool testForwardSens(
    const Ptr<IntegratorBase<Scalar> > &fwdSensIntegrator
    );

  //@}

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED

  /** \name Public friend functions */
  //@{

  ///
  friend RCP< ForwardSensitivityStepperTester<Scalar> >
  forwardSensitivityStepperTester<>();

  //@}

#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED
  

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED
private:
#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED

  ForwardSensitivityStepperTester(); // Note defined and not to be called

private:

  ScalarMag errorTol_;

};


} // namespace Rythmos


#endif //Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DECL_H
