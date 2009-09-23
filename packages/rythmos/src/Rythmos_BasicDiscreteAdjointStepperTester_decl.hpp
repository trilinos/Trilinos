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

#ifndef Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DECL_H
#define Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DECL_H


#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Rythmos {


template<class Scalar> class BasicDiscreteAdjointStepperTester;


namespace BasicDiscreteAdjointStepperTesterUtils {

const std::string ErrorTol_name = "Error Tol";
const double ErrorTol_default = 1e-6;

} // namespace BasicDiscreteAdjointStepperTesterUtils


/** \brief Nonmember constructor.
 *
 * \relates BasicDiscreteAdjointStepperTester
 */
template<class Scalar>
RCP<BasicDiscreteAdjointStepperTester<Scalar> >
basicDiscreteAdjointStepperTester();


/** \brief Nonmember constructor.
 *
 * \relates BasicDiscreteAdjointStepperTester
 */
template<class Scalar>
RCP<BasicDiscreteAdjointStepperTester<Scalar> >
basicDiscreteAdjointStepperTester(const RCP<ParameterList> &paramList);


/** \brief Concrete testing class for basic adjoint calculation.
 *
 * This testing class performs the most basic test of an adjoint computation
 * for a nonlinear model that you can possibly check.  The basic response
 * problem is:

 \verbatim

    f(x_dot_, x_, t) = 0, for t <: [t_0, t_f]
              x(t_0) = x_init + B*p
          x_dot(t_0) = x_dot_int


    d_hat(p) = h(x(t_f,p)) = 0.5 * x^T * x

 \endverbatim

 * This formulation assumes that the mass matrix d(f)/d(x_dot) is full rank
 * which will be needed to compute the adjoint initial condition..
 *
 * The intial condition vectors x_init and x_dot_init are taken from the
 * orginal forward problem's intial condition as is t_0.  The time t_f is
 * taken from an initalized integrator.
 *
 * The multi-vector B can be chosen by the user or can be computed
 * automatically internally.  If B is not choses by the user, it will be
 * computed automatically as a single column with random numbers.
 *
 * The forward sensitivity equations (with S = d(x)/d(p)) that are solved with
 * the reduced response sensitivity are then:

 \verbatim

   d(f)/d(x_dot) * S_dot + d(f)/d(x) * S = 0, for t <: [t_0, t_f]
                                  S(t_0) = B
                              S_dot(t_0) = 0

   d(d_hat)/d(p)^T = S^T * x, at t = t_f

 \endverbatim

 * The adjoint equations that are solved for the reduced sensitivity are then:

 \verbatim

   d(f)/d(x_dot)^T * lambda_dot - d(f)/d(x)^T * lambda = 0, for t <: [t_0, t_f]
                              d(f)/d(x_dot)^T * lambda = x, at t = t_f

   d(d_hat)/d(p)^T = B^T * d(f)/d(x_dot)^T * lambda, at t = t_0

 \endverbatim

 * Note that if d(f)/d(x_dot) is full rank, then the adjoint initial condition
 * at t_f reduces to:

 \verbatim

   lambda(t_f) = d(f)/d(x_dot)^{-T} * x(t_f)

 \endverbatim

 * which is the form of the initial condition used in this test (nice and
 * simple).
 *
 * NOTE: However, if this is a general DAE where d(f)/d(x_dot) is rank
 * deficient, then the adjoint initial value calcuation at t_f gets more
 * complicated and this testing class can not handle those cases.
 */
template<class Scalar> 
class BasicDiscreteAdjointStepperTester
  : virtual public Teuchos::VerboseObject<BasicDiscreteAdjointStepperTester<Scalar> >,
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

  /** \brief Test the the AdjointStepper object for a given forward
   * simulation.
   *
   * \param adjointModel [in] The basic adjoint model ready to be used to
   * integrate the adjoint.  On output, this stepper will have been used to
   * integate the adjoint.
   *
   * \param forwardIntegrator [in/out] The basic forward integrator ready to
   * integrate the forward problem.  This integrator algorithm will be cloned
   * to integrate the forward sensitivities and the adjoint.  This integator
   * should be set up to take fixed time steps.  There is no need for adaptive
   * time steps for a test like this.  On output, this integrator will have
   * been run to the output time.
   *
   * NOTE: This function is declared non-const since it can technically change
   * the parameter list as the fuctions are performed.
   */
  bool testAdjointStepper(
    Thyra::ModelEvaluator<Scalar>& adjointModel,
    const Ptr<IntegratorBase<Scalar> >& forwardIntegrator
    );

  //@}

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED

  /** \name Public friend functions */
  //@{

  ///
  friend RCP< BasicDiscreteAdjointStepperTester<Scalar> >
  basicDiscreteAdjointStepperTester<>();

  //@}

#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED
  

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED
private:
#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED

  BasicDiscreteAdjointStepperTester(); // Note defined and not to be called

private:

  ScalarMag errorTol_;

};


} // namespace Rythmos


#endif //Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DECL_H
