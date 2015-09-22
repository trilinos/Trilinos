/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"


namespace Thyra {


using Teuchos::null;
using Teuchos::as;
using Teuchos::getParameter;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DirectionalFiniteDiffCalculator, defaultConstruct, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  const RCP<DirectionalFiniteDiffCalculator<Scalar> > fdCalc =
    directionalFiniteDiffCalculator<Scalar>();

  TEST_EQUALITY_CONST( fdCalc->fd_method_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO );
  TEST_EQUALITY_CONST( fdCalc->fd_step_select_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size(), as<ScalarMag>(-1.0) );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size_min(), as<ScalarMag>(-1.0) );

  const RCP<const ParameterList> validPL = fdCalc->getValidParameters();
  TEST_ASSERT(nonnull(validPL));
  TEST_EQUALITY_CONST(getParameter<std::string>(*validPL, "FD Method"),
    "order-one");
  TEST_EQUALITY_CONST(getParameter<std::string>(*validPL, "FD Step Select Type"),
    "Absolute");
  TEST_EQUALITY_CONST(getParameter<double>(*validPL, "FD Step Length"),
    as<double>(-1.0));

  // ToDo: Add more tests!

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( DirectionalFiniteDiffCalculator, defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DirectionalFiniteDiffCalculator, plConstruct, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  const RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("FD Method", "order-two-central");
  pl->set("FD Step Select Type", "Relative");
  pl->set("FD Step Length", as<double>(1e-5));

  const RCP<DirectionalFiniteDiffCalculator<Scalar> > fdCalc =
    directionalFiniteDiffCalculator<Scalar>(pl);

  TEST_EQUALITY_CONST( fdCalc->fd_method_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL );
  TEST_EQUALITY_CONST( fdCalc->fd_step_select_type(),
    DirectionalFiniteDiffCalculatorTypes::FD_STEP_RELATIVE );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size(), as<ScalarMag>(1e-5) );
  TEST_EQUALITY_CONST( fdCalc->fd_step_size_min(), as<ScalarMag>(-1.0) );

  // ToDo: Add more tests!

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( DirectionalFiniteDiffCalculator, plConstruct )


} // namespace Thyra
