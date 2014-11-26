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


#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TestForException.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::null;
using Teuchos::is_null;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::fancyOStream;


//
// Unit Tests
//


const Ordinal m = 5;
const Ordinal n = 3;


TEUCHOS_UNIT_TEST( DefaultAddedLinearOp, defaultConstruct )
{

  typedef double Scalar;

  const RCP<DefaultAddedLinearOp<Scalar> > M = defaultAddedLinearOp<Scalar>();

  TEST_ASSERT(is_null(M->range()));
  TEST_ASSERT(is_null(M->domain()));
  
  const std::string M_description = M->description();
  TEST_EQUALITY_CONST(M_description, "Thyra::DefaultAddedLinearOp<double>{numOps=0,rangeDim=0,domainDim=0}");

  {
    std::ostringstream describe_msg;
    describe_msg << "'";
    M->describe(*fancyOStream(rcpFromRef(describe_msg)), Teuchos::VERB_LOW);
    describe_msg << "'";
    
    std::ostringstream expected_msg;
    expected_msg
      << "' " << M_description << "\n'";
    
    TEST_EQUALITY_CONST( describe_msg.str(), expected_msg.str() );
  }


  {
    std::ostringstream describe_msg;
    describe_msg << "'";
    M->describe(*fancyOStream(rcpFromRef(describe_msg)), Teuchos::VERB_EXTREME);
    describe_msg << "'";
    
    std::ostringstream expected_msg;
    expected_msg
      << "' " << M_description << "\n"
      << "  Constituent LinearOpBase objects for M = Op[0]*...*Op[numOps-1]:\n'";

    TEST_EQUALITY_CONST( describe_msg.str(), expected_msg.str() );
  }
  
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( DefaultAddedLinearOp, addConst )
{

  typedef double Scalar;

  const RCP<const VectorSpaceBase<Scalar> > vs = defaultSpmdVectorSpace<Scalar>(m);
  const RCP<const MultiVectorBase<Scalar> > A = createMembers(vs, n, "A");

  TEST_THROW(
    const RCP<const LinearOpBase<Scalar> > M = add<Scalar>(A, adjoint<Scalar>(A)),
    Exceptions::IncompatibleVectorSpaces
    );
  
}


#endif // TEUCHOS_DEBUG


} // namespace Thyra
