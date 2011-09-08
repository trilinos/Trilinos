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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"

#include "OperatorSolveHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::inOutArg;


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorLinearOpWithSolve, defaultConstruct,
  Scalar )
{
  const RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> > dmvlows =
    multiVectorLinearOpWithSolve<Scalar>();
  TEST_ASSERT(nonnull(dmvlows));
  TEST_EQUALITY_CONST(dmvlows->range(), null);
  TEST_EQUALITY_CONST(dmvlows->domain(), null);
  out << "dmvlows = " << *dmvlows;
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorLinearOpWithSolve,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorLinearOpWithSolve, basic,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  const Ordinal dim = 4;
  const int numBlocks = 3;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<const MultiVectorBase<Scalar> > M =
    createNonsingularMultiVector(vs);

  const RCP<const LinearOpWithSolveFactoryBase<Scalar> > lowsf = 
    defaultSerialDenseLinearOpWithSolveFactory<Scalar>();

  const RCP<LinearOpWithSolveBase<Scalar> > Minv = 
    linearOpWithSolve<Scalar>(*lowsf, M);
      
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > dmvpvs =
    multiVectorProductVectorSpace<Scalar>(vs, numBlocks);

  const RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> > dmvlows =
    multiVectorLinearOpWithSolve<Scalar>(Minv, dmvpvs, dmvpvs);

  TEST_ASSERT(nonnull(dmvlows));
  TEST_EQUALITY(dmvlows->range(), dmvpvs);
  TEST_EQUALITY(dmvlows->domain(), dmvpvs);
  out << "dmvlows = " << *dmvlows;

  Thyra::LinearOpTester<Scalar> linearOpTester;
  const bool checkOpResult = linearOpTester.check(*dmvlows, inOutArg(out));
  TEST_ASSERT(checkOpResult);

  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  const bool checkSolveResult = linearOpWithSolveTester.check(*dmvlows, &out);
  TEST_ASSERT(checkSolveResult);

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorLinearOpWithSolve,
  basic )


} // namespace Thyra
