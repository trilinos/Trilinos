// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

using Teuchos::fancyOStream;
using Teuchos::is_null;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

//
// Unit Tests
//

TEUCHOS_UNIT_TEST(DefaultAddedLinearOp, defaultConstruct) {
  typedef double Scalar;

  const RCP<DefaultAddedLinearOp<Scalar> > M = defaultAddedLinearOp<Scalar>();

  TEST_ASSERT(is_null(M->range()));
  TEST_ASSERT(is_null(M->domain()));
  TEST_ASSERT(isFullyUninitialized(*M));
  TEST_ASSERT(!isPartiallyInitialized(*M));
  TEST_ASSERT(!isFullyInitialized(*M));

#if defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)

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

    TEST_EQUALITY_CONST(describe_msg.str(), expected_msg.str());
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

    TEST_EQUALITY_CONST(describe_msg.str(), expected_msg.str());
  }

#endif  // defined(HAVE_GCC_ABI_DEMANGLE) && defined(HAVE_TEUCHOS_DEMANGLE)
}

#ifdef TEUCHOS_DEBUG

const Ordinal m = 5;
const Ordinal n = 3;

TEUCHOS_UNIT_TEST(DefaultAddedLinearOp, addConst) {
  typedef double Scalar;

  const RCP<const VectorSpaceBase<Scalar> > vs = defaultSpmdVectorSpace<Scalar>(m);
  const RCP<const MultiVectorBase<Scalar> > A  = createMembers(vs, n, "A");

  TEST_THROW(
      const RCP<const LinearOpBase<Scalar> > M = add<Scalar>(A, adjoint<Scalar>(A)),
      Exceptions::IncompatibleVectorSpaces);
}

#endif  // TEUCHOS_DEBUG

}  // namespace Thyra
