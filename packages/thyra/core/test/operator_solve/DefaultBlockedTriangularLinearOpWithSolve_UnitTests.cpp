// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
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
using Teuchos::inOutArg;


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  defaultConstruct, Scalar )
{
  const RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows =
    defaultBlockedTriangularLinearOpWithSolve<Scalar>();
  TEST_ASSERT(nonnull(dbtlows));
  TEST_EQUALITY_CONST(dbtlows->range(), null);
  TEST_EQUALITY_CONST(dbtlows->domain(), null);
  out << "dbtlows = " << *dbtlows;
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  basic, Scalar )
{

  // typedef Teuchos::ScalarTraits<Scalar> ST; // unused

  const Ordinal dim = 4;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<const MultiVectorBase<Scalar> > M =
    createNonsingularMultiVector(vs);

  const RCP<const LinearOpWithSolveBase<Scalar> > lows =
    linearOpWithSolve<Scalar>(
      *defaultSerialDenseLinearOpWithSolveFactory<Scalar>(), M );

  const int numBlocks = 3;

  const RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows =
    defaultBlockedTriangularLinearOpWithSolve<Scalar>();
  dbtlows->beginBlockFill(numBlocks, numBlocks);

  for (int block_i = 0; block_i < numBlocks; ++block_i) {
    dbtlows->setLOWSBlock(block_i, block_i, lows);
  }

  dbtlows->endBlockFill();

  out << "dbtlows = " << *dbtlows;

  Thyra::LinearOpTester<Scalar> linearOpTester;
  const bool checkOpResult = linearOpTester.check(*dbtlows, inOutArg(out));
  TEST_ASSERT(checkOpResult);

  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  linearOpWithSolveTester.check_adjoint_default(true);
  const bool checkSolveResult = linearOpWithSolveTester.check(*dbtlows, &out);
  TEST_ASSERT(checkSolveResult);

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  basic )


} // namespace Thyra
