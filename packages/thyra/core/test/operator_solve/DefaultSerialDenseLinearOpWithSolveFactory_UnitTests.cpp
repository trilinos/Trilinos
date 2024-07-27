// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"

#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"

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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSerialDenseLinearOpWithSolveFactory,
  LinearOpWithSolveFactoryExamples, Scalar )
{

  const Ordinal dim = 4;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<const MultiVectorBase<Scalar> > A =
    createNonsingularMultiVector(vs);

  const RCP<const LinearOpWithSolveFactoryBase<Scalar> > lowsFactory =
    defaultSerialDenseLinearOpWithSolveFactory<Scalar>();

  nonExternallyPreconditionedLinearSolveUseCases<Scalar>(*A, *lowsFactory,
    true, out);

  // Just get it to compile!
  bool runPrecTest = false;
  if (runPrecTest) {
    const Ptr<const PreconditionerFactoryBase<Scalar> > pfb_ptr;
    externallyPreconditionedLinearSolveUseCases<Scalar>(
      *A, *lowsFactory, *pfb_ptr, false, true, out);
  }

}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSerialDenseLinearOpWithSolveFactory,
  LinearOpWithSolveFactoryExamples )


} // namespace Thyra
