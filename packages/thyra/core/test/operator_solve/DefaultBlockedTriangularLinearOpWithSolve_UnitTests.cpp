// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultBlockedTriangularLinearOpWithSolve.hpp"
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
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


/** \brief Block triangular operator M with numBlocks x numBlocks blocks of
 * dimension 4: dense nonsingular diagonal blocks and random blocks below
 * (lower) or above (!lower) the diagonal, and the
 * DefaultBlockedTriangularLinearOpWithSolve of M (dense solves of the
 * diagonal blocks).
 */
template<class Scalar>
void createBlockTriangular(
  const bool lower, const int numBlocks,
  RCP<const LinearOpBase<Scalar> > &M,
  RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > &dbtlows
  )
{
  const Ordinal dim = 4;
  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<DefaultBlockedLinearOp<Scalar> > blo =
    defaultBlockedLinearOp<Scalar>();
  blo->beginBlockFill(numBlocks, numBlocks);
  dbtlows = defaultBlockedTriangularLinearOpWithSolve<Scalar>();
  dbtlows->beginBlockFill(numBlocks, numBlocks);

  seed_randomize<Scalar>(23);
  for (int i = 0; i < numBlocks; ++i) {
    const RCP<const MultiVectorBase<Scalar> > M_ii =
      createNonsingularMultiVector(vs);
    blo->setBlock(i, i, M_ii);
    dbtlows->setLOWSBlock(i, i,
      linearOpWithSolve<Scalar>(
        *defaultSerialDenseLinearOpWithSolveFactory<Scalar>(), M_ii));
    for (int j = 0; j < numBlocks; ++j) {
      if (lower ? j < i : j > i) {
        const RCP<MultiVectorBase<Scalar> > M_ij = createMembers(vs, vs);
        randomize(as<Scalar>(-1.0), as<Scalar>(1.0), M_ij.ptr());
        blo->setBlock(i, j, M_ij);
      }
    }
  }
  blo->endBlockFill();
  dbtlows->endBlockFill();
  dbtlows->setBlocks(blo);
  M = blo;
}


/** \brief Check that the apply of the DefaultBlockedTriangularLinearOpWithSolve
 * of a block triangular operator M is the one of M, and that its solves are
 * the inverse of M (also for the adjoints).
 */
template<class Scalar>
void checkBlockTriangular(
  const bool lower, Teuchos::FancyOStream &out, bool &success
  )
{
  RCP<const LinearOpBase<Scalar> > M;
  RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows;
  createBlockTriangular<Scalar>(lower, 3, M, dbtlows);

  out << "dbtlows = " << *dbtlows;

  Thyra::LinearOpTester<Scalar> linearOpTester;
  const bool checkOpResult = linearOpTester.check(*dbtlows, inOutArg(out));
  TEST_ASSERT(checkOpResult);
  const bool compareOpResult =
    linearOpTester.compare(*dbtlows, *M, inOutArg(out));
  TEST_ASSERT(compareOpResult);

  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  linearOpWithSolveTester.check_adjoint_default(true);
  const bool checkSolveResult = linearOpWithSolveTester.check(*dbtlows, &out);
  TEST_ASSERT(checkSolveResult);
}


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  blockLowerTriangular, Scalar )
{
  checkBlockTriangular<Scalar>(true, out, success);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  blockLowerTriangular )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  blockUpperTriangular, Scalar )
{
  checkBlockTriangular<Scalar>(false, out, success);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  blockUpperTriangular )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  notBlockTriangular, Scalar )
{
  // Blocks below and above the diagonal are rejected
  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(2);
  const RCP<const MultiVectorBase<Scalar> > M_d =
    createNonsingularMultiVector(vs);
  const RCP<DefaultBlockedLinearOp<Scalar> > blo =
    defaultBlockedLinearOp<Scalar>();
  blo->beginBlockFill(2, 2);
  blo->setBlock(0, 0, M_d);
  blo->setBlock(0, 1, M_d);
  blo->setBlock(1, 0, M_d);
  blo->setBlock(1, 1, M_d);
  blo->endBlockFill();

  const RCP<const LinearOpWithSolveBase<Scalar> > lows =
    linearOpWithSolve<Scalar>(
      *defaultSerialDenseLinearOpWithSolveFactory<Scalar>(), M_d );
  const RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows =
    defaultBlockedTriangularLinearOpWithSolve<Scalar>();
  dbtlows->beginBlockFill(2, 2);
  dbtlows->setLOWSBlock(0, 0, lows);
  dbtlows->setLOWSBlock(1, 1, lows);
  dbtlows->endBlockFill();
  TEST_THROW(dbtlows->setBlocks(blo), std::logic_error);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  notBlockTriangular )


/** \brief The block diagonal of the 3x3 block operator M of
 * createBlockTriangular() as a blocked operator. */
template<class Scalar>
RCP<const LinearOpBase<Scalar> > blockDiagonal(
  const LinearOpBase<Scalar> &M
  )
{
  const PhysicallyBlockedLinearOpBase<Scalar> &blo =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(M);
  const int numBlocks = blo.productRange()->numBlocks();
  const RCP<DefaultBlockedLinearOp<Scalar> > diag =
    defaultBlockedLinearOp<Scalar>();
  diag->beginBlockFill(numBlocks, numBlocks);
  for (int k = 0; k < numBlocks; ++k)
    diag->setBlock(k, k, blo.getBlock(k,k));
  diag->endBlockFill();
  return diag;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  offDiagonalBlockAccess, Scalar )
{
  // The off-diagonal blocks are those of the blocked operator of setBlocks()
  RCP<const LinearOpBase<Scalar> > M;
  RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows;
  createBlockTriangular<Scalar>(true, 3, M, dbtlows);
  const PhysicallyBlockedLinearOpBase<Scalar> &blo =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(*M);

  TEST_ASSERT(dbtlows->blockExists(1, 0));
  TEST_ASSERT(dbtlows->blockExists(2, 1));
  TEST_ASSERT(!dbtlows->blockExists(0, 1));
  TEST_EQUALITY(dbtlows->getBlock(2, 1).get(), blo.getBlock(2, 1).get());
  TEST_ASSERT(is_null(dbtlows->getBlock(1, 2)));
  TEST_EQUALITY(dbtlows->getBlock(1, 1).get(), dbtlows->getLOWSBlock(1, 1).get());
  // The blocks were given as const
  TEST_THROW(dbtlows->getNonconstBlock(1, 0), Teuchos::NonconstAccessError);
  // and are not set in the block fill
  dbtlows->uninitialize();
  dbtlows->beginBlockFill(3, 3);
  TEST_ASSERT(!dbtlows->acceptsBlock(1, 0));
  TEST_THROW(dbtlows->setBlock(1, 0, blo.getBlock(1, 0)), std::logic_error);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  offDiagonalBlockAccess )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  newBlockFillRemovesOffDiagonalBlocks, Scalar )
{
  // A new block fill, after uninitialize() or not, without setBlocks() gives
  // the block diagonal operator
  for (const bool callUninitialize : {true, false}) {
    out << "callUninitialize = " << callUninitialize << "\n";
    RCP<const LinearOpBase<Scalar> > M;
    RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows;
    createBlockTriangular<Scalar>(true, 3, M, dbtlows);
    Teuchos::Array<RCP<const LinearOpWithSolveBase<Scalar> > > diag;
    for (int k = 0; k < 3; ++k)
      diag.push_back(dbtlows->getLOWSBlock(k, k));

    if (callUninitialize)
      dbtlows->uninitialize();
    dbtlows->beginBlockFill(3, 3);
    for (int k = 0; k < 3; ++k)
      dbtlows->setLOWSBlock(k, k, diag[k]);
    dbtlows->endBlockFill();

    TEST_ASSERT(!dbtlows->blockExists(1, 0));
    Thyra::LinearOpTester<Scalar> linearOpTester;
    const bool compareOpResult =
      linearOpTester.compare(*dbtlows, *blockDiagonal(*M), inOutArg(out));
    TEST_ASSERT(compareOpResult);
  }
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  newBlockFillRemovesOffDiagonalBlocks )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolve,
  rejectedBlocks, Scalar )
{
  // Blocked operators with another number of blocks, or not block triangular,
  // are rejected, and leave the operator (its apply and its solves) unchanged
  RCP<const LinearOpBase<Scalar> > M, M2, M3;
  RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows, dbtlows2, dbtlows3;
  createBlockTriangular<Scalar>(true, 3, M, dbtlows);
  createBlockTriangular<Scalar>(true, 2, M2, dbtlows2);
  const RCP<const PhysicallyBlockedLinearOpBase<Scalar> > blo2 =
    Teuchos::rcp_dynamic_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(M2, true);
  TEST_THROW(dbtlows->setBlocks(blo2), std::logic_error);

  // An upper triangular operator with the lower triangle of M: not triangular
  createBlockTriangular<Scalar>(false, 3, M3, dbtlows3);
  const PhysicallyBlockedLinearOpBase<Scalar> &blo =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(*M);
  const PhysicallyBlockedLinearOpBase<Scalar> &blo3 =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(*M3);
  const RCP<DefaultBlockedLinearOp<Scalar> > full =
    defaultBlockedLinearOp<Scalar>();
  full->beginBlockFill(3, 3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      full->setBlock(i, j, i >= j ? blo.getBlock(i,j) : blo3.getBlock(i,j));
  full->endBlockFill();
  TEST_THROW(dbtlows->setBlocks(full), std::logic_error);

  Thyra::LinearOpTester<Scalar> linearOpTester;
  const bool compareOpResult =
    linearOpTester.compare(*dbtlows, *M, inOutArg(out));
  TEST_ASSERT(compareOpResult);
  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  linearOpWithSolveTester.check_adjoint_default(true);
  const bool checkSolveResult = linearOpWithSolveTester.check(*dbtlows, &out);
  TEST_ASSERT(checkSolveResult);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolve,
  rejectedBlocks )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedTriangularLinearOpWithSolveFactory,
  reinitializeRejected, Scalar )
{
  // The factory solves a block triangular operator, and rejects, when
  // reinitializing, an operator that is not block triangular or has another
  // number of blocks, before changing the diagonal blocks
  RCP<const LinearOpBase<Scalar> > M, M2, M3;
  RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> > dbtlows;
  createBlockTriangular<Scalar>(true, 3, M, dbtlows);
  createBlockTriangular<Scalar>(true, 2, M2, dbtlows);
  createBlockTriangular<Scalar>(false, 3, M3, dbtlows);

  const RCP<LinearOpWithSolveFactoryBase<Scalar> > denseLowsf =
    defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > lowsf =
    defaultBlockedTriangularLinearOpWithSolveFactory<Scalar>(denseLowsf);
  const RCP<LinearOpWithSolveBase<Scalar> > lows =
    linearOpWithSolve<Scalar>(*lowsf, M);

  // Not block triangular, with other diagonal blocks (those of M3)
  const PhysicallyBlockedLinearOpBase<Scalar> &blo =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(*M);
  const PhysicallyBlockedLinearOpBase<Scalar> &blo3 =
    Teuchos::dyn_cast<const PhysicallyBlockedLinearOpBase<Scalar> >(*M3);
  const RCP<DefaultBlockedLinearOp<Scalar> > full =
    defaultBlockedLinearOp<Scalar>();
  full->beginBlockFill(3, 3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      full->setBlock(i, j, i > j ? blo.getBlock(i,j) : blo3.getBlock(i,j));
  full->endBlockFill();
  TEST_THROW(Thyra::initializeOp<Scalar>(*lowsf, full, lows.ptr()),
    std::logic_error);
  TEST_THROW(Thyra::initializeOp<Scalar>(*lowsf, M2, lows.ptr()),
    std::logic_error);

  // lows is still the solver of M
  Thyra::LinearOpTester<Scalar> linearOpTester;
  const bool compareOpResult = linearOpTester.compare(*lows, *M, inOutArg(out));
  TEST_ASSERT(compareOpResult);
  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  const bool checkSolveResult = linearOpWithSolveTester.check(*lows, &out);
  TEST_ASSERT(checkSolveResult);

  // and can be reinitialized with a block triangular operator
  Thyra::initializeOp<Scalar>(*lowsf, M3, lows.ptr());
  const bool compareOp3Result = linearOpTester.compare(*lows, *M3, inOutArg(out));
  TEST_ASSERT(compareOp3Result);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedTriangularLinearOpWithSolveFactory,
  reinitializeRejected )


} // namespace Thyra
