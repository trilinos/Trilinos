//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_ScaledIdentityLinearOpWithSolve.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Thyra::RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM;
using Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM;
using Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM;
using Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM;

typedef Thyra::VectorSpaceBase<double> VectorSpace;
typedef Thyra::VectorBase<double> Vector;
typedef Thyra::ScaledIdentityLinearOpWithSolve<double> ScaledIdentity;

static const double tol = 1.0e-14;

// ************************************************************
// Row statistics of s*I: |s|, and 1/|s| for the inverses (s != 0).
TEUCHOS_UNIT_TEST(ScaledIdentityLinearOpWithSolve, RowStats)
{
  RCP<const VectorSpace> space = Thyra::defaultSpmdVectorSpace<double>(3);
  RCP<ScaledIdentity> op       = Thyra::scaledIdentity(space, -2.0);
  RCP<Vector> v                = Thyra::createMember(space);

  op->getRowStat(ROW_STAT_ROW_SUM, v.ptr());
  TEST_FLOATING_EQUALITY(Thyra::min(*v), 2.0, tol);
  TEST_FLOATING_EQUALITY(Thyra::max(*v), 2.0, tol);
  op->getRowStat(ROW_STAT_COL_SUM, v.ptr());
  TEST_FLOATING_EQUALITY(Thyra::min(*v), 2.0, tol);
  TEST_FLOATING_EQUALITY(Thyra::max(*v), 2.0, tol);
  op->getRowStat(ROW_STAT_INV_ROW_SUM, v.ptr());
  TEST_FLOATING_EQUALITY(Thyra::min(*v), 0.5, tol);
  TEST_FLOATING_EQUALITY(Thyra::max(*v), 0.5, tol);
  op->getRowStat(ROW_STAT_INV_COL_SUM, v.ptr());
  TEST_FLOATING_EQUALITY(Thyra::min(*v), 0.5, tol);
  TEST_FLOATING_EQUALITY(Thyra::max(*v), 0.5, tol);

  op->setScale(0.0);
  TEST_ASSERT(op->rowStatIsSupported(ROW_STAT_ROW_SUM));
  TEST_ASSERT(!op->rowStatIsSupported(ROW_STAT_INV_ROW_SUM));
  TEST_ASSERT(!op->rowStatIsSupported(ROW_STAT_INV_COL_SUM));
}

// ************************************************************
// Left and right scaling by uniform vectors: diag(d) s I = (d s) I.
TEUCHOS_UNIT_TEST(ScaledIdentityLinearOpWithSolve, UniformScaling)
{
  RCP<const VectorSpace> space = Thyra::defaultSpmdVectorSpace<double>(3);
  RCP<ScaledIdentity> op       = Thyra::scaledIdentity(space, -2.0);
  TEST_ASSERT(op->supportsScaleLeft());
  TEST_ASSERT(op->supportsScaleRight());

  RCP<Vector> d = Thyra::createMember(space);
  Thyra::put_scalar(3.0, d.ptr());
  op->scaleLeft(*d);
  TEST_FLOATING_EQUALITY(op->scale(), -6.0, tol);
  Thyra::put_scalar(0.5, d.ptr());
  op->scaleRight(*d);
  TEST_FLOATING_EQUALITY(op->scale(), -3.0, tol);

  // Apply and solve with the scaled operator
  RCP<Vector> x = Thyra::createMember(space);
  for (int i = 0; i < 3; ++i) Thyra::set_ele(i, i + 1.0, x.ptr());
  RCP<Vector> y = Thyra::createMember(space);
  Thyra::apply<double>(*op, Thyra::NOTRANS, *x, y.ptr());
  RCP<Vector> r = y->clone_v();
  Thyra::Vp_StV(r.ptr(), 3.0, *x);  // y + 3 x = 0
  TEST_COMPARE(Thyra::norm_inf(*r), <=, tol);
  RCP<Vector> z = Thyra::createMember(space);
  op->solve(Thyra::NOTRANS, *y, z.ptr());
  Thyra::Vp_StV(z.ptr(), -1.0, *x);  // z - x = 0
  TEST_COMPARE(Thyra::norm_inf(*z), <=, tol);

  // A nonuniform scaling is not a scaled identity
  Thyra::set_ele(1, 4.0, d.ptr());
  TEST_THROW(op->scaleLeft(*d), std::logic_error);
  TEST_THROW(op->scaleRight(*d), std::logic_error);
  TEST_FLOATING_EQUALITY(op->scale(), -3.0, tol);
}

// ************************************************************
// Row sum scaling (as done by NOX with "Function Scaling" = "Row Sum") of
// the operator of AdjointAuxSensitivityModelEvaluator, [A 0; 0 alpha I], with
// A a MultiVectorLinearOp and null off-diagonal blocks: after the scaling by
// the inverse absolute row sums, all the absolute row sums are 1.
TEUCHOS_UNIT_TEST(ScaledIdentityLinearOpWithSolve, BlockedRowSumScaling)
{
  RCP<const VectorSpace> space           = Thyra::defaultSpmdVectorSpace<double>(2);
  RCP<Thyra::MultiVectorBase<double> > A = Thyra::createMembers(space, 2);
  {
    Thyra::DetachedMultiVectorView<double> Av(*A);
    Av(0, 0) = 2.0;
    Av(0, 1) = -1.0;
    Av(1, 0) = 0.5;
    Av(1, 1) = 4.0;
  }
  const int num_adjoint                   = 2;
  RCP<Thyra::LinearOpBase<double> > mv_op = Thyra::nonconstMultiVectorLinearOp<double>(A, num_adjoint);
  RCP<const VectorSpace> z_space          = Thyra::defaultSpmdVectorSpace<double>(3);
  const double alpha                      = -4.0;
  RCP<ScaledIdentity> si_op               = Thyra::scaledIdentity(z_space, alpha);
  RCP<Thyra::LinearOpBase<double> > null_op;
  RCP<Thyra::LinearOpBase<double> > W =
      Thyra::nonconstBlock2x2<double>(mv_op, null_op, null_op, si_op);

  RCP<const Thyra::RowStatLinearOpBase<double> > W_rs =
      rcp_dynamic_cast<const Thyra::RowStatLinearOpBase<double> >(W, true);
  RCP<Thyra::ScaledLinearOpBase<double> > W_sc =
      rcp_dynamic_cast<Thyra::ScaledLinearOpBase<double> >(W, true);
  TEST_ASSERT(W_rs->rowStatIsSupported(ROW_STAT_INV_ROW_SUM));
  TEST_ASSERT(W_sc->supportsScaleLeft());

  RCP<Vector> inv_row_sum = Thyra::createMember(W->range());
  W_rs->getRowStat(ROW_STAT_INV_ROW_SUM, inv_row_sum.ptr());
  W_sc->scaleLeft(*inv_row_sum);

  RCP<Vector> row_sum = Thyra::createMember(W->range());
  W_rs->getRowStat(ROW_STAT_ROW_SUM, row_sum.ptr());
  TEST_FLOATING_EQUALITY(Thyra::min(*row_sum), 1.0, tol);
  TEST_FLOATING_EQUALITY(Thyra::max(*row_sum), 1.0, tol);
  TEST_FLOATING_EQUALITY(si_op->scale(), -1.0, tol);
}

}  // namespace Tempus_Unit_Test
