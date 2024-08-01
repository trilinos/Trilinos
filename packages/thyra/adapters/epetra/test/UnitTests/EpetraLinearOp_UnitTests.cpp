// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_UnitTestHelpers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "EpetraLinearOpTestHelpers.hpp"

#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

namespace {



} // namespace


namespace Thyra {


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( EpetraLinearOp, ScaledLinearOpBase )
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<double> ST;

  // Set up an EpetraLinearOp

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;
  const RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows, numCols);
  const RCP<LinearOpBase<double> > epetraOp = nonconstEpetraLinearOp(epetraCrsM);

  const RCP<ScaledLinearOpBase<double> > scaledOp =
    rcp_dynamic_cast<ScaledLinearOpBase<double> >(epetraOp, true);

  // Get the original mat-vec

  const double two = 2.0;

  const RCP<VectorBase<double> > rhs_vec =
    createMember<double>(epetraOp->domain());
  assign<double>(rhs_vec.ptr(), two);

  const RCP<VectorBase<double> > lhs_orig_vec =
    createMember<double>(epetraOp->range());

  apply<double>(*epetraOp, Thyra::NOTRANS, *rhs_vec, lhs_orig_vec.ptr());

  if (g_dumpAll) {
    out << "epetraOp = " << *epetraOp;
    out << "rhs_vec = " << *rhs_vec;
    out << "lhs_orig_vec = " << *lhs_orig_vec;
  }

  // Scale the op from the left (row scaling)

  const double three = 3.0;

  const RCP<VectorBase<double> > row_scaling =
    createMember<double>(epetraOp->range());
  assign<double>(row_scaling.ptr(), three);

  scaledOp->scaleLeft(*row_scaling);

  if (g_dumpAll) {
    out << "row_scaling = " << *row_scaling;
    out << "epetraOp left scaled = " << *epetraOp;
  }

  // Test that resulting left scaling

  const RCP<VectorBase<double> > lhs_left_scaled_vec =
    createMember<double>(epetraOp->range());

  apply<double>(*epetraOp, NOTRANS, *rhs_vec, lhs_left_scaled_vec.ptr());

  if (g_dumpAll) {
    out << "lhs_left_scaled_vec = " << *lhs_left_scaled_vec;
  }

  TEST_FLOATING_EQUALITY(
    three * sum<double>(*lhs_orig_vec),
    sum<double>(*lhs_left_scaled_vec),
    as<double>(10.0 * ST::eps())
    );

  // Left scale the matrix back

  const RCP<VectorBase<double> > inv_row_scaling =
    createMember<double>(epetraOp->range());
  reciprocal<double>(*row_scaling, inv_row_scaling.ptr());

  scaledOp->scaleLeft(*inv_row_scaling);

  if (g_dumpAll) {
    out << "inv_row_scaling = " << *row_scaling;
    out << "epetraOp left scaled back to orig = " << *epetraOp;
  }

  const RCP<VectorBase<double> > lhs_orig2_vec =
    createMember<double>(epetraOp->range());

  apply<double>(*epetraOp, NOTRANS, *rhs_vec, lhs_orig2_vec.ptr());

  if (g_dumpAll) {
    out << "lhs_orig2_vec = " << *lhs_orig2_vec;
  }

  TEST_FLOATING_EQUALITY(
    sum<double>(*lhs_orig_vec),
    sum<double>(*lhs_orig2_vec),
    as<double>(10.0 * ST::eps())
    );
  // NOTE: Above, it would ask for exact binary match except if one uses
  // threading it will not match exactly!

  // Scale the op from the right (col scaling)

  const double four = 4.0;

  const RCP<VectorBase<double> > col_scaling =
    createMember<double>(epetraOp->domain());
  assign<double>(col_scaling.ptr(), four);

  scaledOp->scaleRight(*col_scaling);

  if (g_dumpAll) {
    out << "col_scaling = " << *col_scaling;
    out << "epetraOp right scaled = " << *epetraOp;
  }

  // Test that resulting right scaling

  const RCP<VectorBase<double> > lhs_right_scaled_vec =
    createMember<double>(epetraOp->range());

  apply<double>(*epetraOp, NOTRANS, *rhs_vec, lhs_right_scaled_vec.ptr());

  if (g_dumpAll) {
    out << "lhs_right_scaled_vec = " << *lhs_right_scaled_vec;
  }

  TEST_FLOATING_EQUALITY(
    four * sum<double>(*lhs_orig_vec),
    sum<double>(*lhs_right_scaled_vec),
    as<double>(10.0 * ST::eps())
    );

  // Right scale the matrix back

  const RCP<VectorBase<double> > inv_col_scaling =
    createMember<double>(epetraOp->domain());
  reciprocal<double>(*col_scaling, inv_col_scaling.ptr());

  scaledOp->scaleRight(*inv_col_scaling);

  if (g_dumpAll) {
    out << "inv_col_scaling = " << *col_scaling;
    out << "epetraOp right scaled back to orig = " << *epetraOp;
  }

  const RCP<VectorBase<double> > lhs_orig3_vec =
    createMember<double>(epetraOp->range());

  apply<double>(*epetraOp, NOTRANS, *rhs_vec, lhs_orig3_vec.ptr());

  if (g_dumpAll) {
    out << "lhs_orig3_vec = " << *lhs_orig3_vec;
  }

  TEST_FLOATING_EQUALITY(
    sum<double>(*lhs_orig_vec),
    sum<double>(*lhs_orig3_vec),
    as<double>(10.0 * ST::eps())
    );
  // NOTE: Above, it would ask for exact binary match except if one uses
  // threading it will not match exactly!

}


TEUCHOS_UNIT_TEST( EpetraLinearOp, RowStatLinearOpBase )
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<double> ST;

  // Set up the EpetraLinearOp

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;
  const RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows, numCols);
  const double two = 2.0;
  epetraCrsM->PutScalar(-two); // put in negative two just to be extra "tricky"
  const RCP<LinearOpBase<double> > epetraOp = nonconstEpetraLinearOp(epetraCrsM);

  const RCP<RowStatLinearOpBase<double> > rowStatOp =
    rcp_dynamic_cast<RowStatLinearOpBase<double> >(epetraOp, true);

  if (g_dumpAll) {
    out << "epetraOp = " << *epetraOp;
  }

  // Get the inverse row sums

  const RCP<VectorBase<double> > inv_row_sums =
    createMember<double>(epetraOp->range());
  const RCP<VectorBase<double> > row_sums =
    createMember<double>(epetraOp->range());

  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  if (g_dumpAll) {
    out << "inv_row_sums = " << *inv_row_sums;
    out << "row_sums = " << *row_sums;
  }

  TEST_FLOATING_EQUALITY(
    sum<double>(*inv_row_sums),
    as<double>((1.0 / (two * numCols)) * numRows),
    as<double>(10.0 * ST::eps())
    );

  TEST_FLOATING_EQUALITY(
    sum<double>(*row_sums),
    as<double>(two * numCols * numRows),
    as<double>(10.0 * ST::eps())
    );
}

RCP<Epetra_CrsMatrix> getMyEpetraMatrix(int numRows, int numCols, double shift=0.0)
{
  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const Epetra_Map rowMap(numRows, 0, *comm);
  const Epetra_Map domainMap(numCols, numCols, 0, *comm);

  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap,domainMap,0));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j)
    columnIndices[j] = j;

  const int numLocalRows = rowMap.NumMyElements();
  for (int i = 0; i < numLocalRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10 + shift;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );
  }

  epetraCrsM->FillComplete(domainMap,rowMap);
  return epetraCrsM;
}

TEUCHOS_UNIT_TEST( EpetraLinearOp, Blocked_ScaledLinearOpBase)
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;
  using Teuchos::rcp_dynamic_cast;
  // typedef ScalarTraits<double> ST; // unused

  // Set up the EpetraLinearOp

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows ;

  out << "numRows = " << numRows << ", numCols = " << numCols << std::endl;

  const RCP<Epetra_CrsMatrix> epetraCrsM00_base = getMyEpetraMatrix(numRows, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM01_base = getMyEpetraMatrix(numRows, numCols);
  const RCP<Epetra_CrsMatrix> epetraCrsM10_base = getMyEpetraMatrix(numCols, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM11_base = getMyEpetraMatrix(numCols, numCols);
  epetraCrsM00_base->PutScalar(2.0);
  epetraCrsM01_base->PutScalar(-8.0);
  epetraCrsM10_base->PutScalar(-9.0);
  epetraCrsM11_base->PutScalar(3.0);
  const RCP<Epetra_CrsMatrix> epetraCrsM00 = getMyEpetraMatrix(numRows, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM01 = getMyEpetraMatrix(numRows, numCols);
  const RCP<Epetra_CrsMatrix> epetraCrsM10 = getMyEpetraMatrix(numCols, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM11 = getMyEpetraMatrix(numCols, numCols);
  epetraCrsM00->PutScalar(2.0);
  epetraCrsM01->PutScalar(-8.0);
  epetraCrsM10->PutScalar(-9.0);
  epetraCrsM11->PutScalar(3.0);

  const RCP<const LinearOpBase<double> > op00_base = epetraLinearOp(epetraCrsM00_base);
  const RCP<const LinearOpBase<double> > op01_base = epetraLinearOp(epetraCrsM01_base);
  const RCP<const LinearOpBase<double> > op10_base = epetraLinearOp(epetraCrsM10_base);
  const RCP<const LinearOpBase<double> > op11_base = epetraLinearOp(epetraCrsM11_base);

  const RCP<LinearOpBase<double> > op00 = nonconstEpetraLinearOp(epetraCrsM00);
  const RCP<LinearOpBase<double> > op01 = nonconstEpetraLinearOp(epetraCrsM01);
  const RCP<LinearOpBase<double> > op10 = nonconstEpetraLinearOp(epetraCrsM10);
  const RCP<LinearOpBase<double> > op11 = nonconstEpetraLinearOp(epetraCrsM11);

  const RCP<const LinearOpBase<double> > blocked_base = block2x2(op00_base,op01_base,op10_base,op11_base);
  const RCP<LinearOpBase<double> > blocked = nonconstBlock2x2(op00,op01,op10,op11);

  const RCP<VectorBase<double> > left_scale  = createMember<double>(blocked_base->range());
  const RCP<VectorBase<double> > right_scale = createMember<double>(blocked_base->domain());

  put_scalar(7.0,left_scale.ptr());
  put_scalar(-4.0,right_scale.ptr());

  rcp_dynamic_cast<ScaledLinearOpBase<double> >(blocked)->scaleLeft(*left_scale);

  {
    LinearOpTester<double> tester;
    tester.set_all_error_tol(1e-10);
    tester.show_all_tests(true);
    tester.dump_all(true);
    tester.num_random_vectors(5);
    const RCP<const LinearOpBase<double> > left_op = Thyra::diagonal(left_scale);
    const RCP<const LinearOpBase<double> > ref_op = multiply(left_op,blocked_base);

    updateSuccess(tester.compare(*ref_op, *blocked, ptrFromRef(out)), success);
  }

  rcp_dynamic_cast<ScaledLinearOpBase<double> >(blocked)->scaleRight(*right_scale);

  {
    LinearOpTester<double> tester;
    tester.set_all_error_tol(1e-10);
    tester.show_all_tests(true);
    tester.dump_all(true);
    tester.num_random_vectors(5);
    const RCP<const LinearOpBase<double> > left_op = Thyra::diagonal(left_scale);
    const RCP<const LinearOpBase<double> > right_op = Thyra::diagonal(right_scale);
    const RCP<const LinearOpBase<double> > ref_op = multiply(left_op,blocked_base,right_op);

    updateSuccess(tester.compare(*ref_op, *blocked, ptrFromRef(out)), success);
  }
}

TEUCHOS_UNIT_TEST( EpetraLinearOp, Blocked_RowStatLinearOpBase )
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<double> ST;

  // Set up the EpetraLinearOp

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  out << "numRows = " << numRows << ", numCols = " << numCols << std::endl;

  const RCP<Epetra_CrsMatrix> epetraCrsM00 = getEpetraMatrix(numRows, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM01 = getEpetraMatrix(numRows, numCols);
  const RCP<Epetra_CrsMatrix> epetraCrsM10 = getEpetraMatrix(numCols, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM11 = getEpetraMatrix(numCols, numCols);
  epetraCrsM00->PutScalar(2.0);
  epetraCrsM01->PutScalar(-8.0);
  epetraCrsM10->PutScalar(-9.0);
  epetraCrsM11->PutScalar(3.0);

  const RCP<const LinearOpBase<double> > op00 = epetraLinearOp(epetraCrsM00);
  const RCP<const LinearOpBase<double> > op01 = epetraLinearOp(epetraCrsM01);
  const RCP<const LinearOpBase<double> > op10 = epetraLinearOp(epetraCrsM10);
  const RCP<const LinearOpBase<double> > op11 = epetraLinearOp(epetraCrsM11);

  const RCP<const LinearOpBase<double> > blocked = block2x2(op00,op01,op10,op11);

  const RCP<const RowStatLinearOpBase<double> > rowStatOp =
    rcp_dynamic_cast<const RowStatLinearOpBase<double> >(blocked, true);

  if (g_dumpAll) {
    out << "epetraOp = " << *blocked;
  }

  // Get the inverse row sums

  const RCP<VectorBase<double> > inv_row_sums =
    createMember<double>(blocked->range());
  const RCP<VectorBase<double> > row_sums =
    createMember<double>(blocked->range());

  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  if (g_dumpAll) {
    out << "inv_row_sums = " << *inv_row_sums;
    out << "row_sums = " << *row_sums;
  }

  TEST_FLOATING_EQUALITY(
    sum<double>(*inv_row_sums),
    as<double>((1.0/(numRows*2.0+numCols*8.0))*numRows + (1.0/(numRows*9.0+numCols*3.0))*numCols),
    as<double>(10.0 * ST::eps())
    );
  TEST_FLOATING_EQUALITY(
    sum<double>(*row_sums),
    as<double>((numRows*2.0+numCols*8.0)*numRows + (numRows*9.0+numCols*3.0)*numCols),
    as<double>(10.0 * ST::eps())
    );
}

TEUCHOS_UNIT_TEST( EpetraLinearOp, Blocked_ScalingWithMultiVectors)
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<double> ST;

  // Set up the EpetraLinearOp

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const RCP<const Teuchos::Comm<Ordinal> > tComm =
      Teuchos::DefaultComm<Ordinal>::getComm();
  const int numLocalRows = 4;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  out << "numRows = " << numRows << ", numCols = " << numCols << std::endl;

  const RCP<Epetra_CrsMatrix> epetraCrsM00 = getMyEpetraMatrix(numRows, numRows);
  const RCP<Epetra_CrsMatrix> epetraCrsM00_base = getMyEpetraMatrix(numRows, numRows);
  epetraCrsM00->PutScalar(2.0);
  epetraCrsM00_base->PutScalar(2.0);

  const RCP<LinearOpBase<double> > op00 = nonconstEpetraLinearOp(epetraCrsM00);
  const RCP<const LinearOpBase<double> > op00_base = epetraLinearOp(epetraCrsM00_base);

  RCP<const Thyra::VectorSpaceBase<double> > vs_0 = op00->range();
  RCP<const Thyra::VectorSpaceBase<double> > vs_1 = Thyra::locallyReplicatedDefaultSpmdVectorSpace<double>(tComm,numCols);

  RCP<Thyra::MultiVectorBase<double> > vec_01  = Thyra::createMembers(vs_0,numCols);
  RCP<Thyra::MultiVectorBase<double> > vec_10t = Thyra::createMembers(op00->domain(),numCols); // tranposed
  RCP<Thyra::MultiVectorBase<double> > vec_01_base  = Thyra::createMembers(vs_0,numCols);
  RCP<Thyra::MultiVectorBase<double> > vec_10t_base = Thyra::createMembers(op00->domain(),numCols); // tranposed
  const RCP<LinearOpBase<double> > op10t = vec_10t;
  const RCP<const LinearOpBase<double> > op10t_base = vec_10t_base;
  assign(vec_01.ptr(),-8.0);
  assign(vec_10t.ptr(),-9.0);
  assign(vec_01_base.ptr(),-8.0);
  assign(vec_10t_base.ptr(),-9.0);

  const RCP<LinearOpBase<double> > op01 = vec_01;
  const RCP<LinearOpBase<double> > op10 = nonconstAdjoint(op10t);
  const RCP<LinearOpBase<double> > op11 = nonconstZero(vec_01->domain(),vec_01->domain());

  const RCP<const LinearOpBase<double> > op01_base = vec_01_base;
  const RCP<const LinearOpBase<double> > op10_base = adjoint(op10t_base);
  const RCP<const LinearOpBase<double> > op11_base = zero(vec_01->domain(),vec_01->domain());

  out << "FIRST" << std::endl;
  const RCP<LinearOpBase<double> > blocked = nonconstBlock2x2(op00,op01,op10,op11);
  out << "SECOND" << Teuchos::describe(*blocked,Teuchos::VERB_EXTREME) << std::endl;
  const RCP<const LinearOpBase<double> > blocked_base = block2x2(op00_base,op01_base,op10_base,op11_base);

  const RCP<const RowStatLinearOpBase<double> > rowStatOp =
    rcp_dynamic_cast<const RowStatLinearOpBase<double> >(blocked, true);

  // Get the inverse row sums

  const RCP<VectorBase<double> > inv_row_sums =
    createMember<double>(blocked->range());
  const RCP<VectorBase<double> > row_sums =
    createMember<double>(blocked->range());

  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  TEST_FLOATING_EQUALITY(
    sum<double>(*inv_row_sums),
    as<double>((1.0/(numRows*2.0+numCols*8.0))*numRows + (1.0/(numRows*9.0+numCols*0.0))*numCols),
    as<double>(10.0 * ST::eps())
    );
  TEST_FLOATING_EQUALITY(
    sum<double>(*row_sums),
    as<double>((numRows*2.0+numCols*8.0)*numRows + (numRows*9.0+numCols*0.0)*numCols),
    as<double>(10.0 * ST::eps())
    );

  {
    const RCP<VectorBase<double> > left_scale  = createMember<double>(blocked_base->range());
    const RCP<VectorBase<double> > right_scale = createMember<double>(blocked_base->domain());

    put_scalar(7.0,left_scale.ptr());
    put_scalar(-4.0,right_scale.ptr());

    rcp_dynamic_cast<ScaledLinearOpBase<double> >(blocked)->scaleLeft(*left_scale);

    {
      LinearOpTester<double> tester;
      tester.set_all_error_tol(1e-10);
      tester.show_all_tests(true);
      tester.dump_all(false);
      tester.num_random_vectors(2);
      const RCP<const LinearOpBase<double> > left_op = Thyra::diagonal(left_scale);
      const RCP<const LinearOpBase<double> > ref_op = multiply(left_op,blocked_base);

      updateSuccess(tester.compare(*ref_op, *blocked, ptrFromRef(out)), success);
    }

    rcp_dynamic_cast<ScaledLinearOpBase<double> >(blocked)->scaleRight(*right_scale);

    {
      LinearOpTester<double> tester;
      tester.set_all_error_tol(1e-10);
      tester.show_all_tests(true);
      tester.dump_all(false);
      tester.num_random_vectors(5);
      const RCP<const LinearOpBase<double> > left_op = Thyra::diagonal(left_scale);
      const RCP<const LinearOpBase<double> > right_op = Thyra::diagonal(right_scale);
      const RCP<const LinearOpBase<double> > ref_op = multiply(left_op,blocked_base,right_op);

      updateSuccess(tester.compare(*ref_op, *blocked, ptrFromRef(out)), success);
    }
  }
}


TEUCHOS_UNIT_TEST( EpetraLinearOp, rectangular )
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;

  const RCP<const Epetra_Comm> comm = getEpetraComm();
  const int numProcs = comm->NumProc();

  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  const RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows, numCols);

  const RCP<const LinearOpBase<double> > epetraOp = epetraLinearOp(epetraCrsM);

  LinearOpTester<double> linearOpTester;
  linearOpTester.check_adjoint(numProcs == 1);
  linearOpTester.show_all_tests(g_show_all_tests);
  linearOpTester.dump_all(g_dumpAll);
  updateSuccess(linearOpTester.check(*epetraOp, inOutArg(out)), success);

  // NOTE: Above, it would seem the Epetra_CrsMatrix::Apply(...) does not work
  // when doing and adjoint where the RowMap has empty processes.

}


TEUCHOS_UNIT_TEST( EpetraLinearOp, blocked_op )
{

  if (Teuchos::GlobalMPISession::getNProc() > 2) {
    out << "Skipping test if numProc > 2 since it fails for some reason\n";
    return;
  }

  using Teuchos::describe;

  // build sub operators
  RCP<const LinearOpBase<double> > A00 =
    epetraLinearOp(getEpetraMatrix(4,4,0));
  RCP<const LinearOpBase<double> > A01 =
    epetraLinearOp(getEpetraMatrix(4,3,1));
  RCP<const LinearOpBase<double> > A02 =
    epetraLinearOp(getEpetraMatrix(4,2,2));
  RCP<const LinearOpBase<double> > A10 =
    epetraLinearOp(getEpetraMatrix(3,4,3));
  RCP<const LinearOpBase<double> > A11 =
    epetraLinearOp(getEpetraMatrix(3,3,4));
  RCP<const LinearOpBase<double> > A12 =
    epetraLinearOp(getEpetraMatrix(3,2,5));
  RCP<const LinearOpBase<double> > A20 =
    epetraLinearOp(getEpetraMatrix(2,4,6));
  RCP<const LinearOpBase<double> > A21 =
    epetraLinearOp(getEpetraMatrix(2,3,8));
  RCP<const LinearOpBase<double> > A22 =
    epetraLinearOp(getEpetraMatrix(2,2,8));

  const Teuchos::EVerbosityLevel verbLevel =
    (g_dumpAll ? Teuchos::VERB_HIGH : Teuchos::VERB_MEDIUM);

  out << "Sub operators built" << std::endl;

  {
     // build composite operator
     RCP<const LinearOpBase<double> > A =
       block2x2<double>(
         block2x2<double>(A00, A01, A10, A11),   block2x1<double>(A02,A12),
         block1x2<double>(A20, A21),             A22
         );

     out << "First composite operator built" << std::endl;

     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);

     randomize(-1.0, 1.0, x.ptr());

     out << "A = \n" << describe(*A, verbLevel) << std::endl;
     out << "x = \n" << describe(*x, verbLevel) << std::endl;
     out << "y = \n" << describe(*y, verbLevel) << std::endl;

     // perform a matrix vector multiply
     apply(*A, NOTRANS, *x, y.ptr());

     out << "First composite operator completed" << std::endl;
  }

  {
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       A11,                          block1x2<double>(A10, A12),
       block2x1<double>(A01, A21),   block2x2<double>(A00, A02, A20, A22)
       );

     out << "Second composite operator built" << std::endl;

     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);

     randomize(-1.0, 1.0, x.ptr());

     out << "A = \n" << describe(*A, verbLevel) << std::endl;
     out << "x = \n" << describe(*x, verbLevel) << std::endl;
     out << "y = \n" << describe(*y, verbLevel) << std::endl;

     // perform a matrix vector multiply
     apply(*A, NOTRANS, *x, y.ptr());

     out << "Second composite operator completed" << std::endl;
  }

  out << "Test complete" << std::endl;

}


} // namespace Thyra
