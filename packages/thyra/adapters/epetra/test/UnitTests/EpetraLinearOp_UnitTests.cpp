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
1//
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

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


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
  epetraCrsM->PutScalar(two);
  const RCP<LinearOpBase<double> > epetraOp = nonconstEpetraLinearOp(epetraCrsM);

  const RCP<RowStatLinearOpBase<double> > rowStatOp =
    rcp_dynamic_cast<RowStatLinearOpBase<double> >(epetraOp, true);

  if (g_dumpAll) {
    out << "epetraOp = " << *epetraOp;
  }

  // Get the inverse row sums

  const RCP<VectorBase<double> > inv_row_sums =
    createMember<double>(epetraOp->range());

  rowStatOp->getRowStat(RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());

  if (g_dumpAll) {
    out << "inv_row_sums = " << *inv_row_sums;
  }

  TEST_FLOATING_EQUALITY(
    sum<double>(*inv_row_sums),
    as<double>((1.0 / (two * numCols)) * numRows),
    as<double>(10.0 * ST::eps())
    );

}


TEUCHOS_UNIT_TEST( EpetraLinearOp, rectangular )
{
  using Teuchos::null;
  using Teuchos::inOutArg;
  using Teuchos::updateSuccess;

  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  const RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows, numCols);

  const RCP<const LinearOpBase<double> > epetraOp = epetraLinearOp(epetraCrsM);

  LinearOpTester<double> linearOpTester;
  linearOpTester.show_all_tests(g_show_all_tests);
  linearOpTester.dump_all(g_dumpAll);
  updateSuccess(linearOpTester.check(*epetraOp, inOutArg(out)), success);
   
}


TEUCHOS_UNIT_TEST( EpetraLinearOp, blocked_op )
{

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
  
  out << "Sub operators built" << std::endl;

  {
     // build composite operator
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       block2x2<double>(A00, A01, A10, A11),   block2x1<double>(A02,A12),
       block1x2<double>(A20, A21),             A22
       );
   
     out << "First composite operator built" << std::endl;
     
     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);
     
     randomize(-1.0, 1.0, x.ptr());
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
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
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
     // perform a matrix vector multiply
     apply(*A, NOTRANS, *x, y.ptr());

     out << "Second composite operator completed" << std::endl;
  }

  out << "Test complete" << std::endl;

}


} // namespace Thyra
