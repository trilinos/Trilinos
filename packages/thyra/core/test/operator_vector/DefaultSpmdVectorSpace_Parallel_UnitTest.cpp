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
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_SpmdLocalDataAccess.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
#include "Thyra_MultiVectorStdOpsTester.hpp"
#include "RTOpPack_SPMD_apply_op_decl.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


int g_localDim = 4;
int g_numCols1 = 2;
int g_numCols2 = 2;
bool g_show_all_tests = false;
bool g_dump_objects = false;
bool g_dumpRTOps = false;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Local dimension of each vector." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "num-cols-1", &g_numCols1, "" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "num-cols-2", &g_numCols2, "" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "show-all-tests", "no-show-all-tests", &g_show_all_tests,
    "Set if all tests are shown or not." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-objects", "no-dump-objects", &g_dump_objects,
    "Set if vector, multivector, etc. objects are printed or not." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-rtops", "no-dump-rtops", &g_dumpRTOps,
    "Set if RTOps are dumped or not." );
}


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::get_extra_data;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::ConstDetachedVectorView;
using Thyra::DetachedVectorView;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::DetachedSpmdVectorView;
using Thyra::assign;
using Thyra::apply;
typedef Thyra::Ordinal Ordinal;


#define PRINT_VAR(varName) \
  out << #varName" = " << (varName) << "\n"


template<class Scalar>
RCP<const DefaultSpmdVectorSpace<Scalar> >
createZeroEleProcVS(const Ordinal localSize)
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  const Ordinal thisLocalSize = comm->getRank()==0 ? 0 : localSize;
  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, thisLocalSize, -1);

  return vs;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcConstruct,
  Scalar )
{
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createZeroEleProcVS<Scalar>(2);
  TEST_EQUALITY_CONST(vs->dim(), as<Ordinal>(2*(vs->getComm()->getSize()-1)));
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcAssignSumVec,
  Scalar )
{
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  ECHO(assign(v.ptr(), as<Scalar>(1.5)));
  TEST_EQUALITY_CONST(sum(*v), as<Scalar>(1.5*g_localDim*(vs->getComm()->getSize()-1)));
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcAssignSumVec )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcGetFullSubVector,
  Scalar )
{
  const Ordinal localSubDim = g_localDim;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(localSubDim);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  ECHO(assign(v.ptr(), as<Scalar>(1.5)));
  RTOpPack::ConstSubVectorView<Scalar> subVec;
  v->acquireDetachedView(Teuchos::Range1D(), &subVec);
  TEST_EQUALITY(subVec.subDim(),
    as<Ordinal>(localSubDim*(vs->getComm()->getSize()-1)));
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcGetFullSubVector )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcPrintVec,
  Scalar )
{
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  ECHO(assign(v.ptr(), as<Scalar>(1.5)));
  out << "v = " << *v;
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcPrintVec )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcAssignSumMultiVec,
  Scalar )
{
  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createZeroEleProcVS<Scalar>(localDim);
  const Ordinal numCols = 3; 
  PRINT_VAR(numCols);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, numCols);
  const Scalar val = 1.5;
  PRINT_VAR(val);
  ECHO(assign(mv.ptr(), as<Scalar>(val)));
  Teuchos::Array<Scalar> sums(numCols);
  Thyra::sums<Scalar>(*mv, sums());
  for (int j = 0; j < numCols; ++j) {
    PRINT_VAR(j);
    TEST_EQUALITY(sums[j],  
      as<Scalar>(val*localDim*(vs->getComm()->getSize()-1)));
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcAssignSumMultiVec )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcPrintMultiVec,
  Scalar )
{
  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const Ordinal numCols = 3;
  PRINT_VAR(numCols);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createZeroEleProcVS<Scalar>(localDim);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, numCols);
  ECHO(assign(mv.ptr(), as<Scalar>(1.5)));
  out << "mv = " << *mv;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcPrintMultiVec )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcSimpleMultiVecAdjointApply,
  Scalar )
{
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const Ordinal numCols1 = g_numCols1;
  const Ordinal numCols2= g_numCols2;
  PRINT_VAR(numCols1);
  PRINT_VAR(numCols2);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(localDim);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, numCols1);
  const Scalar val1 = 2.0;
  PRINT_VAR(val1);
  ECHO(assign<Scalar>(mv.ptr(), val1));
  out << "mv = " << *mv;
  ECHO(const RCP<MultiVectorBase<Scalar> > Y = createMembers<Scalar>(mv->domain(), numCols2));
  ECHO(const RCP<MultiVectorBase<Scalar> > X = createMembers<Scalar>(mv->range(), numCols2));
  const Scalar val2 = 3.0;
  PRINT_VAR(val2);
  ECHO(assign<Scalar>(X.ptr(), val2));
  out << "X = " << *X;
  ECHO(apply<Scalar>(*mv, Thyra::CONJTRANS, *X, Y.ptr()));
  out << "Y = " << *Y;
  RTOpPack::ConstSubMultiVectorView<Scalar> Y_smvv =
    Thyra::getLocalSubMultiVectorView<Scalar>(Y);
  for (int i = 0; i < numCols1; ++i) {
    for (int j = 0; j < numCols2; ++j) {
      out << "i = " << i << ", j = " << j << ": ";
      TEST_EQUALITY(Y_smvv(i,j), as<Scalar>(val1 * val2 * vs->dim()));
    }
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcSimpleMultiVecAdjointApply )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, emptyProcVectorSpaceTester,
  Scalar )
{

  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(localDim);

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType  ScalarMag;
  Scalar tol = 1.0e-12;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(g_show_all_tests);
  vectorSpaceTester.dump_all(g_dump_objects);

  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  vectorStdOpsTester.error_tol(tol);

  Thyra::MultiVectorStdOpsTester<Scalar> multiVectorStdOpsTester;
  multiVectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  multiVectorStdOpsTester.error_tol(tol);

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  out << "\nTesting the VectorSpaceBase interface of vs ...\n";
  TEST_ASSERT(vectorSpaceTester.check(*vs, &out));

  out << "\nTesting standard vector ops for vs ...\n";
  TEST_ASSERT(vectorStdOpsTester.checkStdOps(*vs, &out));

  out << "\nTesting standard multivector ops for vs ...\n";
  TEST_ASSERT(multiVectorStdOpsTester.checkStdOps(*vs, &out));

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  emptyProcVectorSpaceTester )


// ToDo:

// Test a locally repricated vectors space ...

// Test a vector space with the wrong size global dim.



} // namespace
