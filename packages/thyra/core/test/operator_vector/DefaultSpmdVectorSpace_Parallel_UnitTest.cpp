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
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

//#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
#  include "RTOpPack_SPMD_apply_op_decl.hpp"
#  include "Thyra_SpmdVectorBase.hpp"
#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
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
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createZeroEleProcVS<Scalar>(2);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
#ifdef RTOPPACK_ENABLE_SHOW_DUMP
//  RTOpPack::show_spmd_apply_op_dump = true;
#endif
  ECHO(assign(v.ptr(), as<Scalar>(1.5)));
  TEST_EQUALITY_CONST(sum(*v), as<Scalar>(1.5*2*(vs->getComm()->getSize()-1)));
#ifdef RTOPPACK_ENABLE_SHOW_DUMP
//  RTOpPack::show_spmd_apply_op_dump = false;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcAssignSumVec )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace_Parallel, emptyProcAssignSumMultiVec,
  Scalar )
{
  const Ordinal localDim = 2;
  PRINT_VAR(localDim);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createZeroEleProcVS<Scalar>(localDim);
  const Ordinal numCols = 3; 
  PRINT_VAR(numCols);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, numCols);
#ifdef RTOPPACK_ENABLE_SHOW_DUMP
//  RTOpPack::show_spmd_apply_op_dump = true;
#endif
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
#ifdef RTOPPACK_ENABLE_SHOW_DUMP
//  RTOpPack::show_spmd_apply_op_dump = false;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace_Parallel,
  emptyProcAssignSumMultiVec )


#if 0


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, emptyProcVectorSpaceTester,
  Scalar )
{
  const Ordinal localDim = 2;
  PRINT_VAR(localDim);
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(localDim);

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType  ScalarMag;
  Scalar tol = 1.0e-12;
  bool showAllTests=true, dumpAll=true;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  vectorStdOpsTester.error_tol(tol);

#ifdef RTOPPACK_ENABLE_SHOW_DUMP
  RTOpPack::show_spmd_apply_op_dump = true;
#endif

  out << "\nTesting the VectorSpaceBase interface of vs ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*vs, &out), out, success);

  out << "\nTesting standard vector ops for vs ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*vs, &out), out, success);

#ifdef RTOPPACK_ENABLE_SHOW_DUMP
  RTOpPack::show_spmd_apply_op_dump = false;
#endif

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  emptyProcVectorSpaceTester )


#endif // if 0


// ToDo:

// Test a locally repricated vectors space ...

// Test a vector space with the wrong size global dim.



} // namespace
