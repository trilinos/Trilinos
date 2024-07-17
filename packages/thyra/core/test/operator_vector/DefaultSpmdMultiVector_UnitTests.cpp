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
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_MultiVectorTester.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inOutArg;
using Teuchos::Range1D;
using Thyra::VectorSpaceBase;
using Thyra::MultiVectorBase;
using Thyra::createMembers;
using Thyra::DefaultSpmdMultiVector;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
typedef Thyra::Ordinal Ordinal;


const int g_localDim = 6; // ToDo: Make variable!
const int g_numCols = 5;  // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> >
createSpmdVectorSpace(const Teuchos_Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
    localDim, -1 );
}


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, defaultConstruct,
  Scalar )
{
  const RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers(*vs, g_numCols);
  Teuchos::Array<Scalar> mv_sums(g_numCols);

#ifndef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  Thyra::assign<Scalar>(mv.ptr(), 1.0);   // Can't touch uninitialized memory!
#endif

  Thyra::sums<Scalar>(*mv, mv_sums());
  out << "sums(mv) = " << mv_sums << "\n";

#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  typedef Teuchos::ScalarTraits<Scalar> ST;
  for (int i = 0; i < g_numCols; ++i) {
    TEST_ASSERT(ST::isnaninf(mv_sums[i]));
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, defaultTester,
  Scalar )
{
  //typedef Teuchos::ScalarTraits<Scalar> ST; // unused
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  Thyra::MultiVectorTester<Scalar> mvTester;
  const bool mvTesterResult = mvTester.checkMultiVector(*vs, inOutArg(out));
  TEST_ASSERT(mvTesterResult);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  defaultTester )


// Make sure the currect public member access protections are in place
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, memberAccess,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<DefaultSpmdMultiVector<Scalar> > spmdMv =
    rcp_dynamic_cast<DefaultSpmdMultiVector<Scalar> >(mv, true);
  TEST_ASSERT(nonnull(spmdMv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  memberAccess )


//
// Make sure that DefaultSpmdMultiVector::subview(...) return a
// DefaultSpmdMultiVector object!
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, constContigSubViewImpl,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<const MultiVectorBase<Scalar> > mvv = mv.getConst()->subView(Range1D(0, 1));
  RCP<const DefaultSpmdMultiVector<Scalar> > spmdMvv =
    rcp_dynamic_cast<const DefaultSpmdMultiVector<Scalar> >(mvv, true);
  TEST_ASSERT(nonnull(spmdMvv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  constContigSubViewImpl )


//
// Make sure that DefaultSpmdMultiVector::subview(...) return a
// DefaultSpmdMultiVector object!
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, nonconstContigSubViewImpl,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<MultiVectorBase<Scalar> > mvv = mv->subView(Range1D(0, 1));
  RCP<DefaultSpmdMultiVector<Scalar> > spmdMvv =
    rcp_dynamic_cast<DefaultSpmdMultiVector<Scalar> >(mvv, true);
  TEST_ASSERT(nonnull(spmdMvv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  nonconstContigSubViewImpl )


//
// Test that a dangling non-const column subviews don't write back their data.
// If the parent object is gone then why write back the data?  This is a
// performance optimization.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, danglingSubViews,
  Scalar )
{
  using Teuchos::tuple;
  //typedef Teuchos::ScalarTraits<Scalar> ST; // unused
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(*vs, g_numCols);
  RCP<MultiVectorBase<Scalar> > mvView = mv->subView(tuple<int>(0, 1));
  mv = null;
#ifdef THYRA_DEBUG
  const int startingNumCopyBack = DefaultSpmdMultiVector<Scalar>::numSkipCopyBack;
#endif
  mvView = null; // Should not write back data since parent mv is gone now.
#ifdef THYRA_DEBUG
  TEST_EQUALITY(DefaultSpmdMultiVector<Scalar>::numSkipCopyBack, startingNumCopyBack+1);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  danglingSubViews )


#ifdef THYRA_DEBUG


//
// Test that invalid column indexes throw the right exception messages.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, invalidSubviews,
  Scalar )
{
  using Teuchos::tuple;
  // typedef Teuchos::ScalarTraits<Scalar> ST; // unused
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<const MultiVectorBase<Scalar> > mv = createMembers(*vs, 1);
  TEST_THROW(mv->subView(tuple<int>(0, 1)), std::invalid_argument);
  TEST_THROW(mv->subView(tuple<int>(-1)), std::out_of_range);
  TEST_THROW(mv->subView(tuple<int>(1)), std::out_of_range);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  invalidSubviews )


#endif // THYRA_DEBUG


} // namespace
