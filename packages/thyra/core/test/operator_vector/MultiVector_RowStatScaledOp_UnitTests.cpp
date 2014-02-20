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
#include "Thyra_VectorStdOps.hpp"
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
using Thyra::VectorBase;
using Thyra::put_scalar;
using Thyra::createMembers;
using Thyra::DefaultSpmdMultiVector;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
typedef Thyra::Ordinal Ordinal;


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiVector_RowStatScaledOp, RowStat,
  Scalar )
{
  std::size_t numVecs = 5;
  Teuchos_Ordinal myRows = 10;
  int numProcs = Teuchos::DefaultComm<Teuchos_Ordinal>::getComm()->getSize();

  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(myRows);

  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numVecs);

  RCP<VectorBase<Scalar> > sv_1 = createMember(mv->range());
  RCP<VectorBase<Scalar> > sv_2 = createMember(mv->range());
  RCP<VectorBase<Scalar> > sv_3 = createMember(mv->domain());
  RCP<VectorBase<Scalar> > sv_4 = createMember(mv->domain());

  for(std::size_t i=0;i<numVecs;i++) {
    Scalar value = -2.0;
    RCP<VectorBase<Scalar> > ptr = mv->col(i);
    put_scalar(std::pow(value,i),ptr.ptr());
  }

  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,sv_1);
  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,sv_2);

  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM,sv_3);
  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM,sv_4);

  RCP<VectorBase<Scalar> > ref_range = Thyra::createMember(mv->range());
  RCP<VectorBase<Scalar> > ref_domain = Thyra::createMember(mv->domain());

  {
    Scalar normVal = Thyra::norm_2(*sv_1);
    Thyra::Vp_S(sv_1.ptr(),Scalar(-31.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_1)/normVal,0.0,1e-12);
  }

  {
    Scalar normVal = Thyra::norm_2(*sv_2);
    Thyra::Vp_S(sv_2.ptr(),Scalar(-1.0/31.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_2)/normVal,0.0,1e-12);
  }

  // build the absolute column sum
  {
    RTOpPack::SubVectorView<Scalar> view;
    ref_domain->acquireDetachedView(Thyra::Range1D(),&view);
    for(std::size_t i=0;i<numVecs;i++)  {
      Scalar value = -2.0;
      view.values()[i] = std::abs(std::pow(value,i))*numProcs*myRows;
    }
  }

  {
    Scalar normVal = Thyra::norm_2(*sv_3);
    Thyra::Vp_V(sv_3.ptr(),*ref_domain,Scalar(-1.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_3)/normVal,0.0,1e-12);
  }

  {
    Scalar normVal = Thyra::norm_2(*sv_4);
    Thyra::reciprocal(*ref_domain,ref_domain.ptr()); // I hope this can be done in place
    Thyra::Vp_V(sv_4.ptr(),*ref_domain,Scalar(-1.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_4)/normVal,0.0,1e-12);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( MultiVector_RowStatScaledOp,
  RowStat )

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiVector_RowStatScaledOp, ScaledOp,
  Scalar )
{
  std::size_t numVecs = 5;
  Teuchos_Ordinal myRows = 10;
  int numProcs = Teuchos::DefaultComm<Teuchos_Ordinal>::getComm()->getSize();

  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(myRows);

  RTOpPack::SubVectorView<Scalar> view; // for use when needed

  // intialize test multivector
  ///////////////////////////////////////////////////
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numVecs);
  for(std::size_t i=0;i<numVecs;i++) {
    Scalar value = -2.0;
    RCP<VectorBase<Scalar> > ptr = mv->col(i);
    put_scalar(std::pow(value,i),ptr.ptr());
  }

  // initialize scaling vectors
  ///////////////////////////////////////////////////
  RCP<VectorBase<Scalar> > scale_1 = Thyra::createMember(mv->range());
  RCP<VectorBase<Scalar> > scale_2 = Thyra::createMember(mv->domain());

  put_scalar(2.0,scale_1.ptr());

  {
    scale_2->acquireDetachedView(Thyra::Range1D(),&view);
    for(std::size_t i=0;i<numVecs;i++) 
      view.values()[i] = i+1.0;
  }

  // initialize exact value vectors
  ///////////////////////////////////////////////////
  RCP<MultiVectorBase<Scalar> > result_1 = Thyra::createMembers(vs,numVecs);
  RCP<MultiVectorBase<Scalar> > result_2 = Thyra::createMembers(vs,numVecs);

  {
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = result_1->col(i);
      put_scalar(2.0*std::pow(value,i),ptr.ptr());
    }
  }

  {
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = result_2->col(i);
      put_scalar((i+1)*2.0*std::pow(value,i),ptr.ptr());
    }
  }

  // test scale left
  mv->scaleLeft(*scale_1);

  {
    std::vector<Scalar> norms(result_1->domain()->dim()); 
    Teuchos::ArrayView<Scalar> av_norms = Teuchos::arrayViewFromVector(norms);
    Thyra::V_VmV(result_1.ptr(),*result_1,*mv);
    Thyra::norms_2(*result_1,av_norms);

    std::vector<Scalar> zeros(norms.size(),0.0); 
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,1e-12);
  }

  // test scale right
  mv->scaleRight(*scale_2);

  {
    std::vector<Scalar> norms(result_2->domain()->dim()); 
    Teuchos::ArrayView<Scalar> av_norms = Teuchos::arrayViewFromVector(norms);
    Thyra::V_VmV(result_2.ptr(),*result_2,*mv);
    Thyra::norms_2(*result_2,av_norms);

    std::vector<Scalar> zeros(norms.size(),0.0); 
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,1e-12);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( MultiVector_RowStatScaledOp,
  ScaledOp )

} // namespace
