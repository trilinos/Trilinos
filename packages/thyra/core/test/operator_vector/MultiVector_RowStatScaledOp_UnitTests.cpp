// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  {
    Scalar powValue = 1.0;
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = mv->col(i);
      put_scalar(powValue,ptr.ptr()); // pow(value,i)
      powValue *= value; 
    }
  }

  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,sv_1.ptr());
  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,sv_2.ptr());

  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM,sv_3.ptr());
  mv->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM,sv_4.ptr());

  RCP<VectorBase<Scalar> > ref_range = Thyra::createMember(mv->range());
  RCP<VectorBase<Scalar> > ref_domain = Thyra::createMember(mv->domain());

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormType;
  NormType tol = 1e-12;
  {
    NormType normVal = Thyra::norm_2(*sv_1);
    Thyra::Vp_S(sv_1.ptr(),Scalar(-31.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_1)/normVal,NormType(0.0),tol);
  }

  {
    NormType normVal = Thyra::norm_2(*sv_2);
    Thyra::Vp_S(sv_2.ptr(),Scalar(-1.0/31.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_2)/normVal,NormType(0.0),tol);
  }

  // build the absolute column sum
  {
    RTOpPack::SubVectorView<Scalar> view;
    ref_domain->acquireDetachedView(Thyra::Range1D(),&view);
    Scalar powValue = 1.0;
    for(std::size_t i=0;i<numVecs;i++)  {
      Scalar value = -2.0;
      view.values()[i] = std::abs(powValue)*numProcs*myRows;
      powValue *= value; 
    }
  }

  {
    NormType normVal = Thyra::norm_2(*sv_3);
    Thyra::Vp_V(sv_3.ptr(),*ref_domain,Scalar(-1.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_3)/normVal,NormType(0.0),tol);
  }

  {
    NormType normVal = Thyra::norm_2(*sv_4);
    Thyra::reciprocal(*ref_domain,ref_domain.ptr()); // I hope this can be done in place
    Thyra::Vp_V(sv_4.ptr(),*ref_domain,Scalar(-1.0));
    TEST_FLOATING_EQUALITY(Thyra::norm_2(*sv_4)/normVal,NormType(0.0),tol);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( MultiVector_RowStatScaledOp,
  RowStat )

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MultiVector_RowStatScaledOp, ScaledOp,
  Scalar )
{
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormType;
  std::size_t numVecs = 5;
  Teuchos_Ordinal myRows = 10;

  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(myRows);

  RTOpPack::SubVectorView<Scalar> view; // for use when needed

  // intialize test multivector
  ///////////////////////////////////////////////////
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numVecs);
  {
    Scalar powValue = 1.0;
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = mv->col(i);
      put_scalar(powValue,ptr.ptr());
      powValue *= value; 
    }
  }

  // initialize scaling vectors
  ///////////////////////////////////////////////////
  RCP<VectorBase<Scalar> > scale_1 = Thyra::createMember(mv->range());
  RCP<VectorBase<Scalar> > scale_2 = Thyra::createMember(mv->domain());

  put_scalar(Scalar(2.0),scale_1.ptr());

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
    Scalar powValue = 1.0;
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = result_1->col(i);
      put_scalar(Scalar(2.0)*powValue,ptr.ptr());
      powValue *= value; 
    }
  }

  {
    Scalar powValue = 1.0;
    for(std::size_t i=0;i<numVecs;i++) {
      Scalar value = -2.0;
      RCP<VectorBase<Scalar> > ptr = result_2->col(i);
      put_scalar(Scalar((i+1)*2.0)*powValue,ptr.ptr());
      powValue *= value; 
    }
  }

  // test scale left
  mv->scaleLeft(*scale_1);

  {
    std::vector<NormType> norms(result_1->domain()->dim()); 
    Teuchos::ArrayView<NormType> av_norms = Teuchos::arrayViewFromVector(norms);
    Thyra::V_VmV(result_1.ptr(),*result_1,*mv);
    Thyra::norms_2(*result_1,av_norms);

    std::vector<NormType> zeros(norms.size(),NormType(0.0)); 
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,NormType(1e-12));
  }

  // test scale right
  mv->scaleRight(*scale_2);

  {
    std::vector<NormType> norms(result_2->domain()->dim()); 
    Teuchos::ArrayView<NormType> av_norms = Teuchos::arrayViewFromVector(norms);
    Thyra::V_VmV(result_2.ptr(),*result_2,*mv);
    Thyra::norms_2(*result_2,av_norms);

    std::vector<NormType> zeros(norms.size(),NormType(0.0)); 
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,NormType(1e-12));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( MultiVector_RowStatScaledOp,
  ScaledOp )

} // namespace
