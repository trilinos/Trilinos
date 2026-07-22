// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Thyra::createMember;
using Thyra::V_S;
using Thyra::defaultSpmdVectorSpace;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::DetachedSpmdVectorView;


const Teuchos_Ordinal g_localDim = 4; // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createSpmdVectorSpace(const Teuchos_Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
    localDim, -1 );
}


//
// ConstDetachedSpmdVectorView unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ConstDetachedSpmdVectorView, construct_null,
  Scalar )
{
  ECHO(ConstDetachedSpmdVectorView<Scalar> dvv(null));
  TEST_EQUALITY_CONST(dvv.globalOffset(), 0);
  TEST_EQUALITY_CONST(dvv.subDim(), 0);
  TEST_EQUALITY_CONST(dvv.values(), null);
  TEST_EQUALITY_CONST(dvv.stride(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ConstDetachedSpmdVectorView,
  construct_null )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ConstDetachedSpmdVectorView, basic, Scalar )
{
  const int procRank = Teuchos::GlobalMPISession::getRank();
  ECHO(const RCP<const VectorSpaceBase<Scalar> >
    vs = createSpmdVectorSpace<Scalar>(g_localDim));
  ECHO(RCP<VectorBase<Scalar> > v = createMember(vs));
  ECHO(V_S(v.ptr(), as<Scalar>(2.0)));
  ECHO(ConstDetachedSpmdVectorView<Scalar> dvv(v));
  TEST_INEQUALITY(dvv.spmdSpace(), null);
  TEST_ASSERT(vs->isCompatible(*dvv.spmdSpace()));
  TEST_EQUALITY_CONST(dvv.globalOffset(), procRank*g_localDim);
  TEST_EQUALITY(dvv.subDim(), g_localDim);
  TEST_ASSERT(!is_null(dvv.values()));
  TEST_EQUALITY(dvv.values().size(), g_localDim);
  TEST_EQUALITY_CONST(dvv.stride(), 1);
  {
    out << "\nTest that dvv[i] == 2.0 ... ";
    bool local_success = true;
    for (Teuchos_Ordinal i = 0; i < dvv.subDim(); ++i) {
      TEST_ARRAY_ELE_EQUALITY( dvv, i, as<Scalar>(2.0) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }
  // ToDo: Write a better test than this because this does not give me a lot
  // of confidence.
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ConstDetachedSpmdVectorView,
  basic )


//
// DetachedSpmdVectorView unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DetachedSpmdVectorView, construct_null,
  Scalar )
{
  ECHO(DetachedSpmdVectorView<Scalar> dvv(null));
  TEST_EQUALITY_CONST(dvv.globalOffset(), 0);
  TEST_EQUALITY_CONST(dvv.subDim(), 0);
  TEST_EQUALITY_CONST(dvv.values(), null);
  TEST_EQUALITY_CONST(dvv.stride(), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DetachedSpmdVectorView, construct_null )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DetachedSpmdVectorView, basic,
  Scalar )
{
  const int procRank = Teuchos::GlobalMPISession::getRank();
  ECHO(const RCP<const VectorSpaceBase<Scalar> >
    vs = createSpmdVectorSpace<Scalar>(g_localDim));
  ECHO(RCP<VectorBase<Scalar> > v = createMember(vs));
  {
    ECHO(DetachedSpmdVectorView<Scalar> dvv(v));
    TEST_INEQUALITY(dvv.spmdSpace(), null);
    TEST_ASSERT(vs->isCompatible(*dvv.spmdSpace()));
    TEST_EQUALITY_CONST(dvv.globalOffset(), procRank*g_localDim);
    TEST_EQUALITY(dvv.subDim(), g_localDim);
    TEST_ASSERT(!is_null(dvv.values()));
    TEST_EQUALITY(dvv.values().size(), g_localDim);
    TEST_EQUALITY_CONST(dvv.stride(), 1);
    for (Teuchos_Ordinal i = 0; i < dvv.subDim(); ++i) {
      dvv[i] = as<Scalar>(2.0); // ToDo: Do something better here!
    }
  }
  {
    ECHO(ConstDetachedSpmdVectorView<Scalar> dvv(v));
    out << "\nTest that dvv[i] == 2.0 ... ";
    bool local_success = true;
    for (Teuchos_Ordinal i = 0; i < dvv.subDim(); ++i) {
      TEST_ARRAY_ELE_EQUALITY( dvv, i, as<Scalar>(2.0) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }
  // ToDo: Write a better test than this because this does not give me a lot
  // of confidence.
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DetachedSpmdVectorView,
  basic )


} // namespace
