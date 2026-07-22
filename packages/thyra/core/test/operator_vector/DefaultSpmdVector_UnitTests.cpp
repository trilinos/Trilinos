// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inoutArg;
using Teuchos::outArg;


int g_localDim = 3;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Local dimension of each vector." );
}


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVector, getMultiVectorLocalData,
  Scalar )
{
  out << "Test that we can grab MV data from Vector ...\n";

  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  const int procRank = Teuchos::DefaultComm<Ordinal>::getComm()->getRank();
  RCP<VectorBase<Scalar> > v = createMember(*vs);
  const ScalarMag tol = 100.0*ScalarTraits<Scalar>::eps();
  const Ordinal globalOffset = procRank * g_localDim;

  out << "Get non-const MV local data and set it ...\n";
  {
    ArrayRCP<Scalar> localValues;
    Ordinal leadingDim = -1;
    rcp_dynamic_cast<SpmdMultiVectorBase<Scalar> >(v,true)->getNonconstLocalData(
      outArg(localValues), outArg(leadingDim));
    TEST_EQUALITY(localValues.size(), g_localDim);
    TEST_EQUALITY(leadingDim, g_localDim);
    for (int i = 0; i < localValues.size(); ++i) {
      localValues[i] = globalOffset + i + 1;
    } 
  }
  const Ordinal n = vs->dim();
  TEST_FLOATING_EQUALITY(sum<Scalar>(*v), as<Scalar>((n*(n+1))/2.0), tol);

  out << "Get const MV local data and check it ...\n";
  {
    ArrayRCP<const Scalar> localValues;
    Ordinal leadingDim = -1;
    rcp_dynamic_cast<const SpmdMultiVectorBase<Scalar> >(v,true)->getLocalData(
      outArg(localValues), outArg(leadingDim));
    TEST_EQUALITY(localValues.size(), g_localDim);
    TEST_EQUALITY(leadingDim, g_localDim);
    for (int i = 0; i < localValues.size(); ++i) {
      TEST_EQUALITY(localValues[i], as<Scalar>(globalOffset + i + 1));
    } 
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVector,
  getMultiVectorLocalData )


} // namespace Thyra
