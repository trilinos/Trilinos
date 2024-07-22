// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;


const int g_localDim = 4; // ToDo: Make variable!


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorProductVectorSpace, defaultConstruct,
  Scalar )
{

  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > vs =
    multiVectorProductVectorSpace<Scalar>();
  TEST_EQUALITY(vs->dim(), as<Ordinal>(-1));
  out << "vs = " << *vs;
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorProductVectorSpace,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorProductVectorSpace, standard,
  Scalar )
{

  const Ordinal numCols = 3;
  RCP<const VectorSpaceBase<Scalar> > spmdVs =
    createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > vs =
    multiVectorProductVectorSpace<Scalar>(spmdVs, numCols);
  TEST_EQUALITY(vs->dim(), (numCols * spmdVs->dim()));
  out << "vs = " << *vs;
  VectorSpaceTester<Scalar> vsTester;
  TEST_ASSERT(vsTester.check(*vs, &out));
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorProductVectorSpace,
  standard )


} // namespace Thyra
