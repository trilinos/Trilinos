// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_Types.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( PrimitiveTypeTraits, basic, Scalar, ConcreteObj )
{
  using Teuchos::Array;
  typedef RTOpPack::PrimitiveTypeTraits<Scalar,ConcreteObj> PTT;
  typedef typename PTT::primitiveType PrimitiveScalar;
  typedef Teuchos::ScalarTraits<ConcreteObj> CST;

  const ConcreteObj v1 = CST::random();

  Array<PrimitiveScalar> primitiveObjs(PTT::numPrimitiveObjs());
  Array<RTOpPack::index_type> indexObjs(PTT::numIndexObjs());
  Array<char> charObjs(PTT::numCharObjs());

  PTT::extractPrimitiveObjs( v1,
    primitiveObjs(), indexObjs(), charObjs() );

  ConcreteObj v2 = CST::zero();

  PTT::loadPrimitiveObjs( primitiveObjs(), indexObjs(), charObjs(),
    Teuchos::outArg(v2) );

  TEST_EQUALITY( v1, v2 );
}


#define UNIT_TEST_INSTANT_SCALAR_CONCRETEOBJ( SCALAR, CONCRETEOBJ ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( PrimitiveTypeTraits, basic, SCALAR, CONCRETEOBJ )


typedef RTOpPack::index_type index_type;


#define UNIT_TEST_INSTANT_SCALAR( SCALAR ) \
  UNIT_TEST_INSTANT_SCALAR_CONCRETEOBJ( SCALAR, SCALAR ) \
  UNIT_TEST_INSTANT_SCALAR_CONCRETEOBJ( SCALAR, index_type )

#ifdef HAVE_TEUCHOS_FLOAT
  UNIT_TEST_INSTANT_SCALAR(float)
#endif

UNIT_TEST_INSTANT_SCALAR(double)

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_FLOAT)
  typedef std::complex<float> ComplexFloat;
  UNIT_TEST_INSTANT_SCALAR( ComplexFloat )
  UNIT_TEST_INSTANT_SCALAR_CONCRETEOBJ( ComplexFloat, float )
#endif

#if defined(HAVE_TEUCHOS_COMPLEX)
  typedef std::complex<double> ComplexDouble;
  UNIT_TEST_INSTANT_SCALAR( ComplexDouble )
  UNIT_TEST_INSTANT_SCALAR_CONCRETEOBJ( ComplexDouble, double )
#endif


} // namespace
