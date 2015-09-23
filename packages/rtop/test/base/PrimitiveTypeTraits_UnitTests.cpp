/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


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
