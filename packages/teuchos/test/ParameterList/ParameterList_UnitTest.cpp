// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {


/*
TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentNames ) {
  // Lists with different names should not be equal
  ParameterList A("Tom");
  ParameterList B("Bob");
  TEST_ASSERT( A != B );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( A != B );
}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentNames ) {
  ParameterList A("Julie");
  ParameterList B("Shannon");
  TEST_ASSERT( !haveSameValues(A,B) );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
}
*/


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityWithEmpty )
{
  // An empty list should not be equal to a full list
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( A == B );
  A.set("Hello","World");
  TEST_ASSERT( A != B );
  B.set("Hello","World");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentSublistNames )
{
  // Sublists with different names should not be equal
  ParameterList A;
  ParameterList B;
  A.sublist("Bob");
  B.sublist("Tom");
  TEST_ASSERT( A != B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentLengths )
{
  ParameterList A;
  ParameterList B;
  A.set("A","a");
  A.set("B","b");
  A.set("C","c");

  B.set("A","a");
  B.set("B","b");

  TEST_ASSERT( A != B );

  B.set("C","c");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesWithEmpty )
{
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( haveSameValues(A,B) );
  A.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
  B.set("Hello","World");
  TEST_ASSERT( haveSameValues(A,B) );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentSublistNames )
{
  ParameterList A;
  ParameterList B;
  A.sublist("Smith").set("People",4);
  B.sublist("Jones").set("People",4);
  TEST_ASSERT( !haveSameValues(A,B) ); // sublist names matter
}


} // namespace Teuchos



