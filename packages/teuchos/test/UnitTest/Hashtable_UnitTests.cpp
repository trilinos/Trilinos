/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Hashtable.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Hashtable, test0, Key, Value )
{
  using Teuchos::as;

  Teuchos::Hashtable<Key,Value> hashtable;

  hashtable.put(as<Key>(1), as<Value>(1));
  hashtable.put(as<Key>(3), as<Value>(9));
  hashtable.put(as<Key>(5), as<Value>(7));

  TEST_EQUALITY( hashtable.size(), 3 );

  TEST_EQUALITY( hashtable.containsKey(as<Key>(3)), true );
  TEST_EQUALITY( hashtable.containsKey(as<Key>(4)), false );

  TEST_EQUALITY( hashtable.get(as<Key>(5)), as<Value>(7) );
}

//
// Instantiations
//


#define UNIT_TEST_GROUP( K, V ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Hashtable, test0, K, V )

UNIT_TEST_GROUP(int, int)
UNIT_TEST_GROUP(int, float)

} // namespace
