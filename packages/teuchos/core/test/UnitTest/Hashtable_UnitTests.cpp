// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
