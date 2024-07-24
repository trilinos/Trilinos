// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConstNonconstObjectContainer.hpp"

#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {

TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, create ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  RCP<std::vector<double> > vec = rcp(new std::vector<double> );
  vectorObj.initialize(vec); // nonconst
  TEST_ASSERT( vectorObj.isConst() == false );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, DefaultConstruct ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  TEST_ASSERT( vectorObj.isConst() == true );
  TEST_ASSERT( vectorObj.getConstObj() == null );
  // This does not throw an exception because the pointer is null
  TEST_ASSERT( vectorObj.getNonconstObj() == null );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, NonconstConstruct ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( vectorObj.isConst() == false );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  RCP<std::vector<double> > vec3 = vectorObj.getNonconstObj();
  TEST_ASSERT( vec == vec3 );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, ConstConstruct) {
  RCP<const std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( vectorObj.isConst() == true );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  RCP<std::vector<double> > vec3;
  TEST_THROW( vec3 = vectorObj.getNonconstObj(), NonconstAccessError );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, NonconstInitialize) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  vectorObj.initialize(vec);
  TEST_ASSERT( vectorObj.isConst() == false );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  RCP<std::vector<double> > vec3 = vectorObj.getNonconstObj();
  TEST_ASSERT( vec == vec3 );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, ConstInitialize) {
  RCP<const std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  vectorObj.initialize(vec);
  TEST_ASSERT( vectorObj.isConst() == true );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  RCP<std::vector<double> > vec3;
  TEST_THROW( vec3 = vectorObj.getNonconstObj(), NonconstAccessError );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, AssignmentFromRCP) {
  RCP<const std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  vectorObj = vec;
  TEST_ASSERT( vectorObj.isConst() == true );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  RCP<std::vector<double> > vec3;
  TEST_THROW( vec3 = vectorObj.getNonconstObj(), NonconstAccessError );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, uninitialize) {
  RCP<const std::vector<double> > vec = rcp(new std::vector<double> );
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( vectorObj.isConst() == true );
  RCP<const std::vector<double> > vec2 = vectorObj.getConstObj();
  TEST_ASSERT( vec == vec2 );
  vectorObj.uninitialize();
  TEST_ASSERT( vectorObj.isConst() == true );
  TEST_ASSERT( vectorObj.getConstObj() == null );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, parens ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  vectorObj.initialize(vec);
  RCP<const std::vector<double> > vec2 = vectorObj();
  TEST_ASSERT( vec == vec2 );
}


// operator->
TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, arrow ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  vec->push_back(25.0);
  vec->push_back(32.0);
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( vectorObj->size() == 2 );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, arrowEmpty ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
#ifdef TEUCHOS_DEBUG
  TEST_THROW( vectorObj->size(), NullReferenceError );
#endif
}


// operator*
TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, dereference ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  vec->push_back(25.0);
  vec->push_back(32.0);
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( (*vectorObj)[0] == 25.0 );
  TEST_ASSERT( (*vectorObj)[1] == 32.0 );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, dereferenceEmpty ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
#ifdef TEUCHOS_DEBUG
  TEST_THROW( (*vectorObj).size(), NullReferenceError );
#endif
}


// implicit cast
//    RCP<const ObjType>  <-
//    ConstNonconstObjectContainer<const ObjType>
TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, castToRCP ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  RCP<const std::vector<double> > vec2(vectorObj);
  TEST_ASSERT( vec == vec2 );
}


// implicit cast
//    RCP<const ObjType>  ->
//    ConstNonconstObjectContainer<const ObjType>
//    This is already done through the constructors on
//    ConstNonconstObjectContainer.
TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, castFromRCP ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( vectorObj.getConstObj() == vec );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, is_null ) {
  ConstNonconstObjectContainer<std::vector<double> > vectorObj;
  TEST_ASSERT( is_null(vectorObj) );
}


TEUCHOS_UNIT_TEST( ConstNonconstObjectContainer, nonnull ) {
  RCP<std::vector<double> > vec = rcp(new std::vector<double>);
  ConstNonconstObjectContainer<std::vector<double> > vectorObj(vec);
  TEST_ASSERT( nonnull(vectorObj) );
}


} // namespace Teuchos



