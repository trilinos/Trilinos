//@HEADER
// ***********************************************************************
//
//                           Teuchos Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER

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


} // namespace Teuchos



