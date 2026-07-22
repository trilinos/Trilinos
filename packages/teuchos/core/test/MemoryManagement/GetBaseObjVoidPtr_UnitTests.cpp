// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_getBaseObjVoidPtr.hpp"
#include "Teuchos_RCP.hpp"

#include "TestClasses.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, polymorphicClasses )
{
  RCP<C> c_ptr(new C);
  RCP<B1> b1_ptr = c_ptr;
  RCP<B2> b2_ptr = c_ptr;
  RCP<A> a_ptr = c_ptr;
  TEST_EQUALITY( getBaseObjVoidPtr(&*c_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*b1_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*b1_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*b2_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*b2_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*a_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*a_ptr), static_cast<const void*>(&*c_ptr) );
}


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, nonPolymorphicClasses )
{
  RCP<E> e_ptr(new E);
  RCP<D> d_ptr = e_ptr;
  TEST_EQUALITY( getBaseObjVoidPtr(&*e_ptr), static_cast<const void*>(&*e_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*d_ptr), static_cast<const void*>(&*e_ptr) );
}


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, nonPolymorphicBuiltInTypes )
{
  RCP<int> i_ptr(new int);
  TEST_EQUALITY( getBaseObjVoidPtr(&*i_ptr), static_cast<const void*>(&*i_ptr) );
}


} // namespace Teuchos namespace
