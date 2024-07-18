// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_getConst.hpp"
#include "TestClasses.hpp"


namespace {


using Teuchos::null;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ptrFromRef;
using Teuchos::rcpFromPtr;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;


TEUCHOS_UNIT_TEST( Ptr, nonnull )
{
  ECHO(A a);
  ECHO(Ptr<A> a_ptr = ptrFromRef(a));
  TEST_EQUALITY_CONST(is_null(a_ptr), false);
  TEST_EQUALITY_CONST(nonnull(a_ptr), true);
  ECHO(a_ptr = null);
  TEST_EQUALITY_CONST(is_null(a_ptr), true);
  TEST_EQUALITY_CONST(nonnull(a_ptr), false);
}


TEUCHOS_UNIT_TEST( Ptr, getConst )
{
  RCP<A> a_rcp(new A);
  Ptr<A> a_ptr = a_rcp.ptr();
  Ptr<const A> ca_ptr = a_ptr.getConst();
  TEST_EQUALITY(a_ptr.getRawPtr(), ca_ptr.getRawPtr());
}


TEUCHOS_UNIT_TEST( Ptr, rcpFromPtr_weakRef )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  ECHO(Ptr<A> a_ptr = a_rcp.ptr());
  ECHO(RCP<A> a_rcp2 = rcpFromPtr(a_ptr));
  TEST_EQUALITY(a_rcp2.getRawPtr(), a_rcp.getRawPtr());
#ifdef TEUCHOS_DEBUG
  TEST_ASSERT(a_rcp2.shares_resource(a_rcp));
#else
  // In an optimized build, the object a_rcp2 has its own RCPNode object that
  // is unrelated to the orgininal a_rcp object.  This cuts down on overhead.
#endif
  ECHO(a_rcp = null);
#ifdef TEUCHOS_DEBUG
  TEST_THROW(a_ptr.getRawPtr(), DanglingReferenceError);
  TEST_THROW(a_rcp2.getRawPtr(), DanglingReferenceError);
#endif

}


TEUCHOS_UNIT_TEST( Ptr, rcpFromPtr_rawRef )
{
  ECHO(A a);
  ECHO(Ptr<A> a_ptr = ptrFromRef(a));
  ECHO(RCP<A> a_rcp2 = rcpFromPtr(a_ptr));
  TEST_EQUALITY(a_rcp2.getRawPtr(), &a);
}


TEUCHOS_UNIT_TEST( Ptr, rcpFromPtr_null )
{
  ECHO(Ptr<A> a_ptr);
  ECHO(RCP<A> a_rcp2 = rcpFromPtr(a_ptr));
  TEST_EQUALITY(a_rcp2, null);
}


} // namespace
