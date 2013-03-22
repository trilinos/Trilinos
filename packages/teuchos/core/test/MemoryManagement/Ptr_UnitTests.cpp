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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

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
