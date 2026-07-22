// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TeuchosCore_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_getBaseObjVoidPtr.hpp"
#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RCPBoostSharedPtrConversions.hpp"
#endif
#ifdef HAVE_TEUCHOSCORE_CXX11
#  include "Teuchos_RCPStdSharedPtrConversions.hpp"
#endif

#include "TestClasses.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::as;
using Teuchos::null;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcpFromUndefRef;
using Teuchos::outArg;
using Teuchos::rcpWithEmbeddedObj;
using Teuchos::getEmbeddedObj;
using Teuchos::getOptionalEmbeddedObj;
using Teuchos::getOptionalNonconstEmbeddedObj;
using Teuchos::set_extra_data;
using Teuchos::get_optional_nonconst_extra_data;
using Teuchos::getConst;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::DuplicateOwningRCPError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::RCPNodeTracer;


TEUCHOS_UNIT_TEST( DeallocNull, free )
{
  Teuchos::DeallocNull<A> d;
  d.free(0);
}


TEUCHOS_UNIT_TEST( RCP, construct_default )
{
  RCP<A> a_rcp;
  TEST_EQUALITY_CONST(a_rcp.get(), 0);
  TEST_EQUALITY_CONST(a_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), false);
}


TEUCHOS_UNIT_TEST( RCP, construct_nonull )
{
  RCP<A> a_rcp(new A);
  TEST_INEQUALITY_CONST(a_rcp.get(), 0);
  TEST_INEQUALITY_CONST(a_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 1);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), true);
}


TEUCHOS_UNIT_TEST( RCP, copy_construct_null )
{
  RCP<A> a_rcp;
  RCP<A> b_rcp(a_rcp);
  TEST_EQUALITY_CONST(a_rcp.get(), 0);
  TEST_EQUALITY_CONST(a_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), false);
  TEST_EQUALITY_CONST(b_rcp.get(), 0);
  TEST_EQUALITY_CONST(b_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(b_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(b_rcp.strong_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.has_ownership(), false);
}


TEUCHOS_UNIT_TEST( RCP, copy_construct_nonnull )
{
  RCP<A> a_rcp(new A);
  RCP<A> b_rcp(a_rcp);
  TEST_EQUALITY(b_rcp.get(), a_rcp.get());
  TEST_EQUALITY(b_rcp.getRawPtr(), a_rcp.getRawPtr());
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(b_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 2);
  TEST_EQUALITY_CONST(b_rcp.strong_count(), 2);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), true);
  TEST_EQUALITY_CONST(b_rcp.has_ownership(), true);
}


TEUCHOS_UNIT_TEST( RCP, move_construct_nonnull )
{
  RCP<A> a_rcp(new A);
  RCP<A> b_rcp(std::move(a_rcp));
  TEST_EQUALITY_CONST(a_rcp.get(), 0);
  TEST_EQUALITY_CONST(a_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), false);
  TEST_INEQUALITY_CONST(b_rcp.get(), 0);
  TEST_INEQUALITY_CONST(b_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(b_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(b_rcp.strong_count(), 1);
  TEST_EQUALITY_CONST(b_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.has_ownership(), true);
}


TEUCHOS_UNIT_TEST( RCP, assign_self_null )
{
  RCP<A> a_rcp;
  a_rcp = a_rcp;
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, copy_assign_self_nonnull )
{
  RCP<A> a_rcp(new A);
  A *a_raw_ptr = a_rcp.getRawPtr();
  a_rcp = a_rcp;
  TEST_ASSERT(nonnull(a_rcp));
  TEST_EQUALITY(a_rcp.getRawPtr(), a_raw_ptr);
}


TEUCHOS_UNIT_TEST( RCP, copy_assign_nonnull )
{
  RCP<A> a_rcp(new A);
  RCP<A> b_rcp;
  b_rcp = a_rcp;
  TEST_EQUALITY(b_rcp.get(), a_rcp.get());
  TEST_EQUALITY(b_rcp.getRawPtr(), a_rcp.getRawPtr());
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(b_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 2);
  TEST_EQUALITY_CONST(b_rcp.strong_count(), 2);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), true);
  TEST_EQUALITY_CONST(b_rcp.has_ownership(), true);
}


TEUCHOS_UNIT_TEST( RCP, move_assign_nonnull )
{
  RCP<A> a_rcp(new A);
  RCP<A> b_rcp;
  b_rcp = std::move(a_rcp);
  TEST_EQUALITY_CONST(a_rcp.get(), 0);
  TEST_EQUALITY_CONST(a_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(a_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(a_rcp.strong_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(a_rcp.has_ownership(), false);
  TEST_INEQUALITY_CONST(b_rcp.get(), 0);
  TEST_INEQUALITY_CONST(b_rcp.getRawPtr(), 0);
  TEST_EQUALITY_CONST(b_rcp.strength(), RCP_STRONG);
  TEST_EQUALITY_CONST(b_rcp.strong_count(), 1);
  TEST_EQUALITY_CONST(b_rcp.weak_count(), 0);
  TEST_EQUALITY_CONST(b_rcp.has_ownership(), true);
}


TEUCHOS_UNIT_TEST( RCP, operator_bool_null )
{
  RCP<A> a_rcp;
  TEST_EQUALITY_CONST(bool(a_rcp), false);
}


TEUCHOS_UNIT_TEST( RCP, operator_bool_nonnull )
{
  RCP<A> a_rcp(new A);
  TEST_EQUALITY_CONST(bool(a_rcp), true);
}


TEUCHOS_UNIT_TEST( RCP, operator_bool_if )
{
  RCP<A> a_rcp(new A);
  RCP<A> b_rcp;
  bool a_is_null;
  bool b_is_null;

  if (a_rcp)
    a_is_null = false;
  else
    a_is_null = true;
  if (b_rcp)
    b_is_null = false;
  else
    b_is_null = true;
 
  TEST_EQUALITY_CONST(a_is_null, false);
  TEST_EQUALITY_CONST(b_is_null, true );
}


TEUCHOS_UNIT_TEST( RCP, getConst )
{
  RCP<A> a_rcp(new A);
  RCP<const A> ca_rcp = a_rcp.getConst();
  TEST_EQUALITY(a_rcp.getRawPtr(), ca_rcp.getRawPtr());
}


TEUCHOS_UNIT_TEST( RCP, explicit_null )
{
  RCP<A> a_rcp(0);
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, explicit_dealloc_null )
{
  RCP<A> a_rcp = rcpWithDealloc(static_cast<A*>(0), Teuchos::DeallocNull<A>(), false);
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, explicit_null_null )
{
  RCP<A> a_rcp(0, null);
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, explicit_null_nonnull )
{
  A *a = new A;
  RCP<A> a_rcp(a, null);
  TEST_ASSERT(nonnull(a_rcp));
  delete a;
}


TEUCHOS_UNIT_TEST( RCP, rcpFromRef_raw_ref )
{
  A a;
  RCP<A> a_rcp = rcpFromRef(a);
  TEST_EQUALITY(a_rcp.getRawPtr(), &a);
  TEST_ASSERT(nonnull(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, rcpFromRef_from_rcp )
{
  RCP<A> a_rcp1 = rcp<A>(new A);
  RCP<A> a_rcp2 = rcpFromRef(*a_rcp1);
  TEST_EQUALITY(a_rcp2.getRawPtr(), a_rcp1.getRawPtr());
  if (RCPNodeTracer::isTracingActiveRCPNodes())
  {
    TEST_EQUALITY_CONST(a_rcp2.strong_count(), 1);
    TEST_EQUALITY_CONST(a_rcp2.weak_count(), 1);
    TEST_EQUALITY_CONST(a_rcp2.has_ownership(), true);
  }
  else {
    TEST_EQUALITY_CONST(a_rcp2.strong_count(), 1);
    TEST_EQUALITY_CONST(a_rcp2.weak_count(), 0);
    TEST_EQUALITY_CONST(a_rcp2.has_ownership(), false);
  }
}


TEUCHOS_UNIT_TEST( RCP, rcpFromUndefRef )
{
  A a;
  RCP<A> a_rcp = rcpFromUndefRef(a);
  TEST_ASSERT(nonnull(a_rcp));
}

/**
 * @test Test @ref Teuchos::make_rcp without constructor argument.
 */
TEUCHOS_UNIT_TEST( RCP, make_rcp_no_constructor_arg )
{
  Teuchos::RCP<A> a_rcp = Teuchos::make_rcp<A>();
  TEST_ASSERT(Teuchos::nonnull(a_rcp));
}

/**
 * @test Test @ref Teuchos::make_rcp with constructor arguments.
 */
TEUCHOS_UNIT_TEST( RCP, make_rcp )
{
  Teuchos::RCP<A> a_rcp = Teuchos::make_rcp<A>(1,2);
  TEST_ASSERT(Teuchos::nonnull(a_rcp));
  TEST_EQUALITY(a_rcp->A_g(),1);
  TEST_EQUALITY(a_rcp->A_f(),2);
}


//
// Test rcpCloneNode(...)
//


TEUCHOS_UNIT_TEST( RCP, rcpCloneNode_null )
{
  ECHO(RCP<RCP<int> > rcp1 = null);
  ECHO(RCP<RCP<int> > rcp2 = rcpCloneNode(rcp1));
  TEST_EQUALITY(rcp2, null);
}


TEUCHOS_UNIT_TEST( RCP, rcpCloneNode_basic )
{

  ECHO(RCP<int> rcp1 = rcp(new int(0)));

  ECHO(RCP<int> rcp2 = rcpCloneNode(rcp1));
  TEST_ASSERT(nonnull(rcp2));
  TEST_EQUALITY(rcp1.strong_count(), 2);
  TEST_EQUALITY(rcp2.strong_count(), 1);

  ECHO(RCP<int> rcp3 = rcp2);
  TEST_EQUALITY(rcp1.strong_count(), 2);
  TEST_EQUALITY(rcp2.strong_count(), 2);
  TEST_EQUALITY(rcp3.strong_count(), 2);

  ECHO(RCP <int> rcp4 = rcp1);
  TEST_EQUALITY(rcp1.strong_count(), 3);
  TEST_EQUALITY(rcp2.strong_count(), 2);
  TEST_EQUALITY(rcp3.strong_count(), 2);

  ECHO(rcp4 = null);
  TEST_EQUALITY(rcp1.strong_count(), 2);
  TEST_EQUALITY(rcp2.strong_count(), 2);
  TEST_EQUALITY(rcp3.strong_count(), 2);
  TEST_EQUALITY(rcp4.strong_count(), 0);

  ECHO(rcp1 = null);
  TEST_EQUALITY(rcp1.strong_count(), 0);
  TEST_EQUALITY(rcp2.strong_count(), 2);
  TEST_EQUALITY(rcp3.strong_count(), 2);
  TEST_EQUALITY(rcp4.strong_count(), 0);

  ECHO(rcp2 = null);
  TEST_EQUALITY(rcp2.strong_count(), 0);
  TEST_EQUALITY(rcp3.strong_count(), 1);

  ECHO(rcp3 = null);
  TEST_EQUALITY(rcp3.strong_count(), 0);

}


//
// Test duplicate owning RCP objects
//


// Test that shows that we can detect trying to create two owning RCPs
// pointing to the same polymorphic object but having different interfaces
// with different addresses.  This happens due to virtual base classes. Only
// works when we have a working getBaseObjVoidPtr(...) function.
TEUCHOS_UNIT_TEST( RCP, duplicate_rcp_owning_polymorphic )
{
  SET_RCPNODE_TRACING();
  ECHO(C *c_ptr = new C);
  ECHO(A *a_ptr = c_ptr);
  ECHO(RCP<C> c_rcp = rcp(c_ptr)); // Okay
#if defined(TEUCHOS_DEBUG) && defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR)
  // With determine they are pointed to the same object!
  TEST_THROW(RCP<A> a_rcp = rcp(a_ptr), DuplicateOwningRCPError);
#else
  // Will not determine they are point to the same object!
  ECHO(RCP<A> a_rcp = rcp(a_ptr));
  TEST_EQUALITY(a_rcp.getRawPtr(), a_ptr);
  ECHO(a_rcp.release()); // Better or we will get a segfault!
#endif
}


// Test that shows that we can detect trying to create two owning RCPs
// pointing to the same polymorphic object with the same type and therefore
// the same address.  This works even if these use virtual base classes.  This
// works even without a working getBaseObjVoidPtr(...) function.
TEUCHOS_UNIT_TEST( RCP, duplicate_rcp_owning_polymorphic_different_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(A *a_ptr1 = new C);
  ECHO(A *a_ptr2 = a_ptr1);
  ECHO(RCP<A> a_rcp1 = rcp(a_ptr1)); // Okay
#if defined(TEUCHOS_DEBUG)
  // With determine they are pointed to the same object!
  TEST_THROW(RCP<A> a_rcp2 = rcp(a_ptr2), DuplicateOwningRCPError);
#else
  // Will not determine they are point to the same object!
  ECHO(RCP<A> a_rcp2 = rcp(a_ptr2));
  TEST_EQUALITY(a_rcp2.getRawPtr(), a_ptr2);
  ECHO(a_rcp2.release()); // Better or we will get a segfault!
#endif
}


// Test that shows that we can always detect trying to create two owning RCPs
// pointing to the same nonpolymorphic object having different interfaces but
// the same address (single non-virtual inheritance).  Works just fine without
// a working getBaseObjVoidPtr(...) function.
TEUCHOS_UNIT_TEST( RCP, duplicate_rcp_owning_nonpolymorphic_same_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(E *e_ptr = new E);
  ECHO(E *d_ptr = e_ptr);
  ECHO(RCP<E> e_rcp = rcp(e_ptr)); // Okay
#if defined(TEUCHOS_DEBUG)
  // With determine they are pointed to the same object even without support
  // for getBaseObjVoidPtr(...) because no dynamic_cast is needed.
  TEST_THROW(RCP<D> d_rcp = rcp(d_ptr), DuplicateOwningRCPError);
#else
  // Will not determine they are point to the same object!
  ECHO(RCP<D> d_rcp = rcp(d_ptr));
  TEST_EQUALITY(d_rcp.getRawPtr(), d_ptr);
  ECHO(d_rcp.release()); // Better or we will get a segfault!
#endif
}


//
// These next tests shows that we can detect when two RCPs are create to the same
// object, one owning and the other non-owning.  When we have a working
// getBaseObjVoidPtr(...) function, the new non-owning RCP will actually be a
// weak RCP that can be used to detect circular dependencies.
//


// rcp


TEUCHOS_UNIT_TEST( RCP, rcp_duplicate_rcp_nonowning_polymorphic_different_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<C> c_rcp(new C));
  ECHO(A &a_ref = *c_rcp);
  ECHO(RCP<A> a_rcp = rcp(&a_ref, false));
  ECHO(c_rcp = null);
#if defined(TEUCHOS_DEBUG) && defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR)
  TEST_THROW(a_rcp->A_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(a_rcp.getRawPtr());
#endif
}


TEUCHOS_UNIT_TEST( RCP, rcp_duplicate_rcp_nonowning_polymorphic_same_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<A> a_rcp1(new C));
  ECHO(A &a_ref = *a_rcp1);
  ECHO(RCP<A> a_rcp2 = rcp(&a_ref, false));
  ECHO(a_rcp1 = null);
#if defined(TEUCHOS_DEBUG)
  TEST_THROW(a_rcp2->A_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(a_rcp2.getRawPtr());
#endif
}


TEUCHOS_UNIT_TEST( RCP, rcp_duplicate_rcp_nonowning_nonpolymorphic )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<E> e_rcp(new E));
  ECHO(D &d_ref = *e_rcp);
  ECHO(RCP<D> d_rcp = rcp(&d_ref, false));
  ECHO(e_rcp = null);
#if defined(TEUCHOS_DEBUG)
  TEST_THROW(d_rcp->D_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(d_rcp.getRawPtr());
#endif
}


// rcpFromRef


TEUCHOS_UNIT_TEST( RCP, rcpFromRef_duplicate_rcp_nonowning_polymorphic_different_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<C> c_rcp(new C));
  ECHO(A &a_ref = *c_rcp);
  ECHO(RCP<A> a_rcp = rcpFromRef(a_ref));
  ECHO(c_rcp = null);
#if defined(TEUCHOS_DEBUG) && defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR)
  TEST_THROW(a_rcp->A_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(a_rcp.getRawPtr());
#endif
}


TEUCHOS_UNIT_TEST( RCP, rcpFromRef_duplicate_rcp_nonowning_polymorphic_same_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<A> a_rcp1(new C));
  ECHO(A &a_ref = *a_rcp1);
  ECHO(RCP<A> a_rcp2 = rcpFromRef(a_ref));
  ECHO(a_rcp1 = null);
#if defined(TEUCHOS_DEBUG)
  TEST_THROW(a_rcp2->A_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(a_rcp2.getRawPtr());
#endif
}


TEUCHOS_UNIT_TEST( RCP, rcpFromRef_duplicate_rcp_nonowning_nonpolymorphic )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<E> e_rcp(new E));
  ECHO(D &d_ref = *e_rcp);
  ECHO(RCP<D> d_rcp = rcpFromRef(d_ref));
  ECHO(e_rcp = null);
#if defined(TEUCHOS_DEBUG)
  TEST_THROW(d_rcp->D_g(), DanglingReferenceError);
#else
  TEST_NOTHROW(d_rcp.getRawPtr());
#endif
}


// rcpFromUndefRef (Can never detect dangling references)


TEUCHOS_UNIT_TEST( RCP, rcpFromUndefRef_duplicate_rcp_nonowning_polymorphic_same_addr )
{
  SET_RCPNODE_TRACING();
  ECHO(RCP<A> a_rcp1(new C));
  ECHO(A &a_ref = *a_rcp1);
  ECHO(RCP<A> a_rcp2 = rcpFromUndefRef(a_ref));
  ECHO(a_rcp1 = null);
  TEST_NOTHROW(a_rcp2.getRawPtr());
}


//
// extra data and embedded objects tests
//


TEUCHOS_UNIT_TEST( RCP, get_optional_nonconst_extra_data )
{
  RCP<A> a_rcp = rcp(new A);
  set_extra_data( as<int>(1), "blob", outArg(a_rcp) );
  TEST_EQUALITY_CONST(*get_optional_nonconst_extra_data<int>(a_rcp, "blob"), as<int>(1));
}


TEUCHOS_UNIT_TEST( RCP, getOptionalEmbeddedObj_null )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  const Ptr<const RCP<C> > c_ptr_rcp_1 =
    getOptionalEmbeddedObj<A, RCP<C> >(a_rcp);
  TEST_EQUALITY_CONST( c_ptr_rcp_1, null );
  const Ptr<RCP<C> > c_ptr_rcp_2 =
    getOptionalNonconstEmbeddedObj<A, RCP<C> >(a_rcp);
  TEST_EQUALITY_CONST( c_ptr_rcp_2, null );
}


TEUCHOS_UNIT_TEST( RCP, getOptionalEmbeddedObj_default )
{

  ECHO(RCP<C> c_rcp = rcp(new C));
  ECHO(RCP<A> a_rcp = rcpWithEmbeddedObj(new A, c_rcp));

  Ptr<const RCP<C> > c_ptr_rcp_1 =
    getOptionalEmbeddedObj<A, RCP<C> >(a_rcp);
  TEST_EQUALITY_CONST( is_null(c_ptr_rcp_1), false );
  TEST_EQUALITY( (*c_ptr_rcp_1).getRawPtr(), c_rcp.getRawPtr() );
  TEST_EQUALITY( (*c_ptr_rcp_1)->C_g(), C_g_return );

  Ptr<RCP<C> > c_ptr_rcp_2 =
    getOptionalNonconstEmbeddedObj<A, RCP<C> >(a_rcp);
  TEST_EQUALITY_CONST( is_null(c_ptr_rcp_2), false );
  TEST_EQUALITY( (*c_ptr_rcp_2).getRawPtr(), c_rcp.getRawPtr() );
  TEST_EQUALITY( (*c_ptr_rcp_2)->C_f(), C_f_return );

}


TEUCHOS_UNIT_TEST( RCP, reset_null )
{
  RCP<A> a_rcp = rcp(new A);
  a_rcp.reset();
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, reset_nonnull )
{
  RCP<A> a_rcp = rcp(new A);
  C* c_rawp = new C;
  a_rcp.reset(c_rawp);
  A* a_rawp = c_rawp;
  TEST_EQUALITY( a_rcp.getRawPtr(), a_rawp );
}


TEUCHOS_UNIT_TEST( RCP, nonnull )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  TEST_EQUALITY_CONST(is_null(a_rcp), false);
  TEST_EQUALITY_CONST(nonnull(a_rcp), true);
  ECHO(a_rcp = null);
  TEST_EQUALITY_CONST(is_null(a_rcp), true);
  TEST_EQUALITY_CONST(nonnull(a_rcp), false);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RCP, weakDelete, T )
{

  ECHO(RCP<T> rcp_strong = rcp(new T));

  TEST_EQUALITY_CONST( rcp_strong.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( rcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.total_count(), 1 );

  ECHO(RCP<T> rcp_weak1 = rcp_strong.create_weak());

  TEST_EQUALITY_CONST( rcp_weak1.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( rcp_weak1.is_null(), false );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp_weak1.weak_count(), 1 );
  TEST_EQUALITY_CONST( rcp_weak1.total_count(), 2 );

  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( rcp_strong.weak_count(), 1 );
  TEST_EQUALITY_CONST( rcp_strong.total_count(), 2 );

  TEST_EQUALITY_CONST( rcp_weak1.shares_resource(rcp_strong), true );

  TEST_EQUALITY( rcp_weak1.get(), rcp_weak1.getRawPtr() );
  TEST_EQUALITY( rcp_weak1.get(), rcp_strong.get() );
  TEST_EQUALITY( rcp_weak1.getRawPtr(), rcp_strong.getRawPtr() );

  ECHO(RCP<T> rcp_weak2 = rcp_weak1);

  TEST_EQUALITY_CONST( rcp_weak2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( rcp_weak2.is_null(), false );
  TEST_EQUALITY_CONST( rcp_weak2.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak2.total_count(), 3 );

  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( rcp_strong.weak_count(), 2 );
  TEST_EQUALITY_CONST( rcp_strong.total_count(), 3 );

  TEST_EQUALITY_CONST( rcp_weak1.shares_resource(rcp_strong), true );
  TEST_EQUALITY_CONST( rcp_weak1.shares_resource(rcp_weak2), true );
  TEST_EQUALITY_CONST( rcp_weak2.shares_resource(rcp_strong), true );

  TEST_EQUALITY( rcp_weak2.get(), rcp_strong.get() );
  TEST_EQUALITY( rcp_weak2.getRawPtr(), rcp_strong.getRawPtr() );

  ECHO(rcp_strong = null); // This deletes the underlying object of type T!

  TEST_EQUALITY_CONST( rcp_strong.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( rcp_strong.is_null(), true );
  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.total_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( rcp_strong.shares_resource(rcp_weak1), false );
  TEST_EQUALITY_CONST( rcp_strong.shares_resource(rcp_weak2), false );

  TEST_EQUALITY_CONST( rcp_weak1.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.weak_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak1.total_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak1.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( rcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak2.total_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( rcp_weak1.shares_resource(rcp_weak2), true );

  ECHO(rcp_weak1.assert_not_null()); // Does not throw!
  ECHO(rcp_weak2.assert_not_null()); // Does not throw!

  TEST_THROW( rcp_weak1.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  TEST_THROW( rcp_weak1.operator->(), DanglingReferenceError );
  TEST_THROW( *rcp_weak1, DanglingReferenceError );
  TEST_THROW( rcp_weak1.create_weak(), DanglingReferenceError );
  TEST_THROW( rcp_weak1.get(), DanglingReferenceError );
  TEST_THROW( rcp_weak1.getRawPtr(), DanglingReferenceError );
  TEST_THROW( rcp_weak1(), DanglingReferenceError );
  TEST_THROW( rcp_weak1.release(), DanglingReferenceError );
#endif // TEUCHOS_DEBUG

  ECHO(rcp_weak1 = null); // Just deicrements weak count!

  TEST_EQUALITY_CONST( rcp_weak1.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( rcp_weak1.is_null(), true );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.weak_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.total_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( rcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak2.weak_count(), 1 );
  TEST_EQUALITY_CONST( rcp_weak2.total_count(), 1 );
  TEST_EQUALITY_CONST( rcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( rcp_weak1.shares_resource(rcp_weak2), false );

  TEST_THROW( rcp_weak2.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  TEST_THROW( rcp_weak2.operator->(), DanglingReferenceError );
  TEST_THROW( *rcp_weak2, DanglingReferenceError );
  TEST_THROW( rcp_weak2.create_weak(), DanglingReferenceError );
  TEST_THROW( rcp_weak2.get(), DanglingReferenceError );
  TEST_THROW( rcp_weak2.getRawPtr(), DanglingReferenceError );
  TEST_THROW( rcp_weak2(), DanglingReferenceError );
  TEST_THROW( rcp_weak2.release(), DanglingReferenceError );
#endif // TEUCHOS_DEBUG

}


TEUCHOS_UNIT_TEST( RCP, weak_strong )
{

  ECHO(RCP<A> rcp1(rcp(new A)));
  TEST_EQUALITY_CONST( rcp1.strength(), RCP_STRONG );

  ECHO(RCP<A> rcp2 = rcp1.create_weak());

  TEST_EQUALITY_CONST( rcp2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( rcp1.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( rcp2.strong_count(), 1 );
  TEST_EQUALITY_CONST( rcp2.weak_count(), 1 );

  ECHO(RCP<A> rcp3 = rcp2.create_strong());

  TEST_EQUALITY_CONST( rcp3.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( rcp1.strong_count(), 2 );
  TEST_EQUALITY_CONST( rcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( rcp2.strong_count(), 2 );
  TEST_EQUALITY_CONST( rcp2.weak_count(), 1 );

  // This will make the underlying object A gets deleted!
  ECHO(rcp1 = null);
  ECHO(rcp3 = null);

  ECHO(rcp2 = null); // Should make the underlying node go away

}


//
// circularReference
//


TEUCHOS_UNIT_TEST( RCP, circularReference_a_then_c )
{
  {
    // Create objects a and c
    ECHO(RCP<A> a = rcp(new A));
    ECHO(RCP<C> c = rcp(new C));

    // Create a circular reference where 'a' owns 'c' strongly but 'c' only
    // owns 'a' weakly.
    ECHO(a->set_C(c));
    ECHO(c->set_A(a.create_weak()));

    TEST_EQUALITY( a->call_C_f(), C_f_return );
    TEST_EQUALITY( c->call_A_g(), A_g_return );

    // Remove 'a' first and then remove 'c'.  Since 'a' is only weakly held by
    // 'c', this will result in 'a' being deleted right away.  In this case,
    // if anyone tries to access 'a' after this debug mode will throw.

    ECHO(a = null);

    // c is now dead!
#ifdef TEUCHOS_DEBUG
    TEST_THROW(c->call_A_g(), DanglingReferenceError);
#endif

    // Finally, when 'c' goes away implicitly, it will do nothing because it
    // only had weak ownership of a. So it will clean up with any trouble.
  }
}


TEUCHOS_UNIT_TEST( RCP, circularReference_c_then_a )
{
  {
    // Create objects a and c
    ECHO(RCP<A> a = rcp(new A));
    ECHO(RCP<C> c = rcp(new C));

    // Create a circular reference where 'a' owns 'c' strongly but 'c' only
    // owns 'a' weakly.
    ECHO(a->set_C(c));
    ECHO(c->set_A(a.create_weak()));

    TEST_EQUALITY( a->call_C_f(), C_f_return );
    TEST_EQUALITY( c->call_A_g(), A_g_return );

    // Remove 'c' first and then remove 'a' implicitly at the end of the
    // block.  Since 'c' is held strongly by 'a' and since we are keeping the
    // strong pointer for 'a' alive, we can call functions on 'a' all we want
    // with no fear of accessing dead memory.

    ECHO(c = null);

    TEST_EQUALITY( a->call_C_f(), C_f_return ); // a is still alive!

    // Finally, when 'a' goes away implicitly, it will take 'c' with it.  In
    // the complex set of nested calls that take place due to the circular
    // reference, everything will get cleaned up correctly.  Also, if any
    // client code where to try to access an object as it is being deleted, an
    // exception will get thrown and no memory error will occur (unless an
    // abort(...) is called when an exception gets thrown from a destructor
    // when an exception is already active).
  }
}


TEUCHOS_UNIT_TEST( RCP, circularReference_self )
{
  {
    // Create one 'c' object
    ECHO(RCP<C> c = rcp(new C));
    // Create a weak circular reference where 'c' points back to itself
    ECHO(c->set_A(c.create_weak()));
    // Now, try to set 'c' to null.
    ECHO(c = null); // All memory should be cleaned up here!
  }
}


TEUCHOS_UNIT_TEST( RCP, danglingPtr1 )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  ECHO(Ptr<A> a_ptr = a_rcp());
  ECHO(A *badPtr = a_rcp.getRawPtr());
  ECHO(a_rcp = null);
#ifdef TEUCHOS_DEBUG
  TEST_THROW( *a_ptr, DanglingReferenceError );
  (void)badPtr;
#else
  TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
#endif
}


TEUCHOS_UNIT_TEST( RCP, danglingPtr2 )
{
  ECHO(Ptr<A> a_ptr);
  ECHO(A *badPtr = 0);
  {
    ECHO(RCP<A> a_rcp = rcp(new A));
    ECHO(badPtr = a_rcp.getRawPtr());
    ECHO(a_ptr = a_rcp.ptr());
    TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW( *a_ptr, DanglingReferenceError );
  (void)badPtr;
#else
  TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
#endif
}


TEUCHOS_UNIT_TEST( RCP, danglingPtr3 )
{
  ECHO(Ptr<A> a_ptr);
  ECHO(A *badPtr = 0);
  {
    ECHO(RCP<A> a_rcp = rcp(new A));
    ECHO(badPtr = a_rcp.getRawPtr());
    ECHO(Ptr<A> a_ptr2(a_rcp.ptr()));
    ECHO(Ptr<A> a_ptr3(a_ptr2));
    ECHO(a_ptr = a_ptr3);
    TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW( *a_ptr, DanglingReferenceError );
  (void)badPtr;
#else
  TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
#endif
}


TEUCHOS_UNIT_TEST( RCP, danglingPtr4 )
{
  ECHO(Ptr<A> a_ptr);
  ECHO(A *badPtr = 0);
  {
    ECHO(RCP<C> c_rcp = rcp(new C));
    ECHO(badPtr = c_rcp.getRawPtr());
    ECHO(Ptr<A> a_ptr2(c_rcp.ptr()));
    ECHO(a_ptr = a_ptr2);
    TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW( *a_ptr, DanglingReferenceError );
  (void)badPtr;
#else
  TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
#endif
}


#ifdef TEUCHOS_DEBUG

/* ToDo: Comment this back in once I have everything working

// Test that the RCPNode tracing machinary can detect if an owning RCPNode is
// being created that would result in a double delete.
TEUCHOS_UNIT_TEST( RCP, multiRcpCreateError )
{
  C *c_ptr = new C;
#if !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
  Teuchos::setTracingActiveRCPNodes(true);
#endif
  RCP<C> c_rcp = rcp(c_ptr); // Okay
  RCP<C> c_rcp2;
  TEST_THROW(c_rcp2 = rcp(c_ptr), DuplicateOwningRCPError);
#if !defined(HAVE_TEUCHOS_DEBUG_RCP_NODE_TRACING)
  Teuchos::setTracingActiveRCPNodes(false);
#endif
  // Clean up memory so no leaks and not double deletes no matter what.
  c_rcp.release();
  c_rcp2.release();
  delete c_ptr;
}

*/

#endif // TEUCHOS_DEBUG


//
// invertObjectOwnership
//


RCP<C> createCAFactory()
{
  RCP<C> c = rcp(new C);
  c->set_A(rcp(new A));
  return c;
}


RCP<A> createACFactory()
{
  RCP<C> c = createCAFactory();
  return Teuchos::rcpWithInvertedObjOwnership(c->get_A(), c);
}


RCP<C> extractCFromA(const RCP<A> &a)
{
  return Teuchos::getInvertedObjOwnershipParent<C>(a);
}


TEUCHOS_UNIT_TEST( RCP, invertObjectOwnership_basic )
{
  RCP<A> a = createACFactory();
  RCP<C> c = extractCFromA(a);
  TEST_EQUALITY_CONST( a.strong_count(), 1 );
  TEST_EQUALITY_CONST( c->get_A().strong_count(), 3 );
  TEST_ASSERT( !a.shares_resource(c->get_A()) );
  TEST_EQUALITY( a.getRawPtr(), c->get_A().getRawPtr() );
  TEST_EQUALITY( a->A_g(), A_g_return );
  TEST_EQUALITY( c->C_g(), C_g_return );
}


// This unit test shows that you can remove the RCP in the C object
// and the A object will still live on.
TEUCHOS_UNIT_TEST( RCP, invertObjectOwnership_remove_A )
{
  RCP<A> a = createACFactory();
  extractCFromA(a)->set_A(null);
  RCP<C> c = extractCFromA(a);
  TEST_EQUALITY_CONST( a.strong_count(), 1 );
  TEST_EQUALITY_CONST( c->get_A(), null );
  TEST_EQUALITY( a->A_g(), A_g_return );
  TEST_EQUALITY( c->C_g(), C_g_return );
}


//
// createRCPWithBadDealloc
//


RCP<A> createRCPWithBadDealloc()
{
  return rcp(new A[1]); // Will use delete but should use delete []!
}


template<typename T>
class DeallocArrayDeleteExtraData {
public:
  static RCP<DeallocArrayDeleteExtraData<T> > create(T *ptr)
    { return rcp(new DeallocArrayDeleteExtraData(ptr)); }
  ~DeallocArrayDeleteExtraData() { delete [] ptr_; }
private:
  T *ptr_;
  DeallocArrayDeleteExtraData(T *ptr) : ptr_(ptr) {}
  // Not defined!
  DeallocArrayDeleteExtraData();
  DeallocArrayDeleteExtraData(const DeallocArrayDeleteExtraData&);
  DeallocArrayDeleteExtraData& operator=(const DeallocArrayDeleteExtraData&);
};


// This unit test shows how you can use extra data to fix a bad deallocation
// policy
TEUCHOS_UNIT_TEST( RCP, Fix_createRCPWithBadDealloc )
{
  using Teuchos::inOutArg;
  using Teuchos::set_extra_data;
  // Create object with bad deallocator
  RCP<A> a = createRCPWithBadDealloc();
  TEST_ASSERT(nonnull(a));
  // Disable default (incorrect) dealloc and set a new deallocation policy as extra data!
  a.release();
  set_extra_data( DeallocArrayDeleteExtraData<A>::create(a.getRawPtr()), "dealloc",
    inOutArg(a));
}


/**
 * @test Test the aliasing constructor.
 */
TEUCHOS_UNIT_TEST( RCP, aliasing_constructor )
{
  Teuchos::RCP<A> a_rcp = Teuchos::make_rcp<A>(1,2);
  Teuchos::RCP<A> b_rcp = Teuchos::make_rcp<A>(2,3);
  Teuchos::RCP<A> c_rcp(a_rcp, b_rcp.get());

  TEST_ASSERT(Teuchos::nonnull(c_rcp));
  TEST_EQUALITY_CONST(c_rcp.get(), b_rcp.get());
  TEST_EQUALITY(c_rcp->A_g(), b_rcp->A_g());
  TEST_EQUALITY(c_rcp->A_f(), b_rcp->A_f());
  TEST_ASSERT( c_rcp.shares_resource(a_rcp) );
}


//
// Test RCP/boost::shared_ptr covnersions
//


#ifdef HAVE_TEUCHOS_BOOST


TEUCHOS_UNIT_TEST( boost_shared_ptr, nonnull_is_null )
{
  using boost::shared_ptr;
  ECHO(shared_ptr<A> a_sptr(new A));
  TEST_EQUALITY_CONST(is_null(a_sptr), false);
  TEST_EQUALITY_CONST(nonnull(a_sptr), true);
  ECHO(a_sptr = shared_ptr<A>());
  TEST_EQUALITY_CONST(is_null(a_sptr), true);
  TEST_EQUALITY_CONST(nonnull(a_sptr), false);
}


#endif // HAVE_TEUCHOS_BOOST


//
// Test RCP/std::shared_ptr covnersions
//


#ifdef HAVE_TEUCHOSCORE_CXX11


TEUCHOS_UNIT_TEST( std_shared_ptr, nonnull_is_null )
{
  using std::shared_ptr;
  ECHO(shared_ptr<A> a_sptr(new A));
  TEST_EQUALITY_CONST(Teuchos::is_null(a_sptr), false);
  TEST_EQUALITY_CONST(Teuchos::nonnull(a_sptr), true);
  ECHO(a_sptr = shared_ptr<A>());
  TEST_EQUALITY_CONST(Teuchos::is_null(a_sptr), true);
  TEST_EQUALITY_CONST(Teuchos::nonnull(a_sptr), false);
}


TEUCHOS_UNIT_TEST( std_shared_ptr, convert_to_RCP_null )
{
  const std::shared_ptr<A> a_sptr1;
  const RCP<A> a_rsptr1 = rcp(a_sptr1);
  TEST_EQUALITY( a_rsptr1.get(), a_sptr1.get() );
  TEST_ASSERT(is_null(a_rsptr1));
  
  const std::shared_ptr<A> a_sptr2 = get_shared_ptr(a_rsptr1);
  TEST_EQUALITY( a_sptr1.get(), a_sptr1.get() );
  TEST_ASSERT(Teuchos::is_null(a_sptr2));
}


TEUCHOS_UNIT_TEST( std_shared_ptr, convert_to_RCP )
{
  const std::shared_ptr<A> a_sptr1(new C());
  const RCP<A> a_rsptr1 = rcp(a_sptr1);
  TEST_EQUALITY( a_rsptr1.get(), a_sptr1.get() );
  TEST_EQUALITY( a_rsptr1.getRawPtr(), a_sptr1.get() );
  TEST_EQUALITY( a_rsptr1.get(), a_rsptr1.getRawPtr() );

  const std::shared_ptr<A> a_sptr2 = get_shared_ptr(a_rsptr1);
  TEST_EQUALITY( a_sptr1.get(), a_sptr1.get() );
  // NOTE: There is no portable way to check that two std::shared_ptr objects
  // share the same underlying shared node :-(
}


TEUCHOS_UNIT_TEST( std_shared_ptr, convert_from_RCP_null )
{
  const RCP<A> a_rcp1;
  const std::shared_ptr<A> a_sptr1(get_shared_ptr(a_rcp1));
  TEST_EQUALITY( a_rcp1.get(), a_sptr1.get() );
  TEST_ASSERT(Teuchos::is_null(a_sptr1));

  const RCP<A> a_rcp2 = rcp(a_sptr1);
  TEST_EQUALITY(a_rcp2.get(), a_sptr1.get());
  TEST_ASSERT(is_null(a_rcp2));
}


TEUCHOS_UNIT_TEST( std_shared_ptr, convert_from_RCP )
{
  const RCP<A> a_rcp1(new C());
  const std::shared_ptr<A> a_sptr1(get_shared_ptr(a_rcp1));
  TEST_EQUALITY( a_rcp1.get(), a_sptr1.get() );
  TEST_EQUALITY( a_rcp1.getRawPtr(), a_sptr1.get() );
  TEST_EQUALITY( a_sptr1.get(), a_rcp1.getRawPtr() );

  const RCP<A> a_rcp2 = rcp(a_sptr1);
  TEST_EQUALITY( a_rcp2.get(), a_sptr1.get() );
  TEST_EQUALITY( a_rcp2.get(), a_rcp1.get() );
  TEST_ASSERT( a_rcp1.shares_resource(a_rcp2) );
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST( std_shared_ptr, convert_from_RCP_lookup_node )
{
  const RCP<A> a_rcp1(new C());
  const std::shared_ptr<A> a_sptr1(get_shared_ptr(a_rcp1));
  TEST_EQUALITY( a_sptr1.get(), a_rcp1.get() );

  const RCP<C> c_rcp1 = rcp(std::dynamic_pointer_cast<C>(a_sptr1));
  TEST_EQUALITY( c_rcp1.get(), a_rcp1.get() );
  TEST_ASSERT( c_rcp1.shares_resource(a_rcp1) );
  // NOTE: The above test shows how Teuchos::RCP machinery will automatically
  // look up the underlying RCPNode object and use it!
}


#endif // TEUCHOS_DEBUG


#endif // HAVE_TEUCHOSCORE_CXX11



//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( RCP, weakDelete, T ) \
  DEBUG_UNIT_TEST_GROUP(T)


UNIT_TEST_GROUP(A)
UNIT_TEST_GROUP(C)
UNIT_TEST_GROUP(D)
UNIT_TEST_GROUP(E)


} // namespace
