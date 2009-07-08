#include "Teuchos_RCP.hpp"
#include "Teuchos_getConst.hpp"
#ifdef HAVE_TEUCHOS_BOOST
#  include "Teuchos_RCPBoostSharedPtrConversions.hpp"
#endif

#include "TestClasses.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::as;
using Teuchos::null;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::rcpWithEmbeddedObj;
using Teuchos::getOptionalEmbeddedObj;
using Teuchos::getOptionalNonconstEmbeddedObj;
using Teuchos::set_extra_data;
using Teuchos::get_optional_nonconst_extra_data;
using Teuchos::getConst;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::RCP_STRENGTH_INVALID;


TEUCHOS_UNIT_TEST( DeallocNull, free )
{
  Teuchos::DeallocNull<A> d;
  d.free(0);
}


TEUCHOS_UNIT_TEST( RCP, assignSelf )
{
  RCP<A> a_rcp;
  a_rcp = a_rcp;
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, explicit_null )
{
  RCP<A> a_rcp(0, false);
  TEST_ASSERT(is_null(a_rcp));
}


TEUCHOS_UNIT_TEST( RCP, explicit_dealloc_null )
{
  RCP<A> a_rcp = rcp<A>(0, Teuchos::DeallocNull<A>(), false);
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


TEUCHOS_UNIT_TEST( RCP, reset )
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


#ifdef HAVE_TEUCHOS_BOOST


TEUCHOS_UNIT_TEST( shared_ptr, nonnull )
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

  TEST_EQUALITY_CONST( rcp_strong.strength(), RCP_STRENGTH_INVALID );
  TEST_EQUALITY_CONST( rcp_strong.is_null(), true );
  TEST_EQUALITY_CONST( rcp_strong.count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.total_count(), 0 );
  TEST_EQUALITY_CONST( rcp_strong.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( rcp_strong.shares_resource(rcp_weak1), false );
  TEST_EQUALITY_CONST( rcp_strong.shares_resource(rcp_weak2), false );

  TEST_EQUALITY_CONST( rcp_weak1.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak1.count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.weak_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak1.total_count(), 2 );
  TEST_EQUALITY_CONST( rcp_weak1.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( rcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak2.count(), 0 );
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
  TEST_THROW( rcp_weak1.ptr(), DanglingReferenceError );
  TEST_THROW( rcp_weak1.release(), DanglingReferenceError );
#endif // TEUCHOS_DEBUG

  ECHO(rcp_weak1 = null); // Just deicrements weak count!

  TEST_EQUALITY_CONST( rcp_weak1.strength(), RCP_STRENGTH_INVALID );
  TEST_EQUALITY_CONST( rcp_weak1.is_null(), true );
  TEST_EQUALITY_CONST( rcp_weak1.count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.weak_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.total_count(), 0 );
  TEST_EQUALITY_CONST( rcp_weak1.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( rcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( rcp_weak2.count(), 0 );
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
  TEST_THROW( rcp_weak2.ptr(), DanglingReferenceError );
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


TEUCHOS_UNIT_TEST( RCP, circularReference_a_then_c )
{

  //TEST_EQUALITY_CONST(Teuchos::numActiveRCPNodes(), 0);

  {

    // Create objects a and c

    ECHO(RCP<A> a = rcp(new A));
    ECHO(RCP<C> c = rcp(new C));

    // Create a circular reference where 'a' owns 'c' strongly but 'c' only
    // owns 'a' weakly.

    ECHO(a->set_C(c));
    ECHO(c->set_A(a.create_weak()));

#ifdef TEUCHOS_DEBUG
    ECHO(c->call_A_on_delete(true));
    // Here, we set 'c' to call 'a' when it is deleted which will result in an
    // exception being thrown in a call to delete.  NOTE: It is *very* bad
    // practice to allow exceptions to be thrown from destructors but I am
    // allowing it so that I can detect such bad bahavior below!
#endif

    TEST_EQUALITY( a->call_C_f(), C_f_return );
    TEST_EQUALITY( c->call_A_g(), A_g_return );

    // Remove 'a' first and then remove 'c'.  Since 'a' is only weakly held by
    // 'c', this will result in 'a' being deleted right away.  In this case,
    // if anyone tries to access 'a' after this (like 'c' in its destructor),
    // then an exception will get thrown in debug mode!

    ECHO(a = null);

    // Now, remove 'c'.  In this case, since 'a' has already been deleted and
    // 'c' is going to try to call 'a' on its way out, this will thrown an
    // exception.

#ifdef TEUCHOS_DEBUG

    TEST_THROW(c = null, DanglingReferenceError);
    // NOTE: Above, operator==(...) exhibits the 'strong' guarantee?

    // Since an exception was thrown, the 'c' object never got deleted.
    // Therefore, we need to disable 'c' calling 'a' on delete and the object
    // will get cleaned up correctly when this function exists (I hope).
    ECHO(c->call_A_on_delete(false));

    ECHO(c = null); // All memory should be cleaned up here!

#endif // TEUCHOS_DEBUG

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

    ECHO(c->call_A_on_delete(false));
    // Here, we set 'c' to not call 'a' when it is deleted.  It turns out that
    // the set of calls to delete and destructors that takes place is very
    // complex and in order to avoid trouble, an object that holds an RCP to
    // another object weakly should *never* try to call any members on the
    // wrapped object as it gets deleted!
    
    TEST_EQUALITY( a->call_C_f(), C_f_return );
    TEST_EQUALITY( c->call_A_g(), A_g_return );

    // Remove 'c' first and then remove 'a' implicitly at the end of the
    // block.  Since 'c' is held strongly by 'a' and since we are keeping the
    // strong pointer for 'a' alive, we can call functions on 'a' all we want
    // with no fear of accessing dead memory.

    ECHO(c = null);

    TEST_EQUALITY( a->call_C_f(), C_f_return ); // C is still alive!

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


TEUCHOS_UNIT_TEST( RCP, danglingPtr )
{
  ECHO(RCP<A> a_rcp = rcp(new A));
  ECHO(Ptr<A> a_ptr = a_rcp.ptr());
  ECHO(A *badPtr = a_rcp.getRawPtr());
  ECHO(a_rcp = null);
#ifdef TEUCHOS_DEBUG
  TEST_THROW( *a_ptr, DanglingReferenceError );
  (void)badPtr;
#else
  TEST_EQUALITY( a_ptr.getRawPtr(), badPtr );
#endif
}


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
