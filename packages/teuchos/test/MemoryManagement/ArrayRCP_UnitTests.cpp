#include "Teuchos_UnitTestHarness.hpp"
#include "Array_UnitTest_helpers.hpp"
#include "Teuchos_ArrayRCP.hpp"


namespace {

using ArrayUnitTestHelpers::n;
using ArrayUnitTestHelpers::generateArray;

using Teuchos::null;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::ArrayView;
using Teuchos::getConst;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RangeError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::RCP_STRENGTH_INVALID;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, weakDelete, T )
{

  ECHO(ArrayRCP<T> arcp_strong = arcp<T>(n));

  TEST_EQUALITY_CONST( arcp_strong.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 1 );

  ECHO(ArrayRCP<T> arcp_weak1 = arcp_strong.create_weak());

  TEST_EQUALITY_CONST( arcp_weak1.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp_weak1.is_null(), false );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 2 );

  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 2 );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_strong), true );

  TEST_EQUALITY( arcp_weak1.get(), arcp_weak1.getRawPtr() );
  TEST_EQUALITY( arcp_weak1.get(), arcp_strong.get() );
  TEST_EQUALITY( arcp_weak1.getRawPtr(), arcp_strong.getRawPtr() );

  ECHO(ArrayRCP<T> arcp_weak2 = arcp_weak1);

  TEST_EQUALITY_CONST( arcp_weak2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp_weak2.is_null(), false );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 3 );

  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), false );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 3 );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_strong), true );
  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), true );
  TEST_EQUALITY_CONST( arcp_weak2.shares_resource(arcp_strong), true );

  TEST_EQUALITY( arcp_weak2.get(), arcp_strong.get() );
  TEST_EQUALITY( arcp_weak2.getRawPtr(), arcp_strong.getRawPtr() );

  ECHO(arcp_strong = null); // This deletes the underlying object of type T!

  TEST_EQUALITY_CONST( arcp_strong.strength(), RCP_STRENGTH_INVALID );
  TEST_EQUALITY_CONST( arcp_strong.is_null(), true );
  TEST_EQUALITY_CONST( arcp_strong.count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.total_count(), 0 );
  TEST_EQUALITY_CONST( arcp_strong.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( arcp_strong.shares_resource(arcp_weak1), false );
  TEST_EQUALITY_CONST( arcp_strong.shares_resource(arcp_weak2), false );

  TEST_EQUALITY_CONST( arcp_weak1.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak1.count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak1.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak2.count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 2 );
  TEST_EQUALITY_CONST( arcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), true );

  ECHO(arcp_weak1.assert_not_null()); // Does not throw!
  ECHO(arcp_weak2.assert_not_null()); // Does not throw!

  TEST_THROW( arcp_weak1.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  TEST_THROW( arcp_weak1.operator->(), DanglingReferenceError );
  TEST_THROW( *arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1.create_weak(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.get(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.getRawPtr(), DanglingReferenceError );
  TEST_THROW( arcp_weak1[0], DanglingReferenceError );
  TEST_THROW( ++arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1++, DanglingReferenceError );
  TEST_THROW( --arcp_weak1, DanglingReferenceError );
  TEST_THROW( arcp_weak1--, DanglingReferenceError );
  TEST_THROW( arcp_weak1+=1, DanglingReferenceError );
  TEST_THROW( arcp_weak1-=1, DanglingReferenceError );
  TEST_THROW( arcp_weak1+1, DanglingReferenceError );
  TEST_THROW( arcp_weak1-1, DanglingReferenceError );
  TEST_THROW( arcp_weak1.getConst(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.persistingView(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1.lowerOffset(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.upperOffset(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.size(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.begin(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.end(), DanglingReferenceError );
  TEST_THROW( arcp_weak1.view(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1(0,n), DanglingReferenceError );
  TEST_THROW( arcp_weak1(), DanglingReferenceError );
  TEST_THROW( {ArrayView<T> av = arcp_weak1;}, DanglingReferenceError );
  TEST_THROW( {ArrayRCP<const T> ap = getConst(arcp_weak1);},
    DanglingReferenceError );
  TEST_THROW( arcp_weak1.release(), DanglingReferenceError );
#endif // TEUCHOS_DEBUG

  ECHO(arcp_weak1 = null); // Just deicrements weak count!

  TEST_EQUALITY_CONST( arcp_weak1.strength(), RCP_STRENGTH_INVALID );
  TEST_EQUALITY_CONST( arcp_weak1.is_null(), true );
  TEST_EQUALITY_CONST( arcp_weak1.count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.weak_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.total_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak1.is_valid_ptr(), true );

  TEST_EQUALITY_CONST( arcp_weak2.has_ownership(), true );
  TEST_EQUALITY_CONST( arcp_weak2.count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.strong_count(), 0 );
  TEST_EQUALITY_CONST( arcp_weak2.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.total_count(), 1 );
  TEST_EQUALITY_CONST( arcp_weak2.is_valid_ptr(), false );

  TEST_EQUALITY_CONST( arcp_weak1.shares_resource(arcp_weak2), false );

  TEST_THROW( arcp_weak2.assert_valid_ptr(), DanglingReferenceError );
#ifdef TEUCHOS_DEBUG
  // ToDo: Fill in
#endif // TEUCHOS_DEBUG

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, danglingArrayView, T )
{
  ArrayView<T> av;
  {
    ArrayRCP<T> arcp1 = arcp<T>(n);
    av = arcp1();
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW( av.size(), DanglingReferenceError );
  TEST_THROW( av.toString(), DanglingReferenceError );
  TEST_THROW( av.getRawPtr(), DanglingReferenceError );
  TEST_THROW( av[0], DanglingReferenceError );
  TEST_THROW( av.front(), DanglingReferenceError );
  TEST_THROW( av.back(), DanglingReferenceError );
  TEST_THROW( av.view(0, n), DanglingReferenceError );
  TEST_THROW( av(0, n), DanglingReferenceError );
  TEST_THROW( av(), DanglingReferenceError );
  TEST_THROW( av.getConst(), DanglingReferenceError );
  TEST_THROW( av.begin(), DanglingReferenceError );
  TEST_THROW( av.end(), DanglingReferenceError );
#endif  
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, outOfBounds, T )
{
  ECHO(ArrayRCP<T> arcp1 = arcp<T>(n));
  TEST_THROW(arcp1(-1,n), RangeError);
  TEST_THROW(arcp1(0,n+1), RangeError);
  TEST_THROW(arcp1(0,-1), RangeError);
  TEST_THROW(arcp1.view(-1,n), RangeError);
  TEST_THROW(arcp1.view(0,n+1), RangeError);
  TEST_THROW(arcp1.view(0,-1), RangeError);
  TEST_THROW(arcp1.persistingView(-1,n), RangeError);
  TEST_THROW(arcp1.persistingView(0,n+1), RangeError);
  TEST_THROW(arcp1.persistingView(0,-1), RangeError);
}


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, outOfBounds, T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, weakDelete, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, danglingArrayView, T ) \
  DEBUG_UNIT_TEST_GROUP(T)


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(double)
UNIT_TEST_GROUP(float)


} // namespace
