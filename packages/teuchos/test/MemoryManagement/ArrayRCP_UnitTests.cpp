#include "Teuchos_UnitTestHarness.hpp"
#include "Array_UnitTest_helpers.hpp"
#include "TestClasses.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_getRawPtr.hpp"

namespace {

using ArrayUnitTestHelpers::n;

using Teuchos::getRawPtr;
using Teuchos::as;
using Teuchos::null;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::arcp_reinterpret_cast;
using Teuchos::ArrayView;
using Teuchos::getConst;
using Teuchos::NullReferenceError;
using Teuchos::DanglingReferenceError;
using Teuchos::RangeError;
using Teuchos::RCP_STRONG;
using Teuchos::RCP_WEAK;
using Teuchos::RCP_STRENGTH_INVALID;
using Teuchos::implicit_ptr_cast;
using Teuchos::getRawPtr;


//
// Templated unit tests
//


TEUCHOS_UNIT_TEST( ArrayRCP, memberPointer )
{
  ArrayRCP<A> a_arcp = arcp<A>(1);
  TEST_EQUALITY_CONST( a_arcp->A_f(), A_f_return );
}


TEUCHOS_UNIT_TEST( ArrayRCP, getConst_null )
{
  const ArrayRCP<A> a1_arcp;
  const ArrayRCP<const A> a2_arcp = a1_arcp.getConst();
  TEST_ASSERT(is_null(a2_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, operator_parenth_ArrayView_null )
{
  const ArrayRCP<A> a_arcp;
  const ArrayView<A> av = a_arcp();
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, operator_parenth_ArrayView_const_null )
{
  const ArrayRCP<const A> a_arcp;
  const ArrayView<const A> av = a_arcp();
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST( ArrayRCP, implicit_ArrayRCP_const )
{
  const ArrayRCP<A> a_arcp;
  const ArrayRCP<const A> ac_arcp = a_arcp;
  TEST_ASSERT(is_null(ac_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, release )
{
  ArrayRCP<A> a_arcp = arcp<A>(1);
  delete [] a_arcp.release();
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_null )
{
  ArrayRCP<A> a_arcp = arcp<A>(0, 0, -1, false);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_dealloc_null )
{
  ArrayRCP<A> a_arcp = arcp<A, Teuchos::DeallocNull<A> >(0, 0, -1,
    Teuchos::DeallocNull<A>(), false);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_vector_null )
{
  const RCP<std::vector<int> > v_rcp;
  const ArrayRCP<int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_const_vector_null )
{
  const RCP<const std::vector<int> > v_rcp;
  const ArrayRCP<const int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_vector_unsized )
{
  const RCP<std::vector<int> > v_rcp = rcp(new std::vector<int>);
  const ArrayRCP<int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, convert_from_const_vector_unsized )
{
  const RCP<const std::vector<int> > v_rcp = rcp(new std::vector<int>);
  const ArrayRCP<const int> a_arcp = arcp(v_rcp);
  TEST_ASSERT(is_null(a_arcp));
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcpWithEmbeddedObj )
{
  const ArrayRCP<const int> a_arcp =
    Teuchos::arcpWithEmbeddedObj<int>(new int[1], 0, 1, as<int>(1), true);
  const int embeddedObj = Teuchos::getEmbeddedObj<int,int>(a_arcp); 
  TEST_EQUALITY_CONST( embeddedObj, as<int>(1) );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, assignSelf, T )
{
  ArrayRCP<T> a_arcp;
  a_arcp = a_arcp;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, nullIterator, T )
{
  typedef ArrayRCP<T> iter_t;
  ArrayRCP<T> arcp1 = Teuchos::NullIteratorTraits<iter_t>::getNull();
  TEST_EQUALITY_CONST(arcp1, Teuchos::null); 
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, implicitConversions, T )
{

  ECHO(ArrayRCP<T> arcp1 = arcp<T>(n));
  ECHO(ArrayRCP<const T> arcp2 = arcp1);
  TEST_ASSERT(arcp1.shares_resource(arcp2));

}


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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, getRawPtr, T )
{
  ArrayRCP<const T> cptr;
  ArrayRCP<T> ptr;
  TEST_EQUALITY_CONST( getRawPtr(cptr), (const T*)NULL );
  TEST_EQUALITY_CONST( getRawPtr(ptr), (T*)NULL );
  cptr = arcp<T>(n);
  ptr  = arcp<T>(n);
  TEST_EQUALITY( getRawPtr(cptr), &cptr[0]);
  TEST_EQUALITY( getRawPtr(ptr),  &ptr[0] );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CPtr, getRawPtr, T )
{
  const T *cptr = NULL;
  T *ptr = NULL;
  TEST_EQUALITY_CONST( getRawPtr(cptr), (const T*)NULL );
  TEST_EQUALITY_CONST( getRawPtr(ptr),  (T*)NULL );
  cptr = new T[n];
  ptr  = new T[n];
  TEST_EQUALITY( getRawPtr(cptr), &cptr[0]);
  TEST_EQUALITY( getRawPtr(ptr),  &ptr[0] );
  delete [] cptr;
  delete [] ptr;
}


#ifdef TEUCHOS_DEBUG


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, arcp_zero, T )
{
  TEST_THROW(ArrayRCP<T> arcp_strong = arcp<T>(0),
    std::out_of_range);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayRCP, arcp_neg, T )
{
  TEST_THROW(ArrayRCP<T> arcp_strong = arcp<T>(-1),
    std::out_of_range);
}


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


#endif // TEUCHOS_DEBUG

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, arcp_zero, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, arcp_neg, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, outOfBounds, T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, assignSelf, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, nullIterator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, implicitConversions, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, weakDelete, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, danglingArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayRCP, getRawPtr, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CPtr, getRawPtr, T) \
  DEBUG_UNIT_TEST_GROUP(T)


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(double)
UNIT_TEST_GROUP(float)


//
// Non templated unit tests
//


TEUCHOS_UNIT_TEST( ArrayRCP, nonnull )
{
  ECHO(ArrayRCP<int> a_arcp = arcp<int>(10));
  TEST_EQUALITY_CONST(is_null(a_arcp), false);
  TEST_EQUALITY_CONST(nonnull(a_arcp), true);
  ECHO(a_arcp = null);
  TEST_EQUALITY_CONST(is_null(a_arcp), true);
  TEST_EQUALITY_CONST(nonnull(a_arcp), false);
}


TEUCHOS_UNIT_TEST( ArrayRCP, weak_strong )
{

  ECHO(ArrayRCP<int> arcp1 = arcp<int>(10));
  TEST_EQUALITY_CONST( arcp1.strength(), RCP_STRONG );

  ECHO(ArrayRCP<int> arcp2 = arcp1.create_weak());

  TEST_EQUALITY_CONST( arcp2.strength(), RCP_WEAK );
  TEST_EQUALITY_CONST( arcp1.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.strong_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.weak_count(), 1 );

  ECHO(ArrayRCP<int> arcp3 = arcp2.create_strong());

  TEST_EQUALITY_CONST( arcp3.strength(), RCP_STRONG );
  TEST_EQUALITY_CONST( arcp1.strong_count(), 2 );
  TEST_EQUALITY_CONST( arcp1.weak_count(), 1 );
  TEST_EQUALITY_CONST( arcp2.strong_count(), 2 );
  TEST_EQUALITY_CONST( arcp2.weak_count(), 1 );

  // This will make the underlying object A gets deleted!
  ECHO(arcp1 = null);
  ECHO(arcp3 = null);

  ECHO(arcp2 = null); // Should make the underlying node go away

}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_null )
{
  ECHO(ArrayRCP<char> arcp_char = null);
  ECHO(ArrayRCP<int> arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY_CONST(arcp_int, null);
}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_char_to_int )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<char> arcp_char = arcp<char>(num_chars));
  ECHO(ArrayRCP<int> arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY(arcp_int.size(), num_ints);
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[0]),
    implicit_ptr_cast<void>(&arcp_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));

  ECHO(arcp_char+=sizeOfInt);
  ECHO(arcp_int = arcp_reinterpret_cast<int>(arcp_char));
  TEST_EQUALITY(arcp_int.size(), num_ints);
  TEST_EQUALITY_CONST( arcp_int.lowerOffset(), -1);
  TEST_EQUALITY( arcp_int.upperOffset(), num_ints-2);
  TEST_EQUALITY( implicit_ptr_cast<void>(&arcp_int[-1]),
    implicit_ptr_cast<void>(&arcp_char[-sizeOfInt])
    );
  TEST_EQUALITY( implicit_ptr_cast<void>((&arcp_int[num_ints-2])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1-sizeOfInt])+1));

}


TEUCHOS_UNIT_TEST( ArrayRCP, arcp_reinterpret_cast_int_to_char )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ECHO(ArrayRCP<int> arcp_int = arcp<int>(num_ints));
  ECHO(ArrayRCP<char> arcp_char = arcp_reinterpret_cast<char>(arcp_int));
  TEST_EQUALITY(arcp_char.size(), num_chars);
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[0]),
    implicit_ptr_cast<void>(&arcp_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));

  ECHO(++arcp_int);
  ECHO(arcp_char = arcp_reinterpret_cast<char>(arcp_int));
  TEST_EQUALITY(as<int>(arcp_char.lowerOffset()), as<int>(-sizeOfInt));
  TEST_EQUALITY(as<int>(arcp_char.upperOffset()), as<int>(num_chars-1-sizeOfInt));
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[-1]),
    implicit_ptr_cast<void>(&arcp_char[-sizeOfInt]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-2])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1-sizeOfInt])+1));

}


TEUCHOS_UNIT_TEST( ArrayRCP, evil_reinterpret_cast )
{
  ECHO(ArrayRCP<ArrayRCP<int> > arcp1 = arcp<ArrayRCP<int> >(n));
  ECHO(ArrayRCP<ArrayRCP<const int> > arcp2 =
    arcp_reinterpret_cast<ArrayRCP<const int> >(arcp1));
  TEST_EQUALITY(arcp2.size(), arcp1.size());
  TEST_EQUALITY(implicit_ptr_cast<const void>(&arcp1[0]),
    implicit_ptr_cast<const void>(&arcp2[0]));
  ECHO(ArrayRCP<const ArrayRCP<const int> > arcp3 = arcp2);
  TEST_EQUALITY(arcp3.size(), arcp1.size());
  TEST_EQUALITY(implicit_ptr_cast<const void>(&arcp1[0]),
    implicit_ptr_cast<const void>(&arcp3[0]));
  out << "arcp3 = " << arcp3 << "\n";
}


} // namespace
