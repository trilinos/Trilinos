#include "Array_UnitTest_helpers.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_ArrayRCP.hpp"


namespace {


using ArrayUnitTestHelpers::n;
using ArrayUnitTestHelpers::generateArray;
using Teuchos::null;
using Teuchos::tuple;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::as;
using Teuchos::getConst;
using Teuchos::DanglingReferenceError;
using Teuchos::fromStringToArray;
using Teuchos::InvalidArrayStringRepresentation;


TEUCHOS_UNIT_TEST( Array, stringToArray )
{

  {
    Array<std::string> arrayVal = fromStringToArray<std::string>("{}");
    Array<std::string> arrayVal_exp;
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    Array<std::string> arrayVal = fromStringToArray<std::string>("{ a, b, c, d }");
    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("a", "b", "c", "d" );
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    Array<std::string> arrayVal = fromStringToArray<std::string>("{ (a), b, c, (d) }");
    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("(a)", "b", "c", "(d)" );
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  // This should work but does not.  I should fix this!
//  {
//    Array<std::string> arrayVal = fromStringToArray<std::string>("{ {a}, 'b', {c }, d }");
//    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("{a}", "'b'", "{c }", "d" );
//    TEST_EQUALITY(arrayVal, arrayVal_exp);
//  }

}


TEUCHOS_UNIT_TEST( Array, stringToArray_invalid )
{
  TEST_THROW(fromStringToArray<std::string>("{ a, b, c"),
    InvalidArrayStringRepresentation);
  TEST_THROW(fromStringToArray<std::string>("a, b, c}"),
    InvalidArrayStringRepresentation);
  TEST_THROW(fromStringToArray<std::string>("a, b, c"),
    InvalidArrayStringRepresentation);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, defaultConstruct, T )
{
  Array<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( as<int>(a2.empty()), true );
  TEST_EQUALITY_CONST( a2.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( getConst(a2).getRawPtr(), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, sizedConstruct, T )
{
  typedef typename Array<T>::size_type size_type;
  Array<T> a(n);
  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( a.length(), n );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_EQUALITY( a.getRawPtr(), &a[0] );
  TEST_EQUALITY( getConst(a).getRawPtr(), &getConst(a)[0] );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, operatorBracket, T )
{
  out << "\nTest that a[i] == i ... ";
  Array<T> a = generateArray<T>(n);
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, constAt, T )
{
  out << "\nTest that a.at(i) == i ...\n";
  Array<T> a = generateArray<T>(n);
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayViewIter_before_block_end, T )
{
  typedef typename ArrayView<T>::iterator iter_t;
  typedef Teuchos::NullIteratorTraits<iter_t> NIT;
  iter_t iter = NIT::getNull();
  {
    ECHO(Array<T> a(n, as<T>(0)));
    ECHO(ArrayView<T> av = a);
    ECHO(iter = av.begin());
    ECHO(av = null);
    TEST_EQUALITY( *iter, a[0] );
    // Above, the iterator to the ArrayView object is still valid even through
    // the ArrayView object is was created from is gone now.  This is just
    // fine since the underlying data is still there in the original Array object.
    iter = NIT::getNull();
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPArray_to_ArrayRCP, T )
{
  const Array<T> a_const = generateArray<T>(n);
  const RCP<Array<T> > a_rcp = rcp(new Array<T>(a_const));
  const ArrayRCP<T> a_arcp = arcp(a_rcp);
  TEST_COMPARE_ARRAYS( a_const(), a_arcp() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPArray_to_ArrayRCP_null, T )
{
  const RCP<Array<T> > a_rcp = null;
  const ArrayRCP<T> a_arcp = arcp(a_rcp);
  TEST_ASSERT( a_arcp == null );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, toVector, T )
{
  const Array<T> a = generateArray<T>(n);
  const std::vector<T> v = a.toVector();
  TEST_COMPARE_ARRAYS( a, v );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, toVector_empty, T )
{
  const Array<T> a;
  const std::vector<T> v = a.toVector();
  TEST_EQUALITY_CONST( as<int>(v.size()), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, view_empty, T )
{
  Array<T> a;
  const ArrayView<T> av = a.view(0, 0);
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, view_const_empty, T )
{
  const Array<T> a;
  const ArrayView<const T> av = a.view(0, 0);
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, implicit_to_ArrayView_empty, T )
{
  Array<T> a;
  const ArrayView<T> av = a();
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, implicit_to_ArrayView_const_empty, T )
{
  const Array<T> a;
  const ArrayView<const T> av = a();
  TEST_ASSERT(is_null(av));
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_implicit, T )
{
  ArrayView<T> av;
  TEST_THROW( { Array<T> a(n); av = a; },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_implicit_const, T )
{
  ArrayView<const T> av;
  TEST_THROW( { Array<T> a(n); av = getConst(a); },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_explicit, T )
{
  ArrayView<T> av;
  TEST_THROW( { Array<T> a(n); av = a(); },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_explicit_const, T )
{
  ArrayView<const T> av;
  TEST_THROW( { Array<T> a(n); av = getConst(a)(); },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_subview, T )
{
  ArrayView<T> av;
  TEST_THROW( { Array<T> a(n); av = a(0,1); },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_subview_const, T )
{
  ArrayView<const T> av;
  TEST_THROW( { Array<T> a(n); av = getConst(a)(0,1); },
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayViewIter, T )
{
  typedef typename ArrayView<T>::iterator iter_t;
  ECHO(Array<T> a(n));
  ECHO(ArrayView<T> av = a);
  ECHO(iter_t iter = av.begin());
  ECHO(av = null);
  TEST_THROW( a.resize(0), DanglingReferenceError );
  // The way that Array::resize() is able to detect that there is still a
  // dangling iterator comes from the way in which all of this is implemented.
  // The reason that this throws is that the RCP<std::vector<T> > object
  // embedded in the Array object gets embedded inside of the dealloc object
  // which is attached to the node.  Therefore, even though the weak
  // ArrayRCP<T> object that was created and embedded in the ArrayView<T>
  // object has gone away, this same weak ArrayRCP<T> object was used to
  // create another weak ArrayRCP<T> object which *is* the iterator object
  // iter above.  What this means is that any reference counted object that
  // gets created for whatever reason based on the underlying
  // RCP<std::vector<T> > object will result in all the Array functions that
  // alter the underlying array memory to throw an exception right away.  I
  // think I really like this "early warning" behavior.  The only disadvantage
  // is that we don't get a very good error message.  It would be nice to find
  // a way so that it was the dangling reference object itself that threw the
  // exception message and was able to provide better debug feedback.
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayViewIter_const, T )
{
  typedef typename ArrayView<const T>::iterator iter_t;
  ECHO(Array<T> a(n));
  ECHO(ArrayView<T> av = a);
  ECHO(iter_t iter = av.begin());
  ECHO(av = null);
  TEST_THROW( a.resize(0), DanglingReferenceError );
  // See comments above.
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, structuralChangeArrayView, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<T> av = a;
  TEST_THROW( a.push_back(a[0]), 
    DanglingReferenceError );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, structuralChangeArrayView_const, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<const T> av = getConst(a);
  TEST_THROW( a.push_back(a[0]), 
    DanglingReferenceError );
}

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK



//
// Instantiations
//

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_implicit, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_implicit_const, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_explicit, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_explicit_const, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_subview, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_subview_const, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter_const, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, structuralChangeArrayView, T ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, structuralChangeArrayView_const, T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, sizedConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, operatorBracket, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constAt, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter_before_block_end, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPArray_to_ArrayRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPArray_to_ArrayRCP_null, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, toVector, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, toVector_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, view_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, view_const_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, implicit_to_ArrayView_empty, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
