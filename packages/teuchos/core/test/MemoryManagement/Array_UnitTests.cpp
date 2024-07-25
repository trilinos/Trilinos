// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
using Teuchos::arcpFromArray;
using Teuchos::as;
using Teuchos::getConst;
using Teuchos::DanglingReferenceError;
using Teuchos::fromStringToArray;
using Teuchos::InvalidArrayStringRepresentation;


TEUCHOS_UNIT_TEST( Utils, trimWhiteSpace_empty )
{
  TEST_EQUALITY(Teuchos::Utils::trimWhiteSpace(""), "");
}



TEUCHOS_UNIT_TEST( Array, TypeNameTraits )
{
  TEST_EQUALITY(Teuchos::TypeNameTraits<Array<double> >::name(),
    std::string("Array(double)"));
}


TEUCHOS_UNIT_TEST( Array, stringToArray )
{

  {
    std::string arrayString="{}";
    std::istringstream arrayStream(arrayString);
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayStreamVal;
    arrayStream >> arrayStreamVal;
    Array<std::string> arrayVal_exp;
    TEST_EQUALITY(arrayVal, arrayVal_exp);
    TEST_EQUALITY(arrayStreamVal, arrayVal_exp);
  }

  {
    std::string arrayString = "{ a, b, c, d }";
    std::istringstream arrayStream(arrayString);
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayStreamVal;
    arrayStream >> arrayStreamVal;
    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("a", "b", "c", "d" );
    TEST_EQUALITY(arrayVal, arrayVal_exp);
    TEST_EQUALITY(arrayStreamVal, arrayVal_exp);
  }

  {
    std::string arrayString = "{ (a), b, c, (d) }";
    std::istringstream arrayStream(arrayString);
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayStreamVal;
    arrayStream >> arrayStreamVal;
    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("(a)", "b", "c", "(d)" );
    TEST_EQUALITY(arrayVal, arrayVal_exp);
    TEST_EQUALITY(arrayStreamVal, arrayVal_exp);
  }

  {
    std::string arrayString = "{ (a ), b, c, (d ) }";
    std::istringstream arrayStream(arrayString);
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayStreamVal;
    arrayStream >> arrayStreamVal;
    Array<std::string> arrayVal_exp = Teuchos::tuple<std::string>("(a )", "b", "c", "(d )" );
    TEST_EQUALITY(arrayVal, arrayVal_exp);
    TEST_EQUALITY(arrayStreamVal, arrayVal_exp);
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


TEUCHOS_UNIT_TEST( Array, stringToArray_string_hyphens )
{

  {
    std::string arrayString="{-}";
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayVal_exp = tuple<std::string>("-");
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    std::string arrayString="{-,-}";
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayVal_exp = tuple<std::string>("-","-");
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    std::string arrayString="{-,1,-}";
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayVal_exp = tuple<std::string>("-","1","-");
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    std::string arrayString="{}";
    Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString);
    Array<std::string> arrayVal_exp;
    TEST_EQUALITY(arrayVal, arrayVal_exp);
  }

  {
    std::string arrayString="{,}";
    TEST_THROW(
      [[maybe_unused]] Array<std::string> arrayVal = fromStringToArray<std::string>(arrayString),
      Teuchos::InvalidArrayStringRepresentation);
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, defaultConstruct, T )
{
  Array<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( as<int>(a2.empty()), true );
  TEST_EQUALITY_CONST( a2.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( a2.getRawPtr(), a2.data() );
  TEST_EQUALITY_CONST( getConst(a2).getRawPtr(), 0 );
  TEST_EQUALITY_CONST( getConst(a2).getRawPtr(), getConst(a2).data() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, sizedConstruct, T )
{
  typedef typename Array<T>::size_type size_type;
  Array<T> a(n);
  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( a.length(), n );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_EQUALITY( a.getRawPtr(), &a[0] );
  TEST_EQUALITY( a.getRawPtr(), a.data() );
  TEST_EQUALITY( getConst(a).getRawPtr(), &getConst(a)[0] );
  TEST_EQUALITY( getConst(a).getRawPtr(), getConst(a).data() );
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPArray_to_ArrayRCP_null, T )
{
  const RCP<Array<T> > a_rcp = null;
  const ArrayRCP<T> a_arcp = arcp(a_rcp);
  TEST_ASSERT( a_arcp == null );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPconstArray_to_ArrayRCP_null, T )
{
  const RCP<const Array<T> > a_rcp = null;
  const ArrayRCP<const T> a_arcp = arcp(a_rcp);
  TEST_ASSERT( a_arcp == null );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPArray_to_ArrayRCP, T )
{
  const Array<T> a_const = generateArray<T>(n);
  const RCP<Array<T> > a_rcp = Teuchos::rcp( new  Array<T>(a_const));
  const ArrayRCP<T> a_arcp = arcp(a_rcp);
  TEST_COMPARE_ARRAYS( a_const(), a_arcp() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, RCPconstArray_to_ArrayRCP, T )
{
  const Array<T> a_const = generateArray<T>(n);
  const RCP<const Array<T> > a_rcp = Teuchos::rcp( new  Array<T>(a_const));
  const ArrayRCP<const T> a_arcp = arcp(a_rcp);
  TEST_COMPARE_ARRAYS( a_const(), a_arcp() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, Array_to_ArrayRCP_null, T )
{
  Array<T> a;
  const ArrayRCP<T> a_arcp = arcpFromArray(a);
  TEST_ASSERT(a_arcp == null);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, constArray_to_ArrayRCP_null, T )
{
  const Array<T> a;
  const ArrayRCP<const T> a_arcp = arcpFromArray(a);
  TEST_ASSERT(a_arcp == null);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, Array_to_ArrayRCP, T )
{
  Array<T> a = generateArray<T>(n);
  const ArrayRCP<T> a_arcp = arcpFromArray(a);
  TEST_COMPARE_ARRAYS( a(), a_arcp() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, constArray_to_ArrayRCP, T )
{
  const Array<T> a = generateArray<T>(n);
  const ArrayRCP<const T> a_arcp = arcpFromArray(a);
  TEST_COMPARE_ARRAYS( a(), a_arcp() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, Array_to_ArrayRCP_dangling, T )
{
  ArrayRCP<T> a_arcp;
  {
    Array<T> a = generateArray<T>(n);
    a_arcp = arcpFromArray(a);
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW(a_arcp[0], DanglingReferenceError);
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, constArray_to_ArrayRCP_dangling, T )
{
  ArrayRCP<const T> a_arcp;
  {
    const Array<T> a = generateArray<T>(n);
    a_arcp = arcpFromArray(a);
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW(a_arcp[0], DanglingReferenceError);
#else
  (void)out; (void)success;
#endif
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, view_empty_func, T )
{
  Array<T> a;
  const ArrayView<T> av = a.view(0, 0);
  TEST_ASSERT(is_null(av));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, view_empty_operator, T )
{
  Array<T> a;
  const ArrayView<T> av = a(0, 0);
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_implicit, T )
{
  ArrayView<T> av;
  { Array<T> a(n); av = a; }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0] = 0, DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_implicit_const, T )
{
  ArrayView<const T> av;
  { Array<T> a(n); av = getConst(a); }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0], DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_explicit, T )
{
  ArrayView<T> av;
  { Array<T> a(n); av = a(); }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0] = 0, DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_explicit_const, T )
{
  ArrayView<const T> av;
  { Array<T> a(n); av = getConst(a)(); }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0], DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_subview, T )
{
  ArrayView<T> av;
  { Array<T> a(n); av = a(0,1); }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0] = 0, DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayView_subview_const, T )
{
  ArrayView<const T> av;
  { Array<T> a(n); av = getConst(a)(0,1); }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0], DanglingReferenceError );
#else
  (void)out; (void)success;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayViewIter, T )
{
  typedef typename ArrayView<T>::iterator iter_t;
  ECHO(Array<T> a(n));
  ECHO(ArrayView<T> av = a);
  ECHO(iter_t iter = av.begin());
  ECHO(av = null);
  ECHO(a.resize(0));
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( *iter = 0, DanglingReferenceError );
#else
  (void)iter; (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, danglingArrayViewIter_const, T )
{
  typedef typename ArrayView<const T>::iterator iter_t;
  ECHO(Array<T> a(n));
  ECHO(ArrayView<T> av = a);
  ECHO(iter_t iter = av.begin());
  ECHO(av = null);
  ECHO(a.resize(0));
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( *iter, DanglingReferenceError );
#else
  (void)iter; (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, erase_empty, T )
{
  ECHO(std::vector<T> v);
  TEST_NOTHROW(v.erase(v.begin(), v.end()));
  ECHO(Array<T> a);
  TEST_NOTHROW(a.erase(a.begin(), a.end()));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, structuralChangeArrayView, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<T> av = a;
  a.push_back(a[0]);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0] = 0, DanglingReferenceError );
#else
  (void)av; (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, structuralChangeArrayView_const, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<const T> av = getConst(a);
  a.push_back(a[0]);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[0], DanglingReferenceError );
#else
  (void)av; (void)out; (void)success;
#endif
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, initializerList, T )
{
  const Array<T> a2{as<T>(1),as<T>(2),as<T>(3),as<T>(4)};
  TEST_EQUALITY( a2.size(), 4);
  for (typename Array<T>::size_type i=0; i < a2.size(); ++i) {
    TEST_EQUALITY( a2[i], as<T>(i+1) );
  }
}


//
// Instantiations
//

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, sizedConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, operatorBracket, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constAt, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter_before_block_end, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPArray_to_ArrayRCP_null, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPconstArray_to_ArrayRCP_null, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPArray_to_ArrayRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, RCPconstArray_to_ArrayRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, Array_to_ArrayRCP_null, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constArray_to_ArrayRCP_null, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, Array_to_ArrayRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constArray_to_ArrayRCP, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, Array_to_ArrayRCP_dangling, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constArray_to_ArrayRCP_dangling, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, toVector, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, toVector_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, view_empty_func, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, view_empty_operator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, implicit_to_ArrayView_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_implicit, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_implicit_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_explicit, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_explicit_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_subview, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayView_subview_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, danglingArrayViewIter_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, erase_empty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, structuralChangeArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, structuralChangeArrayView_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, initializerList, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
