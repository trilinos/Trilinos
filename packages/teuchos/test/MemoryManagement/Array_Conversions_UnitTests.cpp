
#include "Teuchos_ArrayConversions.hpp"
#include "Array_UnitTest_helpers.hpp"
#include "Array_Conversions_UnitTest_helpers.hpp"

namespace {
using ArrayUnitTestHelpers::n;
using ArrayConversionsUnitTestHelpers::generateArray;
using ArrayConversionsUnitTestHelpers::TestArrayViewInput;
using ArrayConversionsUnitTestHelpers::TestArrayViewOutput;
using Teuchos::arrayPtrConv;
using Teuchos::Array;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::as;

// Verify generateArray works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, generateArray, T )
{
  Array<RCP<T> > a_in = generateArray<T>(n);
  TEST_EQUALITY_CONST( as<Teuchos_Ordinal>(a_in.size()), n );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_in[i]), false );
    TEST_EQUALITY_CONST( *a_in[i], as<T>(i) );
  }
}

// Verify TestArrayViewInput works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, TestArrayViewInput, T )  
{
  typedef Teuchos::ScalarTraits<T> ST;
  Array<RCP<T> > a_data = generateArray<T>(n);
  Array<Ptr<const T> > a_in(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    a_in[i] = a_data[i].ptr();
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_in[i]), false );
  }
  T a_out = TestArrayViewInput<T>(a_in);
  TEST_EQUALITY_CONST( a_out, as<T>(n*(n-1)/2) );
}

// Verify TestArrayViewOutput works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, TestArrayViewOutput, T )  
{
  typedef Teuchos::ScalarTraits<T> ST;
  Array<RCP<T> > a_data = generateArray<T>(n);
  Array<Ptr<T> > a_out;
  a_out.reserve(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    *a_data[i] = ST::zero();
    a_out.push_back( a_data[i].ptr() );
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_out[i]), false );
  }
  TestArrayViewOutput<T>(a_out);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], as<T>(i) );
  }

}

// Verify arrayPtrConv works correctly on const objects
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConvConst, T ) 
{
  Array<RCP<T> > a_in = generateArray<T>(n);
  Array<Ptr<const T> > a_out = arrayPtrConv<const T>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( a_out[i].get(), a_in[i].get() );
    TEST_EQUALITY( *a_out[i], *a_in[i] );
  }
}

// Verify arrayPtrConv works correctly on non-const objects
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConvNonConst, T ) 
{
  Array<RCP<T> > a_in = generateArray<T>(n);
  Array<Ptr<T> > a_out = arrayPtrConv<T>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( a_out[i].get(), a_in[i].get() );
    TEST_EQUALITY( *a_out[i], *a_in[i] );
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, passConstObject, T ) 
{
  Array<RCP<T> > a_in = generateArray<T>(n);
  T a = TestArrayViewInput<T>(arrayPtrConv<const T>(a_in));
  T a_exact = as<T>(n*(n-1)/2);
  TEST_EQUALITY( a, a_exact );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, passNonConstObject, T ) 
{
  typedef Teuchos::ScalarTraits<T> ST;
  Array<RCP<T> > a_out = generateArray<T>(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    *a_out[i] = ST::zero();
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], ST::zero() );
  }
  TestArrayViewOutput<T>(arrayPtrConv<T>(a_out));
  TEST_EQUALITY_CONST( as<Teuchos_Ordinal>(a_out.size()), n );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], as<T>(i) );
  }
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, generateArray, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, TestArrayViewInput, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, TestArrayViewOutput, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConvConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConvNonConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, passConstObject, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, passNonConstObject, T ) \
  DEBUG_UNIT_TEST_GROUP( T )


UNIT_TEST_GROUP(Teuchos_Ordinal)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace

