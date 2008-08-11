#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace {


int n = 4;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "n", &n, "Number of elements in the array" );
}


template<class T>
Teuchos::Array<T> generateArray(const int n)
{
  Teuchos::Array<T> a(n);
  for( int i = 0; i < n; ++i )
    a[i] = i; // tests non-const operator[](i)
  return a;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, defaultConstruct, T )
{
  using Teuchos::Array;
  using Teuchos::as;
  Array<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( as<int>(a2.empty()), true );
  TEST_EQUALITY_CONST( a2.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( getConst(a2).getRawPtr(), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, sizedConstruct, T )
{
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::getConst;
  typedef typename Teuchos::Array<T>::size_type size_type;
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
  using Teuchos::Array;
  using Teuchos::as;
  out << "\nTest that a[i] == i ... ";
  Array<T> a = generateArray<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Array, constAt, T )
{
  using Teuchos::Array;
  using Teuchos::as;
  out << "\nTest that a.at(i) == i ...\n";
  Array<T> a = generateArray<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


//
// Instantiations
//


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, sizedConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, operatorBracket, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Array, constAt, T )


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
