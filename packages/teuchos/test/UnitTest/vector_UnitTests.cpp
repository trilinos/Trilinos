#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace {


int n = 4;


class UnitTestSetup {
public:
  UnitTestSetup()
    {
      Teuchos::UnitTestRepository::getCLP().setOption(
        "n", &n, "Number of elements in the vectors" );
    }
} unitTestSetup;


template<class T>
std::vector<T> generatevector(const int n)
{
  std::vector<T> a(n);
  for( int i = 0; i < n; ++i )
    a[i] = i; // tests non-const operator[](i)
  return a;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, defaultConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  vector<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( as<int>(a2.empty()), true );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, defaultConstruct, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, defaultConstruct, float )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, defaultConstruct, double )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, sizedConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  using Teuchos::getConst;
  typedef typename std::vector<T>::size_type size_type;
  vector<T> a(n);
  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, sizedConstruct, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, sizedConstruct, float )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, sizedConstruct, double )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, operatorBracket, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a[i] == i ... ";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
  }
  if (local_success) out << "passed\n";
  else success = false;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, operatorBracket, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, operatorBracket, float )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, operatorBracket, double )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, constAt, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a.at(i) == i ...\n";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
  }
  if (local_success) out << "passed\n";
  else success = false;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, constAt, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, constAt, float )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, constAt, double )


} // namespace
