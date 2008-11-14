#include "Array_UnitTest_helpers.hpp"


namespace {


using ArrayUnitTestHelpers::n;
using ArrayUnitTestHelpers::generateArray;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::ArrayView;
using Teuchos::DanglingReferenceError;
using Teuchos::as;
using Teuchos::null;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, assignmentOperator, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<T> av1;
  av1 = a;
  ArrayView<T> av2;
  av2 = av1;
  TEST_EQUALITY( av1.getRawPtr(), a.getRawPtr() );
  TEST_EQUALITY( av1.size(), as<int>(a.size()) );
  TEST_EQUALITY( av1.getRawPtr(), av2.getRawPtr() );
  TEST_EQUALITY( av1.size(), av2.size() );
  TEST_COMPARE_ARRAYS( av1, a );
  TEST_COMPARE_ARRAYS( av1, av2 );
  av1 = null;
  TEST_EQUALITY_CONST( av1.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( av1.size(), 0 );
  av2 = null;
  TEST_EQUALITY_CONST( av2.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( av2.size(), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, iterators, T )
{
  typedef typename ArrayView<T>::iterator iter_t;
  typedef Teuchos::ScalarTraits<T> ST;
  ECHO(Array<T> a = generateArray<T>(n));
  ECHO(ArrayView<T> av = a);
  ECHO(const iter_t av_begin = av.begin());
  ECHO(const iter_t av_end = av.end());
#ifdef TEUCHOS_DEBUG
  TEST_ASSERT(av_begin.shares_resource(av_end));
#endif
  ECHO(std::fill(av_begin, av_end, ST::random()));
  ECHO(Array<T> a2 = generateArray<T>(n));
  ECHO(ArrayView<T> av2 = a2);
  ECHO(std::copy(av.begin(), av.end(), av2.begin()));
  TEST_COMPARE_ARRAYS(a, a2);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, danglingView_std_vector, T )
{
  ArrayView<T> av;
  T* badPtr = 0;
  {
    std::vector<T> v(n);
    av = v;
    badPtr = &v[0];
  }
  // Access the raw pointer but it now points to invalid memory!
  TEST_EQUALITY(av.getRawPtr(), badPtr);
  // Above, we have no way to detect that the underlying std::vector object
  // has gone away.  This is the whole point of needing Teuchos::Array and
  // having an integrated set of utility classes that all work together!
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, danglingView_rcp_std_vector, T )
{
  ArrayView<T> av;
  {
    ArrayRCP<T> ap = arcp(rcp(new std::vector<T>(n)));
    av = ap;
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW(av.getRawPtr(), DanglingReferenceError);
#endif
  // Above, because we wrapped the initial std::vector in an RCP object, we
  // can sucessfully detect when the object goes away in debug mode!
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


//
// Instantiations
//



#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, assignmentOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, iterators, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, danglingView_std_vector, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, danglingView_rcp_std_vector, T ) \
  DEBUG_UNIT_TEST_GROUP( T )


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
