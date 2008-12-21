#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace ArrayUnitTestHelpers {


extern int n;


template<class T>
Teuchos::Array<T> generateArray(const int n_in)
{
  using Teuchos::as;
  Teuchos::Array<T> a(n_in);
  for( int i = 0; i < n_in; ++i )
    a[i] = as<T>(i); // tests non-const operator[](i)
  return a;
}


} // namespace ArrayUnitTestHelpers
