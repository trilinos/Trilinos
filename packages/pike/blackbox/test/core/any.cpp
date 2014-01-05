#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_Any.hpp"
#include "Teuchos_Array.hpp"
#include <iostream>
#include <vector>

using namespace std;
//using namespace Teuchos;
using namespace pike;

namespace pike_test {
  
  TEUCHOS_UNIT_TEST(any, basic)
  {
    int i = 5;
    
    any a;
    a = i;
    any b = 5;
    any c = i;

    cout << "\nMy Type is: " << a.type().name() << endl;    

    int intVal = any_cast<int>(a);
    TEST_EQUALITY(intVal,i);

    vector<any> v;
    
    Teuchos::Array<double> t(5);
    Teuchos::ArrayView<const double> tv = t.view(0,t.size());
    v.push_back(a);
    v.push_back(i);
    v.push_back(t);
    v.push_back(tv);

    Teuchos::ArrayView<const double>& check = 
      any_cast<Teuchos::ArrayView<const double>& >(v[3]);

    TEST_EQUALITY(check.size(),tv.size());

    // Test additions to boost any
    
  }

}
