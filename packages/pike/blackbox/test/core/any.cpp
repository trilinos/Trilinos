#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_Any.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
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
    intVal = a.as<int>();
    TEST_EQUALITY(intVal,i);
    Teuchos::ArrayView<const double> check2 = 
      v[3].as<Teuchos::ArrayView<const double> >();
    TEST_EQUALITY(check2.size(),tv.size());

    Teuchos::ArrayView<const double>& check3 = 
      v[3].as<Teuchos::ArrayView<const double>& >();
    TEST_EQUALITY(check3.size(),tv.size());
  }

  TEUCHOS_UNIT_TEST(any, const_rcp)
  {
    std::vector<Teuchos::RCP<any> > v(3);
    for (std::vector<Teuchos::RCP<any> >::iterator i=v.begin(); i != v.end(); ++i)
      (*i) = Teuchos::rcp(new any);

    *v[0] = 1.0;
    *v[1] = 5.0;
    Teuchos::Array<double> a(5,2.0);
    *v[2] = a.view(0,5);

    Teuchos::RCP<const any> v0 = v[0];
    Teuchos::RCP<const any> v1 = v[1];
    Teuchos::RCP<const any> v2 = v[2];

    const double r0 = v0->as<double>();
    const double r1 = v1->as<double>();
    const Teuchos::ArrayView<double> r2 = v2->as<Teuchos::ArrayView<double> >();
    
    const double tol = Teuchos::ScalarTraits<double>::eps();
    TEST_FLOATING_EQUALITY(r0,1.0,tol);
    TEST_FLOATING_EQUALITY(r1,5.0,tol);
    TEST_FLOATING_EQUALITY(r2[0],2.0,tol);
  }

}
