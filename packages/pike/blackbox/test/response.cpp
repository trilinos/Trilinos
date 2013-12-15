#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_Response_Scalar.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <cmath>

namespace pike {

  TEUCHOS_UNIT_TEST(response, scalar)
  {
    Teuchos::RCP<pike::ScalarResponse<double> > r = 
      pike::scalarResponse<double>("q");

    double value = 2.0;
    r->set(value);
    double tol = 100.0 * Teuchos::ScalarTraits<double>::eps();
    TEST_ASSERT(std::fabs(value - r->get()) < tol);

    Teuchos::RCP<pike::Response> r_base = r;

    double test_value;
    pike::getScalarResponse(test_value,*r_base);
    TEST_FLOATING_EQUALITY(value,test_value, tol);

    test_value = pike::getScalarResponse<double>(*r_base);
    TEST_FLOATING_EQUALITY(value,test_value, tol);
  }

}
