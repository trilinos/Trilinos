#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "Panzer_ZeroSensitivities.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(zero_sensitivites, default)
  {

    double x = 1.0;
    TEST_FLOATING_EQUALITY(x,1.0,1e-12);
    panzer::zeroSensitivities(x);
    TEST_FLOATING_EQUALITY(x,1.0,1e-12);

    Sacado::Fad::DFad<double> y;
    y.val() = 1.0;
    y.resize(2);
    y.fastAccessDx(0) = 2.0;
    y.fastAccessDx(1) = 3.0;
    TEST_FLOATING_EQUALITY(y.val(),1.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(0),2.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(1),3.0,1e-12);
    panzer::zeroSensitivities(y);
    TEST_FLOATING_EQUALITY(y.val(),1.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(0),0.0,1e-12);
    TEST_FLOATING_EQUALITY(y.fastAccessDx(1),0.0,1e-12);
  }

}
