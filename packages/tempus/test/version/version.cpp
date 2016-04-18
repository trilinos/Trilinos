#include <Teuchos_UnitTestHarness.hpp>
#include "Tempus_Version.hpp"

namespace tempus_test {

  TEUCHOS_UNIT_TEST(version, default)
  {
    tempus::version();
  }

}
