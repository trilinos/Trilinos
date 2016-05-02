#include <Teuchos_UnitTestHarness.hpp>
#include "Tempus_Version.hpp"

namespace Tempus_Test {

  TEUCHOS_UNIT_TEST(version, default)
  {
    Tempus::version();
  }

}
