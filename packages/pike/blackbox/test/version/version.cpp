#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_Version.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(version, default)
  {
    pike::version();
  }

}
