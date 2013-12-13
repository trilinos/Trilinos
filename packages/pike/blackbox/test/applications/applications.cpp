#include "Teuchos_UnitTestHarness.hpp"
#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Version.hpp"

namespace cthulhu {

  TEUCHOS_UNIT_TEST(version, default)
  {
    cthulhu::version();
  }

}
