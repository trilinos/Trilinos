#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_Version.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(version, default)
  {
    panzer::version();
  }

}
