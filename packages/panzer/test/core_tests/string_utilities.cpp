#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_String_Utilities.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(object_builders, StringTokenizerMultiple)
  {
    std::string names =  "fluid,solid";
    std::vector<std::string> tokens;
    panzer::StringTokenizer(tokens, names);
    
    TEST_EQUALITY(tokens.size(), 2);
    TEST_EQUALITY(tokens[0], "fluid");
    TEST_EQUALITY(tokens[1], "solid");
  }

  TEUCHOS_UNIT_TEST(object_builders, StringTokenizerSingle)
  {
    std::string names =  "fluid";
    std::vector<std::string> tokens;
    panzer::StringTokenizer(tokens, names);
    
    TEST_EQUALITY(tokens.size(), 1);
    TEST_EQUALITY(tokens[0], "fluid");
  }

}
