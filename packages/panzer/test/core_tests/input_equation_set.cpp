#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_InputEquationSet.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(input_equation_set, basic)
  {
    panzer::InputEquationSet ies_2;
    ies_2.name = "Continuity";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = 6;
    ies_2.model_factory = "rf";
    ies_2.prefix = "ION_";
    ies_2.params.set("Tau_C", "Codina");

    Teuchos::ParameterList my_list;
    my_list.set("Name", "Continuity");
    my_list.set("Basis", "Q1");
    my_list.set("Integration Order", 1);
    my_list.set("Model ID", 6);
    my_list.set("Model Factory", "rf");
    my_list.set("Prefix", "ION_");
    my_list.sublist("Options").set("Tau_C", "Codina");
    
    panzer::InputEquationSet ies_3(my_list);
    
    TEST_EQUALITY(ies_2.name, ies_3.name);
    TEST_EQUALITY(ies_2.basis, ies_3.basis);
    TEST_EQUALITY(ies_2.integration_order, ies_3.integration_order);
    TEST_EQUALITY(ies_2.model_id, ies_3.model_id);
    TEST_EQUALITY(ies_2.model_factory, ies_3.model_factory);
    TEST_EQUALITY(ies_2.prefix, ies_3.prefix);
    TEST_EQUALITY(ies_2.params, ies_3.params);
  }
}
