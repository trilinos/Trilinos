#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <iostream>

namespace panzer {

  TEUCHOS_UNIT_TEST(parameter_list_acceptance_test,order_preserving)
  {
    using namespace Teuchos;

    ParameterList p;
    updateParametersFromXmlFile("parameter_list_acceptance_test.xml", ptrFromRef(p));

    TEST_ASSERT(p.isSublist("Boundary Conditions"));

    ParameterList bcs = p.sublist("Boundary Conditions");

    typedef ParameterList::ConstIterator pl_it;

    int index = 0;
    for (pl_it bc=bcs.begin(); bc != bcs.end(); ++bc,++index) {
      ParameterList* sublist=0;
      TEST_EQUALITY(bc->second.getValue(sublist).get<int>("order"),index);

    }

    out << p << std::endl;
  }

  TEUCHOS_UNIT_TEST(parameter_list_acceptance_test,repeated_sublist_name)
  {
    using namespace Teuchos;

    ParameterList p;
    updateParametersFromXmlFile("parameter_list_acceptance_test2.xml", ptrFromRef(p));

    out << p << std::endl;
  }
}
