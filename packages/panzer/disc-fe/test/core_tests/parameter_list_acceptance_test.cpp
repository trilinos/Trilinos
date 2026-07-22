// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    TEST_THROW(updateParametersFromXmlFile("parameter_list_acceptance_test2.xml", ptrFromRef(p)),std::logic_error);

    out << p << std::endl;
  }
}
