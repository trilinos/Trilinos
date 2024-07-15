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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#ifdef PANZER_HAVE_PERCEPT
#include "percept/PerceptMesh.hpp"
#endif

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/CreateAdjacentEntities.hpp"

namespace panzer_stk {

// Test to make sure the default option for mesh rebalancing doesn't throw
TEUCHOS_UNIT_TEST(tSTKRebalancing, default)
{
  using namespace Teuchos;
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEST_EQUALITY(Teuchos::DefaultComm<int>::getComm()->getSize(),4);

  {
    out << "\nCreating pamgen mesh" << std::endl;
    RCP<ParameterList> p = parameterList();
    const std::string input_file_name = "pamgen_test.gen";
    p->set("File Name",input_file_name);
    p->set("File Type","Pamgen");
    p->set("Rebalancing","default");

    RCP<STK_ExodusReaderFactory> pamgenFactory = rcp(new STK_ExodusReaderFactory());
    TEST_NOTHROW(pamgenFactory->setParameterList(p));

    RCP<STK_Interface> mesh = pamgenFactory->buildMesh(MPI_COMM_WORLD);
    TEST_ASSERT(mesh!=Teuchos::null);
  }
}

// Test to make sure the none option for mesh rebalancing doesn't throw
TEUCHOS_UNIT_TEST(tSTKRebalancing, none)
{
  using namespace Teuchos;
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEST_EQUALITY(Teuchos::DefaultComm<int>::getComm()->getSize(),4);

  {
    out << "\nCreating pamgen mesh" << std::endl;
    RCP<ParameterList> p = parameterList();
    const std::string input_file_name = "pamgen_test.gen";
    p->set("File Name",input_file_name);
    p->set("File Type","Pamgen");
    p->set("Rebalancing","none");

    RCP<STK_ExodusReaderFactory> pamgenFactory = rcp(new STK_ExodusReaderFactory());
    TEST_NOTHROW(pamgenFactory->setParameterList(p));

    RCP<STK_Interface> mesh = pamgenFactory->buildMesh(MPI_COMM_WORLD);
    TEST_ASSERT(mesh!=Teuchos::null);
  }
}

} // namespace panzer_stk
