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

#include "Shards_BasicTopologies.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/CreateAdjacentEntities.hpp"

int main(int argc, char *argv[])
{
  using namespace Teuchos;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  Kokkos::initialize(argc, argv);

  std::string input_file_name = "periodic_cylinder.pam";
  std::string output_file_name = "periodic_cylinder.gen";
  
  Teuchos::CommandLineProcessor clp;
  clp.setOption("i", &input_file_name, "pamgen input file");
  clp.setOption("o", &output_file_name, "exodus output file");
  
  auto parse_return = clp.parse(argc,argv,&std::cerr);

  TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL,
                             std::runtime_error,
                             "Failed to parse command line. Use --help to find correct options.");

  RCP<ParameterList> p = parameterList();
  p->set("File Name",input_file_name);
  p->set("File Type","Pamgen");
  
  panzer_stk::STK_ExodusReaderFactory pamgenFactory;
  pamgenFactory.setParameterList(p);
  
  auto mesh = pamgenFactory.buildMesh(MPI_COMM_WORLD);
  auto metaData = mesh->getMetaData();
  auto bulkData = mesh->getBulkData();

  mesh->writeToExodus(output_file_name);
  
  return 0;
}
