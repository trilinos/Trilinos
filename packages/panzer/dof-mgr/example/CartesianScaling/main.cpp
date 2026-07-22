// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include "Kokkos_Core.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "PanzerCore_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include "CartesianConnManager.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

template <typename IntrepidBasisType>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern();

int main(int argc,char * argv[])
{
  using CCM = panzer::unit_test::CartesianConnManager;
  using panzer::DOFManager;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);

  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  int np   = comm.getSize(); // number of processors

  // timings output
  std::string timingsFile = "timings.yaml";

  // mesh description
  int nx = 10, ny = 7, nz = 4;
  int px = np, py = 1, pz = 1;
  int bx =  1, by = 2, bz = 1;

  // parse command line arguments
  Teuchos::CommandLineProcessor clp;
  clp.setOption("nx",&nx);
  clp.setOption("ny",&ny);
  clp.setOption("nz",&nz);
  clp.setOption("px",&px);
  clp.setOption("py",&py);
  clp.setOption("pz",&pz);
  clp.setOption("timings-file",&timingsFile);
  auto cmdResult = clp.parse(argc,argv);
  if(cmdResult!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    clp.printHelpMessage(argv[0],std::cout);
    return -1;
  } 

  // build velocity, temperature and pressure fields
  RCP<const panzer::FieldPattern> pattern_U = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::Device::execution_space,double,double>>();
  RCP<const panzer::FieldPattern> pattern_P = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device::execution_space,double,double>>();
  RCP<const panzer::FieldPattern> pattern_T = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device::execution_space,double,double>>();
  RCP<const panzer::FieldPattern> pattern_B = buildFieldPattern<Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::Device::execution_space,double,double>>();
  RCP<const panzer::FieldPattern> pattern_E = buildFieldPattern<Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::Device::execution_space,double,double>>();

  // repeatedly construct DOFManager timing the buildGlobalUnknowns
  for(int repeats=0;repeats<100;repeats++) {
  
    // build the topology
    RCP<CCM> connManager = rcp(new CCM);
    connManager->initialize(comm,
                            Teuchos::as<panzer::GlobalOrdinal>(nx),
                            Teuchos::as<panzer::GlobalOrdinal>(ny),
                            Teuchos::as<panzer::GlobalOrdinal>(nz),
                            px,py,pz,bx,by,bz);
  
    // build the dof manager, and assocaite with the topology
    RCP<DOFManager> dofManager = rcp(new DOFManager);
    dofManager->setConnManager(connManager,*comm.getRawMpiComm());
  
    // add velocity (U) and PRESSURE fields to the MHD element block
    dofManager->addField("eblock-0_0_0","UX",pattern_U);
    dofManager->addField("eblock-0_0_0","UY",pattern_U);
    dofManager->addField("eblock-0_0_0","UZ",pattern_U);
    dofManager->addField("eblock-0_0_0","PRESSURE",pattern_P);
    dofManager->addField("eblock-0_0_0","B",pattern_B);
    dofManager->addField("eblock-0_0_0","E",pattern_E);
  
    // add velocity (U) fields to the solid element block
    dofManager->addField("eblock-0_1_0","UX",pattern_U);
    dofManager->addField("eblock-0_1_0","UY",pattern_U);
    dofManager->addField("eblock-0_1_0","UZ",pattern_U);
  
    // try to get them all synced up
    comm.barrier();

    {
      PANZER_FUNC_TIME_MONITOR("panzer::ScalingTest::buildGlobalUnknowns");
  
      dofManager->buildGlobalUnknowns();
    }
  }

  Teuchos::TimeMonitor::summarize(std::cout,false,true,false);

  if ( timingsFile != "" ){
    std::ofstream fout(timingsFile.c_str());
    Teuchos::RCP<Teuchos::ParameterList> reportParams = parameterList(* (Teuchos::TimeMonitor::getValidReportParameters()));
    reportParams->set("Report format", "YAML");
    reportParams->set("YAML style", "spacious");
    Teuchos::TimeMonitor::report(fout,reportParams);
  }
 
  // this confirms the application passes
  std::cout << "Scaling test completed" << std::endl;

  return 0;
}

template <typename IntrepidBasisType>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // build a field pattern from a single basis
  RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > basis = rcp(new IntrepidBasisType);
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}
