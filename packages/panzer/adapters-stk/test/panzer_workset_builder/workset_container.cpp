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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_DOFManager.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(workset_container, identifiers)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    std::string sideset = "left";
    int workset_size = 10;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",8);
    pl->set("Y Elements",8);
    pl->set("Z Elements",8);

    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const RCP<panzer::ConnManager> 
      conn_manager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dof_manager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> field_pattern = rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      dof_manager->addField(element_block, "T", field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////
    
    RCP<panzer_stk::WorksetFactory> wkstFactory 
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer);

    // I'm surprised the needs are required
    { 
      BasisDescriptor basis_desc(1,"HGrad");
      WorksetNeeds needs;
      needs.cellData = CellData(workset_size,mesh->getCellTopology(element_block));
      needs.addBasis(basis_desc);
      wkstContainer->setNeeds(element_block,needs);
    }

    wkstContainer->setFactory(wkstFactory);
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    // check volume worksets, both for deviation from 0 (required) and number
    // of unique identifiers
    {
      std::vector<Workset> & worksets = *wkstContainer->getWorksets(blockDescriptor(element_block));

      std::set<std::size_t> identifiers;
      for(auto wkst : worksets) {
        // check a unique identifier 
        TEST_ASSERT(wkst.getIdentifier()!=0u);
 
        identifiers.insert(wkst.getIdentifier());
      }

      // check uniqueness of identifiers
      TEST_EQUALITY(worksets.size(),identifiers.size());
    }

    // check side set worksets, both for deviation from 0 (required) and number
    // of unique identifiers
    {

      Teuchos::RCP<std::map<unsigned,Workset> > rcp_worksets = wkstContainer->getSideWorksets(sidesetDescriptor(element_block,sideset));

      // because this is a boundary workset, sometimes these things are not relevant
      if(rcp_worksets!=Teuchos::null) {
        std::map<unsigned,Workset> & worksets = *rcp_worksets;

        TEST_ASSERT(worksets.size()>0);

        std::set<std::size_t> identifiers;
        for(auto wkst : worksets) {
          out << "IS THIS ZERO ???? " << wkst.second.getIdentifier() << " " << wkst.first << std::endl; 

          // check a unique identifier 
          TEST_ASSERT(wkst.second.getIdentifier()!=0u);
 
          identifiers.insert(wkst.second.getIdentifier());
        }

        // check uniqueness of identifiers
        TEST_EQUALITY(worksets.size(),identifiers.size());
      }
    }

    // check side set worksets, both for deviation from 0 (required) and number
    // of unique identifiers
    {
      Teuchos::RCP<std::map<unsigned,Workset> > rcp_worksets = wkstContainer->getSideWorksets(sidesetVolumeDescriptor(element_block,sideset));


      if(rcp_worksets!=Teuchos::null) {
        std::map<unsigned,Workset> & worksets = *rcp_worksets;

        TEST_ASSERT(worksets.size()>0);

        std::set<std::size_t> identifiers;
        for(auto wkst : worksets) {
          // check a unique identifier 
          TEST_ASSERT(wkst.second.getIdentifier()!=0u);
 
          identifiers.insert(wkst.second.getIdentifier());
        }

        // check uniqueness of identifiers
        TEST_EQUALITY(worksets.size(),identifiers.size());
      }
    }
  }
}
