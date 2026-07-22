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

#include "Panzer_PureBasis.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

#include "user_app_STKClosureModel_Factory_TemplateBuilder.hpp"

#include <vector>
#include <map>
#include <string>

namespace panzer {

  TEUCHOS_UNIT_TEST(initial_condition_control, control)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set<int>("X Elements",2);
      pl->set<int>("Y Elements",2);
      pl->set<int>("X Blocks",2);
      pl->set<int>("Y Blocks",1);

      panzer_stk::SquareQuadMeshFactory mesh_factory;
      mesh_factory.setParameterList(pl);
      mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
      mesh->writeToExodus("initial_condition_control.exo");
    }

    RCP<const shards::CellTopology> ct = mesh->getCellTopology("eblock-0_0");

    panzer::CellData cellData(4,ct);

    RCP<panzer::IntegrationRule> int_rule = rcp(new panzer::IntegrationRule(2,cellData));

    ICFieldDescriptor densDesc;
    densDesc.fieldName  = "DENSITY";
    densDesc.basisName  = "Const";
    densDesc.basisOrder = 0;

    ICFieldDescriptor condDesc;
    condDesc.fieldName  = "CONDUCTIVITY";
    condDesc.basisName  = "HGrad";
    condDesc.basisOrder = 1;

    RCP<panzer::PureBasis> const_basis = rcp(new panzer::PureBasis(densDesc.basisName,densDesc.basisOrder,cellData));
    RCP<panzer::PureBasis> hgrad_basis = rcp(new panzer::PureBasis(condDesc.basisName,condDesc.basisOrder,cellData));

    RCP<const panzer::FieldPattern> constFP = rcp(new panzer::Intrepid2FieldPattern(const_basis->getIntrepid2Basis()));
    RCP<const panzer::FieldPattern> hgradFP = rcp(new panzer::Intrepid2FieldPattern(hgrad_basis->getIntrepid2Basis()));

    // setup DOF manager
    /////////////////////////////////////////////
    RCP<panzer::ConnManager> conn_manager
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dofManager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));
    dofManager->addField(densDesc.fieldName, constFP);
    dofManager->addField(condDesc.fieldName, hgradFP);
    dofManager->buildGlobalUnknowns();

    Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > elof
          = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup worksets
    /////////////////////////////////////////////

    std::map<std::string,panzer::WorksetNeeds> needs;
    needs["eblock-0_0"].cellData = cellData;
    needs["eblock-0_0"].int_rules.push_back(int_rule);
    needs["eblock-0_0"].bases          = { const_basis,    hgrad_basis};
    needs["eblock-0_0"].rep_field_name = {   densDesc.fieldName, condDesc.fieldName};

    needs["eblock-1_0"].cellData = cellData;
    needs["eblock-1_0"].int_rules.push_back(int_rule);
    needs["eblock-1_0"].bases          = { const_basis,    hgrad_basis};
    needs["eblock-1_0"].rep_field_name = {   densDesc.fieldName, condDesc.fieldName};

    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
        = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer              // attach it to a workset container (uses lazy evaluation)
        = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,needs));

    // setup field manager builder
    /////////////////////////////////////////////

    // Add in the application specific closure model factory
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    user_app::STKModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList user_data("User Data");

    Teuchos::ParameterList ic_closure_models("Initial Conditions");
    ic_closure_models.sublist("eblock-0_0").sublist(densDesc.fieldName).set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-0_0").sublist(condDesc.fieldName).set<double>("Value",9.0);
    ic_closure_models.sublist("eblock-1_0").sublist(densDesc.fieldName).set<double>("Value",3.0);
    ic_closure_models.sublist("eblock-1_0").sublist(condDesc.fieldName).set<double>("Value",9.0);

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
    block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

    std::map<std::string,std::vector<ICFieldDescriptor> > block_ids_to_fields;
    block_ids_to_fields["eblock-0_0"] = {densDesc,condDesc};
    block_ids_to_fields["eblock-1_0"] = {densDesc,condDesc};

    int workset_size = 4;

    Teuchos::RCP<panzer::LinearObjContainer> loc  = lof->buildLinearObjContainer();
    lof->initializeContainer(panzer::LinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::EpetraLinearObjContainer> eloc = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
    Teuchos::RCP<Thyra::VectorBase<double> > vec = eloc->get_x_th();

    // this is the Function under test
    panzer::setupControlInitialCondition(block_ids_to_cell_topo,
                                         block_ids_to_fields,
                                         *wkstContainer,
                                         *lof,cm_factory,ic_closure_models,user_data,
                                         workset_size,
                                         0.0, // t0
                                         vec);

    Teuchos::RCP<Epetra_Vector> x = eloc->get_x();
    out << x->GlobalLength() << " " << x->MyLength() << std::endl;
    for (int i=0; i < x->MyLength(); ++i) {
      double v = (*x)[i];
      TEST_ASSERT(v==3.0 || v==9.0);
    }

  }

}
