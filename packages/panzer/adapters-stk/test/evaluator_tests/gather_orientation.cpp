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

#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GatherOrientation.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Epetra_MpiComm.h"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

  TEUCHOS_UNIT_TEST(gather_orientation, gather_constr)
  {

    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    Teuchos::RCP<panzer::PureBasis> basis_qedge1 = buildBasis(workset_size,"QEdge1");

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb);

    const int default_int_order = 1;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                            physicsBlock->getWorksetNeeds()); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::DOFManager> dofManager = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));
    dofManager->addField(fieldName_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_qedge1->getIntrepid2Basis())));
    dofManager->setOrientationsRequired(true);
    dofManager->buildGlobalUnknowns();

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<PHX::FieldTag> evalField_q1, evalField_qedge1;
    {
       Teuchos::RCP<std::vector<std::string> > dofNames = Teuchos::rcp(new std::vector<std::string>);
       dofNames->push_back(fieldName_q1);

       Teuchos::ParameterList pl;
       pl.set("Indexer Names",dofNames);
       pl.set("DOF Names",dofNames);
       pl.set("Basis",basis_q1);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator
         = Teuchos::rcp(new panzer::GatherOrientation<panzer::Traits::Residual,panzer::Traits,int,panzer::GlobalOrdinal>(dofManager,pl));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);
       evalField_q1 = evaluator->evaluatedFields()[0];

       TEST_EQUALITY(evalField_q1->name(),basis_q1->name()+" Orientation");
       TEST_EQUALITY(evalField_q1->dataLayout().extent(0),basis_q1->functional->extent(0));
       TEST_EQUALITY(evalField_q1->dataLayout().extent(1),basis_q1->functional->extent(1));

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }
    {
       Teuchos::RCP<std::vector<std::string> > dofNames = Teuchos::rcp(new std::vector<std::string>);
       dofNames->push_back(fieldName_qedge1);

       Teuchos::ParameterList pl;
       pl.set("Indexer Names",dofNames);
       pl.set("DOF Names",dofNames);
       pl.set("Basis",basis_qedge1);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
         = Teuchos::rcp(new panzer::GatherOrientation<panzer::Traits::Residual,panzer::Traits,int,panzer::GlobalOrdinal>(dofManager,pl));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);
       evalField_qedge1 = evaluator->evaluatedFields()[0];

       TEST_EQUALITY(evalField_qedge1->name(),basis_qedge1->name()+" Orientation");
       TEST_EQUALITY(evalField_qedge1->dataLayout().extent(0),basis_qedge1->functional->extent(0));
       TEST_EQUALITY(evalField_qedge1->dataLayout().extent(1),basis_qedge1->functional->extent(1));

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);

    // <cell,basis>
    PHX::MDField<panzer::Traits::Residual::ScalarT> 
       fieldData_q1(evalField_q1->name(),basis_q1->functional);
    // <cell,basis>
    PHX::MDField<panzer::Traits::Residual::ScalarT> 
       fieldData_qedge1(evalField_qedge1->name(),basis_qedge1->functional);

    fm.getFieldData<panzer::Traits::Residual>(fieldData_q1);
    fm.getFieldData<panzer::Traits::Residual>(fieldData_qedge1);

    for(int i=0;i<static_cast<int>(fieldData_q1.size());i++) {
       TEST_EQUALITY(fieldData_q1[i],1);
    }

    for(int i=0;i<fieldData_qedge1.extent_int(0);i++) {
       TEST_EQUALITY(fieldData_qedge1(i,0), 1);
       TEST_EQUALITY(fieldData_qedge1(i,1), 1);
       TEST_EQUALITY(fieldData_qedge1(i,2),-1);
       TEST_EQUALITY(fieldData_qedge1(i,3),-1);
    }
  }

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName)
  { 
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);
     return Teuchos::rcp(new panzer::PureBasis(basisName,1,cellData)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",elemX);
    pl->set("Y Elements",elemY);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

    return mesh;
  }

  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb)
  {
    // Physics block
    ipb->setName("test physics");
    {
      Teuchos::ParameterList& p = ipb->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = ipb->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HCurl");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    
  }

}
