// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GatherOrientation.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Epetra_MpiComm.h"

#include <cstdio> // for get char
#include <vector>
#include <string>

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(panzer::InputPhysicsBlock& ipb);
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

    panzer::InputPhysicsBlock ipb;
    testInitialization(ipb);
    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,"eblock-0_0",ipb,workset_size); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::DOFManagerFEI<int,int> > dofManager = Teuchos::rcp(new panzer::DOFManagerFEI<int,int>(conn_manager,MPI_COMM_WORLD));
    dofManager->addField(fieldName_q1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_q1->getIntrepidBasis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_qedge1->getIntrepidBasis())));
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
          = Teuchos::rcp(new panzer::GatherOrientation<panzer::Traits::Residual,panzer::Traits,int,int>(dofManager,pl));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);
       evalField_q1 = evaluator->evaluatedFields()[0];

       TEST_EQUALITY(evalField_q1->name(),fieldName_q1+" Orientation");
       TEST_EQUALITY(evalField_q1->dataLayout().dimension(0),basis_q1->functional->dimension(0));
       TEST_EQUALITY(evalField_q1->dataLayout().dimension(1),basis_q1->functional->dimension(1));

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
          = Teuchos::rcp(new panzer::GatherOrientation<panzer::Traits::Residual,panzer::Traits,int,int>(dofManager,pl));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);
       evalField_qedge1 = evaluator->evaluatedFields()[0];

       TEST_EQUALITY(evalField_qedge1->name(),fieldName_qedge1+" Orientation");
       TEST_EQUALITY(evalField_qedge1->dataLayout().dimension(0),basis_qedge1->functional->dimension(0));
       TEST_EQUALITY(evalField_qedge1->dataLayout().dimension(1),basis_qedge1->functional->dimension(1));

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.ghostedLinContainer = Teuchos::null;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);

    PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
       fieldData_q1(evalField_q1->name(),basis_q1->functional);
    PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
       fieldData_qedge1(evalField_qedge1->name(),basis_qedge1->functional);

    fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData_q1);
    fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData_qedge1);

    for(int i=0;i<fieldData_q1.size();i++) {
       TEST_EQUALITY(fieldData_q1[i],1);
    }

    for(int i=0;i<fieldData_qedge1.dimension(0);i++) {
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

     panzer::CellData cellData(worksetSize,2,topo);
     return Teuchos::rcp(new panzer::PureBasis(basisName,cellData)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
    typedef panzer_stk::STK_Interface::VectorFieldType CoordinateField;

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

  void testInitialization(panzer::InputPhysicsBlock& ipb)
  {
    panzer::InputEquationSet ies;
    ies.name = "Energy";
    ies.basis = "Q1";
    ies.integration_order = 1;
    ies.model_id = "solid";
    ies.prefix = "";

    panzer::InputEquationSet iesb;
    iesb.name = "Energy";
    iesb.basis = "QEdge1";
    iesb.integration_order = 1;
    iesb.model_id = "solid";
    iesb.prefix = "";

    ipb.physics_block_id = "1";
    ipb.eq_sets.push_back(ies);
    ipb.eq_sets.push_back(iesb);
  }

}
