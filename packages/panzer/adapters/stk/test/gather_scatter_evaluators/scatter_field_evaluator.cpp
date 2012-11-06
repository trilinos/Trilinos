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

#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_GatherFields.hpp"
#include "Panzer_STK_ScatterFields.hpp"
#include "Panzer_STK_ScatterCellAvgQuantity.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

#include "Epetra_MpiComm.h"

#include <cstdio> // for get char
#include <vector>
#include <string>

namespace panzer {
  Teuchos::RCP<panzer::BasisIRLayout> linBasis;

  //! Interpolates basis DOF values to IP DOF values
  PHX_EVALUATOR_CLASS(XCoordinate)
     PHX::MDField<ScalarT,Cell,NODE> xcoord;
  PHX_EVALUATOR_CLASS_END

  PHX_EVALUATOR_CTOR(XCoordinate,p)
  {
     xcoord = PHX::MDField<ScalarT,Cell,NODE>("x-coord",p.get<Teuchos::RCP<PHX::DataLayout> >("Data Layout"));
     this->addEvaluatedField(xcoord);
  }

  PHX_POST_REGISTRATION_SETUP(XCoordinate,sd,fm)
  { this->utils.setFieldData(xcoord,fm); }

  PHX_EVALUATE_FIELDS(XCoordinate,workset)
  { 
     std::size_t numcells = workset.num_cells;

     for(std::size_t n=0;n<numcells;n++) {
        for(std::size_t v=0;v<4;v++) {
           xcoord(n,v) = workset.cell_vertex_coordinates(n,v,0);
        }
     }
  }

  Teuchos::RCP<panzer::BasisIRLayout> buildLinearBasis(std::size_t worksetSize);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY,bool solution);
  void testInitialzation(panzer::InputPhysicsBlock& ipb,std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(gs_evaluators, gather_constr)
  {
    const std::size_t workset_size = 20;
    linBasis = buildLinearBasis(workset_size);

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(20,20,true);

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::ParameterList pl;
    pl.set("Data Layout",linBasis->functional);
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm = 
      Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    fm->registerEvaluator<panzer::Traits::Residual>(Teuchos::rcp(new XCoordinate<panzer::Traits::Residual,panzer::Traits>(pl)));

    Teuchos::RCP<std::vector<std::string> > fieldNames
          = Teuchos::rcp(new std::vector<std::string>);
    fieldNames->push_back("x-coord");

    pl.set("Mesh",mesh);
    pl.set("Basis",linBasis);
    pl.set("Field Names",fieldNames);
    pl.set("Scatter Name", "xcoord-scatter-residual");
    Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer_stk::ScatterFields<panzer::Traits::Residual,panzer::Traits>(pl));
    fm->registerEvaluator<panzer::Traits::Residual>(eval);
    fm->requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,panzer::InputPhysicsBlock> eb_id_to_ipb;
    eb_id_to_ipb["eblock-0_0"] = ipb;
    
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,eb_id_to_ipb, workset_size);

    panzer::Traits::SetupData sd;
    sd.worksets_ = volume_worksets["eblock-0_0"];
    fm->postRegistrationSetupForType<panzer::Traits::Residual>(sd);
    fm->writeGraphvizFile<panzer::Traits::Residual>("resi-eval-graph.dot");

    std::vector<panzer::Workset> & worksets = *volume_worksets["eblock-0_0"];
    GlobalEvaluationDataContainer preEvalData;
    fm->preEvaluate<panzer::Traits::Residual>(preEvalData);
    for(std::size_t ws=0;ws<worksets.size();ws++) {
       fm->evaluateFields<panzer::Traits::Residual>(worksets[ws]);
    }
    fm->postEvaluate<panzer::Traits::Residual>(0);

    if(mesh->isWritable()) 
       mesh->writeToExodus("x-coord.exo");
  }

  TEUCHOS_UNIT_TEST(gs_evaluators, cell_field)
  {
    const std::size_t workset_size = 20;
    linBasis = buildLinearBasis(workset_size);

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(20,20,false);

    RCP<Epetra_Comm> Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    panzer::CellData cellData(workset_size,2,topo);
    Teuchos::RCP<panzer::IntegrationRule> intRule = Teuchos::rcp(new panzer::IntegrationRule(1,cellData));

    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm = 
      Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    {
       Teuchos::ParameterList pl;
       pl.set("Data Layout",linBasis->functional);
       fm->registerEvaluator<panzer::Traits::Residual>(Teuchos::rcp(new XCoordinate<panzer::Traits::Residual,panzer::Traits>(pl)));
    }

    {
       Teuchos::ParameterList pl;
       pl.set("Data Layout",intRule->dl_scalar);
       fm->registerEvaluator<panzer::Traits::Residual>(Teuchos::rcp(new XCoordinate<panzer::Traits::Residual,panzer::Traits>(pl)));
    }

    {
       Teuchos::RCP<std::vector<std::string> > fieldNames
             = Teuchos::rcp(new std::vector<std::string>);
       fieldNames->push_back("x-coord");
    
       Teuchos::ParameterList pl;
       pl.set("Mesh",mesh);
       pl.set("IR",intRule);
       pl.set("Field Names",fieldNames);
       pl.set("Scatter Name", "xcoord-scatter-cell-residual");
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
             = Teuchos::rcp(new panzer_stk::ScatterCellAvgQuantity<panzer::Traits::Residual,panzer::Traits>(pl));
       fm->registerEvaluator<panzer::Traits::Residual>(eval);
       fm->requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
    }

    {
       Teuchos::RCP<std::vector<std::string> > fieldNames
             = Teuchos::rcp(new std::vector<std::string>);
       fieldNames->push_back("x-coord");

       Teuchos::ParameterList pl;
       pl.set("Mesh",mesh);
       pl.set("Basis",linBasis);
       pl.set("Field Names",fieldNames);
       pl.set("Scatter Name", "xcoord-scatter-residual");
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
             = Teuchos::rcp(new panzer_stk::ScatterFields<panzer::Traits::Residual,panzer::Traits>(pl));
       fm->registerEvaluator<panzer::Traits::Residual>(eval);
       fm->requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    std::map<std::string,panzer::InputPhysicsBlock> eb_id_to_ipb;
    eb_id_to_ipb["eblock-0_0"] = ipb;
    
    std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
      volume_worksets = panzer_stk::buildWorksets(*mesh,eb_id_to_ipb, workset_size);

    panzer::Traits::SetupData sd;
    sd.worksets_ = volume_worksets["eblock-0_0"];
    fm->postRegistrationSetupForType<panzer::Traits::Residual>(sd);
    fm->writeGraphvizFile<panzer::Traits::Residual>("resi-eval-graph.dot");

    std::vector<panzer::Workset> & worksets = *volume_worksets["eblock-0_0"];
    GlobalEvaluationDataContainer preEvalData;
    fm->preEvaluate<panzer::Traits::Residual>(preEvalData);
    for(std::size_t ws=0;ws<worksets.size();ws++) {
       fm->evaluateFields<panzer::Traits::Residual>(worksets[ws]);
    }
    fm->postEvaluate<panzer::Traits::Residual>(0);

    if(mesh->isWritable()) 
       mesh->writeToExodus("x-coord-cell.exo");
  }

  Teuchos::RCP<panzer::BasisIRLayout> buildLinearBasis(std::size_t worksetSize)
  {
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,2,topo);
     panzer::IntegrationRule intRule(1,cellData);

     return Teuchos::rcp(new panzer::BasisIRLayout("Q1",intRule)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY,bool solution)
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
 
    // add in some fields
    mesh->addSolutionField("x-coord","eblock-0_0");
    if(!solution)
       mesh->addCellField("x-coord","eblock-0_0");

    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

    return mesh;
  }

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    panzer::InputEquationSet ies_1;
    ies_1.name = "Energy";
    ies_1.basis = "Q2";
    ies_1.integration_order = 1;
    ies_1.model_id = "solid";
    ies_1.prefix = "";
    ies_1.params.set<int>("junk", 1);

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "solid";
    ies_2.prefix = "ION_";
    ies_2.params.set<int>("junk", 1);

    ipb.physics_block_id = "4";
    ipb.eq_sets.push_back(ies_1);
    ipb.eq_sets.push_back(ies_2);


    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }    
    {
      std::size_t bc_id = 1;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }   
    {
      std::size_t bc_id = 2;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "TEMPERATURE";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
  }
}
