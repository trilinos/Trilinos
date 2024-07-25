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

#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_PureBasis.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
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
#include "Phalanx_MDField.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
  Teuchos::RCP<panzer::PureBasis> linBasis;

  //! Interpolates basis DOF values to IP DOF values
  template<typename EvalT, typename Traits>
  class XCoordinate
    :
    public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {
    public:

      XCoordinate(
        const Teuchos::ParameterList& p);

      void
      evaluateFields(
        typename Traits::EvalData d);

    private:

      using ScalarT = typename EvalT::ScalarT;
     PHX::MDField<ScalarT,Cell,NODE> xcoord;
     int nodes;
  }; // end of class XCoordinate


  template<typename EvalT, typename Traits>
  XCoordinate<EvalT, Traits>::
  XCoordinate(
    const Teuchos::ParameterList& p)
  {
     nodes = 4;
     if(p.isType<int>("Nodes"))
        nodes = p.get<int>("Nodes");

     xcoord = PHX::MDField<ScalarT,Cell,NODE>("x-coord",p.get<Teuchos::RCP<PHX::DataLayout> >("Data Layout"));
     this->addEvaluatedField(xcoord);
  }

  template<typename EvalT, typename Traits>
  void
  XCoordinate<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
     std::size_t numcells = workset.num_cells;
     int l_nodes = nodes;
     auto xcoord_v = PHX::as_view(xcoord);
     auto cnc = PHX::as_view(this->wda(workset).cell_node_coordinates);

     Kokkos::parallel_for(numcells, KOKKOS_LAMBDA (int n) {
        for(int v=0;v<l_nodes;v++) {
	  xcoord_v(n,v) = cnc(n,v,0);
        }
       });
     Kokkos::fence();
  }

  Teuchos::RCP<panzer::PureBasis> buildLinearBasis(std::size_t worksetSize);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY,bool solution);
  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(scatter_field_evaluators, gather_constr)
  {

    const std::size_t workset_size = 20;
    linBasis = buildLinearBasis(workset_size);
    Kokkos::push_finalize_hook( [=] {
      linBasis = Teuchos::RCP<panzer::PureBasis>();
    });

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(20,20,true);

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm
        = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
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

    std::string scatterName = "xcoord-scatter-residual";
    Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer_stk::ScatterFields<panzer::Traits::Residual,panzer::Traits>(scatterName,mesh,linBasis,*fieldNames));
    fm->registerEvaluator<panzer::Traits::Residual>(eval);
    fm->requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
      user_app::BCFactory bc_factory;
      const int default_integration_order = 1;

      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
				 eqset_factory,
				 gd,
				 false,
				 physicsBlocks);
    }


    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(9+4);
    fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer::PhysicsBlock> physics_block_one = panzer::findPhysicsBlock("eblock-0_0",physicsBlocks);

    Teuchos::RCP<std::vector<panzer::Workset> > volume_worksets = panzer_stk::buildWorksets(*mesh,physics_block_one->elementBlockID(),
                                                                                            physics_block_one->getWorksetNeeds());

    panzer::Traits::SD sd;
    sd.worksets_ = volume_worksets;
    fm->postRegistrationSetupForType<panzer::Traits::Residual>(sd);
    fm->writeGraphvizFile<panzer::Traits::Residual>("resi-eval-graph.dot");

    std::vector<panzer::Workset> & worksets = *volume_worksets;
    panzer::Traits::PED preEvalData;
    fm->preEvaluate<panzer::Traits::Residual>(preEvalData);
    for(std::size_t ws=0;ws<worksets.size();ws++) {
       fm->evaluateFields<panzer::Traits::Residual>(worksets[ws]);
    }
    fm->postEvaluate<panzer::Traits::Residual>(0);

    if(mesh->isWritable())
       mesh->writeToExodus("x-coord.exo");
  }

  TEUCHOS_UNIT_TEST(scatter_field_evaluators, cell_field)
  {

    const std::size_t workset_size = 5;
    linBasis = buildLinearBasis(workset_size);
    Kokkos::push_finalize_hook( [=] {
      linBasis = Teuchos::RCP<panzer::PureBasis>();
    });

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(5,5,false);

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm
        = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

     Teuchos::RCP<shards::CellTopology> topo =
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
    panzer::CellData cellData(workset_size,topo);
    Teuchos::RCP<panzer::IntegrationRule> intRule = Teuchos::rcp(new panzer::IntegrationRule(1,cellData));

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
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
       pl.set("Nodes",1);
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

       std::string scatterName = "xcoord-scatter-residual";
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
             = Teuchos::rcp(new panzer_stk::ScatterFields<panzer::Traits::Residual,panzer::Traits>(scatterName,mesh,linBasis,*fieldNames));
       fm->registerEvaluator<panzer::Traits::Residual>(eval);
       fm->requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
    }

    // build physics blocks
    //////////////////////////////////////////////////////////////
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
      user_app::BCFactory bc_factory;
      const int default_integration_order = 1;

      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");

      Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
				 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
				 eqset_factory,
				 gd,
				 false,
				 physicsBlocks);
    }

    // register jacobian size
    //////////////////////////////////////////////////////////////
    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(9+4);
    fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    // build worksets
    //////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::PhysicsBlock> physics_block_one = panzer::findPhysicsBlock("eblock-0_0",physicsBlocks);
    Teuchos::RCP<std::vector<panzer::Workset> > volume_worksets = panzer_stk::buildWorksets(*mesh,physics_block_one->elementBlockID(),
                                                                                            physics_block_one->getWorksetNeeds());


    panzer::Traits::SD sd;
    sd.worksets_ = volume_worksets;
    fm->postRegistrationSetupForType<panzer::Traits::Residual>(sd);
    fm->writeGraphvizFile<panzer::Traits::Residual>("resi-eval-graph.dot");

    std::vector<panzer::Workset> & worksets = *volume_worksets;
    panzer::Traits::PED preEvalData;
    fm->preEvaluate<panzer::Traits::Residual>(preEvalData);
    for(std::size_t ws=0;ws<worksets.size();ws++) {
       fm->evaluateFields<panzer::Traits::Residual>(worksets[ws]);
    }
    fm->postEvaluate<panzer::Traits::Residual>(0);

    if(mesh->isWritable())
       mesh->writeToExodus("x-coord-cell.exo");
  }

  Teuchos::RCP<panzer::PureBasis> buildLinearBasis(std::size_t worksetSize)
  {
     Teuchos::RCP<shards::CellTopology> topo =
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);

     return Teuchos::rcp(new panzer::PureBasis("HGrad",1,cellData));
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY,bool solution)
  {
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

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
    }

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
