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
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_STKClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Panzer_ResponseContainer.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "TestEvaluators.hpp"

#include "Epetra_MpiComm.h"

#include <vector>
#include <map>
#include <string>

#include "Panzer_TypeAssocMap.hpp"

#include "Sacado_mpl_vector.hpp"

namespace panzer {

  void testInitialzation(panzer::InputPhysicsBlock& ipb,
			 std::vector<panzer::BC>& bcs);

  struct Builder {
    template <typename T>
    std::string build() const { return "other"; }
  };
  template < > std::string Builder::build<int>() const { return "Sint"; }
  template < > std::string Builder::build<short>() const { return "Sshort"; }
  template < > std::string Builder::build<char>() const { return "Schar"; }

  TEUCHOS_UNIT_TEST(type_assoc_map, test)
  {
    typedef Sacado::mpl::vector<char,short> VecType;
    TypeAssocMap<VecType,std::string> tMap;

    Builder builder;
    tMap.buildObjects(builder);
  
    TEST_EQUALITY(tMap.get<char>(),"Schar");
    TEST_EQUALITY(tMap.get<short>(),"Sshort");

    tMap.set<short>("not char");
    TEST_EQUALITY(tMap.get<short>(),"not char");
  }

  struct RespFactoryFunc_Builder {
    MPI_Comm comm;

    template <typename T>
    Teuchos::RCP<ResponseEvaluatorFactoryBase> build() const
    { return Teuchos::rcp(new ResponseEvaluatorFactory_Functional<T>(comm)); }
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(response_library_stk, test, EvalT)
  {
    using Teuchos::RCP;

  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif

    panzer_stk::SquareQuadMeshFactory mesh_factory;
    user_app::MyFactory eqset_factory;
    user_app::BCFactory bc_factory;
    const std::size_t workset_size = 20;

    panzer::FieldManagerBuilder fmb;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup physic blocks
    /////////////////////////////////////////////
    panzer::InputPhysicsBlock ipb;
    std::vector<panzer::BC> bcs;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    {
       std::map<std::string,panzer::InputPhysicsBlock> 
             physics_id_to_input_physics_blocks;

       testInitialzation(ipb, bcs);

       std::map<std::string,std::string> block_ids_to_physics_ids;
       block_ids_to_physics_ids["eblock-0_0"] = "test physics";
       block_ids_to_physics_ids["eblock-1_0"] = "test physics";

       std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
       block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
       block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
    
       physics_id_to_input_physics_blocks["test physics"] = ipb;

       Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

       panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                  block_ids_to_cell_topo,
                                  physics_id_to_input_physics_blocks,
                                  2,workset_size,
                                  eqset_factory,
				  gd,
		    	          false,
                                  physics_blocks);
    }

    // setup worksets
    /////////////////////////////////////////////
 
     std::vector<std::string> validEBlocks;
     mesh->getElementBlockNames(validEBlocks);

    // build WorksetContainer
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory 
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physics_blocks,workset_size));
 
    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager 
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    Teuchos::RCP<const panzer::UniqueGlobalIndexerFactory<int,int,int,int> > indexerFactory
          = Teuchos::rcp(new panzer::DOFManagerFactory<int,int>);
    const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
          = indexerFactory->buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physics_blocks,conn_manager);

    // and linear object factory
    Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > elof 
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(comm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup field manager builder
    /////////////////////////////////////////////
      
    // Add in the application specific closure model factory
    user_app::STKModelFactory_TemplateBuilder cm_builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory; 
    cm_factory.buildObjects(cm_builder);

    double iValue = -2.3;
    double tValue = 82.9;

    Teuchos::ParameterList closure_models("Closure Models");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("FIELD_A").set<double>("Value",tValue);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("FIELD_B").set<double>("Value",iValue);

    Teuchos::ParameterList user_data("User Data");
    user_data.sublist("Panzer Data").set("Mesh", mesh);
    user_data.sublist("Panzer Data").set("DOF Manager", dofManager);
    user_data.sublist("Panzer Data").set("Linear Object Factory", lof);
    user_data.set<int>("Workset Size",workset_size);

    // setup and evaluate ResponseLibrary
    ///////////////////////////////////////////////////
 
    out << "Adding responses" << std::endl;

    RCP<ResponseLibrary<Traits> > rLibrary 
          = Teuchos::rcp(new ResponseLibrary<Traits>(wkstContainer,dofManager,lof));
          // = Teuchos::rcp(new ResponseLibrary<Traits>());

    RespFactoryFunc_Builder builder;
    builder.comm = MPI_COMM_WORLD;
    std::vector<std::string> blocks(1);
    blocks[0] = "eblock-0_0";
    rLibrary->addResponse("FIELD_A",blocks,builder);
    blocks[0] = "eblock-1_0";
    rLibrary->addResponse("FIELD_B",blocks,builder);

    Teuchos::RCP<ResponseBase> tResp = rLibrary->getResponse<panzer::Traits::Residual>("FIELD_A");
    Teuchos::RCP<ResponseBase> iResp = rLibrary->getResponse<panzer::Traits::Residual>("FIELD_B");

    TEST_ASSERT(tResp!=Teuchos::null);
    TEST_ASSERT(iResp!=Teuchos::null);

    TEST_EQUALITY(tResp->getName(),"FIELD_A");
    TEST_EQUALITY(iResp->getName(),"FIELD_B");

    TEST_EQUALITY(tResp->getLookupName(),"RESPONSE_FIELD_A");
    TEST_EQUALITY(iResp->getLookupName(),"RESPONSE_FIELD_B");

    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(tResp,true));
    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(iResp,true));

    RCP<Epetra_Vector> eVec;
    RCP<Thyra::VectorBase<double> > tVec;
    {
      RCP<const Epetra_Map> map = Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(iResp)->getMap();
      RCP<const Thyra::VectorSpaceBase<double> > vs = Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(tResp)->getVectorSpace();

      eVec = Teuchos::rcp(new Epetra_Vector(*map)); 
      tVec = Thyra::createMember<double>(vs);
      
      Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(iResp)->setVector(eVec);
      Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(tResp)->setVector(tVec);

      // test epetra or thyra only logic
      TEST_THROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(iResp)->setVector(tVec),std::logic_error);
      TEST_THROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(tResp)->setVector(eVec),std::logic_error);
    }

    std::vector<Teuchos::RCP<ResponseBase> > v;
    v.push_back(Teuchos::null);

    rLibrary->getResponses<panzer::Traits::Residual>(v);
    TEST_EQUALITY(v.size(),2);

    TEST_ASSERT(v[0]->getName()=="FIELD_A" || v[0]->getName()=="FIELD_B");
    TEST_ASSERT(v[1]->getName()=="FIELD_A" || v[1]->getName()=="FIELD_B");

    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(v[0],true));
    TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(v[1],true));

    TEST_ASSERT(!rLibrary->responseEvaluatorsBuilt());

    rLibrary->buildResponseEvaluators(physics_blocks,
  				      cm_factory,
                                      closure_models,
  				      user_data,true);

    TEST_ASSERT(rLibrary->responseEvaluatorsBuilt());

    Teuchos::RCP<panzer::LinearObjContainer> loc = elof->buildLinearObjContainer();
    elof->initializeContainer(panzer::EpetraLinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::LinearObjContainer> gloc = elof->buildGhostedLinearObjContainer();
    elof->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X,*gloc);

    panzer::AssemblyEngineInArgs ae_inargs(gloc,loc);
    ae_inargs.addGlobalEvaluationData(iResp->getLookupName(),iResp);
    ae_inargs.addGlobalEvaluationData(tResp->getLookupName(),tResp);

    rLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);

    Teuchos::ArrayRCP<double> tData;
    Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVec)->getNonconstLocalData(Teuchos::outArg(tData));

    TEST_FLOATING_EQUALITY((*eVec)[0],0.5*iValue,1e-14);
    TEST_FLOATING_EQUALITY(tData[0],0.5*tValue,1e-14);
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

    panzer::InputEquationSet ies_2;
    ies_2.name = "Energy";
    ies_2.basis = "Q1";
    ies_2.integration_order = 1;
    ies_2.model_id = "ion solid";
    ies_2.prefix = "ION_";

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

  typedef Traits::Residual ResidualType;
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(response_library_stk,test,ResidualType)

  #ifdef HAVE_STOKHOS
  typedef Traits::SGResidual SGResidualType;
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(response_library_stk,test,SGResidualType)
  #endif

}
