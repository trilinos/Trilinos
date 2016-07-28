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

#include "Phalanx_KokkosUtilities.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_STKClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "TestEvaluators.hpp"

#include <vector>
#include <map>
#include <string>

#include "Panzer_TypeAssocMap.hpp"

#include "Sacado_mpl_vector.hpp"

using Teuchos::RCP;

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  std::pair<RCP<ResponseLibrary<Traits> >,RCP<LinearObjFactory<panzer::Traits> > > buildResponseLibrary(
                                                           std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physics_blocks,
                                                           panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                                           Teuchos::ParameterList & closure_models,
                                                           Teuchos::ParameterList & user_data);

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
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
    Teuchos::RCP<const panzer::UniqueGlobalIndexer<int,int> > globalIndexer;

    template <typename T>
    Teuchos::RCP<ResponseEvaluatorFactoryBase> build() const
    { return Teuchos::rcp(new ResponseEvaluatorFactory_Functional<T,int,int>(comm,1,true,"",linearObjFactory)); }
  };

  TEUCHOS_UNIT_TEST(response_library2, test)
  {

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models("Closure Models");
    Teuchos::ParameterList user_data("User Data");

    // setup and evaluate ResponseLibrary
    ///////////////////////////////////////////////////
 
    out << "Adding responses" << std::endl;

    std::pair< RCP<ResponseLibrary<Traits> >, RCP<panzer::LinearObjFactory<panzer::Traits> > > data 
          = buildResponseLibrary(physics_blocks,cm_factory,closure_models,user_data);
    RCP<ResponseLibrary<Traits> > rLibrary = data.first;
    RCP<panzer::LinearObjFactory<panzer::Traits> > lof = data.second;

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
      // eVec->PutScalar(0.0);

      // TEST_EQUALITY(eVec->MyLength(),1);

      tVec = Thyra::createMember<double>(vs);
      // Thyra::assign(tVec.ptr(),0.0);
      
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

    Teuchos::RCP<panzer::LinearObjContainer> loc = lof->buildLinearObjContainer();
    lof->initializeContainer(panzer::LinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::LinearObjContainer> gloc = lof->buildGhostedLinearObjContainer();
    lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*gloc);

    panzer::AssemblyEngineInArgs ae_inargs(gloc,loc);

    rLibrary->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
    rLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);

    Teuchos::ArrayRCP<double> tData;
    Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVec)->getNonconstLocalData(Teuchos::outArg(tData));

    double iValue = -2.3;
    double tValue = 82.9;

    TEST_FLOATING_EQUALITY((*eVec)[0],0.5*iValue,1e-14);
    TEST_FLOATING_EQUALITY(tData[0],0.5*tValue,1e-14);
  }

  TEUCHOS_UNIT_TEST(response_library2, test_surface)
  {

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models("Closure Models");
    Teuchos::ParameterList user_data("User Data");

    // setup and evaluate ResponseLibrary
    ///////////////////////////////////////////////////
 
    out << "Adding responses" << std::endl;

    std::pair< RCP<ResponseLibrary<Traits> >, RCP<panzer::LinearObjFactory<panzer::Traits> > > data 
          = buildResponseLibrary(physics_blocks,cm_factory,closure_models,user_data);
    RCP<ResponseLibrary<Traits> > rLibrary = data.first;
    RCP<panzer::LinearObjFactory<panzer::Traits> > lof = data.second;
    RCP<const panzer::UniqueGlobalIndexer<int,int> > globalIndexer 
        = user_data.sublist("Panzer Data").get<RCP<panzer::UniqueGlobalIndexer<int,int> > >("DOF Manager");

    RespFactoryFunc_Builder builder;
    builder.comm = MPI_COMM_WORLD;
    builder.linearObjFactory = lof;
    builder.globalIndexer = globalIndexer;
    std::vector<std::string> blocks(1);
    blocks[0] = "eblock-0_0";
    rLibrary->addResponse("FIELD_A",blocks,builder);

    builder.linearObjFactory = Teuchos::null;
    builder.globalIndexer = Teuchos::null;
    std::vector<std::pair<std::string,std::string> > sidesets;
    sidesets.push_back(std::make_pair("bottom","eblock-0_0")); // 0.5
    sidesets.push_back(std::make_pair("top","eblock-0_0"));    // 0.5
    sidesets.push_back(std::make_pair("right","eblock-1_0"));    // 1.0
    rLibrary->addResponse("FIELD_B",sidesets,builder);

    Teuchos::RCP<ResponseBase> blkResp = rLibrary->getResponse<panzer::Traits::Residual>("FIELD_A");
    Teuchos::RCP<ResponseBase> ssResp = rLibrary->getResponse<panzer::Traits::Residual>("FIELD_B");

    Teuchos::RCP<ResponseBase> blkRespJac = rLibrary->getResponse<panzer::Traits::Jacobian>("FIELD_A");
    TEST_ASSERT(blkRespJac!=Teuchos::null);

    // no response should be build for this one
    TEST_ASSERT(rLibrary->getResponse<panzer::Traits::Jacobian>("FIELD_B")==Teuchos::null);
 
    RCP<Epetra_Vector> eVec, eVec2;
    {
      RCP<const Epetra_Map> map = Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(ssResp)->getMap();

      eVec = Teuchos::rcp(new Epetra_Vector(*map)); 
      eVec2 = Teuchos::rcp(new Epetra_Vector(*map)); 
      
      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(blkResp)->setVector(eVec));
      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Residual> >(ssResp)->setVector(eVec2));

      RCP<Epetra_MultiVector> dVec = Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Jacobian> >(blkRespJac,true)->buildEpetraDerivative();
      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Functional<panzer::Traits::Jacobian> >(blkRespJac,true)->setDerivative(dVec));
    }

    rLibrary->buildResponseEvaluators(physics_blocks,
  				      cm_factory,
                                      closure_models,
  				      user_data,true);

    TEST_ASSERT(rLibrary->responseEvaluatorsBuilt());

    Teuchos::RCP<panzer::LinearObjContainer> loc = lof->buildLinearObjContainer();
    lof->initializeContainer(panzer::LinearObjContainer::X,*loc);
    Teuchos::RCP<panzer::LinearObjContainer> gloc = lof->buildGhostedLinearObjContainer();
    lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*gloc);


    {
      panzer::AssemblyEngineInArgs ae_inargs(gloc,loc);
      rLibrary->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
      rLibrary->evaluate<panzer::Traits::Residual>(ae_inargs);
    }

    // evaluate derivatives
    {
      panzer::AssemblyEngineInArgs ae_inargs(gloc,loc);
      rLibrary->addResponsesToInArgs<panzer::Traits::Jacobian>(ae_inargs);
      rLibrary->evaluate<panzer::Traits::Jacobian>(ae_inargs);
    }

    double iValue = -2.3;
    double tValue = 82.9;

    TEST_FLOATING_EQUALITY((*eVec)[0],0.5*tValue,1e-14);
    TEST_FLOATING_EQUALITY((*eVec2)[0],2.0*iValue,1e-14);
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
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
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

  std::pair<RCP<ResponseLibrary<Traits> >,RCP<LinearObjFactory<panzer::Traits> > > buildResponseLibrary(
                                                           std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physics_blocks,
                                                           panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                                           Teuchos::ParameterList & closure_models,
                                                           Teuchos::ParameterList & user_data)
  {
    using Teuchos::RCP;

  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif

    panzer_stk::SquareQuadMeshFactory mesh_factory;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
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
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    {
       testInitialzation(ipb, bcs);

       std::map<std::string,std::string> block_ids_to_physics_ids;
       block_ids_to_physics_ids["eblock-0_0"] = "test physics";
       block_ids_to_physics_ids["eblock-1_0"] = "test physics";

       std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
       block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
       block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
    
       Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();

      int default_integration_order = 1;
      
       panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                  block_ids_to_cell_topo,
				  ipb,
				  default_integration_order,
				  workset_size,
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
           = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));

    Teuchos::RCP<const panzer::UniqueGlobalIndexerFactory<int,int,int,int> > indexerFactory
          = Teuchos::rcp(new panzer::DOFManagerFactory<int,int>);
    const Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
          = indexerFactory->buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physics_blocks,conn_manager);

    // and linear object factory
    Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > elof 
          = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof = elof;

    // setup field manager builder
    /////////////////////////////////////////////
      
    // Add in the application specific closure model factory
    user_app::STKModelFactory_TemplateBuilder cm_builder;
    cm_factory.buildObjects(cm_builder);

    double iValue = -2.3;
    double tValue = 82.9;

    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("FIELD_A").set<double>("Value",tValue);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("FIELD_B").set<double>("Value",iValue);

    user_data.sublist("Panzer Data").set("Mesh", mesh);
    user_data.sublist("Panzer Data").set("DOF Manager", dofManager);
    user_data.sublist("Panzer Data").set("Linear Object Factory", lof);
    user_data.set<int>("Workset Size",workset_size);

    RCP<ResponseLibrary<Traits> > rLibrary 
          = Teuchos::rcp(new ResponseLibrary<Traits>(wkstContainer,dofManager,lof));

    return std::make_pair(rLibrary,lof);
  }

}
