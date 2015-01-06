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

#include "Kokkos_DefaultNode.hpp"
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
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_Response_Residual.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Thyra_LinearOpTester.hpp"

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  struct AssemblyPieces {
    RCP<panzer::FieldManagerBuilder> fmb;  
    RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    RCP<panzer::GlobalData> gd;
    RCP<panzer::LinearObjFactory<panzer::Traits> > lof;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager;
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer;
    Teuchos::ParameterList user_data;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    Teuchos::RCP<panzer::EquationSetFactory> eqset_factory;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    Teuchos::ParameterList closure_models;
    std::vector<panzer::BC> bcs;
  };

  void buildAssemblyPieces(bool parameter_on,
                           AssemblyPieces & ap);

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout=false);

  // Test that the response library can build the correct residual and jacobian
  TEUCHOS_UNIT_TEST(response_residual, basic)
  {
    typedef panzer::Traits::RealType RealType;
    typedef Thyra::VectorBase<RealType> VectorType;
    typedef Thyra::LinearOpBase<RealType> OperatorType;

    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    bool parameter_on = true;
    AssemblyPieces ap;
    buildAssemblyPieces(parameter_on,ap);

    // Use the model evaluator to setup a bunch of comparisons
    ///////////////////////////////////////////////////////////////////////
    RCP<VectorType> x;
    RCP<VectorType> x_dot;
    RCP<VectorType> f_me;
    RCP<OperatorType> J_me;
    double alpha = 1.3, beta = 0.2;
    {
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef panzer::ModelEvaluator<double> PME;

      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      bool build_transient_support = true;

      RCP<PME> me = Teuchos::rcp(new PME(ap.fmb,ap.rLibrary,ap.lof,p_names,Teuchos::null,ap.gd,build_transient_support,0.0));

      x = Thyra::createMember(*me->get_x_space());
      x_dot = Thyra::createMember(*me->get_x_space());

      Thyra::randomize(-1.0,1.0,x.ptr());
      Thyra::randomize(-1.0,1.0,x_dot.ptr());
      Thyra::assign(x_dot.ptr(),0.0);

      InArgs in_args = me->createInArgs();
      in_args.set_x(x);
      in_args.set_x_dot(x_dot);
      in_args.set_alpha(alpha);
      in_args.set_beta(beta);
      
      OutArgs out_args = me->createOutArgs();
      f_me = Thyra::createMember(*me->get_f_space());
      J_me = me->create_W_op();
      out_args.set_f(f_me);
      out_args.set_W_op(J_me);

      me->evalModel(in_args, out_args);
    }

    bool residualType = true;
    RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(ap.wkstContainer,ap.dofManager,ap.lof,residualType)); 
 
    // verify the residual type
    TEST_ASSERT(rLibrary->isResidualType());

    // test that responses are the right type
    {
      RCP<ResponseBase> response_residual = rLibrary->getResponse<Traits::Residual>("RESIDUAL");
      RCP<ResponseBase> response_jacobian = rLibrary->getResponse<Traits::Jacobian>("RESIDUAL");

      TEST_ASSERT(response_residual!=Teuchos::null);
      TEST_ASSERT(response_jacobian!=Teuchos::null);

      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Residual<Traits::Residual> >(response_residual,true));
      TEST_NOTHROW(Teuchos::rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(response_jacobian,true));
    }

    // build residual response evaluators
    {
      user_app::BCFactory bc_factory;

      TEST_ASSERT(!rLibrary->responseEvaluatorsBuilt());

      rLibrary->buildResidualResponseEvaluators(ap.physicsBlocks,*ap.eqset_factory,ap.bcs,bc_factory,
                                                ap.cm_factory,ap.closure_models,ap.user_data);

      TEST_ASSERT(rLibrary->responseEvaluatorsBuilt());
    }

    RCP<Response_Residual<Traits::Residual> > response_residual = 
      rcp_dynamic_cast<Response_Residual<Traits::Residual> >(rLibrary->getResponse<Traits::Residual>("RESIDUAL"));
    RCP<Response_Residual<Traits::Jacobian> > response_jacobian = 
      rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(rLibrary->getResponse<Traits::Jacobian>("RESIDUAL"));

    // evaluate residual responses
    {
      // set solution vectors
      RCP<LinearObjContainer> loc = ap.lof->buildLinearObjContainer();
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_x_th(x);
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_dxdt_th(x_dot);

      // allocate fill vectors
      RCP<LinearObjContainer> gloc = ap.lof->buildGhostedLinearObjContainer();
      ap.lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::DxDt,*gloc);

      // setup output arguments for the residual response 
      response_residual->setResidual(response_residual->allocateResidualVector());

      // setup in arguments
      AssemblyEngineInArgs ae_inargs(gloc,loc);
      ae_inargs.alpha = alpha;
      ae_inargs.beta = beta;
      ae_inargs.evaluate_transient_terms = true;
      rLibrary->addResponsesToInArgs<Traits::Residual>(ae_inargs);

      // evaluate
      rLibrary->evaluate<Traits::Residual>(ae_inargs);

      TEST_ASSERT(testEqualityOfVectorValues(*response_residual->getResidual(),*f_me,1e-16));
    }

    // evaluate jacobian responses
    {
      // set solution vectors
      RCP<LinearObjContainer> loc = ap.lof->buildLinearObjContainer();
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_x_th(x);
      rcp_dynamic_cast<ThyraObjContainer<RealType> >(loc)->set_dxdt_th(x_dot);

      // allocate fill vectors
      RCP<LinearObjContainer> gloc = ap.lof->buildGhostedLinearObjContainer();
      ap.lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::DxDt,*gloc);

      // setup output arguments for the residual response 
      response_jacobian->setJacobian(response_jacobian->allocateJacobian());

      // setup in arguments
      AssemblyEngineInArgs ae_inargs(gloc,loc);
      ae_inargs.alpha = alpha;
      ae_inargs.beta = beta;
      ae_inargs.evaluate_transient_terms = true;
      rLibrary->addResponsesToInArgs<Traits::Jacobian>(ae_inargs);

      // evaluate
      rLibrary->evaluate<Traits::Jacobian>(ae_inargs);

      Thyra::LinearOpTester<double> tester;
      tester.show_all_tests(true);
      tester.set_all_error_tol(1e-16);      
      tester.num_random_vectors(20);      

      Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
      const bool op_cmp = tester.compare( *J_me, *response_jacobian->getJacobian(), Teuchos::ptrFromRef(fout));
      TEST_ASSERT(op_cmp);
    }
  }

  bool testEqualityOfVectorValues(const Thyra::VectorBase<double> & a, 
                                  const Thyra::VectorBase<double> & b, 
                                  double tolerance, bool write_to_cout)
  {  
    bool is_equal = true;

    TEUCHOS_ASSERT(a.space()->dim() == b.space()->dim());

    Teuchos::ArrayRCP<const double> a_data,b_data;
    dynamic_cast<const Thyra::SpmdVectorBase<double> &>(a).getLocalData(Teuchos::ptrFromRef(a_data));
    dynamic_cast<const Thyra::SpmdVectorBase<double> &>(b).getLocalData(Teuchos::ptrFromRef(b_data));

    for (int i = 0; i < a_data.size(); ++i) {
      
      std::string output = "    equal!: ";
      
      if (std::fabs(a_data[i] - b_data[i]) > tolerance) {
	is_equal = false;
	output = "NOT equal!: ";
      }
      
      if (write_to_cout)
	std::cout << output << a_data[i] << " - " << b_data[i] << " = " << (a_data[i] - b_data[i]) << std::endl;
    }

    int globalSuccess = -1;
    Teuchos::RCP<const Teuchos::Comm<Teuchos::Ordinal> > comm = Teuchos::DefaultComm<Teuchos::Ordinal>::getComm();
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, is_equal ? 0 : 1, Teuchos::outArg(globalSuccess) );
    return (globalSuccess==0);
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
      p.set("Basis Order",1);
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
  
  void buildAssemblyPieces(bool parameter_on,
                           AssemblyPieces & ap)
  {
    using Teuchos::RCP;
  
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",6);
    pl->set("Y Elements",4);
    
    panzer_stk_classic::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> & bcs = ap.bcs;
    testInitialzation(ipb, bcs);
    

    ap.fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    ap.eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    ap.gd = panzer::createGlobalData();
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");
      
      int default_integration_order = 1;
      
      bool build_transient_support = true;
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
                                 ap.eqset_factory,
				 ap.gd,
			         build_transient_support,
                                 ap.physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    // build WorksetContainer
    Teuchos::RCP<panzer_stk_classic::WorksetFactory> wkstFactory 
       = Teuchos::rcp(new panzer_stk_classic::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,ap.physicsBlocks,workset_size));
    ap.wkstContainer = wkstContainer;

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<int>(mesh));

    panzer::DOFManagerFactory<int,int> globalIndexerFactory;
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager 
         = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),ap.physicsBlocks,conn_manager);
    ap.dofManager = dofManager;

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
        // = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64>(Comm.getConst(),dofManager));
        = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpiComm,dofManager));
    ap.lof = linObjFactory;
    linObjFactory = ap.lof;

    ap.rLibrary = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,ap.dofManager,linObjFactory)); 

    // setup field manager build
    /////////////////////////////////////////////////////////////
 
    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    ap.cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    if(parameter_on)
       closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<std::string>("Type","Parameter");
    closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("DENSITY").set<double>("Value",1.0);
    closure_models.sublist("solid").sublist("HEAT_CAPACITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").set<double>("Value",1.0);
    closure_models.sublist("ion solid").sublist("ION_HEAT_CAPACITY").set<double>("Value",1.0);
    ap.closure_models = closure_models;

    ap.user_data = Teuchos::ParameterList("User Data");

    ap.fmb->setWorksetContainer(wkstContainer);
    ap.fmb->setupVolumeFieldManagers(ap.physicsBlocks,ap.cm_factory,closure_models,*linObjFactory,ap.user_data);
    ap.fmb->setupBCFieldManagers(bcs,ap.physicsBlocks,*ap.eqset_factory,ap.cm_factory,bc_factory,closure_models,*linObjFactory,ap.user_data);
  }


}
