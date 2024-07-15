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

#include "Panzer_NodeType.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ExplicitModelEvaluator.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include <cstdio> // for get char
#include <fstream>

namespace panzer {

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  void buildAssemblyPieces(bool parameter_on,
                           Teuchos::RCP<panzer::FieldManagerBuilder> & fmb,
                           Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & rLibrary,
                           Teuchos::RCP<panzer::GlobalData> & gd,
                           Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof);

  TEUCHOS_UNIT_TEST(explicit_model_evaluator, basic)
  {
    using Teuchos::RCP;


    bool parameter_on = true;
    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > rLibrary;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lof;
    Teuchos::RCP<panzer::GlobalData> gd;

    buildAssemblyPieces(parameter_on,fmb,rLibrary,gd,lof);

    // Test a transient me
    {
      typedef Thyra::ModelEvaluatorBase MEB;
      typedef Thyra::ModelEvaluatorBase::InArgs<double> InArgs;
      typedef Thyra::ModelEvaluatorBase::OutArgs<double> OutArgs;
      typedef Thyra::VectorBase<double> VectorType;
      typedef Thyra::LinearOpBase<double> OperatorType;
      typedef panzer::ModelEvaluator<double> PME;
      typedef panzer::ExplicitModelEvaluator<double> ExpPME;

      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
      std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
      bool build_transient_support = true;

      Stratimikos::DefaultLinearSolverBuilder builder;
      Teuchos::RCP<Teuchos::ParameterList> validList = Teuchos::rcp(new Teuchos::ParameterList(*builder.getValidParameters()));
      builder.setParameterList(validList);
      RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = builder.createLinearSolveStrategy("Amesos");

      RCP<PME> me = Teuchos::rcp(new PME(fmb,rLibrary,lof,p_names,p_values,lowsFactory,gd,build_transient_support,0.0));
      RCP<ExpPME> exp_me = Teuchos::rcp(new ExpPME(me,true,false)); // constant mass, use lumped

      RCP<VectorType> exp_f, f;

      // explicit evaluation
      {
        // set the nominal values
        InArgs nom_vals = exp_me->getNominalValues();
        TEST_ASSERT(nom_vals.supports(MEB::IN_ARG_x));
        TEST_ASSERT(nom_vals.supports(MEB::IN_ARG_x_dot)); // this is supported for stabilization purposes
        TEST_ASSERT(nom_vals.supports(MEB::IN_ARG_alpha)); // alpha and beta support needed for outputting responses
        TEST_ASSERT(nom_vals.supports(MEB::IN_ARG_beta));

        // create in args
        InArgs in_args = exp_me->createInArgs();
        TEST_ASSERT(in_args.supports(MEB::IN_ARG_x));
        TEST_ASSERT(in_args.supports(MEB::IN_ARG_x_dot)); // this is supported for stabilization purposes
        TEST_ASSERT(in_args.supports(MEB::IN_ARG_alpha)); // alpha and beta support needed for outputting responses
        TEST_ASSERT(in_args.supports(MEB::IN_ARG_beta));
        InArgs nomValues = exp_me->getNominalValues();
        RCP<VectorType> x = Thyra::createMember(*exp_me->get_x_space());
        RCP<VectorType> x_dot = Thyra::createMember(*exp_me->get_x_space());
        Thyra::assign(x_dot.ptr(),0.0);
        Thyra::assign(x.ptr(),5.0);
        in_args.set_x(x);
        in_args.set_x_dot(x_dot);

        // create out args
        OutArgs out_args = exp_me->createOutArgs();
        TEST_ASSERT(out_args.supports(MEB::OUT_ARG_f));
        TEST_ASSERT(!out_args.supports(MEB::OUT_ARG_W_op));
        TEST_ASSERT(!out_args.supports(MEB::OUT_ARG_W));
        exp_f = Thyra::createMember(*exp_me->get_f_space());
        out_args.set_f(exp_f);

        exp_me->evalModel(in_args, out_args);
      }

      // implicit evaluation
      RCP<OperatorType> mass = me->create_W_op();
      {
        // create in args
        InArgs in_args = me->createInArgs();
        InArgs nomValues = me->getNominalValues();
        RCP<VectorType> x = Thyra::createMember(*me->get_x_space());
        RCP<VectorType> x_dot = Thyra::createMember(*me->get_x_space());
        Thyra::assign(x_dot.ptr(),0.0);
        Thyra::assign(x.ptr(),5.0);
        in_args.set_x(x);
        in_args.set_x_dot(x_dot);

        // create out args
        OutArgs out_args = me->createOutArgs();
        f = Thyra::createMember(*me->get_f_space());
        out_args.set_f(f);

        me->evalModel(in_args, out_args);

        in_args.set_x(x);
        in_args.set_x_dot(x_dot);
        in_args.set_alpha(1.0);
        in_args.set_beta(0.0);

        out_args.set_f(Teuchos::null);
        out_args.set_W_op(mass);

        me->setOneTimeDirichletBeta(1.0);
        me->evalModel(in_args, out_args);
      }

      Teuchos::RCP<Thyra::VectorBase<double> > mass_exp_f = Thyra::createMember(*exp_me->get_f_space());
      Thyra::apply(*mass,Thyra::NOTRANS,*exp_f,mass_exp_f.ptr());

      out << "f_i = \n" << Teuchos::describe(*f,Teuchos::VERB_EXTREME) << std::endl;
      out << "f_e = \n" << Teuchos::describe(*mass_exp_f,Teuchos::VERB_EXTREME) << std::endl;

      // it should be that exp_f = -f b/c x_dot=0 in the evaluation
      Thyra::Vp_StV(mass_exp_f.ptr(),1.0,*f);

      out << "Error = " << Thyra::norm_2(*mass_exp_f) << std::endl;
      TEST_ASSERT(Thyra::norm_2(*mass_exp_f)<=1e-16);
    }
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
      p.set("Integration Order",2);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",2);
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
                           Teuchos::RCP<panzer::FieldManagerBuilder> & fmb,
                           Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & rLibrary,
                           Teuchos::RCP<panzer::GlobalData> & gd,
                           Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof
                           )
  {
    using Teuchos::RCP;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",4);
    pl->set("Y Elements",4);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Comm);

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList("Physics Blocks");
    std::vector<panzer::BC> bcs;
    testInitialzation(ipb, bcs);

    fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);

    // build physics blocks
    //////////////////////////////////////////////////////////////
    const std::size_t workset_size = 20;
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory;
    gd = panzer::createGlobalData();
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      std::map<std::string,std::string> block_ids_to_physics_ids;
      block_ids_to_physics_ids["eblock-0_0"] = "test physics";
      block_ids_to_physics_ids["eblock-1_0"] = "test physics";

      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      block_ids_to_cell_topo["eblock-0_0"] = mesh->getCellTopology("eblock-0_0");
      block_ids_to_cell_topo["eblock-1_0"] = mesh->getCellTopology("eblock-1_0");

      const int default_integration_order = 1;

      bool build_transient_support = true;
      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
				 ipb,
				 default_integration_order,
				 workset_size,
                                 eqset_factory,
				 gd,
			         build_transient_support,
                                 physicsBlocks);
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    // build WorksetContainer
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    for(size_t i=0;i<physicsBlocks.size();i++)
      wkstContainer->setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer->setWorksetSize(workset_size);

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    panzer::DOFManagerFactory globalIndexerFactory;
    RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory
        = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpiComm,dofManager,false));
    lof = linObjFactory;

    rLibrary = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,dofManager,linObjFactory));

    // setup field manager build
    /////////////////////////////////////////////////////////////

    // Add in the application specific closure model factory
    user_app::MyModelFactory_TemplateBuilder cm_builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    cm_factory.buildObjects(cm_builder);

    Teuchos::ParameterList closure_models("Closure Models");
    if (parameter_on)
      closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").
        set<std::string>("Type", "Parameter");
    else
      closure_models.sublist("solid").sublist("SOURCE_TEMPERATURE").
        set<double>("Value", 1.0);

    closure_models.sublist("solid").sublist("DENSITY").
      set<double>("Value", 1.0);
    closure_models.sublist("solid").sublist("HEAT_CAPACITY").
      set<double>("Value", 1.0);
    closure_models.sublist("ion solid").sublist("SOURCE_ION_TEMPERATURE").
      set<double>("Value", 1.0);
    closure_models.sublist("ion solid").sublist("ION_DENSITY").
      set<double>("Value", 1.0);
    closure_models.sublist("ion solid").sublist("ION_HEAT_CAPACITY").
      set<double>("Value", 1.0);

    Teuchos::ParameterList user_data("User Data");

    fmb->setWorksetContainer(wkstContainer);
    fmb->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory,user_data);

    if (parameter_on)
      gd->pl->setRealValueForAllTypes(std::string("SOURCE_TEMPERATURE"), 1.0);
  }


}
