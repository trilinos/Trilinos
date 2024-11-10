// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

#include "Panzer_NodeType.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_ModelEvaluator.hpp"

#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupLOWSFactory.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

#include "NOX_Thyra.H"

#include "Step02_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Step02_EquationSetFactory.hpp"
#include "Step02_BCStrategy_Factory.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>

Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                          const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
                          const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
                          const Teuchos::RCP<panzer::GlobalIndexer> & globalIndexer,
                          const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                          const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                          const Teuchos::ParameterList & closure_model_pl,
                          const Teuchos::ParameterList & user_data);

void writeToExodus(const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                   const panzer::ModelEvaluator<double> & model,
                   panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
                   panzer_stk::STK_Interface & mesh);

int main(int argc, char *argv[])
{
  typedef panzer::ModelEvaluator<double> PME;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  PHX::InitializeKokkosDevice();

  int status = 0;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));
  if (mpiSession.getNProc() > 1) {
    out->setShowProcRank(true);
    out->setOutputToRootOnly(0);
  }

  try {

    Teuchos::RCP<Teuchos::Time> total_time =
      Teuchos::TimeMonitor::getNewTimer("User App: Total Time");

    Teuchos::TimeMonitor timer(*total_time);

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm
        = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    // Parse the command line arguments
    std::string input_file_name = "input.xml";
    {
      Teuchos::CommandLineProcessor clp;

      clp.setOption("i", &input_file_name, "User_App input xml filename");

      Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
         clp.parse(argc,argv,&std::cerr);

      TEUCHOS_TEST_FOR_EXCEPTION(parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL,
                            std::runtime_error, "Failed to parse command line!");
    }

    // Parse the input file and broadcast to other processes
    Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
    Teuchos::updateParametersFromXmlFileAndBroadcast(input_file_name, input_params.ptr(), *comm);

    RCP<Teuchos::ParameterList> mesh_pl             = rcp(new Teuchos::ParameterList(input_params->sublist("Mesh")));
    RCP<Teuchos::ParameterList> physics_blocks_pl   = rcp(new Teuchos::ParameterList(input_params->sublist("Physics Blocks")));
    RCP<Teuchos::ParameterList> lin_solver_pl       = rcp(new Teuchos::ParameterList(input_params->sublist("Linear Solver")));
    Teuchos::ParameterList & block_to_physics_pl    = input_params->sublist("Block ID to Physics ID Mapping");
    Teuchos::ParameterList & bcs_pl                 = input_params->sublist("Boundary Conditions");
    Teuchos::ParameterList & closure_models_pl      = input_params->sublist("Closure Models");
    Teuchos::ParameterList & user_data_pl           = input_params->sublist("User Data");

    user_data_pl.set<RCP<const Teuchos::Comm<int> > >("Comm", comm);

    RCP<panzer::GlobalData> globalData = panzer::createGlobalData();
    RCP<user_app::EquationSetFactory> eqset_factory = Teuchos::rcp(new user_app::EquationSetFactory);
    user_app::BCStrategyFactory bc_factory;

    user_app::ClosureModelFactory_TemplateBuilder cm_builder;
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    cm_factory.buildObjects(cm_builder);

    // read in mesh database, build un committed data
    ////////////////////////////////////////////////////////////////
    RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareQuadMeshFactory);
    mesh_factory->setParameterList(mesh_pl);

    RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);

    // read in physics blocks
    ////////////////////////////////////////////////////////////

    std::map<std::string,std::string> block_ids_to_physics_ids;
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, block_to_physics_pl);

    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    for(auto itr=block_ids_to_physics_ids.begin();itr!=block_ids_to_physics_ids.end();itr++)
      block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;

    // setup some defaults
    int workset_size = 20;
    int default_integration_order = 2;
    bool build_transient_support = false;
    std::vector<std::string> tangentParamNames;

    panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                               block_ids_to_cell_topo,
                               physics_blocks_pl,
                               default_integration_order,
                               workset_size,
                               eqset_factory,
                               globalData,
                               build_transient_support,
                               physicsBlocks,
                               tangentParamNames);

   // Add fields to the mesh data base (this is a peculiarity of how STK classic requires the
   // fields to be setup)
   //////////////////////////////////////////////////////////////////////////////////////////
   for(std::size_t i=0;i<physicsBlocks.size();i++) {
      RCP<panzer::PhysicsBlock> pb = physicsBlocks[i]; // we are assuming only one physics block

      const std::vector<panzer::StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end());

      // build string for modifiying vectors
      std::vector<std::string> dimenStr(3);
      dimenStr[0] = "X"; dimenStr[1] = "Y"; dimenStr[2] = "Z";

      // add basis to DOF manager: block specific
      std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
         Teuchos::RCP<const panzer::PureBasis> basis = fieldItr->second;
         if(basis->getElementSpace()==panzer::PureBasis::HGRAD)
            mesh->addSolutionField(fieldItr->first,pb->elementBlockID());
         else if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
            for(int i=0;i<basis->dimension();i++)
               mesh->addCellField(fieldItr->first+dimenStr[i],pb->elementBlockID());
         }
      }

      mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
    }

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // build the state dof manager and LOF
    RCP<panzer::GlobalIndexer> dofManager;
    RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    {
      panzer::DOFManagerFactory globalIndexerFactory;
      dofManager = globalIndexerFactory.buildGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
      linObjFactory = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(comm,dofManager));
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
    wkstContainer->setGlobalIndexer(dofManager);
    wkstContainer->setWorksetSize(workset_size);

    // build linear solver
    /////////////////////////////////////////////////////////////

    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory
        = panzer_stk::buildLOWSFactory(false, dofManager, conn_manager,
                                               Teuchos::as<int>(mesh->getDimension()),
                                               comm, lin_solver_pl,Teuchos::null);

    // build and setup model evaluatorlinear solver
    /////////////////////////////////////////////////////////////

    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs,bcs_pl,globalData);

    RCP<PME> physics = Teuchos::rcp(new PME(linObjFactory,lowsFactory,globalData,build_transient_support,0.0));
    physics->setupModel(wkstContainer,physicsBlocks,bcs,
                   *eqset_factory,
                   bc_factory,
                   cm_factory,
                   cm_factory,
                   closure_models_pl,
                   user_data_pl,false,"");

    // setup a response library to write to the mesh
    /////////////////////////////////////////////////////////////

    RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
        = buildSTKIOResponseLibrary(physicsBlocks,linObjFactory,wkstContainer,dofManager,cm_factory,mesh,
                                    closure_models_pl,user_data_pl);


    // Allocate vectors and matrix for linear solve
    /////////////////////////////////////////////////////////////
    RCP<Thyra::VectorBase<double> > solution_vec = Thyra::createMember(physics->get_x_space());
    Thyra::assign(solution_vec.ptr(),0.0); // some random initializization

    RCP<Thyra::VectorBase<double> > residual = Thyra::createMember(physics->get_f_space());
    RCP<Thyra::LinearOpWithSolveBase<double> > jacobian = physics->create_W();

    // do the assembly, this is where the evaluators are called and the graph is execueted.
    /////////////////////////////////////////////////////////////

    {
      Thyra::ModelEvaluatorBase::InArgs<double> inArgs = physics->createInArgs();
      inArgs.set_x(solution_vec);

      Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = physics->createOutArgs();
      outArgs.set_f(residual);
      outArgs.set_W(jacobian);

      // construct the residual and jacobian
      physics->evalModel(inArgs,outArgs);
    }

    // do a linear solve
    /////////////////////////////////////////////////////////////

    jacobian->solve(Thyra::NOTRANS,*residual,solution_vec.ptr());
    Thyra::scale(-1.0,solution_vec.ptr());

    // write to an exodus file
    /////////////////////////////////////////////////////////////
    writeToExodus(solution_vec,*physics,*stkIOResponseLibrary,*mesh);

  }
  catch (std::exception& e) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << e.what() << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (std::string& msg) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << msg << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }
  catch (...) {
    *out << "*********** Caught Exception: Begin Error Report ***********" << std::endl;
    *out << "Caught UNKOWN exception" << std::endl;
    *out << "************ Caught Exception: End Error Report ************" << std::endl;
    status = -1;
  }

  // Teuchos::TimeMonitor::summarize(*out,false,true,false);

  if (status == 0)
    *out << "panzer::MainDriver run completed." << std::endl;

  PHX::FinalizeKokkosDevice();

  return status;
}

Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                          const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
                          const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
                          const Teuchos::RCP<panzer::GlobalIndexer> & globalIndexer,
                          const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                          const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                          const Teuchos::ParameterList & closure_model_pl,
                          const Teuchos::ParameterList & user_data)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
      = rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,globalIndexer,linObjFactory));

  // get a vector of all the element blocks
  std::vector<std::string> eBlocks;
  mesh->getElementBlockNames(eBlocks);

  panzer_stk::RespFactorySolnWriter_Builder builder;
  builder.mesh = mesh;

  stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);

  std::map<std::string,std::vector<std::string> > nodalFields,cellFields;

  // this automatically adds in the nodal fields
  panzer_stk::IOClosureModelFactory_TemplateBuilder<panzer::Traits> io_cm_builder(cm_factory,mesh,
                                                                                          nodalFields,
                                                                                          cellFields);
  panzer::ClosureModelFactory_TemplateManager<panzer::Traits> io_cm_factory;
  io_cm_factory.buildObjects(io_cm_builder);

  stkIOResponseLibrary->buildResponseEvaluators(physicsBlocks,
                                    io_cm_factory,
                                    closure_model_pl,
                                    user_data);

  return stkIOResponseLibrary;
}

void writeToExodus(const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
                   const panzer::ModelEvaluator<double> & model,
                   panzer::ResponseLibrary<panzer::Traits> & stkIOResponseLibrary,
                   panzer_stk::STK_Interface & mesh)
{
  // fill STK mesh objects
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = model.createInArgs();
  inArgs.set_x(x);

  panzer::AssemblyEngineInArgs respInput;
  model.setupAssemblyInArgs(inArgs,respInput);

  stkIOResponseLibrary.addResponsesToInArgs<panzer::Traits::Residual>(respInput);
  stkIOResponseLibrary.evaluate<panzer::Traits::Residual>(respInput);

  mesh.writeToExodus("output.exo");
}


