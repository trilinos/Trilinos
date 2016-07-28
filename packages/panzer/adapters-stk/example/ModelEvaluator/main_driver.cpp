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

#include "Phalanx_KokkosUtilities.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
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

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include <Ioss_SerializeIO.h>

#include <string>
#include <iostream>

Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >
buildSTKIOResponseLibrary(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                          const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linObjFactory,
                          const Teuchos::RCP<panzer::WorksetContainer> & wkstContainer,
                          const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & globalIndexer,
                          const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                          const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                          const Teuchos::ParameterList & closure_model_pl,
                          const Teuchos::ParameterList & user_data);

void writeToExodus(double time_stamp,
                   const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
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
    std::string input_file_name = "user_app.xml";
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
    Teuchos::ParameterList & nonlinsolver_pl        = input_params->sublist("Nonlinear Solver");

    user_data_pl.set<RCP<const Teuchos::Comm<int> > >("Comm", comm);

    RCP<panzer::GlobalData> globalData = panzer::createGlobalData();
    RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    user_app::BCFactory bc_factory; 

    user_app::MyModelFactory_TemplateBuilder cm_builder;
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
 
      mesh->setupTransientExodusFile("output.exo");
    }

    // build worksets
    //////////////////////////////////////////////////////////////
    
    // build WorksetContainer
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory 
       = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer             // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

    // build DOF Manager
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const Teuchos::RCP<panzer::ConnManager<int,int> > 
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));

    // build the state dof manager and LOF
    RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager;
    RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    {
      panzer::DOFManagerFactory<int,int> globalIndexerFactory;
      dofManager = globalIndexerFactory.buildUniqueGlobalIndexer(Teuchos::opaqueWrapper(MPI_COMM_WORLD),physicsBlocks,conn_manager);
      linObjFactory = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(comm,dofManager));
    }

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

  
    // setup the nonlinear solver
    /////////////////////////////////////////////////////////////
    
    RCP<Thyra::NOXNonlinearSolver> noxSolver_obj = rcp(new Thyra::NOXNonlinearSolver);
    noxSolver_obj->setParameterList(rcp(new Teuchos::ParameterList(nonlinsolver_pl)));

    // do a nonlinear solve
    /////////////////////////////////////////////////////////////
    
    RCP<Thyra::VectorBase<double> > solution_vec = Thyra::createMember(physics->get_x_space());
    Thyra::assign(solution_vec.ptr(),0.0);

    {
      // set the model to use and the default parameters
      noxSolver_obj->setModel(physics);
      noxSolver_obj->setBasePoint(physics->createInArgs());

      Thyra::SolveCriteria<double> solve_criteria; // this object is ignored
      Thyra::assign(solution_vec.ptr(),0.0);
      const Thyra::SolveStatus<double> solve_status = noxSolver_obj->solve(&*solution_vec,&solve_criteria,NULL);

      TEUCHOS_TEST_FOR_EXCEPTION(
        solve_status.solveStatus != Thyra::SOLVE_STATUS_CONVERGED,
        std::runtime_error,
        "Nonlinear solver failed to converge");
    }

    // write to an exodus file
    /////////////////////////////////////////////////////////////
    writeToExodus(0,solution_vec,*physics,*stkIOResponseLibrary,*mesh);
     
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
                          const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & globalIndexer,
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

void writeToExodus(double time_stamp,
                   const Teuchos::RCP<const Thyra::VectorBase<double> > & x,
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

  mesh.writeToExodus(time_stamp);
}


