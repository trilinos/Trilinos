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

#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_DOFManagerFEI.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_ExplicitModelEvaluator.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_PamgenReaderFactory.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_MultiBlockMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Panzer_STK_RythmosObserverFactory.hpp"
#include "Panzer_STK_ParameterListCallback.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

#include <vector>
#include <iostream>
#include <fstream>

// Piro solver objects
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"
#include "Piro_LOCASolver.hpp"
#include "Piro_RythmosSolver.hpp"

#include "Epetra_MpiComm.h"

#include "EpetraExt_VectorOut.h"

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_CrsMatrix.hpp"

#ifdef HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

#ifdef HAVE_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueluTpetraHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"
#endif

#ifdef HAVE_IFPACK2
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

namespace panzer_stk_classic {

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
  {
    paramList->validateParametersAndSetDefaults(*this->getValidParameters());

    // add in some addtional defaults that are hard to validate externally (this is because of the "disableRecursiveValidation" calls)

    if(!paramList->sublist("Initial Conditions").isType<bool>("Zero Initial Conditions"))
      paramList->sublist("Initial Conditions").set<bool>("Zero Initial Conditions",false);
        
    paramList->sublist("Initial Conditions").sublist("Vector File").validateParametersAndSetDefaults(
      getValidParameters()->sublist("Initial Conditions").sublist("Vector File"));

    this->setMyParamList(paramList);
  }

  template<typename ScalarT>
  Teuchos::RCP<const Teuchos::ParameterList> ModelEvaluatorFactory<ScalarT>::getValidParameters() const
  {
    static Teuchos::RCP<const Teuchos::ParameterList> validPL;
    if (is_null(validPL)) {
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList());

      pl->sublist("Physics Blocks").disableRecursiveValidation();
      pl->sublist("Closure Models").disableRecursiveValidation();
      pl->sublist("Boundary Conditions").disableRecursiveValidation();
      pl->sublist("Solution Control").disableRecursiveValidation();
      pl->set<bool>("Use Discrete Adjoint",false);

      pl->sublist("Mesh").disableRecursiveValidation();

      pl->sublist("Initial Conditions").set<bool>("Zero Initial Conditions",false);
      pl->sublist("Initial Conditions").sublist("Transient Parameters").disableRecursiveValidation();
      pl->sublist("Initial Conditions").sublist("Vector File");
      pl->sublist("Initial Conditions").sublist("Vector File").set("File Name","");
      pl->sublist("Initial Conditions").sublist("Vector File").set<bool>("Enabled",false);
      pl->sublist("Initial Conditions").disableRecursiveValidation();

      pl->sublist("Output").set("File Name","panzer.exo");
      pl->sublist("Output").set("Write to Exodus",true);
      pl->sublist("Output").sublist("Cell Average Quantities").disableRecursiveValidation();
      pl->sublist("Output").sublist("Cell Quantities").disableRecursiveValidation();
      pl->sublist("Output").sublist("Cell Average Vectors").disableRecursiveValidation();
      pl->sublist("Output").sublist("Nodal Quantities").disableRecursiveValidation();
      pl->sublist("Output").sublist("Allocate Nodal Quantities").disableRecursiveValidation();

      // Assembly sublist
      {
	Teuchos::ParameterList& p = pl->sublist("Assembly");
	p.set<int>("Workset Size", 1);
	p.set<int>("Default Integration Order",-1);
	p.set<std::string>("Field Order","");
	p.set<bool>("Use DOFManager FEI",false);
	p.set<bool>("Load Balance DOFs",false);
	p.set<bool>("Use Tpetra",false);
	p.set<bool>("Use Epetra ME",true);
	p.set<bool>("Lump Explicit Mass",false);
	p.set<Teuchos::RCP<const panzer::EquationSetFactory> >("Equation Set Factory", Teuchos::null);
	p.set<Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >("Closure Model Factory", Teuchos::null);
	p.set<Teuchos::RCP<const panzer::BCStrategyFactory> >("BC Factory",Teuchos::null);
        p.set<std::string>("Excluded Blocks","");
      }

      pl->sublist("Block ID to Physics ID Mapping").disableRecursiveValidation();
      pl->sublist("Options").disableRecursiveValidation();
      pl->sublist("Active Parameters").disableRecursiveValidation();
      pl->sublist("User Data").disableRecursiveValidation();
      pl->sublist("User Data").sublist("Panzer Data").disableRecursiveValidation();

      validPL = pl;
    }
    return validPL;
  }

  template<typename ScalarT>
  void  ModelEvaluatorFactory<ScalarT>::buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                                     const Teuchos::RCP<panzer::GlobalData>& global_data,
                                                     const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                                                     const panzer::BCStrategyFactory & bc_factory,
                                                     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & user_cm_factory)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->getParameterList()), std::runtime_error,
		       "ParameterList must be set before objects can be built!");

    TEUCHOS_ASSERT(nonnull(comm));
    TEUCHOS_ASSERT(nonnull(global_data));
    TEUCHOS_ASSERT(nonnull(global_data->os));
    TEUCHOS_ASSERT(nonnull(global_data->pl));

    Teuchos::FancyOStream& fout = *global_data->os;

    // for convience cast to an MPI comm
    const Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);

    // this function will need to be broken up eventually and probably
    // have parts moved back into panzer.  Just need to get something
    // running.

    Teuchos::ParameterList& p = *this->getNonconstParameterList();

    // "parse" parameter list
    Teuchos::ParameterList & mesh_params     = p.sublist("Mesh");
    Teuchos::ParameterList & assembly_params = p.sublist("Assembly");
    Teuchos::ParameterList & solncntl_params = p.sublist("Solution Control");
    Teuchos::ParameterList & output_list = p.sublist("Output");

    Teuchos::ParameterList & user_data_params = p.sublist("User Data");
    Teuchos::ParameterList & panzer_data_params = user_data_params.sublist("Panzer Data");

    // Build mesh factory and uncommitted mesh
    Teuchos::RCP<panzer_stk_classic::STK_MeshFactory> mesh_factory = this->buildSTKMeshFactory(mesh_params);
    Teuchos::RCP<panzer_stk_classic::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*(mpi_comm->getRawMpiComm()));
    m_mesh = mesh;

    m_eqset_factory = eqset_factory;

    // setup physical mappings and boundary conditions
    std::map<std::string,std::string> block_ids_to_physics_ids;
    panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, p.sublist("Block ID to Physics ID Mapping"));

    // build cell ( block id -> cell topology ) mapping
    std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
    for(std::map<std::string,std::string>::const_iterator itr=block_ids_to_physics_ids.begin();
        itr!=block_ids_to_physics_ids.end();++itr) {
       block_ids_to_cell_topo[itr->first] = mesh->getCellTopology(itr->first);
       TEUCHOS_ASSERT(block_ids_to_cell_topo[itr->first]!=Teuchos::null);
    }

    Teuchos::RCP<Teuchos::ParameterList> physics_block_plist = Teuchos::sublist(this->getMyNonconstParamList(),"Physics Blocks");

    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs, p.sublist("Boundary Conditions"));

    // extract assembly information
    std::size_t workset_size = Teuchos::as<std::size_t>(assembly_params.get<int>("Workset Size"));
    std::string field_order  = assembly_params.get<std::string>("Field Order"); // control nodal ordering of unknown
                                                                                   // global IDs in linear system
    bool use_dofmanager_fei  = assembly_params.get<bool>("Use DOFManager FEI"); // use FEI if true, otherwise use internal dof manager
    bool use_load_balance = assembly_params.get<bool>("Load Balance DOFs");
    bool useTpetra = assembly_params.get<bool>("Use Tpetra");
    bool useThyraME = !assembly_params.get<bool>("Use Epetra ME");

    // this is weird...we are accessing the solution control to determine if things are transient
    // it is backwards!
    bool is_transient  = solncntl_params.get<std::string>("Piro Solver") == "Rythmos" ? true : false;
    // for pseudo-transient, we need to enable transient solver support to get time derivatives into fill
    if (solncntl_params.get<std::string>("Piro Solver") == "NOX") {
      if (solncntl_params.sublist("NOX").get<std::string>("Nonlinear Solver") == "Pseudo-Transient")
	is_transient = true;
    }
    // for eigenvalues, we need to enable transient solver support to
    // get time derivatives into generalized eigenvale problem
    if (solncntl_params.get<std::string>("Piro Solver") == "LOCA") {
      if (solncntl_params.sublist("LOCA").sublist("Stepper").get<bool>("Compute Eigenvalues"))
	is_transient = true;
    }
    m_is_transient = is_transient;

    bool useDiscreteAdjoint = p.get<bool>("Use Discrete Adjoint");

    // build physics blocks

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                               block_ids_to_cell_topo,
			       physics_block_plist,
			       assembly_params.get<int>("Default Integration Order"),
			       workset_size,
			       eqset_factory,
			       global_data,
			       is_transient,
			       physicsBlocks);
    m_physics_blocks = physicsBlocks; // hold onto physics blocks for safe keeping

    panzer_stk_classic::IOClosureModelFactory_TemplateBuilder<panzer::Traits> io_cm_builder(user_cm_factory,mesh,output_list);
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    cm_factory.buildObjects(io_cm_builder);

    // register cell averaged scalar fields
    Teuchos::ParameterList & cellAvgQuants = output_list.sublist("Cell Average Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=cellAvgQuants.begin();
        itr!=cellAvgQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh->addCellField(tokens[i],blockId);
    }

    // register cell averaged components of vector fields 
    // just allocate space for the fields here. The actual calculation and writing 
    // are done by panzer_stk_classic::ScatterCellAvgVector.
    Teuchos::ParameterList & cellAvgVectors = output_list.sublist("Cell Average Vectors");
    for(Teuchos::ParameterList::ConstIterator itr = cellAvgVectors.begin();
        itr != cellAvgVectors.end(); ++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i = 0; i < tokens.size(); i++) {
          std::string d_mod[3] = {"X","Y","Z"};
          for(std::size_t d = 0; d < mesh->getDimension(); d++) 
              mesh->addCellField(tokens[i]+d_mod[d],blockId);  
       }   
    }

    // register cell quantities
    Teuchos::ParameterList & cellQuants = output_list.sublist("Cell Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=cellQuants.begin();
        itr!=cellQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh->addCellField(tokens[i],blockId);
    }

    // register ndoal quantities
    Teuchos::ParameterList & nodalQuants = output_list.sublist("Nodal Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=nodalQuants.begin();
        itr!=nodalQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh->addSolutionField(tokens[i],blockId);
    }

    Teuchos::ParameterList & allocNodalQuants = output_list.sublist("Allocate Nodal Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=allocNodalQuants.begin();
        itr!=allocNodalQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh->addSolutionField(tokens[i],blockId);
    }

    // finish building mesh, set required field variables and mesh bulk data
    ////////////////////////////////////////////////////////////////////////
    try {
       // this throws some exceptions, catch them as neccessary
       this->finalizeMeshConstruction(*mesh_factory,physicsBlocks,*mpi_comm,*mesh);
    } catch(const panzer_stk_classic::STK_Interface::ElementBlockException & ebexp) {
       fout << "*****************************************\n\n";
       fout << "Element block exception, could not finalize the mesh, printing block and sideset information:\n";
       fout.pushTab(3);
       mesh->printMetaData(fout);
       fout.popTab();
       fout << std::endl;

       throw ebexp;
    } catch(const panzer_stk_classic::STK_Interface::SidesetException & ssexp) {
       fout << "*****************************************\n\n";
       fout << "Sideset exception, could not finalize the mesh, printing block and sideset information:\n";
       fout.pushTab(3);
       mesh->printMetaData(fout);
       fout.popTab();
       fout << std::endl;

       throw ssexp;
    }

    mesh->print(fout);
    if(p.sublist("Output").get<bool>("Write to Exodus"))
      mesh->setupTransientExodusFile(p.sublist("Output").get<std::string>("File Name"));

    // build DOF Manager
    /////////////////////////////////////////////////////////////

    // build the connection manager
    Teuchos::RCP<panzer::ConnManagerBase<int> > conn_manager;
    if(useTpetra)
      conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<panzer::Ordinal64>(mesh));
    else
      conn_manager = Teuchos::rcp(new panzer_stk_classic::STKConnManager<int>(mesh));

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    Teuchos::RCP<panzer::UniqueGlobalIndexerBase> globalIndexer;

    std::string loadBalanceString = ""; // what is the load balancing information
    bool blockedAssembly = false;

    if(panzer::BlockedDOFManagerFactory<int,int>::requiresBlocking(field_order) && !useTpetra) {
       const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager_int
         = Teuchos::rcp_dynamic_cast<panzer::ConnManager<int,int> >(conn_manager,true);

       // use a blocked DOF manager
       blockedAssembly = true;

       panzer::BlockedDOFManagerFactory<int,int> globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(use_dofmanager_fei);

       Teuchos::RCP<panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > dofManager
         = globalIndexerFactory.buildUniqueGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager_int,field_order);
       globalIndexer = dofManager;

       Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > bloLinObjFactory
        = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpi_comm,
                                                          Teuchos::rcp_dynamic_cast<panzer::BlockedDOFManager<int,int> >(dofManager)));

       // parse any explicitly excluded pairs or blocks
       const std::string excludedBlocks = assembly_params.get<std::string>("Excluded Blocks");
       std::vector<std::string> stringPairs;
       panzer::StringTokenizer(stringPairs,excludedBlocks,";",true);
       for(std::size_t i=0;i<stringPairs.size();i++) {
          std::vector<std::string> sPair;
          std::vector<int> iPair;
          panzer::StringTokenizer(sPair,stringPairs[i],",",true);
          panzer::TokensToInts(iPair,sPair);

          TEUCHOS_TEST_FOR_EXCEPTION(iPair.size()!=2,std::logic_error,
                        "Input Error: The correct format for \"Excluded Blocks\" parameter in \"Assembly\" sub list is:\n"
                        "   <int>,<int>; <int>,<int>; ...; <int>,<int>\n"
                        "Failure on string pair " << stringPairs[i] << "!");

          bloLinObjFactory->addExcludedPair(iPair[0],iPair[1]);
       }

       linObjFactory = bloLinObjFactory;

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
    }
    else if(panzer::BlockedDOFManagerFactory<int,panzer::Ordinal64>::requiresBlocking(field_order) && useTpetra) {
       const Teuchos::RCP<panzer::ConnManager<int,panzer::Ordinal64> > conn_manager_long
         = Teuchos::rcp_dynamic_cast<panzer::ConnManager<int,panzer::Ordinal64> >(conn_manager,true);

       // use a blocked DOF manager
       blockedAssembly = true;

       TEUCHOS_ASSERT(!use_dofmanager_fei);
       panzer::BlockedDOFManagerFactory<int,panzer::Ordinal64> globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(false);

       Teuchos::RCP<panzer::UniqueGlobalIndexer<int,std::pair<int,panzer::Ordinal64> > > dofManager
         = globalIndexerFactory.buildUniqueGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager_long,field_order);
       globalIndexer = dofManager;

       Teuchos::RCP<panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64> > bloLinObjFactory
        = Teuchos::rcp(new panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64>(mpi_comm,
                                                          Teuchos::rcp_dynamic_cast<panzer::BlockedDOFManager<int,panzer::Ordinal64> >(dofManager)));

       // parse any explicitly excluded pairs or blocks
       const std::string excludedBlocks = assembly_params.get<std::string>("Excluded Blocks");
       std::vector<std::string> stringPairs;
       panzer::StringTokenizer(stringPairs,excludedBlocks,";",true);
       for(std::size_t i=0;i<stringPairs.size();i++) {
          std::vector<std::string> sPair;
          std::vector<int> iPair;
          panzer::StringTokenizer(sPair,stringPairs[i],",",true);
          panzer::TokensToInts(iPair,sPair);

          TEUCHOS_TEST_FOR_EXCEPTION(iPair.size()!=2,std::logic_error,
                        "Input Error: The correct format for \"Excluded Blocks\" parameter in \"Assembly\" sub list is:\n"
                        "   <int>,<int>; <int>,<int>; ...; <int>,<int>\n"
                        "Failure on string pair " << stringPairs[i] << "!");

          bloLinObjFactory->addExcludedPair(iPair[0],iPair[1]);
       }

       linObjFactory = bloLinObjFactory;

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
    }
    else if(useTpetra) {
       const Teuchos::RCP<panzer::ConnManager<int,panzer::Ordinal64> > conn_manager_long
         = Teuchos::rcp_dynamic_cast<panzer::ConnManager<int,panzer::Ordinal64> >(conn_manager,true);

       // use a flat DOF manager

       TEUCHOS_ASSERT(!use_dofmanager_fei);
       panzer::DOFManagerFactory<int,panzer::Ordinal64> globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(false);
       globalIndexerFactory.setUseTieBreak(use_load_balance);
       Teuchos::RCP<panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> > dofManager
         = globalIndexerFactory.buildUniqueGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager_long,field_order);
       globalIndexer = dofManager;

       TEUCHOS_ASSERT(!useDiscreteAdjoint); // safety check
       linObjFactory = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64>(mpi_comm,dofManager));

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
    }
    else {
       const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager_int
         = Teuchos::rcp_dynamic_cast<panzer::ConnManager<int,int> >(conn_manager,true);

       // use a flat DOF manager

       panzer::DOFManagerFactory<int,int> globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(use_dofmanager_fei);
       globalIndexerFactory.setUseTieBreak(use_load_balance);
       Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dofManager
         = globalIndexerFactory.buildUniqueGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager_int,field_order);
       globalIndexer = dofManager;

       linObjFactory = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,int>(mpi_comm,dofManager,useDiscreteAdjoint));

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
    }

    TEUCHOS_ASSERT(globalIndexer!=Teuchos::null);
    TEUCHOS_ASSERT(linObjFactory!=Teuchos::null);
    m_global_indexer = globalIndexer;
    m_lin_obj_factory = linObjFactory;

    // print out load balancing information
    fout << "Degree of freedom load balancing: " << loadBalanceString << std::endl;

    // build worksets
    //////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk_classic::WorksetFactory> wkstFactory
       = Teuchos::rcp(new panzer_stk_classic::WorksetFactory(mesh)); // build STK workset factory
    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,physicsBlocks,workset_size));

    // set the global indexer so the orientations are evaluated
    wkstContainer->setGlobalIndexer(globalIndexer);

    // Add mesh objects to user data to make available to user ctors
    /////////////////////////////////////////////////////////////

    panzer_data_params.set("STK Mesh", mesh);
    panzer_data_params.set("DOF Manager", globalIndexer);
    panzer_data_params.set("Linear Object Factory", linObjFactory);

    // setup field manager build
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    {
      bool write_dot_files = false;
      std::string prefix = "Panzer_AssemblyGraph_";
      write_dot_files = p.sublist("Options").get("Write Volume Assembly Graphs",write_dot_files);
      prefix = p.sublist("Options").get("Volume Assembly Graph Prefix",prefix);

      fmb = buildFieldManagerBuilder(wkstContainer,physicsBlocks,bcs,*eqset_factory,bc_factory,cm_factory,
                                     user_cm_factory,p.sublist("Closure Models"),*linObjFactory,user_data_params,
                                     write_dot_files,prefix);
    }

    // build response library
    /////////////////////////////////////////////////////////////

    m_response_library = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wkstContainer,globalIndexer,linObjFactory));

    {
       bool write_dot_files = false;
       std::string prefix = "Panzer_ResponseGraph_";
       write_dot_files = p.sublist("Options").get("Write Volume Response Graphs",write_dot_files);
       prefix = p.sublist("Options").get("Volume Response Graph Prefix",prefix);

       Teuchos::ParameterList user_data(p.sublist("User Data"));
       user_data.set<int>("Workset Size",workset_size);

    }

    // build solvers
    /////////////////////////////////////////////////////////////

    // Setup active parameters
    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    if (p.isSublist("Active Parameters")) {
      Teuchos::ParameterList& active_params = p.sublist("Active Parameters");

      int num_param_vecs = active_params.get<int>("Number of Parameter Vectors",0);
      p_names.resize(num_param_vecs);
      for (int i=0; i<num_param_vecs; i++) {
	std::stringstream ss;
	ss << "Parameter Vector " << i;
	Teuchos::ParameterList& pList = active_params.sublist(ss.str());
	int numParameters = pList.get<int>("Number");
	TEUCHOS_TEST_FOR_EXCEPTION(numParameters == 0,
				   Teuchos::Exceptions::InvalidParameter,
				   std::endl << "Error!  panzer::ModelEvaluator::ModelEvaluator():  " <<
				   "Parameter vector " << i << " has zero parameters!" << std::endl);
	p_names[i] =
	  Teuchos::rcp(new Teuchos::Array<std::string>(numParameters));
	for (int j=0; j<numParameters; j++) {
	  std::stringstream ss2;
	  ss2 << "Parameter " << j;
	  (*p_names[i])[j] = pList.get<std::string>(ss2.str());
	}
      }
    }


    // Setup solver factory
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
          buildLOWSFactory(blockedAssembly,globalIndexer,conn_manager,mesh,mpi_comm);

    // Setup physics model evaluator
    /////////////////////////////////////////////////////////////

    double t_init = 0.0;
    if(is_transient)
      t_init = this->getInitialTime(p.sublist("Initial Conditions").sublist("Transient Parameters"), *mesh);

    if(blockedAssembly || useTpetra) // override the user request
      useThyraME = true;

    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me
        = buildPhysicsModelEvaluator(useThyraME, // blockedAssembly || useTpetra, // this determines if a Thyra or Epetra ME is used
                                     fmb,
                                     m_response_library,
                                     linObjFactory,
                                     p_names,
                                     lowsFactory,
                                     global_data,
                                     is_transient,
                                     t_init);

    // Setup initial conditions
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::LinearObjContainer> loc = linObjFactory->buildLinearObjContainer();

    if(p.sublist("Initial Conditions").get<bool>("Zero Initial Conditions")) {
      // zero out the x vector
      Thyra::ModelEvaluatorBase::InArgs<double> nomValues = thyra_me->getNominalValues();

      Thyra::assign(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x()).ptr(),0.0); 
    }
    else if(!p.sublist("Initial Conditions").sublist("Vector File").get<bool>("Enabled")) {
      // read from exodus, or compute using field managers

      bool write_dot_files = false;
      std::string prefix = "Panzer_AssemblyGraph_";
      write_dot_files = p.sublist("Options").get("Write Volume Assembly Graphs",write_dot_files);
      prefix = p.sublist("Options").get("Volume Assembly Graph Prefix",prefix);

      std::map<std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > phx_ic_field_managers;
      panzer::setupInitialConditionFieldManagers(*wkstContainer,
                                                 physicsBlocks,
                                                 user_cm_factory,
                                                 p.sublist("Initial Conditions"),
                                                 *linObjFactory,
                                                 p.sublist("User Data"),
                                                 write_dot_files,
                                                 prefix,
                                                 phx_ic_field_managers);

      // set the vector to be filled
      Teuchos::RCP<panzer::ThyraObjContainer<double> > tloc = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc);
      Thyra::ModelEvaluatorBase::InArgs<double> nomValues = thyra_me->getNominalValues();
      tloc->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x()));

      panzer::evaluateInitialCondition(*wkstContainer, phx_ic_field_managers, loc, 0.0);
   }
   else {
      const std::string & vectorFile = p.sublist("Initial Conditions").sublist("Vector File").get<std::string>("File Name");
      TEUCHOS_TEST_FOR_EXCEPTION(vectorFile=="",std::runtime_error,
                                 "If \"Read From Vector File\" is true, then parameter \"Vector File\" cannot be the empty string.");
 
      // set the vector to be filled
      Teuchos::RCP<panzer::ThyraObjContainer<double> > tloc = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc);
      Thyra::ModelEvaluatorBase::InArgs<double> nomValues = thyra_me->getNominalValues();
      tloc->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x()));
      
      // read the vector
      linObjFactory->readVector(vectorFile,*loc,panzer::LinearObjContainer::X);
   }

   // Write the IC vector into the STK mesh: use response library
   //////////////////////////////////////////////////////////////////////////

   {
      Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > solnWriter
          = initializeSolnWriterResponseLibrary(wkstContainer,globalIndexer,linObjFactory,mesh);

      {
         Teuchos::ParameterList user_data(p.sublist("User Data"));
         user_data.set<int>("Workset Size",workset_size);

         finalizeSolnWriterResponseLibrary(*solnWriter,physicsBlocks,user_cm_factory,p.sublist("Closure Models"),workset_size,user_data);
      }

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = loc;
      ae_inargs.ghostedContainer_ = linObjFactory->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      linObjFactory->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      // do import
      linObjFactory->globalToGhostContainer(*ae_inargs.container_,*ae_inargs.ghostedContainer_,panzer::LinearObjContainer::X);

      // fill STK mesh objects
      solnWriter->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
      solnWriter->evaluate<panzer::Traits::Residual>(ae_inargs);

    }

    // see if field coordinates are required, if so reset the workset container
    // and set the coordinates to be associated with a field in the mesh
    useDynamicCoordinates_ = false;
    for(std::size_t p=0;p<physicsBlocks.size();p++) {
      if(physicsBlocks[p]->getCoordinateDOFs().size()>0) {
         mesh->setUseFieldCoordinates(true);
         useDynamicCoordinates_ = true;
         wkstContainer->clear(); // this serves to refresh the worksets 
                                 // and put in new coordinates
         break;
      }
    }

    m_physics_me = thyra_me;
    m_global_data = global_data;
  }

  //! build STK mesh from a mesh parameter list
  template<typename ScalarT>
  Teuchos::RCP<panzer_stk_classic::STK_MeshFactory> ModelEvaluatorFactory<ScalarT>::buildSTKMeshFactory(const Teuchos::ParameterList & mesh_params) const
  {
    Teuchos::RCP<panzer_stk_classic::STK_MeshFactory> mesh_factory;

    // first contruct the mesh factory
    if (mesh_params.get<std::string>("Source") ==  "Exodus File") {
      mesh_factory = Teuchos::rcp(new panzer_stk_classic::STK_ExodusReaderFactory());
      mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Exodus File"))));
    }
    else if (mesh_params.get<std::string>("Source") ==  "Pamgen Mesh") {
      mesh_factory = Teuchos::rcp(new panzer_stk_classic::STK_PamgenReaderFactory());
      mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Pamgen Mesh"))));
    }
    else if (mesh_params.get<std::string>("Source") ==  "Inline Mesh") {

      int dimension = mesh_params.sublist("Inline Mesh").get<int>("Mesh Dimension");
      std::string typeStr = "";
      if(mesh_params.sublist("Inline Mesh").isParameter("Type"))
         typeStr = mesh_params.sublist("Inline Mesh").get<std::string>("Type");

      if (dimension == 1) {
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::LineMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 2 && typeStr=="Tri") {
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::SquareTriMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 2) {
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::SquareQuadMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 3 && typeStr=="Tet") {
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::CubeTetMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension == 3) {
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::CubeHexMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension==4) { // not really "dimension==4" simply a flag to try this other mesh for testing
	mesh_factory = Teuchos::rcp(new panzer_stk_classic::MultiBlockMeshFactory);
	Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
	*in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
	mesh_factory->setParameterList(in_mesh);
      }
    }

    // get rebalancing parameters
    if(mesh_params.isSublist("Rebalance")) {
      const Teuchos::ParameterList & rebalance = mesh_params.sublist("Rebalance");

      // check to see if its enabled
      bool enabled = false;
      if(rebalance.isType<bool>("Enabled"))
        enabled = rebalance.get<bool>("Enabled");

      // we can also use a list description of what to load balance
      Teuchos::RCP<Teuchos::ParameterList> rebalanceCycles;
      if(enabled && rebalance.isSublist("Cycles"))
        rebalanceCycles = Teuchos::rcp(new Teuchos::ParameterList(rebalance.sublist("Cycles")));

      // setup rebalancing as neccessary
      mesh_factory->enableRebalance(enabled,rebalanceCycles);
    }

    return mesh_factory;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::finalizeMeshConstruction(const STK_MeshFactory & mesh_factory,
                                                                const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                                                const Teuchos::MpiComm<int> mpi_comm,
                                                                STK_Interface & mesh) const
  {
    // finish building mesh, set required field variables and mesh bulk data
    {
      std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
      for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
        // what is the block weight for this element block?
        double blockWeight = 0.0;

	Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
	const std::vector<panzer::StrPureBasisPair> & blockFields = pb->getProvidedDOFs();
	const std::vector<std::vector<std::string> > & coordinateDOFs = pb->getCoordinateDOFs(); 
          // these are treated specially

	// insert all fields into a set
	std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> fieldNames;
	fieldNames.insert(blockFields.begin(),blockFields.end());

        // Now we will set up the coordinate fields (make sure to remove
        // the DOF fields)
        {
          std::set<std::string> fields_to_remove;

          // add mesh coordinate fields, setup their removal from fieldNames 
          // set to prevent duplication
          for(std::size_t i=0;i<coordinateDOFs.size();i++) {
            mesh.addMeshCoordFields(pb->elementBlockID(),coordinateDOFs[i],"DISPL");
            for(std::size_t j=0;j<coordinateDOFs[i].size();j++)
              fields_to_remove.insert(coordinateDOFs[i][j]);
          }

          // remove the already added coordinate fields
	  std::set<std::string>::const_iterator rmItr;
  	  for (rmItr=fields_to_remove.begin();rmItr!=fields_to_remove.end();++rmItr) 
            fieldNames.erase(fieldNames.find(panzer::StrPureBasisPair(*rmItr,Teuchos::null)));
        }

	// add basis to DOF manager: block specific
	std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator fieldItr;
	for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {

          if(fieldItr->second->isScalarBasis() &&  
             fieldItr->second->getElementSpace()==panzer::PureBasis::CONST) {
             mesh.addCellField(fieldItr->first,pb->elementBlockID());
          }
          else if(fieldItr->second->isScalarBasis()) {
             mesh.addSolutionField(fieldItr->first,pb->elementBlockID());
          }
          else if(fieldItr->second->isVectorBasis()) {
            std::string d_mod[3] = {"X","Y","Z"};
            for(int d=0;d<fieldItr->second->dimension();d++) 
              mesh.addCellField(fieldItr->first+d_mod[d],pb->elementBlockID());
          }
          else { TEUCHOS_ASSERT(false); }

          blockWeight += double(fieldItr->second->cardinality());
        }

        // set the compute block weight (this is the sum of the cardinality of all basis
        // functions on this block
        mesh.setBlockWeight(pb->elementBlockID(),blockWeight);
      }

      mesh_factory.completeMeshConstruction(mesh,*(mpi_comm.getRawMpiComm()));
    }
  }


  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluatorFactory<ScalarT>::getPhysicsModelEvaluator()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(m_physics_me), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    return  m_physics_me;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::setNOXObserverFactory(const Teuchos::RCP<const panzer_stk_classic::NOXObserverFactory>& nox_observer_factory)
  {
    m_nox_observer_factory = nox_observer_factory;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::setRythmosObserverFactory(const Teuchos::RCP<const panzer_stk_classic::RythmosObserverFactory>& rythmos_observer_factory)
  {
    m_rythmos_observer_factory = rythmos_observer_factory;
  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluatorFactory<ScalarT>::getResponseOnlyModelEvaluator()
  {
    if(m_rome_me==Teuchos::null)
      m_rome_me = buildResponseOnlyModelEvaluator(m_physics_me,m_global_data);

    return m_rome_me;
  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluatorFactory<ScalarT>::
  buildResponseOnlyModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & thyra_me,
 		                  const Teuchos::RCP<panzer::GlobalData>& global_data,
                                  const Teuchos::RCP<Piro::RythmosSolver<ScalarT> > rythmosSolver,
                                  const Teuchos::Ptr<const panzer_stk_classic::NOXObserverFactory> & in_nox_observer_factory,
                                  const Teuchos::Ptr<const panzer_stk_classic::RythmosObserverFactory> & in_rythmos_observer_factory
                                  )
  {
    using Teuchos::is_null;
    using Teuchos::Ptr;

    TEUCHOS_TEST_FOR_EXCEPTION(is_null(m_lin_obj_factory), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(m_global_indexer), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(m_mesh), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");
    Teuchos::Ptr<const panzer_stk_classic::NOXObserverFactory> nox_observer_factory 
        = is_null(in_nox_observer_factory) ? m_nox_observer_factory.ptr() : in_nox_observer_factory;
    Teuchos::Ptr<const panzer_stk_classic::RythmosObserverFactory> rythmos_observer_factory
        = is_null(in_rythmos_observer_factory) ? m_rythmos_observer_factory.ptr() : in_rythmos_observer_factory;

    Teuchos::ParameterList& p = *this->getNonconstParameterList();
    Teuchos::ParameterList & solncntl_params = p.sublist("Solution Control");
    Teuchos::RCP<Teuchos::ParameterList> piro_params = Teuchos::rcp(new Teuchos::ParameterList(solncntl_params));
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > piro;

    std::string solver = solncntl_params.get<std::string>("Piro Solver");
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me_db
       = Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double> >(thyra_me);
    if ( (solver=="NOX") || (solver == "LOCA") ) {

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(nox_observer_factory), std::runtime_error,
				 "No NOX obersver built!  Please call setNOXObserverFactory() member function if you plan to use a NOX solver.");

      Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = nox_observer_factory->buildNOXObserver(m_mesh,m_global_indexer,m_lin_obj_factory);
      piro_params->sublist("NOX").sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);

      if (solver=="NOX")
	piro = Teuchos::rcp(new Piro::NOXSolver<double>(piro_params,
							Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double> >(thyra_me_db)));
      else if (solver == "LOCA")
	piro = Teuchos::rcp(new Piro::LOCASolver<double>(piro_params,
							 Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double> >(thyra_me_db),
							 Teuchos::null));
      TEUCHOS_ASSERT(nonnull(piro));

      // override printing to use panzer ostream
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Output Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Error Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<int>("Output Processor",global_data->os->getOutputToRootOnly());
    }
    else if (solver=="Rythmos") {

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(rythmos_observer_factory), std::runtime_error,
				 "No NOX obersver built!  Please call setrythmosObserverFactory() member function if you plan to use a Rythmos solver.");

      // install the nox observer
      if(rythmos_observer_factory->useNOXObserver()) {
        Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = nox_observer_factory->buildNOXObserver(m_mesh,m_global_indexer,m_lin_obj_factory);
        piro_params->sublist("NOX").sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);
      }

      // override printing to use panzer ostream
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Output Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Error Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<int>("Output Processor",global_data->os->getOutputToRootOnly());

      // use the user specfied rythmos solver if they pass one in
      Teuchos::RCP<Piro::RythmosSolver<double> > piro_rythmos;
      if(rythmosSolver==Teuchos::null)
        piro_rythmos = Teuchos::rcp(new Piro::RythmosSolver<double>());
      else
        piro_rythmos = rythmosSolver;

      // if you are using explicit RK, make sure to wrap the ME in an explicit model evaluator decorator
      Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > rythmos_me = thyra_me;
      const std::string stepper_type = piro_params->sublist("Rythmos").get<std::string>("Stepper Type");
      if(stepper_type=="Explicit RK" || stepper_type=="Forward Euler") {
        const Teuchos::ParameterList & assembly_params = p.sublist("Assembly");
        bool lumpExplicitMass = assembly_params.get<bool>("Lump Explicit Mass");
        rythmos_me = Teuchos::rcp(new panzer::ExplicitModelEvaluator<ScalarT>(thyra_me,!useDynamicCoordinates_,lumpExplicitMass)); 
      }

      piro_rythmos->initialize(piro_params, rythmos_me, rythmos_observer_factory->buildRythmosObserver(m_mesh,m_global_indexer,m_lin_obj_factory));

      piro = piro_rythmos;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error: Unknown Piro Solver : " << solver);
    }
    return piro;
  }

  template<typename ScalarT>
  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > ModelEvaluatorFactory<ScalarT>::getResponseLibrary()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(m_response_library), std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");

    return m_response_library;
  }

  template<typename ScalarT>
  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & ModelEvaluatorFactory<ScalarT>::getPhysicsBlocks() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(m_physics_blocks.size()==0, std::runtime_error,
		       "Objects are not built yet!  Please call buildObjects() member function.");

    return m_physics_blocks;
  }

  template<typename ScalarT>
  Teuchos::RCP<panzer::FieldManagerBuilder> 
  ModelEvaluatorFactory<ScalarT>::
  buildFieldManagerBuilder(const Teuchos::RCP<panzer::WorksetContainer> & wc,
                           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                           const std::vector<panzer::BC> & bcs,
                           const panzer::EquationSetFactory & eqset_factory,
                           const panzer::BCStrategyFactory& bc_factory,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& volume_cm_factory,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& bc_cm_factory,
                           const Teuchos::ParameterList& closure_models,
                           const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
                           const Teuchos::ParameterList& user_data,
                           bool writeGraph,const std::string & graphPrefix) const
  {
    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
    fmb->setWorksetContainer(wc);
    fmb->setupVolumeFieldManagers(physicsBlocks,volume_cm_factory,closure_models,lo_factory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,bc_cm_factory,bc_factory,closure_models,lo_factory,user_data);

    // Print Phalanx DAGs
    if (writeGraph){
      fmb->writeVolumeGraphvizDependencyFiles(graphPrefix, physicsBlocks);
      fmb->writeBCGraphvizDependencyFiles(graphPrefix+"BC_");
    }
    return fmb;
  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<double> > 
  ModelEvaluatorFactory<ScalarT>::
  cloneWithNewPhysicsBlocks(const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                            const Teuchos::RCP<Teuchos::ParameterList> & physics_block_plist,
                            const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                            const panzer::BCStrategyFactory & bc_factory,
                            const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & user_cm_factory,
                            bool is_transient,bool is_explicit,
                            const Teuchos::Ptr<const Teuchos::ParameterList> & bc_list) const
  {
    typedef panzer::ModelEvaluator<ScalarT> PanzerME;

    const Teuchos::ParameterList& p = *this->getParameterList();

    // build PhysicsBlocks
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
      const Teuchos::ParameterList & assembly_params = p.sublist("Assembly");
 
      // setup physical mappings and boundary conditions
      std::map<std::string,std::string> block_ids_to_physics_ids;
      panzer::buildBlockIdToPhysicsIdMap(block_ids_to_physics_ids, p.sublist("Block ID to Physics ID Mapping"));

      // build cell ( block id -> cell topology ) mapping
      std::map<std::string,Teuchos::RCP<const shards::CellTopology> > block_ids_to_cell_topo;
      for(std::map<std::string,std::string>::const_iterator itr=block_ids_to_physics_ids.begin();
          itr!=block_ids_to_physics_ids.end();++itr) {
         block_ids_to_cell_topo[itr->first] = m_mesh->getCellTopology(itr->first);
         TEUCHOS_ASSERT(block_ids_to_cell_topo[itr->first]!=Teuchos::null);
      }

      std::size_t workset_size = Teuchos::as<std::size_t>(assembly_params.get<int>("Workset Size"));

      panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                 block_ids_to_cell_topo,
			         physics_block_plist,
			         assembly_params.get<int>("Default Integration Order"),
			         workset_size,
	  		         eqset_factory,
			         m_global_data,
			         is_transient,
			         physicsBlocks);
    }

    // build FMB
    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    {
      const Teuchos::ParameterList & user_data_params = p.sublist("User Data");

      bool write_dot_files = false;
      std::string prefix = "Cloned_";

      std::vector<panzer::BC> bcs;
      if(bc_list==Teuchos::null) {
        panzer::buildBCs(bcs, p.sublist("Boundary Conditions"));
      }
      else { 
        panzer::buildBCs(bcs, *bc_list);
      }
      
      fmb = buildFieldManagerBuilder(Teuchos::rcp_const_cast<panzer::WorksetContainer>(m_response_library->getWorksetContainer()),
                                     physicsBlocks,
                                     bcs,
                                     *eqset_factory,
                                     bc_factory,
                                     user_cm_factory,
                                     user_cm_factory,
                                     p.sublist("Closure Models"),
                                     *m_lin_obj_factory,
                                     user_data_params,
                                     write_dot_files,prefix);
    }

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > response_library 
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(m_response_library->getWorksetContainer(),
                                                                   m_response_library->getGlobalIndexer(),
                                                                   m_response_library->getLinearObjFactory()));

    // using the FMB, build the model evaluator
    {
      // get nominal input values, make sure they match with internal me
      Thyra::ModelEvaluatorBase::InArgs<ScalarT> nomVals = m_physics_me->getNominalValues();
  
      // determine if this is a Epetra or Thyra ME
      Teuchos::RCP<Thyra::EpetraModelEvaluator> ep_thyra_me = Teuchos::rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(m_physics_me);
      Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(m_physics_me);
      bool useThyra = true;
      if(ep_thyra_me!=Teuchos::null)
        useThyra = false;
  
      // get parameter names
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(m_physics_me->Np());
      for(std::size_t i=0;i<p_names.size();i++) 
        p_names[i] = Teuchos::rcp(new Teuchos::Array<std::string>(*m_physics_me->get_p_names(i)));
  
      Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me
          = buildPhysicsModelEvaluator(useThyra,
                                       fmb,
                                       response_library,
                                       m_lin_obj_factory,
                                       p_names,
                                       solverFactory,
                                       m_global_data,
                                       is_transient,
                                       nomVals.get_t());
  
      // set the nominal values...does this work???
      thyra_me->getNominalValues() = nomVals;
  
      // build an explicit model evaluator
      if(is_explicit) {
        const Teuchos::ParameterList & assembly_params = p.sublist("Assembly");
        bool lumpExplicitMass = assembly_params.get<bool>("Lump Explicit Mass");
        thyra_me = Teuchos::rcp(new panzer::ExplicitModelEvaluator<ScalarT>(thyra_me,!useDynamicCoordinates_,lumpExplicitMass)); 
      }
  
      return thyra_me;
    }
  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > 
  ModelEvaluatorFactory<ScalarT>::
  buildPhysicsModelEvaluator(bool buildThyraME,
                             const Teuchos::RCP<panzer::FieldManagerBuilder> & fmb,
                             const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & rLibrary,
      	                     const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                             const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > & p_names,
                             const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                             const Teuchos::RCP<panzer::GlobalData> & global_data,
                             bool is_transient,double t_init) const
  {
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me;
    if(!buildThyraME) {
      Teuchos::RCP<panzer::ModelEvaluator_Epetra> ep_me 
          = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,rLibrary,lof, p_names, global_data, is_transient));

      if (is_transient)
        ep_me->set_t_init(t_init);
      
      // Build Thyra Model Evaluator
      thyra_me = Thyra::epetraModelEvaluator(ep_me,solverFactory);
    }
    else {
      thyra_me = Teuchos::rcp(new panzer::ModelEvaluator<double>
                  (fmb,rLibrary,lof,p_names,solverFactory,global_data,is_transient,t_init));
    }

    return thyra_me;
  }

  template<typename ScalarT>
  bool ModelEvaluatorFactory<ScalarT>::determineCoordinateField(const panzer::UniqueGlobalIndexerBase & globalIndexer,std::string & fieldName) const
  {
    std::vector<std::string> elementBlocks;
    globalIndexer.getElementBlockIds(elementBlocks);

    // grab fields for first block
    std::set<int> runningFields;
    {
      const std::vector<int> & fields = globalIndexer.getBlockFieldNumbers(elementBlocks[0]);
      runningFields.insert(fields.begin(),fields.end());
    }

    // grab fields for first block
    for(std::size_t i=1;i<elementBlocks.size();i++) {
      const std::vector<int> & fields = globalIndexer.getBlockFieldNumbers(elementBlocks[i]);

      std::set<int> currentFields(runningFields);
      runningFields.clear();
      std::set_intersection(fields.begin(),fields.end(),
                            currentFields.begin(),currentFields.end(),
                            std::inserter(runningFields,runningFields.begin()));
    }

    if(runningFields.size()<1)
      return false;

    fieldName = globalIndexer.getFieldString(*runningFields.begin());
    return true;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::fillFieldPatternMap(const panzer::UniqueGlobalIndexerBase & globalIndexer,
                                                                  const std::string & fieldName,
                                                                  std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const
  {
    using Teuchos::Ptr;
    using Teuchos::ptrFromRef;
    using Teuchos::ptr_dynamic_cast;
    using panzer::DOFManager;
    using panzer::DOFManagerFEI;

    // first standard dof manager
    {
      Ptr<const DOFManager<int,int> > dofManager = ptr_dynamic_cast<const DOFManager<int,int> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
    {
      Ptr<const DOFManager<int,panzer::Ordinal64> > dofManager = ptr_dynamic_cast<const DOFManager<int,panzer::Ordinal64> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }

    // now FEI dof manager
    {
      Ptr<const DOFManagerFEI<int,int> > dofManager = ptr_dynamic_cast<const DOFManagerFEI<int,int> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
    {
      Ptr<const DOFManagerFEI<int,panzer::Ordinal64> > dofManager = ptr_dynamic_cast<const DOFManagerFEI<int,panzer::Ordinal64> >(ptrFromRef(globalIndexer));

      if(dofManager!=Teuchos::null) {
        fillFieldPatternMap(*dofManager,fieldName,fieldPatterns);
        return;
      }
    }
  }

  template<typename ScalarT>
  template<typename GO>
  void ModelEvaluatorFactory<ScalarT>::fillFieldPatternMap(const panzer::DOFManagerFEI<int,GO> & globalIndexer,
                                                                  const std::string & fieldName,
                                                                  std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const
  {
     std::vector<std::string> elementBlocks;
     globalIndexer.getElementBlockIds(elementBlocks);

     for(std::size_t e=0;e<elementBlocks.size();e++) {
        std::string blockId = elementBlocks[e];

        if(globalIndexer.fieldInBlock(fieldName,blockId))
           fieldPatterns[blockId] =
              Teuchos::rcp_dynamic_cast<const panzer::IntrepidFieldPattern>(globalIndexer.getFieldPattern(blockId,fieldName),true);
     }
  }

  template<typename ScalarT>
  template<typename GO>
  void ModelEvaluatorFactory<ScalarT>::fillFieldPatternMap(const panzer::DOFManager<int,GO> & globalIndexer,
                                                                  const std::string & fieldName,
                                                                  std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const
  {
     std::vector<std::string> elementBlocks;
     globalIndexer.getElementBlockIds(elementBlocks);

     for(std::size_t e=0;e<elementBlocks.size();e++) {
        std::string blockId = elementBlocks[e];

        if(globalIndexer.fieldInBlock(fieldName,blockId))
           fieldPatterns[blockId] =
              Teuchos::rcp_dynamic_cast<const panzer::IntrepidFieldPattern>(globalIndexer.getFieldPattern(blockId,fieldName),true);
     }
  }

  template<typename ScalarT>
  double ModelEvaluatorFactory<ScalarT>::
  getInitialTime(Teuchos::ParameterList& p,
		 const panzer_stk_classic::STK_Interface & mesh) const
  {
    Teuchos::ParameterList validPL;
    {
      Teuchos::setStringToIntegralParameter<int>(
      "Start Time Type",
      "From Input File",
      "Set the start time",
      Teuchos::tuple<std::string>("From Input File","From Exodus File"),
      &validPL
      );

      validPL.set<double>("Start Time",0.0);
    }

    p.validateParametersAndSetDefaults(validPL);

    std::string t_init_type = p.get<std::string>("Start Time Type");
    double t_init = 10.0;

    if (t_init_type == "From Input File")
      t_init = p.get<double>("Start Time");

    if (t_init_type == "From Exodus File")
      t_init = mesh.getInitialStateTime();

    return t_init;
  }

  // Setup STK response library for writing out the solution fields
  ////////////////////////////////////////////////////////////////////////
  template<typename ScalarT>
  Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > ModelEvaluatorFactory<ScalarT>::
  initializeSolnWriterResponseLibrary(const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                      const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & ugi,
                                      const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                                      const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh) const
  {
     Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wc,ugi,lof));

     std::vector<std::string> eBlocks;
     mesh->getElementBlockNames(eBlocks);

     panzer_stk_classic::RespFactorySolnWriter_Builder builder;
     builder.mesh = mesh;
     stkIOResponseLibrary->addResponse("Main Field Output",eBlocks,builder);

     return stkIOResponseLibrary;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  finalizeSolnWriterResponseLibrary(panzer::ResponseLibrary<panzer::Traits> & rl,
                                    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                    const Teuchos::ParameterList & closure_models,
                                    int workset_size, Teuchos::ParameterList & user_data) const
  {
     user_data.set<int>("Workset Size",workset_size);
     rl.buildResponseEvaluators(physicsBlocks, cm_factory, closure_models, user_data);
  }

  template<typename ScalarT>
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > ModelEvaluatorFactory<ScalarT>::
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer::ConnManagerBase<int> > & conn_manager,
                   const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm)
  {
    RCP<panzer_stk_classic::STKConnManager<panzer::Ordinal64> > long_conn = Teuchos::rcp_dynamic_cast<panzer_stk_classic::STKConnManager<panzer::Ordinal64> >(conn_manager);
    if(long_conn!=Teuchos::null)
      return buildLOWSFactory(blockedAssembly,globalIndexer,long_conn,mesh,mpi_comm);

    RCP<panzer_stk_classic::STKConnManager<int> > int_conn = Teuchos::rcp_dynamic_cast<panzer_stk_classic::STKConnManager<int> >(conn_manager);
    if(int_conn!=Teuchos::null)
      return buildLOWSFactory(blockedAssembly,globalIndexer,int_conn,mesh,mpi_comm);

    // should never reach this
    TEUCHOS_ASSERT(false);
    return Teuchos::null;
  }

  template<typename ScalarT>
  template<typename GO>
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > ModelEvaluatorFactory<ScalarT>::
  buildLOWSFactory(bool blockedAssembly,
                   const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                   const Teuchos::RCP<panzer_stk_classic::STKConnManager<GO> > & stkConn_manager,
                   const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm)
  {
    Teuchos::ParameterList& p = *this->getNonconstParameterList();
    Teuchos::ParameterList & solncntl_params = p.sublist("Solution Control");

    // Build stratimikos solver (note that this is a hard coded path to linear solver options in nox list!)
    Teuchos::RCP<Teuchos::ParameterList> strat_params = Teuchos::rcp(new Teuchos::ParameterList);
    {
      *strat_params = solncntl_params.sublist("NOX").sublist("Direction").
	sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos");
    }

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    #ifdef HAVE_TEKO
    if(!blockedAssembly) {

       std::string fieldName;

       // try to set request handler from member variable
       Teuchos::RCP<Teko::RequestHandler> reqHandler = m_req_handler;
       if(m_req_handler==Teuchos::null) {
          reqHandler = Teuchos::rcp(new Teko::RequestHandler);
          m_req_handler = reqHandler;
       }

       // add in the coordinate parameter list callback handler
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > fieldPatterns;
          fillFieldPatternMap(*globalIndexer,fieldName,fieldPatterns);

          Teuchos::RCP<panzer_stk_classic::ParameterListCallback<int,GO> > callback = Teuchos::rcp(new
                panzer_stk_classic::ParameterListCallback<int,GO>(fieldName,fieldPatterns,stkConn_manager,
                Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,GO> >(globalIndexer)));
          reqHandler->addRequestCallback(callback);

          bool writeCoordinates = p.sublist("Options").get("Write Coordinates",false);
          if(writeCoordinates) {
             // force parameterlistcallback to build coordinates
             callback->preRequest(Teko::RequestMesg(Teuchos::rcp(new Teuchos::ParameterList())));

             // extract coordinate vectors
             const std::vector<double> & xcoords = callback->getXCoordsVector();
             const std::vector<double> & ycoords = callback->getYCoordsVector();
             const std::vector<double> & zcoords = callback->getZCoordsVector();

             // use epetra to write coordinates to matrix market files
             Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm()); // this is OK access to RawMpiComm becase its declared on the stack?
                                                                 // and all users of this object are on the stack (within scope of mpi_comm
             Epetra_Map map(-1,xcoords.size(),0,ep_comm);

             Teuchos::RCP<Epetra_Vector> vec;
             switch(mesh->getDimension()) {
             case 3:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("zcoords.mm",*vec);
             case 2:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
                EpetraExt::VectorToMatrixMarketFile("ycoords.mm",*vec);
             case 1:
                vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
                EpetraExt::VectorToMatrixMarketFile("xcoords.mm",*vec);
                break;
             default:
                TEUCHOS_ASSERT(false);
             }
          }

          #ifdef HAVE_MUELU
          if(Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> >(globalIndexer)!=Teuchos::null) {
             if(!writeCoordinates)
                callback->preRequest(Teko::RequestMesg(Teuchos::rcp(new Teuchos::ParameterList())));

             typedef Tpetra::Map<int,panzer::Ordinal64,KokkosClassic::DefaultNode::DefaultNodeType> Map;
             typedef Tpetra::MultiVector<double,int,panzer::Ordinal64,KokkosClassic::DefaultNode::DefaultNodeType> MV;

             // extract coordinate vectors and modify strat_params to include coordinate vectors
             unsigned dim = mesh->getDimension();
             Teuchos::RCP<MV> coords;
             for(unsigned d=0;d<dim;d++) {
               const std::vector<double> & coord = callback->getCoordsVector(d);

               // no coords vector has been build yet, build one
               if(coords==Teuchos::null) {
                 if(globalIndexer->getNumFields()==1) {
                   Teuchos::RCP<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> > ugi
                       = Teuchos::rcp_dynamic_cast<const panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> >(globalIndexer);
                   std::vector<panzer::Ordinal64> ownedIndices;
                   ugi->getOwnedIndices(ownedIndices);
                   Teuchos::RCP<Map> coords_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::Ordinal64>::invalid(),ownedIndices,0,mpi_comm));
                   coords = Teuchos::rcp(new MV(coords_map,dim));
                 }
                 else {
                   Teuchos::RCP<Map> coords_map = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::Ordinal64>::invalid(),coord.size(),0,mpi_comm));
                   coords = Teuchos::rcp(new MV(coords_map,dim));
                 }
               }

               // sanity check the size
               TEUCHOS_ASSERT(coords->getLocalLength()==coord.size());

               // fill appropriate coords vector
               Teuchos::ArrayRCP<double> dest = coords->getDataNonConst(d);
               for(std::size_t i=0;i<coord.size();i++)
                 dest[i] = coord[i];
             }

             // inject coordinates into parameter list
             Teuchos::ParameterList & muelu_params = strat_params->sublist("Preconditioner Types").sublist("MueLu-Tpetra");
             muelu_params.set<Teuchos::RCP<MV> >("Coordinates",coords);
          }
          #endif
       }
       // else write_out_the_mesg("Warning: No unique field determines the coordinates, coordinates unavailable!")

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler);
    }
    else {
       // try to set request handler from member variable
       Teuchos::RCP<Teko::RequestHandler> reqHandler = m_req_handler;
       if(m_req_handler==Teuchos::null) {
          reqHandler = Teuchos::rcp(new Teko::RequestHandler);
          m_req_handler = reqHandler;
       }

       std::string fieldName;
       if(determineCoordinateField(*globalIndexer,fieldName)) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);
          Teuchos::RCP<panzer_stk_classic::ParameterListCallbackBlocked<int,GO> > callback =
                Teuchos::rcp(new panzer_stk_classic::ParameterListCallbackBlocked<int,GO>(stkConn_manager,blkDofs));
          reqHandler->addRequestCallback(callback);
       }

       Teko::addTekoToStratimikosBuilder(linearSolverBuilder,reqHandler);

       bool writeCoordinates = p.sublist("Options").get("Write Coordinates",false);
       if(writeCoordinates) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);

          // loop over blocks
          const std::vector<Teuchos::RCP<panzer::UniqueGlobalIndexer<int,GO> > > & dofVec
             = blkDofs->getFieldDOFManagers();
          for(std::size_t i=0;i<dofVec.size();i++) {
            std::string fieldName;

            // add in the coordinate parameter list callback handler
            TEUCHOS_ASSERT(determineCoordinateField(*dofVec[i],fieldName));

            std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > fieldPatterns;
            fillFieldPatternMap(*dofVec[i],fieldName,fieldPatterns);
            panzer_stk_classic::ParameterListCallback<int,GO> plCall(fieldName,fieldPatterns,stkConn_manager,dofVec[i]);
            plCall.buildArrayToVector();
            plCall.buildCoordinates();

            // extract coordinate vectors
            const std::vector<double> & xcoords = plCall.getXCoordsVector();
            const std::vector<double> & ycoords = plCall.getYCoordsVector();
            const std::vector<double> & zcoords = plCall.getZCoordsVector();

            // use epetra to write coordinates to matrix market files
            Epetra_MpiComm ep_comm(*mpi_comm->getRawMpiComm()); // this is OK access to RawMpiComm becase its declared on the stack?
                                                                // and all users of this object are on the stack (within scope of mpi_comm
            Epetra_Map map(-1,xcoords.size(),0,ep_comm);

            Teuchos::RCP<Epetra_Vector> vec;
            switch(mesh->getDimension()) {
            case 3:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&zcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_zcoords.mm").c_str(),*vec);
            case 2:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&ycoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_ycoords.mm").c_str(),*vec);
            case 1:
               vec = Teuchos::rcp(new Epetra_Vector(Copy,map,const_cast<double *>(&xcoords[0])));
               EpetraExt::VectorToMatrixMarketFile((fieldName+"_xcoords.mm").c_str(),*vec);
               break;
            default:
               TEUCHOS_ASSERT(false);
            }
          }
       }

       bool writeTopo = p.sublist("Options").get("Write Topology",false);
       if(writeTopo) {
          Teuchos::RCP<const panzer::BlockedDOFManager<int,GO> > blkDofs =
             Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager<int,GO> >(globalIndexer);

          writeTopology(*blkDofs);
       }
    }
    #endif

    #ifdef HAVE_MUELU
    {
      Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder); // Register MueLu as a Stratimikos preconditioner strategy for Epetra
      Stratimikos::enableMueLuTpetra<int,panzer::Ordinal64,KokkosClassic::DefaultNode::DefaultNodeType>(linearSolverBuilder,"MueLu-Tpetra");
    }
    #endif // MUELU
    #ifdef HAVE_IFPACK2
    {
      typedef Thyra::PreconditionerFactoryBase<double> Base;
      typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double, int, panzer::Ordinal64,KokkosClassic::DefaultNode::DefaultNodeType> > Impl;

      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    }
    #endif // MUELU

    linearSolverBuilder.setParameterList(strat_params);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);

    return lowsFactory;
  }

  template<typename ScalarT>
  template<typename GO>
  void ModelEvaluatorFactory<ScalarT>::
  writeTopology(const panzer::BlockedDOFManager<int,GO> & blkDofs) const
  {
    using Teuchos::RCP;

    // loop over each field block
    const std::vector<RCP<panzer::UniqueGlobalIndexer<int,GO> > > & blk_dofMngrs = blkDofs.getFieldDOFManagers();
    for(std::size_t b=0;b<blk_dofMngrs.size();b++) {
      RCP<panzer::DOFManagerFEI<int,GO> > dofMngr = Teuchos::rcp_dynamic_cast<panzer::DOFManagerFEI<int,GO> >(blk_dofMngrs[b],true);

      std::vector<std::string> eBlocks;
      dofMngr->getElementBlockIds(eBlocks);

      // build file name
      std::stringstream fileName;
      fileName << "elements_" << b;
      std::ofstream file(fileName.str().c_str());

      // loop over each element block, write out topology
      for(std::size_t e=0;e<eBlocks.size();e++)
        writeTopology(*dofMngr,eBlocks[e],file);
    }
  }

  template<typename ScalarT>
  template <typename GO>
  void ModelEvaluatorFactory<ScalarT>::
  writeTopology(const panzer::DOFManagerFEI<int,GO> & dofs,const std::string & block,std::ostream & os) const
  {
    std::vector<std::string> fields(dofs.getElementBlockGIDCount(block));

    const std::set<int> & fieldIds = dofs.getFields(block);
    for(std::set<int>::const_iterator itr=fieldIds.begin();itr!=fieldIds.end();++itr) {
      std::string field = dofs.getFieldString(*itr);

      // get the layout of each field
      const std::vector<int> & fieldOffsets = dofs.getGIDFieldOffsets(block,*itr);
      for(std::size_t f=0;f<fieldOffsets.size();f++)
        fields[fieldOffsets[f]] = field;

    }

    // print the layout of the full pattern
    os << "#" << std::endl;
    os << "# Element Block \"" << block << "\"" << std::endl;
    os << "#   field pattern = [ " << fields[0];
    for(std::size_t f=1;f<fields.size();f++)
      os << ", " << fields[f];
    os << " ]" << std::endl;
    os << "#" << std::endl;

    const std::vector<int> & elements = dofs.getElementBlock(block);
    for(std::size_t e=0;e<elements.size();e++) {
      std::vector<GO> gids;
      dofs.getElementGIDs(elements[e],gids,block);

      // output gids belonging to this element
      os << "[ " << gids[0];
      for(std::size_t g=1;g<gids.size();g++)
        os << ", " << gids[g];
      os << " ]" << std::endl;
    }
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  buildResponses(const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                 const bool write_graphviz_file,
                 const std::string& graphviz_file_prefix)
  {
    typedef panzer::ModelEvaluator<double> PanzerME;

    Teuchos::ParameterList & p = *this->getNonconstParameterList();
    Teuchos::ParameterList & user_data = p.sublist("User Data");
    Teuchos::ParameterList & closure_models = p.sublist("Closure Models");

    // uninitialize the thyra model evaluator, its respone counts are wrong!
    Teuchos::RCP<Thyra::EpetraModelEvaluator> thyra_me = Teuchos::rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(m_physics_me);
    Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(m_physics_me);

    if(thyra_me!=Teuchos::null && panzer_me==Teuchos::null) {
      Teuchos::RCP<const EpetraExt::ModelEvaluator> const_ep_me;
      Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFactory;
      thyra_me->uninitialize(&const_ep_me,&solveFactory); // this seems dangerous!

      // I don't need no const-ness!
      Teuchos::RCP<EpetraExt::ModelEvaluator> ep_me = Teuchos::rcp_const_cast<EpetraExt::ModelEvaluator>(const_ep_me);
      Teuchos::RCP<panzer::ModelEvaluator_Epetra> ep_panzer_me = Teuchos::rcp_dynamic_cast<panzer::ModelEvaluator_Epetra>(ep_me);

      ep_panzer_me->buildResponses(m_physics_blocks,*m_eqset_factory,cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix);

      // reinitialize the thyra model evaluator, now with the correct responses
      thyra_me->initialize(ep_me,solveFactory);

      return;
    }
    else if(panzer_me!=Teuchos::null && thyra_me==Teuchos::null) {
      panzer_me->buildResponses(m_physics_blocks,*m_eqset_factory,cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix);

      return;
    }

    TEUCHOS_ASSERT(false);
  }
}

#endif
