// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_T_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_InitialCondition_Builder.hpp"
#include "Panzer_ModelEvaluator.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_ExplicitModelEvaluator.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_MultiBlockMeshFactory.hpp"
#include "Panzer_STK_CustomMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_NOXObserverFactory.hpp"
#ifdef PANZER_HAVE_TEMPUS
#include "Panzer_STK_TempusObserverFactory.hpp"
#endif
#include "Panzer_STK_ParameterListCallback.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"
#include "Panzer_STK_IOClosureModel_Factory_TemplateBuilder.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_SetupLOWSFactory.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib> // for std::getenv

// Piro solver objects
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"
#include "Piro_LOCASolver.hpp"
#ifdef PANZER_HAVE_TEMPUS
#include "Piro_TempusSolverForwardOnly.hpp"
#endif

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#include <Panzer_NodeType.hpp>

namespace panzer_stk {

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
        p.set<std::string>("Auxiliary Field Order","");
        p.set<bool>("Use DOFManager FEI",false);
        p.set<bool>("Load Balance DOFs",false);
        p.set<bool>("Use Tpetra",false);
        p.set<bool>("Use Epetra ME",true);
        p.set<bool>("Lump Explicit Mass",false);
        p.set<bool>("Constant Mass Matrix",true);
        p.set<bool>("Apply Mass Matrix Inverse in Explicit Evaluator",true);
        p.set<bool>("Use Conservative IMEX",false);
        p.set<bool>("Compute Real Time Derivative",false);
        p.set<bool>("Use Time Derivative in Explicit Model",false);
        p.set<bool>("Compute Time Derivative at Time Step",false);
        p.set<Teuchos::RCP<const panzer::EquationSetFactory> >("Equation Set Factory", Teuchos::null);
        p.set<Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >("Closure Model Factory", Teuchos::null);
        p.set<Teuchos::RCP<const panzer::BCStrategyFactory> >("BC Factory",Teuchos::null);
        p.set<std::string>("Excluded Blocks","");
        p.sublist("ALE").disableRecursiveValidation();
      }

      pl->sublist("Block ID to Physics ID Mapping").disableRecursiveValidation();
      pl->sublist("Options").disableRecursiveValidation();
      pl->sublist("Active Parameters").disableRecursiveValidation();
      pl->sublist("Controls").disableRecursiveValidation();
      pl->sublist("ALE").disableRecursiveValidation(); // this sucks
      pl->sublist("User Data").disableRecursiveValidation();
      pl->sublist("User Data").sublist("Panzer Data").disableRecursiveValidation();

      validPL = pl;
    }
    return validPL;
  }

  namespace {
    bool hasInterfaceCondition(const std::vector<panzer::BC>& bcs)
    {
      for (std::vector<panzer::BC>::const_iterator bcit = bcs.begin(); bcit != bcs.end(); ++bcit)
        if (bcit->bcType() == panzer::BCT_Interface)
          return true;
      return false;
    }

    Teuchos::RCP<STKConnManager>
    getSTKConnManager(const Teuchos::RCP<panzer::ConnManager>& conn_mgr)
    {
      const Teuchos::RCP<STKConnManager> stk_conn_mgr =
        Teuchos::rcp_dynamic_cast<STKConnManager>(conn_mgr);
      TEUCHOS_TEST_FOR_EXCEPTION(stk_conn_mgr.is_null(), std::logic_error,
                                 "There are interface conditions, but the connection manager"
                                 " does not support the necessary connections.");
      return stk_conn_mgr;
    }

    void buildInterfaceConnections(const std::vector<panzer::BC>& bcs,
                                   const Teuchos::RCP<panzer::ConnManager>& conn_mgr)
    {
      const Teuchos::RCP<STKConnManager> stk_conn_mgr = getSTKConnManager(conn_mgr);
      for (std::vector<panzer::BC>::const_iterator bcit = bcs.begin(); bcit != bcs.end(); ++bcit)
        if (bcit->bcType() == panzer::BCT_Interface)
          stk_conn_mgr->associateElementsInSideset(bcit->sidesetID());
    }

    void checkInterfaceConnections(const Teuchos::RCP<panzer::ConnManager>& conn_mgr,
                                   const Teuchos::RCP<Teuchos::Comm<int> >& comm)
    {
      const Teuchos::RCP<STKConnManager> stk_conn_mgr = getSTKConnManager(conn_mgr);
      std::vector<std::string> sidesets = stk_conn_mgr->checkAssociateElementsInSidesets(*comm);
      if ( ! sidesets.empty()) {
        std::stringstream ss;
        ss << "Sideset IDs";
        for (std::size_t i = 0; i < sidesets.size(); ++i)
          ss << " " << sidesets[i];
        ss << " did not yield associations, but these sidesets correspond to BCT_Interface BCs.";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, ss.str());
      }
    }
  } // namespace

  template<typename ScalarT>
  void  ModelEvaluatorFactory<ScalarT>::buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                                     const Teuchos::RCP<panzer::GlobalData>& global_data,
                                                     const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                                                     const panzer::BCStrategyFactory & bc_factory,
                                                     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & user_cm_factory,
                                                     bool meConstructionOn)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->getParameterList()), std::runtime_error,
                       "ParameterList must be set before objects can be built!");

    TEUCHOS_ASSERT(nonnull(comm));
    TEUCHOS_ASSERT(nonnull(global_data));
    TEUCHOS_ASSERT(nonnull(global_data->os));
    TEUCHOS_ASSERT(nonnull(global_data->pl));

    // begin at the beginning...
    m_global_data = global_data;

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    // Parse input file, setup parameters
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

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

    Teuchos::RCP<Teuchos::ParameterList> physics_block_plist = Teuchos::sublist(this->getMyNonconstParamList(),"Physics Blocks");

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
    bool is_transient  = (solncntl_params.get<std::string>("Piro Solver") == "Tempus") ? true : false;
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

    useDiscreteAdjoint = p.get<bool>("Use Discrete Adjoint");

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    // Do stuff
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    Teuchos::FancyOStream& fout = *global_data->os;

    // for convience cast to an MPI comm
    const Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);

    // Build mesh factory and uncommitted mesh
    ////////////////////////////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory = this->buildSTKMeshFactory(mesh_params);
    Teuchos::RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*(mpi_comm->getRawMpiComm()));
    m_mesh = mesh;

    m_eqset_factory = eqset_factory;

    // setup the physcs blocks
    ////////////////////////////////////////////////////////////////////////////////////////

    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    {
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

      // build physics blocks

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
    }

    // add fields automatically written through the closure model
    ////////////////////////////////////////////////////////////////////////////////////////
    addUserFieldsToMesh(*mesh,output_list);

    // finish building mesh, set required field variables and mesh bulk data
    ////////////////////////////////////////////////////////////////////////////////////////

    try {
       // this throws some exceptions, catch them as neccessary
       this->finalizeMeshConstruction(*mesh_factory,physicsBlocks,*mpi_comm,*mesh);
    } catch(const panzer_stk::STK_Interface::ElementBlockException & ebexp) {
       fout << "*****************************************\n\n";
       fout << "Element block exception, could not finalize the mesh, printing block and sideset information:\n";
       fout.pushTab(3);
       mesh->printMetaData(fout);
       fout.popTab();
       fout << std::endl;

       throw ebexp;
    } catch(const panzer_stk::STK_Interface::SidesetException & ssexp) {
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
      mesh->setupExodusFile(p.sublist("Output").get<std::string>("File Name"));

    // build a workset factory that depends on STK
    ////////////////////////////////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer_stk::WorksetFactory> wkstFactory;
    if(m_user_wkst_factory==Teuchos::null)
       wkstFactory = Teuchos::rcp(new panzer_stk::WorksetFactory()); // build STK workset factory
    else
       wkstFactory = m_user_wkst_factory;

     // set workset factory mesh
     wkstFactory->setMesh(mesh);

    // handle boundary and interface conditions
    ////////////////////////////////////////////////////////////////////////////////////////
    std::vector<panzer::BC> bcs;
    panzer::buildBCs(bcs, p.sublist("Boundary Conditions"), global_data);

    // build the connection manager
    ////////////////////////////////////////////////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    m_conn_manager = conn_manager;

    // build DOF Manager
    ////////////////////////////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    Teuchos::RCP<panzer::GlobalIndexer> globalIndexer;

    std::string loadBalanceString = ""; // what is the load balancing information
    bool blockedAssembly = false;

    const bool has_interface_condition = hasInterfaceCondition(bcs);

    if(panzer::BlockedDOFManagerFactory::requiresBlocking(field_order) && !useTpetra) {

#ifdef PANZER_HAVE_EPETRA_STACK
       // Can't yet handle interface conditions for this system
       TEUCHOS_TEST_FOR_EXCEPTION(has_interface_condition,
                                  Teuchos::Exceptions::InvalidParameter,
                                  "ERROR: Blocked Epetra systems cannot handle interface conditions.");

       // use a blocked DOF manager
       blockedAssembly = true;

       panzer::BlockedDOFManagerFactory globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(use_dofmanager_fei);

       Teuchos::RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager,field_order);
       globalIndexer = dofManager;

       Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > bloLinObjFactory
        = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpi_comm,
                                                          Teuchos::rcp_dynamic_cast<panzer::BlockedDOFManager>(dofManager)));

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
#else
       TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: buildObjects() - Epetra support is NOT enabled in this build!");
#endif
    }
    else if(panzer::BlockedDOFManagerFactory::requiresBlocking(field_order) && useTpetra) {

       // Can't yet handle interface conditions for this system
       TEUCHOS_TEST_FOR_EXCEPTION(has_interface_condition,
                                  Teuchos::Exceptions::InvalidParameter,
                                  "ERROR: Blocked Tpetra system cannot handle interface conditions.");

       // use a blocked DOF manager
       blockedAssembly = true;

       TEUCHOS_ASSERT(!use_dofmanager_fei);
       panzer::BlockedDOFManagerFactory globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(false);

       Teuchos::RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager,field_order);
       globalIndexer = dofManager;

       Teuchos::RCP<panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> > bloLinObjFactory
        = Teuchos::rcp(new panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(mpi_comm,
                                                          Teuchos::rcp_dynamic_cast<panzer::BlockedDOFManager>(dofManager)));

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

       if (has_interface_condition)
         buildInterfaceConnections(bcs, conn_manager);

       // use a flat DOF manager

       TEUCHOS_ASSERT(!use_dofmanager_fei);
       panzer::DOFManagerFactory globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(false);
       globalIndexerFactory.setUseTieBreak(use_load_balance);
       Teuchos::RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager,field_order);
       globalIndexer = dofManager;

       if (has_interface_condition)
         checkInterfaceConnections(conn_manager, dofManager->getComm());

       TEUCHOS_ASSERT(!useDiscreteAdjoint); // safety check
       linObjFactory = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(mpi_comm,dofManager));

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
    }
    else {

#ifdef PANZER_HAVE_EPETRA_STACK
       if (has_interface_condition)
         buildInterfaceConnections(bcs, conn_manager);

       // use a flat DOF manager
       panzer::DOFManagerFactory globalIndexerFactory;
       globalIndexerFactory.setUseDOFManagerFEI(use_dofmanager_fei);
       globalIndexerFactory.setUseTieBreak(use_load_balance);
       globalIndexerFactory.setUseNeighbors(has_interface_condition);
       Teuchos::RCP<panzer::GlobalIndexer> dofManager
         = globalIndexerFactory.buildGlobalIndexer(mpi_comm->getRawMpiComm(),physicsBlocks,conn_manager,
                                                         field_order);
       globalIndexer = dofManager;

       if (has_interface_condition)
         checkInterfaceConnections(conn_manager, dofManager->getComm());

       linObjFactory = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(mpi_comm,dofManager,useDiscreteAdjoint));

       // build load balancing string for informative output
       loadBalanceString = printUGILoadBalancingInformation(*dofManager);
#else
       TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: buildObjects() - Epetra support is NOT enabled in this build!");
#endif
    }

    TEUCHOS_ASSERT(globalIndexer!=Teuchos::null);
    TEUCHOS_ASSERT(linObjFactory!=Teuchos::null);
    m_global_indexer = globalIndexer;
    m_lin_obj_factory = linObjFactory;
    m_blockedAssembly = blockedAssembly;

    // print out load balancing information
    fout << "Degree of freedom load balancing: " << loadBalanceString << std::endl;

    // build worksets
    //////////////////////////////////////////////////////////////

    // build up needs array for workset container
    std::map<std::string,panzer::WorksetNeeds> needs;
    for(std::size_t i=0;i<physicsBlocks.size();i++)
      needs[physicsBlocks[i]->elementBlockID()] = physicsBlocks[i]->getWorksetNeeds();

    Teuchos::RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory,needs));

    wkstContainer->setWorksetSize(workset_size);
    wkstContainer->setGlobalIndexer(globalIndexer); // set the global indexer so the orientations are evaluated

    m_wkstContainer = wkstContainer;

    // find max number of worksets
    std::size_t max_wksets = 0;
    for(std::size_t pb=0;pb<physicsBlocks.size();pb++) {
      const panzer::WorksetDescriptor wd = panzer::blockDescriptor(physicsBlocks[pb]->elementBlockID());
      Teuchos::RCP< std::vector<panzer::Workset> >works = wkstContainer->getWorksets(wd);
      max_wksets = std::max(max_wksets,works->size());
    }
    user_data_params.set<std::size_t>("Max Worksets",max_wksets);
    wkstContainer->clear();

    // Setup lagrangian type coordinates
    /////////////////////////////////////////////////////////////

    // see if field coordinates are required, if so reset the workset container
    // and set the coordinates to be associated with a field in the mesh
    useDynamicCoordinates_ = false;
    for(std::size_t pb=0;pb<physicsBlocks.size();pb++) {
      if(physicsBlocks[pb]->getCoordinateDOFs().size()>0) {
         mesh->setUseFieldCoordinates(true);
         useDynamicCoordinates_ = true;
         wkstContainer->clear(); // this serves to refresh the worksets
                                 // and put in new coordinates
         break;
      }
    }

    // Add mesh objects to user data to make available to user ctors
    /////////////////////////////////////////////////////////////

    panzer_data_params.set("STK Mesh", mesh);
    panzer_data_params.set("DOF Manager", globalIndexer);
    panzer_data_params.set("Linear Object Factory", linObjFactory);

    // If user requested it, short circuit model construction
    ////////////////////////////////////////////////////////////////////////////////////////

    if(!meConstructionOn)
      return;

    // Setup active parameters
    /////////////////////////////////////////////////////////////

    std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names;
    std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values;
    if (p.isSublist("Active Parameters")) {
      Teuchos::ParameterList& active_params = p.sublist("Active Parameters");

      int num_param_vecs = active_params.get<int>("Number of Parameter Vectors",0);
      p_names.resize(num_param_vecs);
      p_values.resize(num_param_vecs);
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
        p_values[i] =
          Teuchos::rcp(new Teuchos::Array<double>(numParameters));
        for (int j=0; j<numParameters; j++) {
          std::stringstream ss2;
          ss2 << "Parameter " << j;
          (*p_names[i])[j] = pList.get<std::string>(ss2.str());
          ss2.str("");

          ss2 << "Initial Value " << j;
          (*p_values[i])[j] = pList.get<double>(ss2.str());

          // this is a band-aid/hack to make sure parameters are registered before they are accessed
          panzer::registerScalarParameter((*p_names[i])[j],*global_data->pl,(*p_values[i])[j]);
        }
      }
    }

    // setup the closure model for automatic writing (during residual/jacobian update)
    ////////////////////////////////////////////////////////////////////////////////////////

    panzer_stk::IOClosureModelFactory_TemplateBuilder<panzer::Traits> io_cm_builder(user_cm_factory,mesh,output_list);
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> cm_factory;
    cm_factory.buildObjects(io_cm_builder);

    // setup field manager build
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    {
      bool write_dot_files = p.sublist("Options").get("Write Volume Assembly Graphs",false);
      std::string dot_file_prefix = p.sublist("Options").get("Volume Assembly Graph Prefix","Panzer_AssemblyGraph");
      bool write_fm_files = p.sublist("Options").get("Write Field Manager Files",false);
      std::string fm_file_prefix = p.sublist("Options").get("Field Manager File Prefix","Panzer_AssemblyGraph");

      // Allow users to override inputs via runtime env
      {
        auto check_write_dag = std::getenv("PANZER_WRITE_DAG");
        if (check_write_dag != nullptr) {
          write_dot_files = true;
          write_fm_files = true;
        }
      }

      fmb = buildFieldManagerBuilder(wkstContainer,physicsBlocks,bcs,*eqset_factory,bc_factory,cm_factory,
                                     user_cm_factory,p.sublist("Closure Models"),*linObjFactory,user_data_params,
                                     write_dot_files,dot_file_prefix,
				     write_fm_files,fm_file_prefix);
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
                                     p_values,
                                     lowsFactory,
                                     global_data,
                                     is_transient,
                                     t_init);

    // Setup initial conditions
    /////////////////////////////////////////////////////////////

    {
      // Create closure model list for use in defining initial conditions
      // For now just remove Global MMS Parameters, if it exists
      const Teuchos::ParameterList& models = p.sublist("Closure Models");
      Teuchos::ParameterList cl_models(models.name());
      for (Teuchos::ParameterList::ConstIterator model_it=models.begin();
           model_it!=models.end(); ++model_it) {
           std::string key = model_it->first;
           if (model_it->first != "Global MMS Parameters")
              cl_models.setEntry(key,model_it->second);
       }
      bool write_dot_files = false;
      std::string prefix = "Panzer_AssemblyGraph_";
      setupInitialConditions(*thyra_me,*wkstContainer,physicsBlocks,user_cm_factory,*linObjFactory,
                             cl_models,
                             p.sublist("Initial Conditions"),
                             p.sublist("User Data"),
                             p.sublist("Options").get("Write Volume Assembly Graphs",write_dot_files),
                             p.sublist("Options").get("Volume Assembly Graph Prefix",prefix));
    }

    // Write the IC vector into the STK mesh: use response library
    //////////////////////////////////////////////////////////////////////////
    writeInitialConditions(*thyra_me,physicsBlocks,wkstContainer,globalIndexer,linObjFactory,mesh,user_cm_factory,
                           p.sublist("Closure Models"),
                           p.sublist("User Data"),workset_size);

    m_physics_me = thyra_me;
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  addUserFieldsToMesh(panzer_stk::STK_Interface & mesh,const Teuchos::ParameterList & output_list) const
  {
    // register cell averaged scalar fields
    const Teuchos::ParameterList & cellAvgQuants = output_list.sublist("Cell Average Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=cellAvgQuants.begin();
        itr!=cellAvgQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh.addCellField(tokens[i],blockId);
    }

    // register cell averaged components of vector fields
    // just allocate space for the fields here. The actual calculation and writing
    // are done by panzer_stk::ScatterCellAvgVector.
    const Teuchos::ParameterList & cellAvgVectors = output_list.sublist("Cell Average Vectors");
    for(Teuchos::ParameterList::ConstIterator itr = cellAvgVectors.begin();
        itr != cellAvgVectors.end(); ++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i = 0; i < tokens.size(); i++) {
          std::string d_mod[3] = {"X","Y","Z"};
          for(std::size_t d = 0; d < mesh.getDimension(); d++)
              mesh.addCellField(tokens[i]+d_mod[d],blockId);
       }
    }

    // register cell quantities
    const Teuchos::ParameterList & cellQuants = output_list.sublist("Cell Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=cellQuants.begin();
        itr!=cellQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh.addCellField(tokens[i],blockId);
    }

    // register ndoal quantities
    const Teuchos::ParameterList & nodalQuants = output_list.sublist("Nodal Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=nodalQuants.begin();
        itr!=nodalQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh.addSolutionField(tokens[i],blockId);
    }

    const Teuchos::ParameterList & allocNodalQuants = output_list.sublist("Allocate Nodal Quantities");
    for(Teuchos::ParameterList::ConstIterator itr=allocNodalQuants.begin();
        itr!=allocNodalQuants.end();++itr) {
       const std::string & blockId = itr->first;
       const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
       std::vector<std::string> tokens;

       // break up comma seperated fields
       panzer::StringTokenizer(tokens,fields,",",true);

       for(std::size_t i=0;i<tokens.size();i++)
          mesh.addSolutionField(tokens[i],blockId);
    }
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  setupInitialConditions(Thyra::ModelEvaluator<ScalarT> & model,
                         panzer::WorksetContainer & wkstContainer,
                         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                         const panzer::LinearObjFactory<panzer::Traits> & lof,
                         const Teuchos::ParameterList & closure_pl,
                         const Teuchos::ParameterList & initial_cond_pl,
                         const Teuchos::ParameterList & user_data_pl,
                         bool write_dot_files,const std::string & dot_file_prefix) const
  {
    using Teuchos::RCP;

    Thyra::ModelEvaluatorBase::InArgs<double> nomValues = model.getNominalValues();
    Teuchos::RCP<Thyra::VectorBase<double> > x_vec = Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(nomValues.get_x());

    if(initial_cond_pl.get<bool>("Zero Initial Conditions")) {
      // zero out the x vector
      Thyra::assign(x_vec.ptr(),0.0);
    }
    else if(!initial_cond_pl.sublist("Vector File").get<bool>("Enabled")) {
      // read from exodus, or compute using field managers

      std::map<std::string, Teuchos::RCP< PHX::FieldManager<panzer::Traits> > > phx_ic_field_managers;

      panzer::setupInitialConditionFieldManagers(wkstContainer,
                                                 physicsBlocks,
                                                 cm_factory,
                                                 closure_pl,
                                                 initial_cond_pl,
                                                 lof,
                                                 user_data_pl,
                                                 write_dot_files,
                                                 dot_file_prefix,
                                                 phx_ic_field_managers);
/*
      panzer::setupInitialConditionFieldManagers(wkstContainer,
                                                 physicsBlocks,
                                                 cm_factory,
                                                 initial_cond_pl,
                                                 lof,
                                                 user_data_pl,
                                                 write_dot_files,
                                                 dot_file_prefix,
                                                 phx_ic_field_managers);
*/

      // set the vector to be filled
      Teuchos::RCP<panzer::LinearObjContainer> loc = lof.buildLinearObjContainer();
      Teuchos::RCP<panzer::ThyraObjContainer<double> > tloc = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc);
      tloc->set_x_th(x_vec);

      panzer::evaluateInitialCondition(wkstContainer, phx_ic_field_managers, loc, lof, 0.0);
    }
    else {
      const std::string & vectorFile = initial_cond_pl.sublist("Vector File").get<std::string>("File Name");
      TEUCHOS_TEST_FOR_EXCEPTION(vectorFile=="",std::runtime_error,
                                 "If \"Read From Vector File\" is true, then parameter \"Vector File\" cannot be the empty string.");

      // set the vector to be filled
      Teuchos::RCP<panzer::LinearObjContainer> loc = lof.buildLinearObjContainer();
      Teuchos::RCP<panzer::ThyraObjContainer<double> > tloc = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc);
      tloc->set_x_th(x_vec);

      // read the vector
      lof.readVector(vectorFile,*loc,panzer::LinearObjContainer::X);
    }
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  writeInitialConditions(const Thyra::ModelEvaluator<ScalarT> & model,
                         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                         const Teuchos::RCP<panzer::WorksetContainer> & wc,
                         const Teuchos::RCP<const panzer::GlobalIndexer> & ugi,
                         const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                         const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                         const Teuchos::ParameterList & closure_model_pl,
                         const Teuchos::ParameterList & user_data_pl,
                         int workset_size) const
  {
    Teuchos::RCP<panzer::LinearObjContainer> loc = lof->buildLinearObjContainer();
    Teuchos::RCP<panzer::ThyraObjContainer<double> > tloc = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc);
    tloc->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(model.getNominalValues().get_x()));

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > solnWriter
        = initializeSolnWriterResponseLibrary(wc,ugi,lof,mesh);

    {
       Teuchos::ParameterList user_data(user_data_pl);
       user_data.set<int>("Workset Size",workset_size);

       finalizeSolnWriterResponseLibrary(*solnWriter,physicsBlocks,cm_factory,closure_model_pl,workset_size,user_data);
    }

    // initialize the assembly container
    panzer::AssemblyEngineInArgs ae_inargs;
    ae_inargs.container_ = loc;
    ae_inargs.ghostedContainer_ = lof->buildGhostedLinearObjContainer();
    ae_inargs.alpha = 0.0;
    ae_inargs.beta = 1.0;
    ae_inargs.evaluate_transient_terms = false;

    // initialize the ghosted container
    lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

    // do import
    lof->globalToGhostContainer(*ae_inargs.container_,*ae_inargs.ghostedContainer_,panzer::LinearObjContainer::X);

    // fill STK mesh objects
    solnWriter->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
    solnWriter->evaluate<panzer::Traits::Residual>(ae_inargs);
  }

  //! build STK mesh from a mesh parameter list
  template<typename ScalarT>
  Teuchos::RCP<panzer_stk::STK_MeshFactory> ModelEvaluatorFactory<ScalarT>::buildSTKMeshFactory(const Teuchos::ParameterList & mesh_params) const
  {
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;

    // first contruct the mesh factory
    if (mesh_params.get<std::string>("Source") ==  "Exodus File") {
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory());
      mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Exodus File"))));
    }
    else if (mesh_params.get<std::string>("Source") ==  "Pamgen Mesh") {
      mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory());
      Teuchos::RCP<Teuchos::ParameterList> pamgenList = Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Pamgen Mesh")));
      pamgenList->set("File Type","Pamgen"); // For backwards compatibility when pamgen had separate factory from exodus
      mesh_factory->setParameterList(pamgenList);
    }
    else if (mesh_params.get<std::string>("Source") ==  "Inline Mesh") {

      int dimension = mesh_params.sublist("Inline Mesh").get<int>("Mesh Dimension");
      std::string typeStr = "";
      if(mesh_params.sublist("Inline Mesh").isParameter("Type"))
         typeStr = mesh_params.sublist("Inline Mesh").get<std::string>("Type");

      if (dimension == 1) {
        mesh_factory = Teuchos::rcp(new panzer_stk::LineMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 2 && typeStr=="Tri") {
        mesh_factory = Teuchos::rcp(new panzer_stk::SquareTriMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 2) {
        mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
      else if (dimension == 3 && typeStr=="Tet") {
        mesh_factory = Teuchos::rcp(new panzer_stk::CubeTetMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension == 3) {
        mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
      else if(dimension==4) { // not really "dimension==4" simply a flag to try this other mesh for testing
        mesh_factory = Teuchos::rcp(new panzer_stk::MultiBlockMeshFactory);
        Teuchos::RCP<Teuchos::ParameterList> in_mesh = Teuchos::rcp(new Teuchos::ParameterList);
        *in_mesh = mesh_params.sublist("Inline Mesh").sublist("Mesh Factory Parameter List");
        mesh_factory->setParameterList(in_mesh);
      }
    }
    else if (mesh_params.get<std::string>("Source") ==  "Custom Mesh") {
      mesh_factory = Teuchos::rcp(new panzer_stk::CustomMeshFactory());
      mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(mesh_params.sublist("Custom Mesh"))));
    }
    else {
      // throw a runtime exception for invalid parameter values
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
  void ModelEvaluatorFactory<ScalarT>::setNOXObserverFactory(const Teuchos::RCP<const panzer_stk::NOXObserverFactory>& nox_observer_factory)
  {
    m_nox_observer_factory = nox_observer_factory;
  }

#ifdef PANZER_HAVE_TEMPUS
  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::setTempusObserverFactory(const Teuchos::RCP<const panzer_stk::TempusObserverFactory>& tempus_observer_factory)
  {
    m_tempus_observer_factory = tempus_observer_factory;
  }
#endif

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::setUserWorksetFactory(Teuchos::RCP<panzer_stk::WorksetFactory>& user_wkst_factory)
  {
    m_user_wkst_factory = user_wkst_factory;
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
#ifdef PANZER_HAVE_TEMPUS
                                  const Teuchos::RCP<Piro::TempusSolverForwardOnly<ScalarT> > tempusSolver,
#endif
                                  const Teuchos::Ptr<const panzer_stk::NOXObserverFactory> & in_nox_observer_factory
#ifdef PANZER_HAVE_TEMPUS
                                  , const Teuchos::Ptr<const panzer_stk::TempusObserverFactory> & in_tempus_observer_factory
#endif
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
    Teuchos::Ptr<const panzer_stk::NOXObserverFactory> nox_observer_factory
        = is_null(in_nox_observer_factory) ? m_nox_observer_factory.ptr() : in_nox_observer_factory;
#ifdef PANZER_HAVE_TEMPUS
    Teuchos::Ptr<const panzer_stk::TempusObserverFactory> tempus_observer_factory
        = is_null(in_tempus_observer_factory) ? m_tempus_observer_factory.ptr() : in_tempus_observer_factory;
#endif

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
#ifdef PANZER_HAVE_TEMPUS
    else if (solver=="Tempus") {

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(tempus_observer_factory), std::runtime_error,
                                 "No Tempus observer built! Please call setTempusObserverFactory() member function if you plan to use a Tempus solver.");

      // install the nox observer
      if(tempus_observer_factory->useNOXObserver()) {
        Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = nox_observer_factory->buildNOXObserver(m_mesh,m_global_indexer,m_lin_obj_factory);
        piro_params->sublist("NOX").sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);
      }

      // override printing to use panzer ostream
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Output Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<Teuchos::RCP<std::ostream> >("Error Stream",global_data->os);
      piro_params->sublist("NOX").sublist("Printing").set<int>("Output Processor",global_data->os->getOutputToRootOnly());

      // use the user specfied tempus solver if they pass one in
      Teuchos::RCP<Piro::TempusSolverForwardOnly<double> > piro_tempus;

      if(tempusSolver==Teuchos::null)
      {
        piro_tempus =
          Teuchos::rcp(new Piro::TempusSolverForwardOnly<double>(piro_params, thyra_me,
                                                                 tempus_observer_factory->buildTempusObserver(m_mesh,m_global_indexer,m_lin_obj_factory)));
      }
      else
      {
        piro_tempus = tempusSolver;
        piro_tempus->initialize(piro_params, thyra_me,
                                tempus_observer_factory->buildTempusObserver(m_mesh,m_global_indexer,m_lin_obj_factory));
      }

      piro = piro_tempus;
    }
#endif
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
                           bool writeGraph,const std::string & graphPrefix,
			   bool write_field_managers,const std::string & field_manager_prefix) const
  {
    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
    fmb->setWorksetContainer(wc);
    fmb->setupVolumeFieldManagers(physicsBlocks,volume_cm_factory,closure_models,lo_factory,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,bc_cm_factory,bc_factory,closure_models,lo_factory,user_data);

    // Print Phalanx DAGs
    if (writeGraph){
      fmb->writeVolumeGraphvizDependencyFiles(graphPrefix, physicsBlocks);
      fmb->writeBCGraphvizDependencyFiles(graphPrefix);
    }
    if (write_field_managers){
      fmb->writeVolumeTextDependencyFiles(graphPrefix, physicsBlocks);
      fmb->writeBCTextDependencyFiles(field_manager_prefix);
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
                            const Teuchos::Ptr<const Teuchos::ParameterList> & bc_list,
                            const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & physics_me_in) const
  {
    typedef panzer::ModelEvaluator<ScalarT> PanzerME;

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > physics_me = physics_me_in==Teuchos::null ? m_physics_me : physics_me_in;

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
        panzer::buildBCs(bcs, p.sublist("Boundary Conditions"), m_global_data);
      }
      else {
        panzer::buildBCs(bcs, *bc_list, m_global_data);
      }

      fmb = buildFieldManagerBuilder(// Teuchos::rcp_const_cast<panzer::WorksetContainer>(
                                     // m_response_library!=Teuchos::null ? m_response_library->getWorksetContainer()
                                     //                                   : m_wkstContainer),
                                     m_wkstContainer,
                                     physicsBlocks,
                                     bcs,
                                     *eqset_factory,
                                     bc_factory,
                                     user_cm_factory,
                                     user_cm_factory,
                                     p.sublist("Closure Models"),
                                     *m_lin_obj_factory,
                                     user_data_params,
                                     write_dot_files,prefix,
				     write_dot_files,prefix);
    }

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > response_library
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(m_wkstContainer,
                                                                   m_global_indexer,
                                                                   m_lin_obj_factory));
        // = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(m_response_library->getWorksetContainer(),
        //                                                            m_response_library->getGlobalIndexer(),
        //                                                            m_response_library->getLinearObjFactory()));

    // using the FMB, build the model evaluator
    {
      // get nominal input values, make sure they match with internal me
      Thyra::ModelEvaluatorBase::InArgs<ScalarT> nomVals = physics_me->getNominalValues();

      // determine if this is a Epetra or Thyra ME
      Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(physics_me);

      bool useThyra = true;
#ifdef PANZER_HAVE_EPETRA_STACK
      Teuchos::RCP<Thyra::EpetraModelEvaluator> ep_thyra_me = Teuchos::rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(physics_me);
      if(ep_thyra_me!=Teuchos::null)
        useThyra = false;
#endif

      // get parameter names
      std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > p_names(physics_me->Np());
      std::vector<Teuchos::RCP<Teuchos::Array<double> > > p_values(physics_me->Np());
      for(std::size_t i=0;i<p_names.size();i++) {
        p_names[i] = Teuchos::rcp(new Teuchos::Array<std::string>(*physics_me->get_p_names(i)));
        p_values[i] = Teuchos::rcp(new Teuchos::Array<double>(p_names[i]->size(),0.0));
      }

      Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me
          = buildPhysicsModelEvaluator(useThyra,
                                       fmb,
                                       response_library,
                                       m_lin_obj_factory,
                                       p_names,
                                       p_values,
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
                             const std::vector<Teuchos::RCP<Teuchos::Array<double> > > & p_values,
                             const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                             const Teuchos::RCP<panzer::GlobalData> & global_data,
                             bool is_transient,double t_init) const
  {
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyra_me;
    if(!buildThyraME) {
#ifdef PANZER_HAVE_EPETRA_STACK
      Teuchos::RCP<panzer::ModelEvaluator_Epetra> ep_me
          = Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb,rLibrary,lof, p_names,p_values, global_data, is_transient));

      if (is_transient)
        ep_me->set_t_init(t_init);

      // Build Thyra Model Evaluator
      thyra_me = Thyra::epetraModelEvaluator(ep_me,solverFactory);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: buildPhysicsModelEvalautor() - Epetra stack is not enabled!");
#endif
    }
    else {
      thyra_me = Teuchos::rcp(new panzer::ModelEvaluator<double>
                  (fmb,rLibrary,lof,p_names,p_values,solverFactory,global_data,is_transient,t_init));
    }

    return thyra_me;
  }

  template<typename ScalarT>
  double ModelEvaluatorFactory<ScalarT>::
  getInitialTime(Teuchos::ParameterList& p,
                 const panzer_stk::STK_Interface & mesh) const
  {
    Teuchos::ParameterList validPL;
    {
      validPL.set<std::string>("Start Time Type", "From Input File", "Set the start time",
        rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("From Input File","From Exodus File"))));

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
                                      const Teuchos::RCP<const panzer::GlobalIndexer> & ugi,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                      const Teuchos::RCP<panzer_stk::STK_Interface> & mesh) const
  {
     Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary
        = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>(wc,ugi,lof));

     std::vector<std::string> eBlocks;
     mesh->getElementBlockNames(eBlocks);

     panzer_stk::RespFactorySolnWriter_Builder builder;
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
                   const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                   const Teuchos::RCP<panzer::ConnManager> & conn_manager,
                   const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                   const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm
                   #ifdef PANZER_HAVE_TEKO
                   , const Teuchos::RCP<Teko::RequestHandler> & reqHandler
                   #endif
                   ) const
  {
    const Teuchos::ParameterList & p = *this->getParameterList();
    const Teuchos::ParameterList & solncntl_params = p.sublist("Solution Control");

    // Build stratimikos solver (note that this is a hard coded path to linear solver options in nox list!)
    Teuchos::RCP<Teuchos::ParameterList> strat_params
       = Teuchos::rcp(new Teuchos::ParameterList(solncntl_params.sublist("NOX").sublist("Direction").
                      sublist("Newton").sublist("Stratimikos Linear Solver").sublist("Stratimikos")));

    bool writeCoordinates = false;
    if(p.sublist("Options").isType<bool>("Write Coordinates"))
      writeCoordinates = p.sublist("Options").get<bool>("Write Coordinates");

    bool writeTopo = false;
    if(p.sublist("Options").isType<bool>("Write Topology"))
      writeTopo = p.sublist("Options").get<bool>("Write Topology");


    return panzer_stk::buildLOWSFactory(
                            blockedAssembly,globalIndexer,conn_manager,
                            Teuchos::as<int>(mesh->getDimension()), mpi_comm, strat_params,
                            #ifdef PANZER_HAVE_TEKO
                            reqHandler,
                            #endif
                            writeCoordinates,
                            writeTopo
                            );
  }

  template<typename ScalarT>
  void ModelEvaluatorFactory<ScalarT>::
  buildResponses(const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                 const bool write_graphviz_file,
                 const std::string& graphviz_file_prefix)
  {
    Teuchos::ParameterList & p = *this->getNonconstParameterList();
    Teuchos::ParameterList & user_data = p.sublist("User Data");
    Teuchos::ParameterList & closure_models = p.sublist("Closure Models");

    // uninitialize the thyra model evaluator, and rebuild the
    // responses to get the correct response counts!

    using PanzerME = panzer::ModelEvaluator<double>;
    Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(m_physics_me);

    if(nonnull(panzer_me)) {
      panzer_me->buildResponses(m_physics_blocks,*m_eqset_factory,cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix);
      return;
    }
#ifdef PANZER_HAVE_EPETRA_STACK
    else {
      Teuchos::RCP<Thyra::EpetraModelEvaluator> epetra_me = Teuchos::rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(m_physics_me);

      if(epetra_me!=Teuchos::null) {
        Teuchos::RCP<const EpetraExt::ModelEvaluator> const_ep_me;
        Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solveFactory;
        epetra_me->uninitialize(&const_ep_me,&solveFactory);

        Teuchos::RCP<EpetraExt::ModelEvaluator> ep_me = Teuchos::rcp_const_cast<EpetraExt::ModelEvaluator>(const_ep_me);
        Teuchos::RCP<panzer::ModelEvaluator_Epetra> ep_panzer_me = Teuchos::rcp_dynamic_cast<panzer::ModelEvaluator_Epetra>(ep_me);
        ep_panzer_me->buildResponses(m_physics_blocks,*m_eqset_factory,cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix);

        // reinitialize the thyra model evaluator, now with the correct responses
        epetra_me->initialize(ep_me,solveFactory);

        return;
      }
    }
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: buildResponses() - could not cast Physics ME to PanzerME!");
  }
}

#endif
