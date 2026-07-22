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

#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_ScatterResidual_BlockedEpetra.hpp"
#include "Panzer_GatherSolution_BlockedEpetra.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_CrsMatrix.h"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

  TEUCHOS_UNIT_TEST(block_assembly, scatter_dirichlet_residual)
  {

   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

    int myRank = tComm->getRank();

    const std::size_t workset_size = 4;
    const std::string fieldName1_q1 = "U";
    const std::string fieldName2_q1 = "V";
    const std::string fieldName_qedge1 = "B";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    Teuchos::RCP<panzer::PureBasis> basis_qedge1 = buildBasis(workset_size,"QEdge1");

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb);

    const int default_int_order = 1;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                            physicsBlock->getWorksetNeeds()); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::BlockedDOFManager> dofManager = Teuchos::rcp(new panzer::BlockedDOFManager(conn_manager,MPI_COMM_WORLD));

    dofManager->addField(fieldName1_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName2_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_qedge1->getIntrepid2Basis())));

    std::vector<std::vector<std::string> > fieldOrder(3);
    fieldOrder[0].push_back(fieldName1_q1);
    fieldOrder[1].push_back(fieldName_qedge1);
    fieldOrder[2].push_back(fieldName2_q1);
    dofManager->setFieldOrder(fieldOrder);

    // dofManager->setOrientationsRequired(true);
    dofManager->buildGlobalUnknowns();

    // setup linear object factory
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > be_lof 
       = Teuchos::rcp(new BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
    Teuchos::RCP<LinearObjFactory<panzer::Traits> > lof = be_lof;
    Teuchos::RCP<LinearObjContainer> dd_loc = be_lof->buildGhostedLinearObjContainer();
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::F,*dd_loc);
    dd_loc->initialize();

    be_lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::F,*loc);
    loc->initialize();

    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_dd_loc = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(dd_loc);
    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_loc = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(loc);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_x());
    Thyra::assign(p_vec->getNonconstVectorBlock(0).ptr(),123.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(1).ptr(),456.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(2).ptr(),789.0+myRank);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    std::string resName = "";
    Teuchos::RCP<std::map<std::string,std::string> > names_map =
       Teuchos::rcp(new std::map<std::string,std::string>);
    names_map->insert(std::make_pair(fieldName1_q1,resName+fieldName1_q1));
    names_map->insert(std::make_pair(fieldName2_q1,resName+fieldName2_q1));
    names_map->insert(std::make_pair(fieldName_qedge1,resName+fieldName_qedge1));

    // evaluators under test
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(resName+fieldName1_q1);
       names->push_back(resName+fieldName2_q1);

       Teuchos::ParameterList pl; 
       pl.set("Scatter Name", "ScatterQ1");
       pl.set("Basis", basis_q1);
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);
       pl.set("Side Subcell Dimension", 1);
       pl.set("Local Side ID", 2);
       pl.set("Check Apply BC", false);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatterDirichlet<panzer::Traits::Residual>(pl);

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(resName+fieldName_qedge1);

       Teuchos::ParameterList pl; 
       pl.set("Scatter Name", "ScatterQEdge1");
       pl.set("Basis", basis_qedge1);
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);
       pl.set("Side Subcell Dimension", 1);
       pl.set("Local Side ID", 2);
       pl.set("Check Apply BC", false);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatterDirichlet<panzer::Traits::Residual>(pl);

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

    // support evaluators
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(fieldName1_q1);
       names->push_back(fieldName2_q1);

       Teuchos::ParameterList pl; 
       pl.set("Basis", basis_q1);
       pl.set("DOF Names",names);
       pl.set("Indexer Names",names);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildGather<panzer::Traits::Residual>(pl);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
    }
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(fieldName_qedge1);

       Teuchos::ParameterList pl; 
       pl.set("Basis", basis_qedge1);
       pl.set("DOF Names",names);
       pl.set("Indexer Names",names);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildGather<panzer::Traits::Residual>(pl);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(12);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    fm.postRegistrationSetup(sd);

    // panzer::Traits::PED ped;
    // ped.dirichletData.ghostedCounter = dd_loc;
    // fm.preEvaluate<panzer::Traits::Residual>(ped);
    panzer::Traits::PED ped;
    ped.gedc->addDataObject("Dirichlet Counter",dd_loc);
    ped.gedc->addDataObject("Solution Gather Container",loc);
    ped.gedc->addDataObject("Residual Scatter Container",loc);
    fm.preEvaluate<panzer::Traits::Residual>(ped);


    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.alpha = 0.0;
    workset.beta = 2.0; // derivatives multiplied by 2
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);

    // test Residual fields
    panzer::index_t dd_count(0);
    Teuchos::ArrayRCP<const double> data, dd_data;
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > f_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_f());
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > dd_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_dd_loc->get_f());

    // check all the residual values. This is kind of crappy test since it simply checks twice the target
    // value and the target. Its this way because you add two entries across elements.

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(0)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 123.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,2*workset.num_cells); // there are 2 nodes on the side and the sides are not shared

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(1)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 456.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,workset.num_cells); // there are 2 nodes on the side and the sides are not shared

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(2)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 789.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,2*workset.num_cells); // there are 2 nodes on the side and the sides are not shared

  }

  TEUCHOS_UNIT_TEST(block_assembly, scatter_dirichlet_jacobian)
  {
    
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

    int myRank = tComm->getRank();

    const std::size_t workset_size = 4;
    const std::string fieldName1_q1 = "U";
    const std::string fieldName2_q1 = "V";
    const std::string fieldName_qedge1 = "B";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    Teuchos::RCP<panzer::PureBasis> basis_qedge1 = buildBasis(workset_size,"QEdge1");

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb);

    const int default_int_order = 1;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                            physicsBlock->getWorksetNeeds()); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager> conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::BlockedDOFManager> dofManager = Teuchos::rcp(new panzer::BlockedDOFManager(conn_manager,MPI_COMM_WORLD));

    dofManager->addField(fieldName1_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName2_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_qedge1->getIntrepid2Basis())));

    std::vector<std::vector<std::string> > fieldOrder(3);
    fieldOrder[0].push_back(fieldName1_q1);
    fieldOrder[1].push_back(fieldName_qedge1);
    fieldOrder[2].push_back(fieldName2_q1);
    dofManager->setFieldOrder(fieldOrder);

    // dofManager->setOrientationsRequired(true);
    dofManager->buildGlobalUnknowns();

    // setup linear object factory
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > be_lof 
       = Teuchos::rcp(new BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),dofManager));
    Teuchos::RCP<LinearObjFactory<panzer::Traits> > lof = be_lof;
    Teuchos::RCP<LinearObjContainer> dd_loc = be_lof->buildGhostedLinearObjContainer();
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::F,*dd_loc);
    dd_loc->initialize();

    be_lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::F | LinearObjContainer::Mat,*loc);
    loc->initialize();

    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_dd_loc = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(dd_loc);
    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_loc = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(loc);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_x());
    Thyra::assign(p_vec->getNonconstVectorBlock(0).ptr(),123.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(1).ptr(),456.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(2).ptr(),789.0+myRank);

    Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blk_A 
       = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(b_loc->get_A());
    double values[] = {123.0+myRank,456.0+myRank,789.0+myRank};
    for(int i=0;i<3;i++) 
       for(int j=0;j<3;j++) 
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(i,j)))->PutScalar(values[i]*values[j]);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    std::string resName = "";
    Teuchos::RCP<std::map<std::string,std::string> > names_map =
       Teuchos::rcp(new std::map<std::string,std::string>);
    names_map->insert(std::make_pair(fieldName1_q1,resName+fieldName1_q1));
    names_map->insert(std::make_pair(fieldName2_q1,resName+fieldName2_q1));
    names_map->insert(std::make_pair(fieldName_qedge1,resName+fieldName_qedge1));

    // evaluators under test
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(resName+fieldName1_q1);
       names->push_back(resName+fieldName2_q1);

       Teuchos::ParameterList pl; 
       pl.set("Scatter Name", "ScatterQ1");
       pl.set("Basis", basis_q1);
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);
       pl.set("Side Subcell Dimension", 1);
       pl.set("Local Side ID", 2);
       pl.set("Check Apply BC", false);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatterDirichlet<panzer::Traits::Jacobian>(pl);

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(resName+fieldName_qedge1);

       Teuchos::ParameterList pl; 
       pl.set("Scatter Name", "ScatterQEdge1");
       pl.set("Basis", basis_qedge1);
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);
       pl.set("Side Subcell Dimension", 1);
       pl.set("Local Side ID", 2);
       pl.set("Check Apply BC", false);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatterDirichlet<panzer::Traits::Jacobian>(pl);

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    // support evaluators
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(fieldName1_q1);
       names->push_back(fieldName2_q1);

       Teuchos::ParameterList pl; 
       pl.set("Basis", basis_q1);
       pl.set("DOF Names",names);
       pl.set("Indexer Names",names);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildGather<panzer::Traits::Jacobian>(pl);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       using Teuchos::RCP;
       using Teuchos::rcp;
       RCP<std::vector<std::string> > names = rcp(new std::vector<std::string>);
       names->push_back(fieldName_qedge1);

       Teuchos::ParameterList pl; 
       pl.set("Basis", basis_qedge1);
       pl.set("DOF Names",names);
       pl.set("Indexer Names",names);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildGather<panzer::Traits::Jacobian>(pl);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(12);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    fm.postRegistrationSetup(sd);

    // panzer::Traits::PED ped;
    // ped.dirichletData.ghostedCounter = dd_loc;
    // fm.preEvaluate<panzer::Traits::Jacobian>(ped);
    panzer::Traits::PED ped;
    ped.gedc->addDataObject("Dirichlet Counter",dd_loc);
    ped.gedc->addDataObject("Solution Gather Container",loc);
    ped.gedc->addDataObject("Residual Scatter Container",loc);
    fm.preEvaluate<panzer::Traits::Jacobian>(ped);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.alpha = 0.0;
    workset.beta = 2.0; // derivatives multiplied by 2
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    // test Residual fields
    panzer::index_t dd_count(0);
    Teuchos::ArrayRCP<const double> data, dd_data;
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > f_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_f());
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > dd_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_dd_loc->get_f());

    // check all the residual values. This is kind of crappy test since it simply checks twice the target
    // value and the target. Its this way because you add two entries across elements.

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(0)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 123.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,2*workset.num_cells); // there are 2 nodes on the side and the sides are not shared

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(1)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 456.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,workset.num_cells); // there are 2 nodes on the side and the sides are not shared

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(data));
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(dd_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(dd_data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(2)->NumMyElements());
    TEST_EQUALITY(data.size(),dd_data.size());
    dd_count = 0;
    for(int i=0;i<data.size();i++) {
 
       double target = 789.0+myRank;
       if(dd_data[i]==0.0)
       {  TEST_EQUALITY(data[i],0.0); }
       else
       {  TEST_EQUALITY(data[i],target); dd_count++; }
    }
    TEST_EQUALITY(dd_count,2*workset.num_cells); // there are 2 nodes on the side and the sides are not shared

  }

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName)
  { 
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);
     return Teuchos::rcp(new panzer::PureBasis(basisName,1,cellData)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",elemX);
    pl->set("Y Elements",elemY);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

    return mesh;
  }

  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb)
  {
    // Physics block
    ipb->setName("test physics");
    {
      Teuchos::ParameterList& p = ipb->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = ipb->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HCurl");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }
    
  }

}
