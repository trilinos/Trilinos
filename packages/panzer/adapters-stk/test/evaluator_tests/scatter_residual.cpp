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

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
  typedef Teuchos::ArrayRCP<const double>::size_type size_type;

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

  TEUCHOS_UNIT_TEST(block_assembly, scatter_solution_residual)
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
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));
    RCP<panzer::BlockedDOFManager<int,int> > dofManager = Teuchos::rcp(new panzer::BlockedDOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

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
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::F,*loc);
    loc->initialize();

    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_loc 
       = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(loc);
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
       pl.set("Basis", basis_q1.getConst());
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatter<panzer::Traits::Residual>(pl);

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
       pl.set("Basis", basis_qedge1.getConst());
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatter<panzer::Traits::Residual>(pl);

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

    panzer::Traits::PED ped;
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
    Teuchos::ArrayRCP<const double> data;
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > f_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_f());

    // check all the residual values. This is kind of crappy test since it simply checks twice the target
    // value and the target. Its this way because you add two entries across elements.
    
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(0)->NumMyElements());
    for(size_type i=0;i<data.size();i++) {
       double target = 123.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target);
    }

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(1)->NumMyElements());
    for(size_type i=0;i<data.size();i++) {
       double target = 456.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target);
    }

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(2)->NumMyElements());
    for(size_type i=0;i<data.size();i++) {
       double target = 789.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target);
    }
  }

  TEUCHOS_UNIT_TEST(block_assembly, scatter_solution_jacobian)
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
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));
    RCP<panzer::BlockedDOFManager<int,int> > dofManager = Teuchos::rcp(new panzer::BlockedDOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

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
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::X | LinearObjContainer::F | LinearObjContainer::Mat,*loc);
    loc->initialize();

    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_loc 
       = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(loc);

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
       pl.set("Basis", basis_q1.getConst());
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatter<panzer::Traits::Jacobian>(pl);

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
       pl.set("Basis", basis_qedge1.getConst());
       pl.set("Dependent Names", names);
       pl.set("Dependent Map", names_map);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator = lof->buildScatter<panzer::Traits::Jacobian>(pl);

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

    panzer::Traits::PED ped;
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

    // test Jacobian fields
    Teuchos::ArrayRCP<const double> data;
    Teuchos::RCP<const Thyra::ProductVectorBase<double> > f_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_f());

    // check all the residual values. This is kind of crappy test since it simply checks twice the target
    // value and the target. Its this way because you add two entries across elements.
    
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(0)->NumMyElements());
    out << std::endl;
    for(size_type i=0;i<data.size();i++) {
       double target = 123.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target);
       out << data[i] << std::endl;
    }

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(1)->NumMyElements());
    out << std::endl;
    for(size_type i=0;i<data.size();i++) {
       double target = 456.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target); 
       out << data[i] << std::endl;
    }

    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_vec->getVectorBlock(2))->getLocalData(Teuchos::ptrFromRef(data));
    TEST_EQUALITY(data.size(),b_loc->getMapForBlock(2)->NumMyElements());
    out << std::endl;
    for(size_type i=0;i<data.size();i++) {
       double target = 789.0+myRank;
       TEST_ASSERT(data[i]==target || data[i]==2.0*target);
       out << data[i] << std::endl;
    }
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
