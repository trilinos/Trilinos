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
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GatherOrientation.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Epetra_MpiComm.h"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

typedef double ScalarT;
typedef int LocalOrdinalT;
typedef int GlobalOrdinalT;
typedef Kokkos::DefaultNode::DefaultNodeType NodeT;

typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
typedef Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> OperatorType;
typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

typedef Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> ThyraLinearOp;

typedef panzer::BlockedTpetraLinearObjContainer<double,int,int> BLOC;
typedef panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,int> BLOFact;

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);

  TEUCHOS_UNIT_TEST(block_assembly, gather_solution)
  {
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

    int myRank = eComm->MyPID();
    int numProcs = eComm->NumProc();

    const std::size_t workset_size = 4/numProcs;
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

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,*physicsBlock); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::BlockedDOFManager<int,int> > dofManager = Teuchos::rcp(new panzer::BlockedDOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

    dofManager->addField(fieldName1_q1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_q1->getIntrepidBasis())));
    dofManager->addField(fieldName2_q1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_q1->getIntrepidBasis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_qedge1->getIntrepidBasis())));

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
       = Teuchos::rcp(new BlockedEpetraLinearObjFactory<panzer::Traits,int>(eComm.getConst(),dofManager));
    Teuchos::RCP<LinearObjFactory<panzer::Traits> > lof = be_lof;
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::X,*loc);

    Teuchos::RCP<BlockedEpetraLinearObjContainer> b_loc 
       = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(loc);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_x());
    Thyra::assign(p_vec->getNonconstVectorBlock(0).ptr(),123.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(1).ptr(),456.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(2).ptr(),789.0+myRank);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<PHX::FieldTag> evalField_q1, evalField_qedge1;
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),2);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),2);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);

    panzer::GlobalEvaluationDataContainer gedc;
    gedc.addDataObject("Solution Gather Container",loc);
    fm.preEvaluate<panzer::Traits::Residual>(gedc);
    fm.preEvaluate<panzer::Traits::Jacobian>(gedc);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.ghostedLinContainer = loc;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 2.0; // derivatives multiplied by 2
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);
    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    // test Residual fields
    {
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData1_q1(fieldName1_q1,basis_q1->functional);
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData2_q1(fieldName2_q1,basis_qedge1->functional);

       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData1_q1);
       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData2_q1);

       TEST_EQUALITY(fieldData1_q1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData1_q1.dimension(1),4);
       TEST_EQUALITY(fieldData2_q1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData2_q1.dimension(1),4);
       TEST_EQUALITY(fieldData1_q1.size(),4*4/numProcs);
       TEST_EQUALITY(fieldData2_q1.size(),4*4/numProcs);
   
       for(int i=0;i<fieldData1_q1.dimension(0);i++) 
          for(int j=0;j<fieldData1_q1.dimension(1);j++) 
             TEST_EQUALITY(fieldData1_q1(i,j),123.0+myRank);

       for(int i=0;i<fieldData2_q1.dimension(0);i++) 
          for(int j=0;j<fieldData2_q1.dimension(1);j++) 
             TEST_EQUALITY(fieldData2_q1(i,j),789.0+myRank);
    }
    {
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData_qedge1(fieldName_qedge1,basis_qedge1->functional);

       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData_qedge1);
 
       TEST_EQUALITY(fieldData_qedge1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData_qedge1.dimension(1),4);
       TEST_EQUALITY(fieldData_qedge1.size(),4*4/numProcs);
   
       for(int i=0;i<fieldData_qedge1.size();i++) 
          TEST_EQUALITY(fieldData_qedge1[i],456.0+myRank);
    }

    // test Jacobian fields
    {
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData1_q1(fieldName1_q1,basis_q1->functional);
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData2_q1(fieldName2_q1,basis_qedge1->functional);
   
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData1_q1);
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData2_q1);
   
       for(int i=0;i<fieldData1_q1.size();i++) {
          TEST_EQUALITY(fieldData1_q1[i],123.0+myRank);
          TEST_EQUALITY(fieldData1_q1[i].availableSize(),12);
       }
       for(int i=0;i<fieldData2_q1.size();i++)  {
          TEST_EQUALITY(fieldData2_q1[i],789.0+myRank);
          TEST_EQUALITY(fieldData2_q1[i].availableSize(),12);
       }
    }
    {
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData_qedge1(fieldName_qedge1,basis_qedge1->functional);
   
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData_qedge1);
   
       for(int i=0;i<fieldData_qedge1.size();i++) {
          TEST_EQUALITY(fieldData_qedge1[i],456.0+myRank);
          TEST_EQUALITY(fieldData_qedge1[i].availableSize(),12);
       }
    }
  }

  TEUCHOS_UNIT_TEST(tpetra_block_assembly, gather_solution)
  {
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

    int myRank = comm->getRank();
    int numProcs = comm->getSize();

    const std::size_t workset_size = 4/numProcs;
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

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,*physicsBlock); 
    TEST_EQUALITY(work_sets->size(),1);

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
    RCP<panzer::BlockedDOFManager<int,int> > dofManager = Teuchos::rcp(new panzer::BlockedDOFManager<int,int>(conn_manager,MPI_COMM_WORLD));

    dofManager->addField(fieldName1_q1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_q1->getIntrepidBasis())));
    dofManager->addField(fieldName2_q1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_q1->getIntrepidBasis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::IntrepidFieldPattern(basis_qedge1->getIntrepidBasis())));

    std::vector<std::vector<std::string> > fieldOrder(3);
    fieldOrder[0].push_back(fieldName1_q1);
    fieldOrder[1].push_back(fieldName_qedge1);
    fieldOrder[2].push_back(fieldName2_q1);
    dofManager->setFieldOrder(fieldOrder);

    // dofManager->setOrientationsRequired(true);
    dofManager->buildGlobalUnknowns();

    // setup linear object factory
    /////////////////////////////////////////////////////////////

    Teuchos::RCP<BLOFact> be_lof 
       = Teuchos::rcp(new BLOFact(comm.getConst(),dofManager));
    Teuchos::RCP<LinearObjFactory<panzer::Traits> > lof = be_lof;
    Teuchos::RCP<LinearObjContainer> loc = be_lof->buildGhostedLinearObjContainer();
    be_lof->initializeGhostedContainer(LinearObjContainer::X,*loc);

    Teuchos::RCP<BLOC> b_loc 
       = Teuchos::rcp_dynamic_cast<BLOC>(loc);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(b_loc->get_x());
    Thyra::assign(p_vec->getNonconstVectorBlock(0).ptr(),123.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(1).ptr(),456.0+myRank);
    Thyra::assign(p_vec->getNonconstVectorBlock(2).ptr(),789.0+myRank);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<PHX::FieldTag> evalField_q1, evalField_qedge1;
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),2);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),2);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
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

       TEST_EQUALITY(evaluator->evaluatedFields().size(),1);

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);

    panzer::GlobalEvaluationDataContainer gedc;
    gedc.addDataObject("Solution Gather Container",loc);
    fm.preEvaluate<panzer::Traits::Residual>(gedc);
    fm.preEvaluate<panzer::Traits::Jacobian>(gedc);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.ghostedLinContainer = loc;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 2.0; // derivatives multiplied by 2
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);
    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    // test Residual fields
    {
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData1_q1(fieldName1_q1,basis_q1->functional);
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData2_q1(fieldName2_q1,basis_qedge1->functional);

       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData1_q1);
       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData2_q1);

       TEST_EQUALITY(fieldData1_q1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData1_q1.dimension(1),4);
       TEST_EQUALITY(fieldData2_q1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData2_q1.dimension(1),4);
       TEST_EQUALITY(fieldData1_q1.size(),4*4/numProcs);
       TEST_EQUALITY(fieldData2_q1.size(),4*4/numProcs);
   
       for(int i=0;i<fieldData1_q1.dimension(0);i++) 
          for(int j=0;j<fieldData1_q1.dimension(1);j++) 
             TEST_EQUALITY(fieldData1_q1(i,j),123.0+myRank);

       for(int i=0;i<fieldData2_q1.dimension(0);i++) 
          for(int j=0;j<fieldData2_q1.dimension(1);j++) 
             TEST_EQUALITY(fieldData2_q1(i,j),789.0+myRank);
    }
    {
       PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData_qedge1(fieldName_qedge1,basis_qedge1->functional);

       fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(fieldData_qedge1);
 
       TEST_EQUALITY(fieldData_qedge1.dimension(0),4/numProcs);
       TEST_EQUALITY(fieldData_qedge1.dimension(1),4);
       TEST_EQUALITY(fieldData_qedge1.size(),4*4/numProcs);
   
       for(int i=0;i<fieldData_qedge1.size();i++) 
          TEST_EQUALITY(fieldData_qedge1[i],456.0+myRank);
    }

    // test Jacobian fields
    {
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData1_q1(fieldName1_q1,basis_q1->functional);
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData2_q1(fieldName2_q1,basis_qedge1->functional);
   
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData1_q1);
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData2_q1);
   
       for(int i=0;i<fieldData1_q1.size();i++) {
          TEST_EQUALITY(fieldData1_q1[i],123.0+myRank);
          TEST_EQUALITY(fieldData1_q1[i].availableSize(),12);
       }
       for(int i=0;i<fieldData2_q1.size();i++)  {
          TEST_EQUALITY(fieldData2_q1[i],789.0+myRank);
          TEST_EQUALITY(fieldData2_q1[i].availableSize(),12);
       }
    }
    {
       PHX::MDField<panzer::Traits::Jacobian::ScalarT,panzer::Cell,panzer::BASIS> 
          fieldData_qedge1(fieldName_qedge1,basis_qedge1->functional);
   
       fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fieldData_qedge1);
   
       for(int i=0;i<fieldData_qedge1.size();i++) {
          TEST_EQUALITY(fieldData_qedge1[i],456.0+myRank);
          TEST_EQUALITY(fieldData_qedge1[i].availableSize(),12);
       }
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
    typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
    typedef panzer_stk::STK_Interface::VectorFieldType CoordinateField;

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
