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

#include <string>
#include <iostream>

#include "Panzer_Traits.hpp"

// for testing gather/scatter construction
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "UnitTest_GlobalIndexer.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

#include "UnitTest_ConnManager.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

Teuchos::RCP<Epetra_CrsMatrix> getSubBlock(int i,int j,Thyra::LinearOpBase<double> & lo)
{
   Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<Epetra_Operator> e_blo = Thyra::get_Epetra_Operator(*blo.getNonconstBlock(i,j));

   return rcp_dynamic_cast<Epetra_CrsMatrix>(e_blo);
}

Teuchos::RCP<const Epetra_CrsMatrix> getSubBlock(int i,int j,const Thyra::LinearOpBase<double> & lo)
{
   const Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<const Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<const Epetra_Operator> e_blo = Thyra::get_Epetra_Operator(*blo.getBlock(i,j));

   return rcp_dynamic_cast<const Epetra_CrsMatrix>(e_blo);
}

template <typename Intrepid2Type>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // build a geometric pattern from a single basis
  RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}

Teuchos::RCP<const panzer::BlockedDOFManager> buildBlockedIndexer(int myRank,int numProc,int numBlocks)
{
  std::string names[] = {"U","V","W","X"};

  Teuchos::RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
  Teuchos::RCP<ConnManager> connManager = rcp(new unit_test::ConnManager(myRank,numProc));
  Teuchos::RCP<panzer::BlockedDOFManager> indexer = rcp(new panzer::BlockedDOFManager());

  indexer->setConnManager(connManager,MPI_COMM_WORLD);

  std::vector<std::vector<std::string> > fieldOrder(numBlocks);
  for(int i=0;i<numBlocks;i++) {
    indexer->addField(names[i],patternC1);

    fieldOrder[i].push_back(names[i]);
  }
  indexer->setFieldOrder(fieldOrder);
  indexer->buildGlobalUnknowns();

  return indexer;
}

RCP<Epetra_MultiVector> getEpetraMultiVector(RCP<Thyra::MultiVectorBase<double> > & vec,const Epetra_Map & eMap)
{
   return Thyra::get_Epetra_MultiVector(eMap,vec);
}

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, gather_scatter_constr)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::GlobalIndexer> indexer 
         = rcp(new panzer::unit_test::GlobalIndexer(myRank,numProc));
 
   // setup factory
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),indexer));

   // build parameter lists for gather and scatters
   //////////////////////////////////////////////////////////////////
   std::size_t numCells = 10;

   Teuchos::ParameterList gatherParams;
   {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      // auxiliary information needed to construct basis object
      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,topo);

      // build DOF names
      RCP<std::vector<std::string> > dofNames = rcp(new std::vector<std::string>);
      dofNames->push_back("ux"); // in practice these probably would not be gathered together!
      dofNames->push_back("p");

      // build basis
      RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,1,cellData));

      // build gather parameter list
      gatherParams.set<RCP<std::vector<std::string> > >("DOF Names",dofNames);
      gatherParams.set<RCP<std::vector<std::string> > >("Indexer Names",dofNames);
      gatherParams.set<RCP<panzer::PureBasis> >("Basis",basis);
   }

   Teuchos::ParameterList scatterParams;
   {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,topo);
   
      // build basis
      RCP<const panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,1,cellData));
   
      std::string scatterName = "Residual_NS";
   
      // build DOF names
      RCP<std::vector<std::string> > evaluatedNames = rcp(new std::vector<std::string>);
      evaluatedNames->push_back("Residual_ux"); // in practice these probably would not be scattered together!
      evaluatedNames->push_back("Residual_p");
   
      // build evaluated map
      RCP<std::map<std::string,std::string> > evaluatedMap = rcp(new std::map<std::string,std::string>);
      evaluatedMap->insert(std::make_pair("Residual_ux","ux")); // in practice these probably would not be scattered together!
      evaluatedMap->insert(std::make_pair("Residual_p","p"));

      // build scatter parameter list
      scatterParams.set<std::string>("Scatter Name",scatterName);
      scatterParams.set<RCP<std::vector<std::string> > >("Dependent Names",evaluatedNames);
      scatterParams.set<RCP<std::map<std::string,std::string> > >("Dependent Map",evaluatedMap);
      scatterParams.set("Basis",basis);
   }

   Teuchos::ParameterList scatterDirichletParams;
   {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,topo);
   
      // build basis
      RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,1,cellData));
   
      std::string scatterName = "Residual_NS";
   
      // build DOF names
      RCP<std::vector<std::string> > evaluatedNames = rcp(new std::vector<std::string>);
      evaluatedNames->push_back("Residual_ux"); // in practice these probably would not be scattered together!
      evaluatedNames->push_back("Residual_p");
   
      // build evaluated map
      RCP<std::map<std::string,std::string> > evaluatedMap = rcp(new std::map<std::string,std::string>);
      evaluatedMap->insert(std::make_pair("Residual_ux","ux")); // in practice these probably would not be scattered together!
      evaluatedMap->insert(std::make_pair("Residual_p","p"));
   
      // build scatter parameter list
      scatterDirichletParams.set<std::string>("Scatter Name",scatterName);
      scatterDirichletParams.set<RCP<std::vector<std::string> > >("Dependent Names",evaluatedNames);
      scatterDirichletParams.set<RCP<std::map<std::string,std::string> > >("Dependent Map",evaluatedMap);
      scatterDirichletParams.set<RCP<panzer::PureBasis> >("Basis",basis);
      scatterDirichletParams.set<int>("Side Subcell Dimension",1);
      scatterDirichletParams.set<int>("Local Side ID",2);
      scatterDirichletParams.set("Check Apply BC",false);
   }

   // Evaluator construction tests
   ///////////////////////////////////////////////////////
   RCP<PHX::Evaluator<panzer::Traits> > evaluator;

   {
      typedef panzer::Traits::Residual EvalType;

      // gather test
      {
         evaluator = la_factory->buildGather<EvalType>(gatherParams);

         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<GatherSolution_Epetra<EvalType,panzer::Traits,int,int> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(gatherSolutionEval!=Teuchos::null);
   
         const std::vector<RCP<PHX::FieldTag> > & fields = gatherSolutionEval->evaluatedFields();
         TEST_EQUALITY(fields.size(),2);
   
         TEST_EQUALITY(fields[0]->name(),"ux");
         TEST_EQUALITY(fields[1]->name(),"p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter test
      {
         evaluator = la_factory->buildScatter<EvalType>(scatterParams);

         out << "SCATTER RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterResidual_Epetra<EvalType,panzer::Traits,int,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter dirichlet test
      {
         evaluator = la_factory->buildScatterDirichlet<EvalType>(scatterDirichletParams);

         out << "SCATTER DIRICHLET RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,int,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

   }

   {
      typedef panzer::Traits::Jacobian EvalType;

      // gather test
      {
         evaluator = la_factory->buildGather<EvalType>(gatherParams);

         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<GatherSolution_Epetra<EvalType,panzer::Traits,int,int> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(gatherSolutionEval!=Teuchos::null);
   
         const std::vector<RCP<PHX::FieldTag> > & fields = gatherSolutionEval->evaluatedFields();
         TEST_EQUALITY(fields.size(),2);
   
         TEST_EQUALITY(fields[0]->name(),"ux");
         TEST_EQUALITY(fields[1]->name(),"p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter test
      {
         evaluator = la_factory->buildScatter<EvalType>(scatterParams);

         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterResidual_Epetra<EvalType,panzer::Traits,int,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1

         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter dirichlet test
      {
         evaluator = la_factory->buildScatterDirichlet<EvalType>(scatterDirichletParams);

         out << "SCATTER DIRICHLET RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,int,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,int,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().extent_int(1),Teuchos::as<int>(4)); // for Q1
      }

   }
}

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, initializeContainer)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef EpetraLinearObjContainer ELOC;

   RCP<panzer::GlobalIndexer> indexer 
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));

   std::vector<panzer::GlobalOrdinal> ownedIndices, ownedAndGhostedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndGhostedIndices(ownedAndGhostedIndices);
 
   // setup factory
   Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm.getConst(),indexer));

   RCP<LinearObjContainer> container = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghostedContainer = la_factory->buildGhostedLinearObjContainer();

   RCP<EpetraLinearObjContainer> eContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(container);
   RCP<EpetraLinearObjContainer> eGhostedContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(ghostedContainer);

   // tests ghosted initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeGhostedContainer(ELOC::X,*ghostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_x()->MyLength(),(int) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::DxDt,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt()->MyLength(),(int) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::F,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f()->MyLength(),(int) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::Mat,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_A()->NumMyRows(),(int) ownedAndGhostedIndices.size());
   
      // jacobian and residual vector output
      la_factory->initializeGhostedContainer(ELOC::F | ELOC::Mat,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeGhostedContainer(ELOC::X | ELOC::DxDt,*ghostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeGhostedContainer(ELOC::X | ELOC::DxDt | ELOC::F | ELOC::Mat,*ghostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
   
      // Epetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeGhostedContainer(ELOC::X,*eGhostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(ELOC::DxDt,*eGhostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(ELOC::F,*eGhostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(ELOC::Mat,*eGhostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      la_factory->initializeGhostedContainer(ELOC::F | ELOC::Mat,*eGhostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeGhostedContainer(ELOC::X | ELOC::DxDt,*eGhostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeGhostedContainer(ELOC::X | ELOC::DxDt | ELOC::F | ELOC::Mat,*eGhostedContainer);
      TEST_ASSERT(eGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
   }

}

TEUCHOS_UNIT_TEST(tBlockedLinearObjFactory, intializeContainer_epetra)
{

   panzer::BlockedEpetraLinearObjContainer container;

   TEST_ASSERT(container.checkCompatibility());
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, epetra_factory_tests)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;

   int numBlocks = 3;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   BlockedEpetraLinearObjFactory<panzer::Traits,int> factory(tComm,blkIndexer);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);

   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // Epetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   }
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, ghostToGlobal)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   out << "Built indexer = " << blkIndexer << std::endl;
   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm,blkIndexer));

   out << "LA Factory = " << la_factory << std::endl;

   Teuchos::RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
   Teuchos::RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();

   out << "LOC = " << global << " " << ghosted << std::endl;

   la_factory->initializeContainer(LinearObjContainer::Mat,*global);
   la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);

   out << "Initialized matrices" << std::endl;

   Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(ghosted)->initializeMatrix(1.0);

   out << "Initialize matrix" << std::endl;

   la_factory->ghostToGlobalContainer(*ghosted,*global,LinearObjContainer::Mat);

   out << "GhostToGlobal" << std::endl;

   Teuchos::RCP<Thyra::LinearOpBase<double> > th_A = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(global)->get_A();
   Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blk_A = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(th_A);
   
   out << "Blk A " << blk_A << std::endl;
   out << "Blk A " << Teuchos::describe(*blk_A,Teuchos::VERB_MEDIUM) << std::endl;
  
   Teuchos::RCP<Epetra_Operator> cA_00 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(0,0));
   out << "00" << std::endl;

   Teuchos::RCP<Epetra_Operator> cA_01 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(0,1));
   out << "01" << std::endl;

   Teuchos::RCP<Epetra_Operator> cA_10 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(1,0));
   out << "10" << std::endl;

   Teuchos::RCP<Epetra_Operator> cA_11 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(1,1));
   out << "11" << std::endl;

   Teuchos::RCP<Epetra_CrsMatrix> A_00 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_00);
   Teuchos::RCP<Epetra_CrsMatrix> A_01 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_01);

   Teuchos::RCP<Epetra_CrsMatrix> A_10 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_10);
   Teuchos::RCP<Epetra_CrsMatrix> A_11 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_11);

   out << "A_00 = \n";
   A_00->Print(out);

   out << "A_01 = \n";
   A_01->Print(out);

   out << "A_10 = \n";
   A_10->Print(out);

   out << "A_11 = \n";
   A_11->Print(out);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, graph_constr)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm,blkIndexer));

   Teuchos::RCP<Epetra_CrsMatrix> A_00 = la_factory->getGhostedEpetraMatrix(0,0); A_00->PutScalar(1.0);
   Teuchos::RCP<Epetra_CrsMatrix> A_01 = la_factory->getGhostedEpetraMatrix(0,1); A_01->PutScalar(1.0);

   Teuchos::RCP<Epetra_CrsMatrix> A_10 = la_factory->getGhostedEpetraMatrix(1,0); A_10->PutScalar(1.0);
   Teuchos::RCP<Epetra_CrsMatrix> A_11 = la_factory->getGhostedEpetraMatrix(1,1); A_11->PutScalar(1.0);

   std::cout << "A_00 = \n";
   A_00->Print(out);

   std::cout << "A_01 = \n";
   A_01->Print(out);

   std::cout << "A_10 = \n";
   A_10->Print(out);

   std::cout << "A_11 = \n";
   A_11->Print(out);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, adjustDirichlet)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 3;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm,blkIndexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   // b_sys->get_A()->PutScalar(-3.0);
   for(int i=0;i<numBlocks;i++)
      for(int j=0;j<numBlocks;j++)
         getSubBlock(i,j,*b_sys->get_A())->PutScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[0] = 1.0; // GID = 0
         data_1[0] = 1.0;
   
         // case 1
         data_0[2] = 1.0; // GID = 2
         data_1[2] = 2.0;
   
         // case 2
         data_1[5] = 2.0; // GID = 5
      }
   }
   else if(myRank==1) {
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[3] = 1.0; // GID = 9
         data_1[3] = 1.0;
   
         // case 1
         data_0[0] = 1.0; // GID =2
         data_1[0] = 2.0;
   
         // case 2
         data_1[6] = 2.0; // GID = 4
      }
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   int numEntries = 0;
   double * values = 0;
   int * indices = 0;

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);      // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->ExtractMyRowView(2,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->ExtractMyRowView(5,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(3,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->ExtractMyRowView(6,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, node_cell)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(tComm,blkIndexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   // b_sys->get_A()->PutScalar(-3.0);
   for(int i=0;i<numBlocks;i++)
      for(int j=0;j<numBlocks;j++)
         getSubBlock(i,j,*b_sys->get_A())->PutScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[0] = 1.0; // GID = 0
      data_1[0] = 1.0;

      // case 1
      data_0[2] = 1.0; // GID = 2
      data_1[2] = 2.0;

      // case 2
      data_1[5] = 2.0; // GID = 5

      {
         x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(1);
         x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(1);

         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         data_1[0] = 2.0;
      }
   }
   else if(myRank==1) {
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[3] = 1.0; // GID = 9
      data_1[3] = 1.0;

      // case 1
      data_0[0] = 1.0; // GID =2
      data_1[0] = 2.0;

      // case 2
      data_1[6] = 2.0; // GID = 4
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   int numEntries = 0;
   double * values = 0;
   int * indices = 0;

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->ExtractMyRowView(2,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->ExtractMyRowView(5,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],0.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(3,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->ExtractMyRowView(6,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, exclusion)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;

   int numBlocks = 3;
   int myRank = tComm->getRank();
   int numProc = tComm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer(myRank,numProc,numBlocks);

   BlockedEpetraLinearObjFactory<panzer::Traits,int> factory(tComm,blkIndexer);
 
   // exclude some pairs
   std::vector<std::pair<int,int> > exPairs;
   exPairs.push_back(std::make_pair(0,2));
   exPairs.push_back(std::make_pair(2,1));
   factory.addExcludedPairs(exPairs);
   factory.addExcludedPair(1,2);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);

   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
      factory.initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);

      RCP<Thyra::BlockedLinearOpBase<double> > blo 
         = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(bContainer->get_A());

      TEST_ASSERT(!blo->getNonconstBlock(0,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(0,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(0,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(1,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(1,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(1,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(2,0).is_null());
      TEST_ASSERT(blo->getNonconstBlock(2,1).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(2,2).is_null());
   }
}

}
