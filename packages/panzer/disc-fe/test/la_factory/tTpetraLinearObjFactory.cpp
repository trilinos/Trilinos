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
#include <set>
#include <algorithm>

#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

// for testing gather/scatter construction
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_GlobalIndexer.hpp"

#ifdef HAVE_MPI
   #include "Teuchos_DefaultMpiComm.hpp"
#else
   #include "NO_SERIAL_BUILD.h"
#endif

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Tpetra::MultiVector<double,int,panzer::GlobalOrdinal> MultiVector;
typedef Tpetra::Vector<double,int,panzer::GlobalOrdinal> Vector;
typedef Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal> CrsMatrix;
typedef Tpetra::CrsGraph<int,panzer::GlobalOrdinal> CrsGraph;
typedef Tpetra::Map<int,panzer::GlobalOrdinal> Map;

namespace panzer {

/*
RCP<MultiVector> getTpetraMultiVector(RCP<Thyra::MultiVectorBase<double> > & vec,const Map & map)
{
   return Thyra::get_Tpetra_MultiVector<double,int,int>(map,vec);
}
*/

TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, gather_scatter_constr)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::Comm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Teuchos::Comm<int> > failure_comm = THIS_,_SERIAL_BUILDS_,_SHOULD_FAIL;
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();

   // panzer::pauseToAttach();

   RCP<panzer::GlobalIndexer> indexer 
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));
 
   // setup factory
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > la_factory
         = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(tComm.getConst(),indexer));

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
         RCP<GatherSolution_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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
         RCP<ScatterResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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
         RCP<ScatterDirichletResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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
         RCP<GatherSolution_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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
         RCP<ScatterResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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
         RCP<ScatterDirichletResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Tpetra<EvalType,panzer::Traits,int,panzer::GlobalOrdinal> >(evaluator);
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

TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, adjustDirichlet)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::Comm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Teuchos::Comm<int> > failure_comm = THIS_,_SERIAL_BUILDS_,_SHOULD_FAIL;
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   typedef TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOC;

   RCP<panzer::GlobalIndexer> indexer 
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));

   // setup factory
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > la_factory
         = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(tComm.getConst(),indexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<LOC> t_0   = rcp_dynamic_cast<LOC>(ghosted_0);
   RCP<LOC> t_1   = rcp_dynamic_cast<LOC>(ghosted_1);
   RCP<LOC> t_sys = rcp_dynamic_cast<LOC>(ghosted_sys);

   Teuchos::ArrayRCP<double> x_0_a = t_0->get_f()->get1dViewNonConst();
   Teuchos::ArrayRCP<double> x_1_a = t_1->get_f()->get1dViewNonConst();
   Teuchos::ArrayRCP<double> f_a = t_sys->get_f()->get1dViewNonConst();

   TEST_ASSERT(!Teuchos::is_null(t_0->get_f()));
   TEST_ASSERT(!Teuchos::is_null(t_1->get_f()));
   TEST_ASSERT(!Teuchos::is_null(t_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(t_sys->get_A()));

   t_sys->get_f()->putScalar(-3.0); // put some garbage in the systems
   t_sys->get_A()->setAllToScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      // case 0
      x_0_a[0] = 1.0; // GID = 0
      x_1_a[0] = 1.0; // GID = 0

      // case 1
      x_0_a[2] = 1.0; // GID = 2
      x_1_a[2] = 2.0; // GID = 2

      // case 2
      x_1_a[5] = 2.0; // GID = 5
   }
   else if(myRank==1) {
      // case 0
      x_0_a[3] = 1.0; // GID = 9
      x_1_a[3] = 1.0; // GID = 9

      // case 1
      x_0_a[0] = 1.0; // GID = 2
      x_1_a[0] = 2.0; // GID = 2

      // case 2
      x_1_a[6] = 2.0; // GID = 4
   }
   else 
      TEUCHOS_ASSERT(false);

   // run test for conditions
   la_factory->beginFill(*ghosted_sys);
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);
   la_factory->endFill(*ghosted_sys);

   std::size_t sz = t_sys->get_A()->getLocalMaxNumRowEntries();
   std::size_t numEntries = 0;
   typename CrsMatrix::nonconst_local_inds_host_view_type indices("indices", sz);
   typename CrsMatrix::nonconst_values_host_view_type values("values", sz);

   if(myRank==0) {   
      TEST_EQUALITY(f_a[0],-3.0);     // case 0
      t_sys->get_A()->getLocalRowCopy(0,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0);

      TEST_EQUALITY(f_a[2],-3.0/2.0); // case 1
      t_sys->get_A()->getLocalRowCopy(2,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0/2.0);

      TEST_EQUALITY(f_a[5],0.0);      // case 2
      t_sys->get_A()->getLocalRowCopy(5,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],0.0);
   }
   else if(myRank==1) {
      TEST_EQUALITY(f_a[3],-3.0);     // case 0
      t_sys->get_A()->getLocalRowCopy(3,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0);

      TEST_EQUALITY(f_a[0],-3.0/2.0); // case 1
      t_sys->get_A()->getLocalRowCopy(0,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0/2.0);

      TEST_EQUALITY(f_a[6],0.0);      // case 2
      t_sys->get_A()->getLocalRowCopy(6,indices,values,numEntries);
      for(std::size_t i=0;i<numEntries;i++) TEST_EQUALITY(values[i],0.0);
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, initializeContainer)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::Comm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Teuchos::Comm<int> > failure_comm = THIS_,_SERIAL_BUILDS_,_SHOULD_FAIL;
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();
 
   typedef TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOC;

   RCP<panzer::GlobalIndexer> indexer 
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));

   std::vector<panzer::GlobalOrdinal> ownedIndices, ownedAndGhostedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndGhostedIndices(ownedAndGhostedIndices);
 
   // setup factory
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > la_factory
         = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(tComm.getConst(),indexer));

   RCP<LinearObjContainer> container = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghostedContainer = la_factory->buildGhostedLinearObjContainer();

   RCP<TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> > tContainer = rcp_dynamic_cast<TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> >(container);
   RCP<TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> > tGhostedContainer = rcp_dynamic_cast<TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> >(ghostedContainer);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeContainer(LOC::X,*container);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tContainer->get_x()->getLocalLength(),(std::size_t) ownedIndices.size());
   
      la_factory->initializeContainer(LOC::DxDt,*container);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tContainer->get_dxdt()->getLocalLength(),(std::size_t) ownedIndices.size());
   
      la_factory->initializeContainer(LOC::F,*container);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tContainer->get_f()->getLocalLength(),(std::size_t) ownedIndices.size());
   
      la_factory->initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_A()->getLocalNumRows(),(std::size_t) ownedIndices.size());
   
      // jacobian and residual vector output
      la_factory->initializeContainer(LOC::F | LOC::Mat,*container);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeContainer(LOC::X | LOC::DxDt,*container);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*container);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
   
      // Tpetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeContainer(LOC::X,*tContainer);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(LOC::DxDt,*tContainer);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(LOC::F,*tContainer);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(LOC::Mat,*tContainer);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      la_factory->initializeContainer(LOC::F | LOC::Mat,*tContainer);
      TEST_EQUALITY(tContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeContainer(LOC::X | LOC::DxDt,*tContainer);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*tContainer);
      TEST_ASSERT(tContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tContainer->get_A()!=Teuchos::null);
   }

   // tests ghosted initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeGhostedContainer(LOC::X,*ghostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tGhostedContainer->get_x()->getLocalLength(),(std::size_t) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(LOC::DxDt,*ghostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tGhostedContainer->get_dxdt()->getLocalLength(),(std::size_t) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(LOC::F,*ghostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
	TEST_EQUALITY(tGhostedContainer->get_f()->getLocalLength(),(std::size_t) ownedAndGhostedIndices.size());
   
      la_factory->initializeGhostedContainer(LOC::Mat,*ghostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_A()->getLocalNumRows(),(std::size_t) ownedAndGhostedIndices.size());
   
      // jacobian and residual vector output
      la_factory->initializeGhostedContainer(LOC::F | LOC::Mat,*ghostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeGhostedContainer(LOC::X | LOC::DxDt,*ghostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeGhostedContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*ghostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
   
      // Tpetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeGhostedContainer(LOC::X,*tGhostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(LOC::DxDt,*tGhostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(LOC::F,*tGhostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeGhostedContainer(LOC::Mat,*tGhostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      la_factory->initializeGhostedContainer(LOC::F | LOC::Mat,*tGhostedContainer);
      TEST_EQUALITY(tGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeGhostedContainer(LOC::X | LOC::DxDt,*tGhostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(tGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(tGhostedContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeGhostedContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*tGhostedContainer);
      TEST_ASSERT(tGhostedContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(tGhostedContainer->get_A()!=Teuchos::null);
   }
}

TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, fe_opt_in)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::Comm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Teuchos::Comm<int> > failure_comm = THIS_,_SERIAL_BUILDS_,_SHOULD_FAIL;
   #endif

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();

   RCP<panzer::GlobalIndexer> indexer
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));

   typedef panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> LOFType;

   // default construction: FE assembly is not enabled, so the FE accessors must throw
   {
      Teuchos::RCP<LOFType> la_factory
            = Teuchos::rcp(new LOFType(tComm.getConst(),indexer));

      TEST_ASSERT(!la_factory->useFEAssembly());
      TEST_THROW(la_factory->getFEGraph(),std::logic_error);
      TEST_THROW(la_factory->getFEMatrix(),std::logic_error);
      TEST_THROW(la_factory->getFEMultiVector(),std::logic_error);
   }

   // opt in to FE assembly
   {
      Teuchos::RCP<LOFType> la_factory
            = Teuchos::rcp(new LOFType(tComm.getConst(),indexer,true));

      TEST_ASSERT(la_factory->useFEAssembly());

      // getFEGraph() returns the graph in the ACTIVE_OWNED state. Tpetra::FECrsGraph/
      // FECrsMatrix have no public accessor to ask "what state are you in" (activeCrsGraph_/
      // activeCrsMatrix_ are private, toggled only by the private switchActiveCrsGraph()/
      // switchActiveCrsMatrix(), invoked from the public beginAssembly()/endAssembly()), so
      // this is known from control flow, not queried: buildFEGraph() (the function backing
      // getFEGraph()) itself calls feGraph->endAssembly() right before returning the graph --
      // that call is what performs the ACTIVE_OWNED_PLUS_SHARED -> ACTIVE_OWNED toggle. So by
      // the time this test ever sees the graph, the toggle has already happened; nothing here
      // does it. The isSameAs() check just below is a deliberate, direct proof of that (rather
      // than just trusting the comment) -- if the graph were still owned+shared, its row map
      // would be the larger ghosted map and this would fail outright.
      RCP<LOFType::FECrsGraphType> feGraph = la_factory->getFEGraph();
      TEST_ASSERT(feGraph!=Teuchos::null);
      TEST_ASSERT(feGraph->isFillComplete());
      TEST_ASSERT(feGraph->getRowMap()->isSameAs(*la_factory->getMap()));

      // the FE graph (built via the V2 FECrsGraph constructor) must have the same
      // sparsity pattern as the "old" owned graph built via buildGraph()/buildGhostedGraph().
      // Both graphs perform the same cross-rank structural merge by this point (the FE graph
      // via FECrsGraph::endFill()'s self-export; the legacy graph via buildGraph()'s explicit
      // doExport(INSERT) from the ghosted graph), so this must be an EXACT match, row for row,
      // GID for GID -- not just matching entry counts.
      RCP<LOFType::CrsGraphType> oldGraph = la_factory->getGraph();
      TEST_EQUALITY(feGraph->getGlobalNumRows(),oldGraph->getGlobalNumRows());
      TEST_EQUALITY(feGraph->getGlobalNumEntries(),oldGraph->getGlobalNumEntries());

      RCP<const Map> ownedMap = feGraph->getRowMap();
      RCP<const Map> feColMap = feGraph->getColMap();
      RCP<const Map> oldColMap = oldGraph->getColMap();
      for(size_t lclRow=0;lclRow<ownedMap->getLocalNumElements();lclRow++) {
         LOFType::CrsGraphType::local_inds_host_view_type feIndices, oldIndices;
         feGraph->getLocalRowView(lclRow,feIndices);
         oldGraph->getLocalRowView(lclRow,oldIndices);

         std::set<panzer::GlobalOrdinal> feGids, oldGids;
         for(std::size_t k=0;k<feIndices.size();k++)
            feGids.insert(feColMap->getGlobalElement(feIndices(k)));
         for(std::size_t k=0;k<oldIndices.size();k++)
            oldGids.insert(oldColMap->getGlobalElement(oldIndices(k)));

         TEST_EQUALITY(feGids.size(),oldGids.size());
         TEST_ASSERT(feGids==oldGids);
      }

      // getFEGraph() is cached -- repeated calls must return the same object
      TEST_EQUALITY(la_factory->getFEGraph().get(),feGraph.get());

      // getFEMatrix() presents the OWNED_PLUS_SHARED (ghosted) view -- but unlike feGraph
      // above, that's not because anything here (or in getFEMatrix()) calls beginAssembly()/
      // endAssembly() -- getFEMatrix() just does `new FECrsMatrixType(getFEGraph())`, no
      // assembly call at all. The starting state is a Tpetra::FECrsMatrix CONSTRUCTOR
      // default: it reads an optional "start owned" ParameterList entry (defaulting to
      // false when, as here, no params are passed), and starts in OWNED_PLUS_SHARED unless
      // that's explicitly set true. So this state comes from construction, not a toggle.
      // The graph-level checks above already verify the FE sparsity pattern matches the
      // non-FE path exactly, so these are just sanity/non-null checks on the thin wrappers
      // (see fe_owned_shared_matches_ghosted for a real content check of this ghosted state).
      RCP<LOFType::FECrsMatrixType> feMatrix = la_factory->getFEMatrix();
      TEST_ASSERT(feMatrix!=Teuchos::null);
      TEST_ASSERT(feMatrix->getLocalNumRows()>0);

      // getFEMultiVector() builds `new FEMultiVectorType(getMap(), feGraph->getImporter(),
      // numVectors)`. FEMultiVector's constructor uses `importer.is_null() ? map :
      // importer->getTargetMap()` as its active map -- since our importer is feGraph's own
      // (domainMap -> colMap) importer, its target map IS feGraph->getColMap(), so feVec's
      // map should be exactly that (isSameAs, not just fitted -- this is a direct
      // construction identity, not merely a compatibility guarantee).
      RCP<LOFType::FEMultiVectorType> feVec = la_factory->getFEMultiVector();
      TEST_ASSERT(feVec!=Teuchos::null);
      TEST_EQUALITY(feVec->getNumVectors(),static_cast<std::size_t>(1));
      TEST_ASSERT(feVec->getLocalLength()>0);
      TEST_ASSERT(feVec->getMap()->isSameAs(*feGraph->getColMap()));

      // The actual "does the multivector have the same GID/row numbering" question: since
      // feGraph's column map is only guaranteed *locally fitted* to the traditional ghosted
      // map (getGhostedMap()) -- not identical to it, per the V2 FECrsGraph column-map
      // discussion in fe_owned_shared_matches_ghosted -- checking isSameAs() against
      // la_factory->getGhostedMap() here would not generally hold. isLocallyFitted() is the
      // correct, weaker check: every local id in the traditional ghosted map must map to the
      // same GID at the same local id in feVec's map, i.e. local assembly code indexing into
      // feVec using the GlobalIndexer's own (ghosted) local ids lands on the correct rows.
      TEST_ASSERT(feVec->getMap()->isLocallyFitted(*la_factory->getGhostedMap()));
   }
}

TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, fe_owned_shared_matches_ghosted)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::Comm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Teuchos::Comm<int> > failure_comm = THIS_,_SERIAL_BUILDS_,_SHOULD_FAIL;
   #endif

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();

   RCP<panzer::GlobalIndexer> indexer
         = rcp(new unit_test::GlobalIndexer(myRank,numProc));

   typedef panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> LOFType;

   Teuchos::RCP<LOFType> la_factory
         = Teuchos::rcp(new LOFType(tComm.getConst(),indexer,true));

   // A freshly built Tpetra::FECrsMatrix starts life presenting the owned+shared
   // ("ghosted") active view -- this is exactly the view local element-loop assembly
   // (beginAssembly()/scatter/endAssembly()) writes into.
   //
   // It is NOT expected to be byte-identical, row for row, to the traditional ghosted
   // Tpetra::CrsMatrix built by the non-FE path (getGhostedTpetraMatrix()/getGhostedGraph()).
   // FECrsGraph::endFill() performs a genuine cross-rank STRUCTURAL merge as part of
   // finishing the graph (FECrsGraph::doOwnedPlusSharedToOwned(), a self-export with
   // Tpetra::ADD): for every row this rank OWNS, it unions in whatever entries ANY other
   // rank locally inserted for that same row (e.g. two elements on different ranks sharing
   // an interface node). The traditional ghosted graph never does this merge -- it only ever
   // reflects what THIS rank's own local elements inserted; cross-rank merging is deferred to
   // a separate, later VALUE-level Export (ghostToGlobalTpetraMatrix). So:
   //   - For rows NOT owned by this rank (pure ghost rows): FECrsGraph never touches them
   //     during the merge (the reverse direction, doOwnedToOwnedPlusShared(), is a documented
   //     no-op), so they must match the traditional ghosted row EXACTLY.
   //   - For rows OWNED by this rank: the FE row is the traditional row PLUS whatever any
   //     other rank's elements also contributed to that row, i.e. a (non-strict) superset.
   // This is by design, and is exactly what makes it correct to use FE's owned+shared matrix
   // as the local assembly target: each rank still only ever scatters its own elements' data
   // using its own ghosted LIDs, and the pre-sized structure is what lets the matrix's own
   // endAssembly() correctly sum in every rank's contribution afterward, without the
   // legacy path's separate explicit ghostToGlobalContainer export step.
   RCP<LOFType::FECrsMatrixType> feMatrix = la_factory->getFEMatrix();
   RCP<const LOFType::CrsGraphType> feMatrixGraph = feMatrix->getCrsGraph();
   RCP<const LOFType::CrsGraphType> oldGhostedGraph = la_factory->getGhostedGraph();

   // Row maps must match exactly: both are built directly from getGhostedMap(), with no
   // Tpetra-internal reordering involved (unlike column maps, row maps are supplied
   // directly rather than derived from the insertion pattern).
   RCP<const Map> feRowMap = feMatrixGraph->getRowMap();
   RCP<const Map> oldRowMap = oldGhostedGraph->getRowMap();
   TEST_ASSERT(feRowMap->isSameAs(*oldRowMap));

   // Column maps need NOT be identical: the FE graph was built with the "V2" FECrsGraph
   // constructor, which only guarantees that its column map is *locally fitted* to the
   // owned+shared domain map supplied at construction (getGhostedColMap()) -- i.e. the
   // traditional ghosted column map's global ids must appear, in the same order, as a
   // leading prefix of the FE column map. Any additional columns beyond that prefix
   // (from the cross-rank merge described above) are appended afterward by Tpetra and are
   // not part of the guarantee. This ordering is the entire reason V2 was chosen over V1: it
   // is what lets local (ghosted) element assembly use the GlobalIndexer's own local ids
   // directly as FE matrix local column indices.
   RCP<const Map> feColMap = feMatrixGraph->getColMap();
   RCP<const Map> oldColMap = oldGhostedGraph->getColMap();
   TEST_ASSERT(feColMap->isLocallyFitted(*oldColMap));

   // Per row: exact match for rows this rank doesn't own, superset for rows it does (see the
   // explanation above). Order within a row is not checked: a row's local storage order is an
   // internal CrsGraph implementation detail, not a documented guarantee -- the per-column
   // *map* ordering checked above (isLocallyFitted) is the property that actually matters for
   // local-index-based assembly, not row-local storage order.
   std::vector<panzer::GlobalOrdinal> ownedIndices;
   indexer->getOwnedIndices(ownedIndices);
   std::set<panzer::GlobalOrdinal> ownedSet(ownedIndices.begin(),ownedIndices.end());

   for(size_t lclRow=0;lclRow<feRowMap->getLocalNumElements();lclRow++) {
      LOFType::CrsGraphType::local_inds_host_view_type feIndices, oldIndices;
      feMatrixGraph->getLocalRowView(lclRow,feIndices);
      oldGhostedGraph->getLocalRowView(lclRow,oldIndices);

      std::set<panzer::GlobalOrdinal> feGids, oldGids;
      for(std::size_t k=0;k<feIndices.size();k++)
         feGids.insert(feColMap->getGlobalElement(feIndices(k)));
      for(std::size_t k=0;k<oldIndices.size();k++)
         oldGids.insert(oldColMap->getGlobalElement(oldIndices(k)));

      const panzer::GlobalOrdinal rowGid = feRowMap->getGlobalElement(lclRow);
      if(ownedSet.count(rowGid)>0) {
         // owned row: FE row must be at least the traditional row's entries (superset)
         TEST_ASSERT(std::includes(feGids.begin(),feGids.end(),oldGids.begin(),oldGids.end()));
      }
      else {
         // ghost-only row: FE never touches it during the merge, so it must match exactly
         TEST_EQUALITY(feGids.size(),oldGids.size());
         TEST_ASSERT(feGids==oldGids);
      }
   }
}

}
