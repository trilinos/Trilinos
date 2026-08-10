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

      // A range-space (residual-side) FE multivector's owned+shared view must be over the
      // ghosted ROW map exactly -- that is the map the scatter evaluators index with the
      // GlobalIndexer's local ids, and the source map getGhostedExport() migrates from.
      // isSameAs, not merely isLocallyFitted: anything weaker would not be a legal source
      // for that export.
      RCP<LOFType::FEMultiVectorType> feVec = la_factory->getFEMultiVector();
      TEST_ASSERT(feVec!=Teuchos::null);
      TEST_EQUALITY(feVec->getNumVectors(),static_cast<std::size_t>(1));
      TEST_ASSERT(feVec->getLocalLength()>0);
      TEST_ASSERT(feVec->getMap()->isSameAs(*la_factory->getGhostedMap()));
      TEST_EQUALITY(feVec->getLocalLength(),la_factory->getGhostedMap()->getLocalNumElements());

      // It must specifically NOT be the FE graph's column map, which the cross-rank clique
      // merge widens beyond the ghosted row map. Building the residual over that map was the
      // original bug this asserts against.
      TEST_ASSERT(!feVec->getMap()->isSameAs(*feGraph->getColMap()));

      // The domain-space counterpart is over the ghosted COLUMN map. With no separate column
      // GlobalIndexer the col maps fall back to the row maps, so here the two agree; the
      // accessors are still distinct so a rectangular/blocked setup gets the right one.
      RCP<LOFType::FEMultiVectorType> feColVec = la_factory->getFEColMultiVector();
      TEST_ASSERT(feColVec!=Teuchos::null);
      TEST_ASSERT(feColVec->getMap()->isSameAs(*la_factory->getGhostedColMap()));
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

// Verifies that an FE matrix placed in a container is actually VISIBLE through the
// container's accessors and DETECTED by the factory's beginFill()/endFill(). This is the
// coverage for the FE branch added to beginFill()/endFill(): those guard on
// `A != null` and then rcp_dynamic_cast<FECrsMatrixType>, so if either step silently
// yielded null the FE state machine would never be driven and ghost-row contributions
// would never migrate to the owned rows -- a silent wrong answer rather than an error.
TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, fe_container_visibility)
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
   typedef panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOCType;

   Teuchos::RCP<LOFType> la_factory
         = Teuchos::rcp(new LOFType(tComm.getConst(),indexer,true));

   RCP<LinearObjContainer> loc = la_factory->buildGhostedLinearObjContainer();
   RCP<LOCType> tloc = rcp_dynamic_cast<LOCType>(loc);
   TEST_ASSERT(tloc!=Teuchos::null);

   // Place FE objects into the container. set_A() takes RCP<CrsMatrixType>; an
   // RCP<FECrsMatrixType> converts implicitly since FECrsMatrix IS-A CrsMatrix, so no
   // container type change was needed for the matrix (unlike the vectors).
   RCP<LOFType::FECrsMatrixType> feMatrix = la_factory->getFEMatrix();
   RCP<LOFType::FEMultiVectorType> feVec = la_factory->getFEMultiVector();
   tloc->set_A(feMatrix);
   tloc->set_f(feVec);

   // The matrix must survive the round trip through the container's CrsMatrix-typed
   // storage: non-null, and still dynamic_cast-able back to FECrsMatrix. If either failed,
   // beginFill()/endFill() would silently fall through to the plain-CrsMatrix branch.
   TEST_ASSERT(tloc->get_A()!=Teuchos::null);
   TEST_ASSERT(Teuchos::rcp_dynamic_cast<LOFType::FECrsMatrixType>(tloc->get_A())!=Teuchos::null);
   TEST_EQUALITY(tloc->get_A().get(),static_cast<LOFType::CrsMatrixType*>(feMatrix.get()));

   // The FEMultiVector is visible through the _mv() accessor...
   TEST_ASSERT(tloc->get_f_mv()!=Teuchos::null);
   TEST_EQUALITY(tloc->get_f_mv().get(),static_cast<LOFType::MultiVectorType*>(feVec.get()));
   // ...and the narrowing get_f() throws loudly rather than silently returning null,
   // since an FEMultiVector is not a Tpetra::Vector.
   TEST_THROW(tloc->get_f(),std::bad_cast);

   // Finally, drive the fill cycle the AssemblyEngine would: beginFill() must put the FE
   // matrix into its owned+shared assembly state, and endFill() must run the owned+shared
   // -> owned migration. Neither should throw, and the FE state machine asserts internally
   // if they are called out of order, so completing this round trip is itself the check.
   // Before assembly the FE matrix presents its owned+shared view, whose row map is exactly
   // getGhostedMap() -- this is what lets the scatter evaluators keep writing through the
   // GlobalIndexer's ghosted LIDs with no change.
   TEST_ASSERT(feMatrix->getRowMap()->isSameAs(*la_factory->getGhostedMap()));

   // The range-space FE multivector is over the ghosted ROW map, so it IS a legal source for
   // getGhostedExport() and is indexed by the same local ids the scatter evaluators use.
   TEST_ASSERT(feVec->getMap()->isSameAs(*la_factory->getGhostedMap()));

   // Drive the fill cycle: beginFill() puts the FE matrix into its owned+shared assembly
   // state, endFill() runs the owned+shared -> owned migration and switches the active view.
   TEST_NOTHROW(la_factory->beginFill(*tloc));
   TEST_NOTHROW(la_factory->endFill(*tloc));
   TEST_ASSERT(feMatrix->getRowMap()->isSameAs(*la_factory->getMap()));
}

// The core Step 3 test: in FE mode the owned and ghosted containers must hold the SAME
// FECrsMatrix, and a full AssemblyEngine-style fill cycle over that shared object must
// produce a correctly cross-rank-summed owned Jacobian -- WITHOUT the traditional manual
// ghost->global matrix export, which endAssembly() replaces.
TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, fe_shared_matrix_assembly_cycle)
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
   typedef panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOCType;
   typedef panzer::LinearObjContainer LOC;

   Teuchos::RCP<LOFType> la_factory
         = Teuchos::rcp(new LOFType(tComm.getConst(),indexer,true));

   // Build both containers exactly as AssemblyEngine/ModelEvaluator would.
   RCP<LinearObjContainer> gLoc = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghLoc = la_factory->buildGhostedLinearObjContainer();
   la_factory->initializeContainer(LOC::Mat,*gLoc);
   la_factory->initializeGhostedContainer(LOC::Mat,*ghLoc);

   RCP<LOCType> gt = rcp_dynamic_cast<LOCType>(gLoc);
   RCP<LOCType> ght = rcp_dynamic_cast<LOCType>(ghLoc);

   // THE key structural property of the FE path: one object, two views. The owned container
   // gets its matrix via getTpetraMatrix() (the same route ModelEvaluator::create_W_op takes
   // to hand NOX its W operator) and the ghosted container via getGhostedTpetraMatrix();
   // in FE mode both must resolve to the identical cached FECrsMatrix.
   TEST_ASSERT(gt->get_A()!=Teuchos::null);
   TEST_ASSERT(ght->get_A()!=Teuchos::null);
   TEST_EQUALITY(gt->get_A().get(),ght->get_A().get());

   RCP<LOFType::FECrsMatrixType> feA
      = Teuchos::rcp_dynamic_cast<LOFType::FECrsMatrixType>(ght->get_A());
   TEST_ASSERT(feA!=Teuchos::null);

   // --- AssemblyEngine-style cycle ---------------------------------------------------
   // 1) begin fill on the ghosted container
   la_factory->beginFill(*ghLoc);

   // 2) "scatter": every rank sums 1.0 into the diagonal of each row it can see through its
   //    GHOSTED map -- i.e. exactly what a local element loop does with ghosted LIDs. Rows
   //    that two ranks both see (the shared interface DOFs) therefore receive 1.0 twice and
   //    must end up at 2.0 only if the cross-rank merge actually happened.
   {
      std::vector<panzer::GlobalOrdinal> ghosted;
      indexer->getOwnedAndGhostedIndices(ghosted);
      Teuchos::Array<panzer::GlobalOrdinal> cols(1);
      Teuchos::Array<double> vals(1); vals[0] = 1.0;
      for(std::size_t i=0;i<ghosted.size();i++) {
         cols[0] = ghosted[i];
         feA->sumIntoGlobalValues(ghosted[i],cols,vals);
      }
   }

   // 3) ghost->global. For the shared FE matrix this must be a NO-OP: the migration is done
   //    in place by endAssembly() below. If this instead attempted the traditional export it
   //    would be exporting the object into itself.
   TEST_NOTHROW(la_factory->ghostToGlobalContainer(*ghLoc,*gLoc,LOC::Mat));

   // 4) the remaining begin/end pair AssemblyEngine issues on the global container, then the
   //    ghosted one. Only one begin/end may actually reach the FE state machine.
   TEST_NOTHROW(la_factory->beginFill(*gLoc));
   TEST_NOTHROW(la_factory->endFill(*gLoc));
   TEST_NOTHROW(la_factory->endFill(*ghLoc));

   // --- verify the assembled OWNED Jacobian -------------------------------------------
   // The matrix now presents its owned view.
   TEST_ASSERT(feA->getRowMap()->isSameAs(*la_factory->getMap()));

   // Which owned rows are also visible (ghosted) on the other rank? In the 2-rank mock those
   // are the shared interface DOFs; each is summed into by BOTH ranks, so its diagonal must
   // be 2.0. Rows only this rank sees must be 1.0. Getting 1.0 everywhere would mean the
   // cross-rank merge silently did not happen.
   std::vector<panzer::GlobalOrdinal> ownedIdx;
   indexer->getOwnedIndices(ownedIdx);
   std::set<panzer::GlobalOrdinal> sharedGids;
   if(myRank==0) { sharedGids.insert(2); sharedGids.insert(3); }
   else          { sharedGids.insert(4); sharedGids.insert(5); }

   RCP<const Map> ownedMap = feA->getRowMap();
   for(std::size_t i=0;i<ownedIdx.size();i++) {
      panzer::GlobalOrdinal gid = ownedIdx[i];
      LOFType::CrsMatrixType::nonconst_global_inds_host_view_type inds("inds",feA->getNumEntriesInGlobalRow(gid));
      LOFType::CrsMatrixType::nonconst_values_host_view_type vals("vals",feA->getNumEntriesInGlobalRow(gid));
      size_t numEnt = 0;
      feA->getGlobalRowCopy(gid,inds,vals,numEnt);

      double diag = 0.0;
      for(size_t k=0;k<numEnt;k++)
         if(inds(k)==gid) diag = vals(k);

      const double expected = (sharedGids.find(gid)!=sharedGids.end()) ? 2.0 : 1.0;
      TEST_FLOATING_EQUALITY(diag,expected,1e-14);
   }
   TEST_ASSERT(ownedMap->getLocalNumElements()>0);
}

// The range-space FE multivector must be a legal, correct source for the traditional
// ghost->global export. This is the point of building it over getGhostedImport() rather than
// the FE graph's (column-map) importer: with the column map its local length exceeds the
// ghosted map's and getGhostedExport() -- whose source map is getGhostedMap() -- cannot
// consume it. Here it is filled through ghosted local ids exactly as a scatter evaluator
// would, then migrated, and the cross-rank sums are checked.
TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, fe_multivector_exports_like_ghosted)
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
   typedef panzer::TpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> LOCType;
   typedef panzer::LinearObjContainer LOC;

   Teuchos::RCP<LOFType> la_factory
         = Teuchos::rcp(new LOFType(tComm.getConst(),indexer,true));

   // ghosted container gets the FE residual; the global container keeps a plain owned vector
   // (which is what NOX supplies in practice, since f_out comes from the Thyra vector space).
   RCP<LinearObjContainer> gLoc = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghLoc = la_factory->buildGhostedLinearObjContainer();
   la_factory->initializeContainer(LOC::F,*gLoc);

   RCP<LOCType> gt = rcp_dynamic_cast<LOCType>(gLoc);
   RCP<LOCType> ght = rcp_dynamic_cast<LOCType>(ghLoc);

   RCP<LOFType::FEMultiVectorType> feF = la_factory->getFEMultiVector();
   ght->set_f(feF);
   feF->putScalar(0.0);
   gt->get_f_mv()->putScalar(0.0);

   // "scatter": write 1.0 at every ghosted local id, exactly as an element loop would.
   {
      auto host = feF->getLocalViewHost(Tpetra::Access::ReadWrite);
      for(size_t i=0;i<feF->getLocalLength();i++)
         host(i,0) = 1.0;
   }

   // migrate ghosted -> owned through the factory's normal path
   TEST_NOTHROW(la_factory->ghostToGlobalContainer(*ghLoc,*gLoc,LOC::F));

   // Owned entries also visible (ghosted) on the other rank were written by both ranks, so
   // they must sum to 2.0; rank-local entries stay at 1.0. Getting 1.0 everywhere would mean
   // the export silently moved nothing.
   std::vector<panzer::GlobalOrdinal> ownedIdx;
   indexer->getOwnedIndices(ownedIdx);
   std::set<panzer::GlobalOrdinal> sharedGids;
   if(myRank==0) { sharedGids.insert(2); sharedGids.insert(3); }
   else          { sharedGids.insert(4); sharedGids.insert(5); }

   RCP<const Map> ownedMap = la_factory->getMap();
   auto ownedHost = gt->get_f_mv()->getLocalViewHost(Tpetra::Access::ReadOnly);
   for(std::size_t i=0;i<ownedIdx.size();i++) {
      panzer::GlobalOrdinal gid = ownedIdx[i];
      int lid = ownedMap->getLocalElement(gid);
      TEST_ASSERT(lid>=0);
      const double expected = (sharedGids.find(gid)!=sharedGids.end()) ? 2.0 : 1.0;
      TEST_FLOATING_EQUALITY(ownedHost(lid,0),expected,1e-14);
   }
}

}
