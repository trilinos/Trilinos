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

#include <string>
#include <iostream>

#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

// for testing gather/scatter construction
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "Thyra_EpetraThyraWrappers.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

RCP<Epetra_MultiVector> getEpetraMultiVector(RCP<Thyra::MultiVectorBase<double> > & vec,const Epetra_Map & eMap)
{
   return Thyra::get_Epetra_MultiVector(eMap,vec);
}

#if 0
TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, vector_constr)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));
 
   // setup factory
   panzer::EpetraLinearObjFactory<panzer::Traits,short> la_factory(eComm.getConst(),indexer);

   // build vectors from factory
   RCP<Thyra::MultiVectorBase<double> > ghostedVec = la_factory.getGhostedVector();
   RCP<Thyra::MultiVectorBase<double> > vec = la_factory.getVector();
   
   // convert to epetra vectors
   RCP<Epetra_Map> eMap = la_factory.getMap();
   RCP<Epetra_Map> eGhostedMap = la_factory.getGhostedMap();
   RCP<Epetra_MultiVector> ptrEVec = getEpetraMultiVector(vec,*eMap);
   RCP<Epetra_MultiVector> ptrEGhostedVec = getEpetraMultiVector(ghostedVec,*eGhostedMap);
   Epetra_Vector & eVec = *(*ptrEVec)(0);
   Epetra_Vector & eGhostedVec = *(*ptrEGhostedVec)(0);

   // check sizes of global epetra vectors
   TEST_EQUALITY(eVec.NumVectors(),1);
   TEST_EQUALITY(eVec.GlobalLength(),12);
   TEST_EQUALITY(eVec.MyLength(),6);

   TEST_EQUALITY(eGhostedVec.NumVectors(),1);
   TEST_EQUALITY(eGhostedVec.GlobalLength(),16);
   TEST_EQUALITY(eGhostedVec.MyLength(),8);

   // fill epetra vectors
   for(int i=0;i<eVec.MyLength();i++) 
      eVec[i] = double(eMap->GID(i))+0.1;
   for(int i=0;i<eGhostedVec.MyLength();i++) 
      eGhostedVec[i] = double(eGhostedMap->GID(i))+0.1;

   // run parallel assembly for global vector
   eVec.PutScalar(-10000.0);
   la_factory.ghostToGlobalVector(*ghostedVec,*vec);

   // check global vector 
   {
      for(int i=0;i<eVec.MyLength();i++) {
         int gid = eMap->GID(i); 
   
         if(gid==2 || gid==3 || gid==4 || gid==5) {
            TEST_FLOATING_EQUALITY(eVec[i],2.0*(double(gid)+0.1),1e-14);
         }
         else {
            TEST_FLOATING_EQUALITY(eVec[i],double(gid)+0.1,1e-14);
         }
      }
   }

   // construct ghosted vector
   eGhostedVec.PutScalar(-10000.0);
   la_factory.globalToGhostVector(*vec,*ghostedVec);

   // check ghosted vector 
   {
      eVec.PutScalar(0.0);
      for(int i=0;i<eGhostedVec.MyLength();i++) {
         int gid = eGhostedMap->GID(i); 
   
         if(gid==2 || gid==3 || gid==4 || gid==5) {
            TEST_FLOATING_EQUALITY(eGhostedVec[i],2.0*(double(gid)+0.1),1e-14);
         }
         else {
            TEST_FLOATING_EQUALITY(eGhostedVec[i],double(gid)+0.1,1e-14);
         }
      }
   }
}
#endif

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, gather_scatter_constr)
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

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<short>(myRank,numProc));
 
   // setup factory
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > la_factory
         = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));

   // build parameter lists for gather and scatters
   //////////////////////////////////////////////////////////////////
   std::size_t numCells = 10;

   Teuchos::ParameterList gatherParams;
   {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      // auxiliary information needed to construct basis object
      int baseCellDim = 2;
      int cubatureDegree = 2;
      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,baseCellDim,topo);

      // build DOF names
      RCP<std::vector<std::string> > dofNames = rcp(new std::vector<std::string>);
      dofNames->push_back("ux"); // in practice these probably would not be gathered together!
      dofNames->push_back("p");

      // build basis
      RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,cellData));

      // build gather parameter list
      gatherParams.set<RCP<std::vector<std::string> > >("DOF Names",dofNames);
      gatherParams.set<RCP<std::vector<std::string> > >("Indexer Names",dofNames);
      gatherParams.set<RCP<panzer::PureBasis> >("Basis",basis);
   }

   Teuchos::ParameterList scatterParams;
   {
      Teuchos::RCP<shards::CellTopology> topo = 
         Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

      int baseCellDim = 2;
      int cubatureDegree = 2;
      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,baseCellDim,topo);
   
      // build basis
      RCP<const panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,cellData));
   
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

      int baseCellDim = 2;
      int cubatureDegree = 2;
      std::string basisType = "Q1";
      panzer::CellData cellData(numCells,baseCellDim,topo);
   
      // build basis
      RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis(basisType,cellData));
   
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
         RCP<GatherSolution_Epetra<EvalType,panzer::Traits,short,int> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(gatherSolutionEval!=Teuchos::null);
   
         const std::vector<RCP<PHX::FieldTag> > & fields = gatherSolutionEval->evaluatedFields();
         TEST_EQUALITY(fields.size(),2);
   
         TEST_EQUALITY(fields[0]->name(),"ux");
         TEST_EQUALITY(fields[1]->name(),"p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter test
      {
         evaluator = la_factory->buildScatter<EvalType>(scatterParams);

         out << "SCATTER RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterResidual_Epetra<EvalType,panzer::Traits,short,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter dirichlet test
      {
         evaluator = la_factory->buildScatterDirichlet<EvalType>(scatterDirichletParams);

         out << "SCATTER DIRICHLET RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,short,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

   }

   {
      typedef panzer::Traits::Jacobian EvalType;

      // gather test
      {
         evaluator = la_factory->buildGather<EvalType>(gatherParams);

         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<GatherSolution_Epetra<EvalType,panzer::Traits,short,int> > gatherSolutionEval 
               = rcp_dynamic_cast<GatherSolution_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(gatherSolutionEval!=Teuchos::null);
   
         const std::vector<RCP<PHX::FieldTag> > & fields = gatherSolutionEval->evaluatedFields();
         TEST_EQUALITY(fields.size(),2);
   
         TEST_EQUALITY(fields[0]->name(),"ux");
         TEST_EQUALITY(fields[1]->name(),"p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter test
      {
         evaluator = la_factory->buildScatter<EvalType>(scatterParams);

         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterResidual_Epetra<EvalType,panzer::Traits,short,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterResidual_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1

         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

      // scatter dirichlet test
      {
         evaluator = la_factory->buildScatterDirichlet<EvalType>(scatterDirichletParams);

         out << "SCATTER DIRICHLET RES NAME: \"" << evaluator->getName() << "\"" << std::endl;
         TEST_ASSERT(evaluator!=Teuchos::null);
         RCP<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,short,int> > scatterResidual 
               = rcp_dynamic_cast<ScatterDirichletResidual_Epetra<EvalType,panzer::Traits,short,int> >(evaluator);
         TEST_ASSERT(scatterResidual!=Teuchos::null);

         const std::vector<RCP<PHX::FieldTag> > & evalFields = scatterResidual->evaluatedFields();
         TEST_EQUALITY(evalFields.size(),1); // this is a dummy holder for the sake of the field manager
   
         const std::vector<RCP<PHX::FieldTag> > & fields = scatterResidual->dependentFields();
         TEST_EQUALITY(fields.size(),2); // these store the residual values
   
         TEST_EQUALITY(fields[0]->name(),"Residual_ux");
         TEST_EQUALITY(fields[1]->name(),"Residual_p");
   
         TEST_EQUALITY(fields[0]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[0]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
   
         TEST_EQUALITY(fields[1]->dataLayout().dimension(0),Teuchos::as<int>(numCells));
         TEST_EQUALITY(fields[1]->dataLayout().dimension(1),Teuchos::as<int>(4)); // for Q1
      }

   }
}

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, adjustDirichlet)
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

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef EpetraLinearObjContainer ELOC;

   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer<short>(myRank,numProc));

   // setup factory
   Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > la_factory
         = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<EpetraLinearObjContainer> e_0   = rcp_dynamic_cast<EpetraLinearObjContainer>(ghosted_0);
   RCP<EpetraLinearObjContainer> e_1   = rcp_dynamic_cast<EpetraLinearObjContainer>(ghosted_1);
   RCP<EpetraLinearObjContainer> e_sys = rcp_dynamic_cast<EpetraLinearObjContainer>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(e_0->get_x()));
   TEST_ASSERT(!Teuchos::is_null(e_1->get_x()));
   TEST_ASSERT(!Teuchos::is_null(e_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(e_sys->get_A()));

   e_sys->get_f()->PutScalar(-3.0); // put some garbage in the systems
   e_sys->get_A()->PutScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      // case 0
      (*(e_0->get_x()))[0] = 1.0; // GID = 0
      (*(e_1->get_x()))[0] = 1.0; // GID = 0

      // case 1
      (*(e_0->get_x()))[2] = 1.0; // GID = 2
      (*(e_1->get_x()))[2] = 2.0; // GID = 2

      // case 2
      (*(e_1->get_x()))[5] = 2.0; // GID = 5
   }
   else if(myRank==1) {
      // case 0
      (*(e_0->get_x()))[3] = 1.0; // GID = 9
      (*(e_1->get_x()))[3] = 1.0; // GID = 9

      // case 1
      (*(e_0->get_x()))[0] = 1.0; // GID = 2
      (*(e_1->get_x()))[0] = 2.0; // GID = 2

      // case 2
      (*(e_1->get_x()))[6] = 2.0; // GID = 4
   }
   else 
      TEUCHOS_ASSERT(false);

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   int numEntries = 0;
   double * values = 0;
   int * indices = 0;

   if(myRank==0) {   
      TEST_EQUALITY((*e_sys->get_f())[0],-3.0);     // case 0
      e_sys->get_A()->ExtractMyRowView(0,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0);

      TEST_EQUALITY((*e_sys->get_f())[2],-3.0/2.0); // case 1
      e_sys->get_A()->ExtractMyRowView(2,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0/2.0);

      TEST_EQUALITY((*e_sys->get_f())[5],0.0);      // case 2
      e_sys->get_A()->ExtractMyRowView(5,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],0.0);
   }
   else if(myRank==1) {
      TEST_EQUALITY((*e_sys->get_f())[3],-3.0);     // case 0
      e_sys->get_A()->ExtractMyRowView(3,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0);

      TEST_EQUALITY((*e_sys->get_f())[0],-3.0/2.0); // case 1
      e_sys->get_A()->ExtractMyRowView(0,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],-3.0/2.0);

      TEST_EQUALITY((*e_sys->get_f())[6],0.0);      // case 2
      e_sys->get_A()->ExtractMyRowView(6,numEntries,values,indices);
      for(int i=0;i<numEntries;i++) TEST_EQUALITY(values[i],0.0);
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, initializeContianer)
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

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef EpetraLinearObjContainer ELOC;

   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer<short>(myRank,numProc));

   std::vector<int> ownedIndices, ownedAndSharedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndSharedIndices(ownedAndSharedIndices);
 
   // setup factory
   Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > la_factory
         = Teuchos::rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));

   RCP<LinearObjContainer> container = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghostedContainer = la_factory->buildGhostedLinearObjContainer();

   RCP<EpetraLinearObjContainer> eContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(container);
   RCP<EpetraLinearObjContainer> eGhostedContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(ghostedContainer);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeContainer(ELOC::X,*container);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_x()->MyLength(),(int) ownedIndices.size());
   
      la_factory->initializeContainer(ELOC::DxDt,*container);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt()->MyLength(),(int) ownedIndices.size());
   
      la_factory->initializeContainer(ELOC::F,*container);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_f()->MyLength(),(int) ownedIndices.size());
   
      la_factory->initializeContainer(ELOC::Mat,*container);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_A()->NumMyRows(),(int) ownedIndices.size());
   
      // jacobian and residual vector output
      la_factory->initializeContainer(ELOC::F | ELOC::Mat,*container);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeContainer(ELOC::X | ELOC::DxDt,*container);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeContainer(ELOC::X | ELOC::DxDt | ELOC::F | ELOC::Mat,*container);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
   
      // Epetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      la_factory->initializeContainer(ELOC::X,*eContainer);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(ELOC::DxDt,*eContainer);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(ELOC::F,*eContainer);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
   
      la_factory->initializeContainer(ELOC::Mat,*eContainer);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      la_factory->initializeContainer(ELOC::F | ELOC::Mat,*eContainer);
      TEST_EQUALITY(eContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      la_factory->initializeContainer(ELOC::X | ELOC::DxDt,*eContainer);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eContainer->get_A(),    Teuchos::null)
   
      // everything
      la_factory->initializeContainer(ELOC::X | ELOC::DxDt | ELOC::F | ELOC::Mat,*eContainer);
      TEST_ASSERT(eContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(eContainer->get_A()!=Teuchos::null);
   }

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
      TEST_EQUALITY(eGhostedContainer->get_x()->MyLength(),(int) ownedAndSharedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::DxDt,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt()->MyLength(),(int) ownedAndSharedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::F,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_A(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f()->MyLength(),(int) ownedAndSharedIndices.size());
   
      la_factory->initializeGhostedContainer(ELOC::Mat,*ghostedContainer);
      TEST_EQUALITY(eGhostedContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(eGhostedContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(eGhostedContainer->get_A()!=Teuchos::null);
      TEST_EQUALITY(eGhostedContainer->get_A()->NumMyRows(),(int) ownedAndSharedIndices.size());
   
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

}
