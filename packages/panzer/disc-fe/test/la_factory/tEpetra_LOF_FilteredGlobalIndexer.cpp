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
#include <Teuchos_DefaultMpiComm.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "Thyra_SpmdVectorSpaceBase.hpp"

#include "Panzer_DOFManager.hpp"
#include "Panzer_Filtered_GlobalIndexer.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"

#include "UnitTest_ConnManager.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

namespace panzer {

template <typename Intrepid2Type>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   Teuchos::RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tEpetra_LOF_FilteredUGI,epetra_lof)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   typedef Thyra::SpmdVectorSpaceBase<double> SpmdSpace;


   // build global (or serial communicator)
   #ifdef HAVE_MPI
      RCP<Epetra_Comm> eComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      RCP<const Teuchos::MpiComm<int> > tComm 
         = rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      PANZER DOES NOT DO SERIAL
   #endif

   // panzer::pauseToAttach();

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager> connManager = rcp(new unit_test::ConnManager(myRank,numProc));
   RCP<DOFManager> dofManager = rcp(new DOFManager); 
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);

   RCP<const panzer::FieldPattern> patternC1 
     = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

   dofManager->addField("T",patternC1); // add it to all three blocks
   dofManager->addField("block_0","Ux",patternC1);
   dofManager->addField("block_0","Uy",patternC1);
   dofManager->addField("block_0","P",patternC1);
   dofManager->addField("block_2","rho",patternC1);

   dofManager->buildGlobalUnknowns();

   // get GIDs on the "bottom" side of local cell 0, use those as the filtering
   // GIDs
   std::vector<panzer::GlobalOrdinal> filtered;
   {
     std::pair<std::vector<int>,std::vector<int> > fieldOffsets
         = dofManager->getGIDFieldOffsets_closure("block_0",dofManager->getFieldNum("Ux"),1,0);
 
     std::vector<panzer::GlobalOrdinal> gids;
     dofManager->getElementGIDs(0,gids);

     filtered.resize(fieldOffsets.first.size());
     for(std::size_t i=0;i<fieldOffsets.first.size();i++)
       filtered[i] = gids[fieldOffsets.first[i]];
   }

   TEST_EQUALITY(filtered.size(),2);

   RCP<Filtered_GlobalIndexer> filtered_ugi = rcp(new Filtered_GlobalIndexer);
   filtered_ugi->initialize(dofManager,filtered);

   out << "check out ownsership" << std::endl;
   {
     std::vector<bool> isOwned;
     dofManager->ownedIndices(filtered,isOwned);

     TEST_EQUALITY(isOwned.size(),2);
     out << "IS OWNED: " << isOwned[0] << " " << isOwned[1] << std::endl;

     filtered_ugi->ownedIndices(filtered,isOwned);

     TEST_EQUALITY(isOwned.size(),2);
     TEST_EQUALITY(isOwned[0],false);
     TEST_EQUALITY(isOwned[1],false);
   }

   out << "test out sizes construction by LOF" << std::endl;
   {
     std::vector<panzer::GlobalOrdinal> indices_f;
     filtered_ugi->getOwnedIndices(indices_f);

     BlockedEpetraLinearObjFactory<panzer::Traits,int> lof(tComm,filtered_ugi);
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof.getThyraDomainSpace(),true)->localSubDim(),Teuchos::as<int>(indices_f.size())); 
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof.getThyraRangeSpace(),true)->localSubDim(),Teuchos::as<int>(indices_f.size())); 

     RCP<Thyra::LinearOpBase<double> > A = lof.getThyraMatrix();
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->range(),true)->localSubDim(),Teuchos::as<int>(indices_f.size())); 
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->domain(),true)->localSubDim(),Teuchos::as<int>(indices_f.size())); 

     // This next chunk of code tests to ensure parallel communication works as
     // expected, in particular that a filtered owned vector can be used with
     // an unfiltered ghosted vector and that the entries in the ghosted vector
     // will be preserved.
     ////////////////////////////////////////////////////////////////////////////////////

     // check that parallel communication works as expected
     RCP<Epetra_Import> importer   = lof.getGhostedImport2(0);
     RCP<Epetra_Map>    ghostedMap = lof.getGhostedMap2(0);
     RCP<Epetra_Map>    ownedMap   = lof.getMap(0);

     TEST_ASSERT(importer   != Teuchos::null);
     TEST_ASSERT(ghostedMap != Teuchos::null);
     TEST_ASSERT(ownedMap   != Teuchos::null);

     Epetra_Vector ghosted_x(*ghostedMap);
     Epetra_Vector owned_x(*ownedMap);

     ghosted_x.PutScalar(0.0);
     owned_x.PutScalar(1.0);

     ghosted_x.Import(owned_x, *importer, Insert);

     // check all the filtered indices remain 0
     int count=0;
     for(std::size_t i=0;i<filtered.size();i++) {
       int lid = ghostedMap->LID(filtered[0]);
       if(lid>=0) {
         count++;
         TEST_EQUALITY(ghosted_x[lid],0.0);
       }
     }

     int sum = 0;
     for(int i=0;i<ghosted_x.MyLength();i++) 
       sum += ghosted_x[i]; 

     // check that ones sum up to the number of ids
     // that were not filtered
     TEST_EQUALITY(sum,ghosted_x.MyLength()-count);

     
     // do a lazy test to ensure you can construct an owned matrix
     RCP<Thyra::LinearOpBase<double> > ownedMatrix  = lof.getThyraMatrix();
     TEST_ASSERT(ownedMatrix != Teuchos::null);
   }


}


}
