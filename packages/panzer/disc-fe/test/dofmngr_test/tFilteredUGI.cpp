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
#include <Teuchos_DefaultMpiComm.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "Phalanx_KokkosUtilities.hpp"

#include "Thyra_SpmdVectorSpaceBase.hpp"

#include "Panzer_DOFManager.hpp"
#include "Panzer_Filtered_UniqueGlobalIndexer.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

#include "UnitTest_ConnManager.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Intrepid2_FieldContainer.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Intrepid2::FieldContainer<double> FieldContainer;

namespace panzer {

template <typename Intrepid2Type>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   Teuchos::RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2Type);
   Teuchos::RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tFilteredUGI,equivalence_test)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager(myRank,numProc));
   RCP<DOFManager<int,int> > dofManager = rcp(new DOFManager<int,int>); 
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);

   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   dofManager->addField("T",patternC1); // add it to all three blocks
   dofManager->addField("block_0","Ux",patternC1);
   dofManager->addField("block_0","Uy",patternC1);
   dofManager->addField("block_0","P",patternC1);
   dofManager->addField("block_2","rho",patternC1);

   dofManager->buildGlobalUnknowns();

   std::vector<int> filtered;
   Filtered_UniqueGlobalIndexer<int,int> filtered_ugi;
   filtered_ugi.initialize(dofManager,filtered);

   // check the GIDs
   {
     std::vector<int> gids,gids_f;

     dofManager->getElementGIDs(0,gids);
     filtered_ugi.getElementGIDs(0,gids_f);

     TEST_EQUALITY(gids.size(),gids_f.size());
     for(std::size_t i=0;i<gids.size();i++) {
       TEST_EQUALITY(gids[i],gids_f[i]);
     }
   }

   // check the LIDs
   {
     const std::vector<int> & lids = dofManager->getElementLIDs(0);
     const std::vector<int> & lids_f = filtered_ugi.getElementLIDs(0);

     TEST_EQUALITY(lids.size(),lids_f.size());
     for(std::size_t i=0;i<lids.size();i++) {
       TEST_EQUALITY(lids[i],lids_f[i]);
     }
   }
   
   // check owned and shared
   {
     std::vector<int> indices, indices_f;

     dofManager->getOwnedAndSharedIndices(indices);
     filtered_ugi.getOwnedAndSharedIndices(indices_f);

     TEST_EQUALITY(indices.size(),indices_f.size());
     for(std::size_t i=0;i<indices.size();i++) {
       TEST_EQUALITY(indices[i],indices_f[i]);
     }
   }

   // check owned 
   {
     std::vector<int> indices, indices_f;

     dofManager->getOwnedIndices(indices);
     filtered_ugi.getOwnedIndices(indices_f);

     TEST_EQUALITY(indices.size(),indices_f.size());
     for(std::size_t i=0;i<indices.size();i++) {
       TEST_EQUALITY(indices[i],indices_f[i]);
     }
   }

}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tFilteredUGI,filtering)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

      RCP<const Teuchos::MpiComm<int> > tComm 
         = rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager(myRank,numProc));
   RCP<DOFManager<int,int> > dofManager = rcp(new DOFManager<int,int>); 
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);

   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   dofManager->addField("T",patternC1); // add it to all three blocks
   dofManager->addField("block_0","Ux",patternC1);
   dofManager->addField("block_0","Uy",patternC1);
   dofManager->addField("block_0","P",patternC1);
   dofManager->addField("block_2","rho",patternC1);

   dofManager->buildGlobalUnknowns();

   std::vector<int> my_filtered;
   dofManager->getElementGIDs(1,my_filtered);

   // now we will compute all the filtered indices in this UGI
   // here we will print some things out for a sanity check
   std::vector<int> all_filtered;
   {
     int mySize = Teuchos::as<int>(my_filtered.size());
     std::vector<int> neighborSizes(numProc,0);
     Teuchos::gatherAll(*tComm,1,&mySize,Teuchos::as<int>(neighborSizes.size()),&neighborSizes[0]);

     out << "MY SZ = " << my_filtered.size() << std::endl;
     out << "SZ = " << neighborSizes[0] << std::endl;
     out << "SZ = " << neighborSizes[1] << std::endl;

     int totalSize = 0;
     for(std::size_t i=0;i<neighborSizes.size();i++)
       totalSize += Teuchos::as<int>(neighborSizes[i]);

     all_filtered.resize(totalSize);
     Teuchos::gatherAll(*tComm,mySize,&my_filtered[0],Teuchos::as<int>(totalSize),&all_filtered[0]);

     out << "MY Filtered = ";
     for(std::size_t i=0;i<my_filtered.size();i++)
       out << my_filtered[i] << " ";
     out << std::endl;

     out << "All Filtered = ";
     for(std::size_t i=0;i<all_filtered.size();i++)
       out << all_filtered[i] << " ";
     out << std::endl;
   }

   Filtered_UniqueGlobalIndexer<int,int> filtered_ugi;
   filtered_ugi.initialize(dofManager,my_filtered);

   // check the GIDs
   {
     std::vector<int> gids,gids_f;

     dofManager->getElementGIDs(0,gids);
     filtered_ugi.getElementGIDs(0,gids_f);

     TEST_EQUALITY(gids.size(),gids_f.size());
     for(std::size_t i=0;i<gids.size();i++) {
       TEST_EQUALITY(gids[i],gids_f[i]);
     }

     dofManager->getElementGIDs(1,gids);
     filtered_ugi.getElementGIDs(1,gids_f);

     TEST_EQUALITY(gids.size(),gids_f.size());
     for(std::size_t i=0;i<gids.size();i++) {
       TEST_EQUALITY(gids[i],gids_f[i]);
     }
   }

   // check the LIDs
   {
     const std::vector<int> & lids = dofManager->getElementLIDs(0);
     const std::vector<int> & lids_f = filtered_ugi.getElementLIDs(0);

     TEST_EQUALITY(lids.size(),lids_f.size());
     for(std::size_t i=0;i<lids.size();i++) {
       TEST_EQUALITY(lids[i],lids_f[i]);
     }
   }
   
   // check owned and shared
   {
     std::vector<int> indices, indices_f;

     dofManager->getOwnedAndSharedIndices(indices);
     filtered_ugi.getOwnedAndSharedIndices(indices_f);

     TEST_EQUALITY(indices.size(),indices_f.size());
     for(std::size_t i=0;i<indices.size();i++) {
       TEST_EQUALITY(indices[i],indices_f[i]);
     }
   }

   // check owned 
   {
     std::vector<int> indices, indices_f,diff;

     dofManager->getOwnedIndices(indices);
     filtered_ugi.getOwnedIndices(indices_f);

     // check that the size didn't grow
     TEST_ASSERT(indices_f.size()<=indices.size());

     // allocate space
     diff.resize(indices.size()); 

     // sort each set of owned indices, then do a set difference and 
     // crop the difference indices
     std::sort(indices.begin(),indices.end());
     std::sort(indices_f.begin(),indices_f.end());
     std::vector<int>::iterator it = std::set_difference(indices.begin(),indices.end(),
                                                         indices_f.begin(),indices_f.end(),diff.begin());
     diff.resize(it-diff.begin());

     // look at the difference, and make sure they are all in the
     // filtering vector
     TEST_ASSERT(diff.size()<=my_filtered.size());
     for(std::size_t i=0;i<diff.size();i++) {
       TEST_ASSERT(std::find(my_filtered.begin(),my_filtered.end(),diff[i])!=my_filtered.end());
     }

     // make sure no index in the owned index vector is in the filtering vector
     for(std::size_t i=0;i<indices_f.size();i++) {
       TEST_ASSERT(std::find(my_filtered.begin(),my_filtered.end(),indices_f[i])==my_filtered.end());
     }
   }

   // test getOwnedAndSharedNotFilteredIndicator
   out << "testing getOwnedAndSharedNotFilteredIndicator" << std::endl;
   {
     std::vector<int> indicator;
     filtered_ugi.getOwnedAndSharedNotFilteredIndicator(indicator);

     std::vector<int> indices;
     filtered_ugi.getOwnedAndSharedIndices(indices);

     TEST_EQUALITY(indices.size(),indicator.size());

     for(std::size_t i=0;i<indicator.size();i++) {
       if(indicator[i]==1) {
         // not filtered, should not be in all_filtered array
         TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices[i])==all_filtered.end());
       }
       else if(indicator[i]==0) {
         // filtered, should be in all_filtered array
         TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices[i])!=all_filtered.end());
       }
       else {
         // protect that neither one nor zero is set (also report a useful error)
         TEST_ASSERT(indicator[i]==0 || indicator[i]==1);
       }
     }
   }

   // test getFilteredOwnedAndSharedIndices
   out << "testing getFilteredOwnedAndSharedIndices" << std::endl;
   {
     std::vector<int> indices_f;
     filtered_ugi.getFilteredOwnedAndSharedIndices(indices_f);

     std::vector<int> indices;
     filtered_ugi.getOwnedAndSharedIndices(indices);

     // check that the size didn't grow
     TEST_ASSERT(indices_f.size()<indices.size());

     // sort each set of owned indices, then do a set difference and 
     // crop the difference indices
     std::vector<int> diff;
     diff.resize(indices.size()); 
     std::sort(indices.begin(),indices.end());
     std::sort(indices_f.begin(),indices_f.end());
     std::vector<int>::iterator it = std::set_difference(indices.begin(),indices.end(),
                                                         indices_f.begin(),indices_f.end(),diff.begin());
     diff.resize(it-diff.begin());

     // look at the difference, and make sure they are all in the
     // filtering vector
     TEST_ASSERT(diff.size()<=all_filtered.size());
     for(std::size_t i=0;i<diff.size();i++) {
       TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),diff[i])!=all_filtered.end());
     }

     // make sure no index in the owned index vector is in the filtering vector
     for(std::size_t i=0;i<indices_f.size();i++) {
       TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices_f[i])==all_filtered.end());
     }
   }

}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tFilteredUGI,epetra_lof)
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

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager(myRank,numProc));
   RCP<DOFManager<int,int> > dofManager = rcp(new DOFManager<int,int>); 
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);

   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   dofManager->addField("T",patternC1); // add it to all three blocks
   dofManager->addField("block_0","Ux",patternC1);
   dofManager->addField("block_0","Uy",patternC1);
   dofManager->addField("block_0","P",patternC1);
   dofManager->addField("block_2","rho",patternC1);

   dofManager->buildGlobalUnknowns();

   // get GIDs on the "bottom" side of local cell 0, use those as the filtering
   // GIDs
   std::vector<int> filtered;
   {
     std::pair<std::vector<int>,std::vector<int> > fieldOffsets
         = dofManager->getGIDFieldOffsets_closure("block_0",dofManager->getFieldNum("Ux"),1,0);
 
     std::vector<int> gids;
     dofManager->getElementGIDs(0,gids);

     filtered.resize(fieldOffsets.first.size());
     for(std::size_t i=0;i<fieldOffsets.first.size();i++)
       filtered[i] = gids[fieldOffsets.first[i]];
   }

   TEST_EQUALITY(filtered.size(),2);

   RCP<Filtered_UniqueGlobalIndexer<int,int> > filtered_ugi = rcp(new Filtered_UniqueGlobalIndexer<int,int>);
   filtered_ugi->initialize(dofManager,filtered);

   // check out ownsership
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

   // test out sizes construction by LOF
   {
     std::vector<int> indices_f;
     filtered_ugi->getOwnedIndices(indices_f);

     EpetraLinearObjFactory<panzer::Traits,int> lof(tComm,filtered_ugi);
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof.getThyraDomainSpace())->localSubDim(),Teuchos::as<int>(indices_f.size())); 
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof.getThyraRangeSpace())->localSubDim(),Teuchos::as<int>(indices_f.size())); 

     RCP<Thyra::LinearOpBase<double> > A = lof.getThyraMatrix();
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->range())->localSubDim(),Teuchos::as<int>(indices_f.size())); 
     TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->domain())->localSubDim(),Teuchos::as<int>(indices_f.size())); 

     // next chunk of code tests to ensure parallel communication works as expected.
     // In particular that a filtered unique vector can be used with an unfiltered 
     // ghosted vector and that the entries in the ghosted vector will be preserved
     ////////////////////////////////////////////////////////////////////////////////////

     // check that parallel communication works as expected
     RCP<Epetra_Import> importer = lof.getGhostedImport();
     RCP<Epetra_Map> ghostedMap = lof.getGhostedMap();
     RCP<Epetra_Map> uniqueMap  = lof.getMap();

     Epetra_Vector ghosted_x(*ghostedMap);
     Epetra_Vector unique_x(*uniqueMap);

     ghosted_x.PutScalar(0.0);
     unique_x.PutScalar(1.0);

     ghosted_x.Import(unique_x,*importer,Insert);

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

     
     // do a lazy test to ensure you can construct a unique matrix
     RCP<Thyra::LinearOpBase<double> > uniqueMatrix  = lof.getThyraMatrix();
     TEST_ASSERT(uniqueMatrix!=Teuchos::null);
   }


}


}
