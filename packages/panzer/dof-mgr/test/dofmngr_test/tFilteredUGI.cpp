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
#include <unordered_set>

#include "Panzer_DOFManager.hpp"
#include "Panzer_Filtered_GlobalIndexer.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_ConnManager.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Kokkos_DynRankView.hpp"

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
TEUCHOS_UNIT_TEST(tFilteredUGI,equivalence_test)
{
   RCP<const Teuchos::MpiComm<int> > tComm 
     = rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = tComm->getRank();
   int numProc = tComm->getSize();   

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

   std::vector<panzer::GlobalOrdinal> filtered;
   Filtered_GlobalIndexer filtered_ugi;
   filtered_ugi.initialize(dofManager,filtered);

   // check the GIDs
   {
     std::vector<panzer::GlobalOrdinal> gids,gids_f;

     dofManager->getElementGIDs(0,gids);
     filtered_ugi.getElementGIDs(0,gids_f);

     TEST_EQUALITY(gids.size(),gids_f.size());
     for(std::size_t i=0;i<gids.size();i++) {
       TEST_EQUALITY(gids[i],gids_f[i]);
     }
   }

   // check the LIDs
   {
     PHX::View<const int*> lids = dofManager->getElementLIDs(0);
     PHX::View<const int*> lids_f = filtered_ugi.getElementLIDs(0);

     auto lids_host = Kokkos::create_mirror_view(lids);
     Kokkos::deep_copy(lids_host,lids);
     auto lids_f_host = Kokkos::create_mirror_view(lids_f);
     Kokkos::deep_copy(lids_f_host,lids_f);

     TEST_EQUALITY(lids.size(),lids_f.size());
     for(std::size_t i=0;i<lids.size();i++) {
       TEST_EQUALITY(lids_host(i),lids_f_host(i));
     }
   }

   // check owned and ghosted
   {
     std::vector<panzer::GlobalOrdinal> indices, indices_f;

     dofManager->getOwnedAndGhostedIndices(indices);
     filtered_ugi.getOwnedAndGhostedIndices(indices_f);

     TEST_EQUALITY(indices.size(),indices_f.size());
     for(std::size_t i=0;i<indices.size();i++) {
       TEST_EQUALITY(indices[i],indices_f[i]);
     }
   }

   // check owned 
   {
     std::vector<panzer::GlobalOrdinal> indices, indices_f;

     dofManager->getOwnedIndices(indices);
     filtered_ugi.getOwnedIndices(indices_f);

     TEST_EQUALITY(indices.size(),indices_f.size());
     for(std::size_t i=0;i<indices.size();i++) {
       TEST_EQUALITY(indices[i],indices_f[i]);
     }
   }

  // check ghosted
  {
    std::vector<panzer::GlobalOrdinal> indices, indices_f;
    dofManager->getGhostedIndices(indices);
    filtered_ugi.getGhostedIndices(indices_f);
    TEST_EQUALITY(indices.size(), indices_f.size());
    for (std::size_t i(0); i < indices.size(); ++i)
      TEST_EQUALITY(indices[i], indices_f[i])
  }
}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tFilteredUGI,filtering)
{
   RCP<const Teuchos::MpiComm<int> > tComm 
     = rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = tComm->getRank(); 
   int numProc = tComm->getSize(); 

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

   std::vector<panzer::GlobalOrdinal> my_filtered;
   dofManager->getElementGIDs(1,my_filtered);

   // now we will compute all the filtered indices in this UGI
   // here we will print some things out for a sanity check
   std::vector<panzer::GlobalOrdinal> all_filtered;
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

   Filtered_GlobalIndexer filtered_ugi;
   filtered_ugi.initialize(dofManager,my_filtered);

   out << "check the GIDs" << std::endl;
   {
     std::vector<panzer::GlobalOrdinal> gids,gids_f;

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

  out << "check owned and ghosted" << std::endl;
  {
    using std::size_t;
    using std::unordered_set;
    using std::vector;
    vector<panzer::GlobalOrdinal> indices, filteredIndices;
    dofManager->getOwnedAndGhostedIndices(indices);
    filtered_ugi.getOwnedAndGhostedIndices(filteredIndices);
    TEST_EQUALITY(indices.size(), filteredIndices.size());
    unordered_set<panzer::GlobalOrdinal> indicesSet;
    for (size_t i(0); i < indices.size(); ++i)
      indicesSet.insert(indices[i]);
    for (size_t i(0); i < filteredIndices.size(); ++i)
      TEST_EQUALITY(indicesSet.count(filteredIndices[i]), 1)
  }

   out << "check owned" << std::endl;
   {
     std::vector<panzer::GlobalOrdinal> indices, indices_f,diff;

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
     std::vector<panzer::GlobalOrdinal>::iterator it = 
       std::set_difference(indices.begin(),indices.end(),indices_f.begin(),indices_f.end(),diff.begin());
     diff.resize(it-diff.begin());

     // look at the difference, and make sure they are all in the
     // filtering vector
     TEST_ASSERT(diff.size()<=all_filtered.size());
     for(std::size_t i=0;i<diff.size();i++) {
       TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),diff[i])!=my_filtered.end());
     }

     // make sure no index in the owned index vector is in the filtering vector
     for(std::size_t i=0;i<indices_f.size();i++) {
       TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices_f[i])==all_filtered.end());
     }
   }

  out << "Check the ghosted indices." << std::endl;
  {
    using std::find;
    using std::set_difference;
    using std::size_t;
    using std::sort;
    using std::vector;
    vector<panzer::GlobalOrdinal> indices, filteredIndices, diff;
    dofManager->getGhostedIndices(indices);
    filtered_ugi.getGhostedIndices(filteredIndices);

    // Check that the size didn't shrink.
    TEST_ASSERT(filteredIndices.size() >= indices.size())

    // Sort each set of ghosted indices, and then do a set difference and 
    // crop the difference indices.
    diff.resize(filteredIndices.size()); 
    sort(indices.begin(), indices.end());
    sort(filteredIndices.begin(), filteredIndices.end());
    vector<panzer::GlobalOrdinal>::iterator it = set_difference(filteredIndices.begin(),
      filteredIndices.end(), indices.begin(), indices.end(), diff.begin());
    diff.resize(it - diff.begin());

    // Look at the difference, and make sure they are all in the filtering
    // vector.
    TEST_ASSERT(diff.size() <= all_filtered.size())
    for (size_t i(0); i < diff.size(); ++i)
      TEST_ASSERT(find(all_filtered.begin(), all_filtered.end(), diff[i]) !=
        all_filtered.end())
  }

   out << "testing getOwnedAndGhostedNotFilteredIndicator" << std::endl;
   {
     std::vector<int> indicator;
     filtered_ugi.getOwnedAndGhostedNotFilteredIndicator(indicator);

     std::vector<panzer::GlobalOrdinal> indices;
     filtered_ugi.getOwnedAndGhostedIndices(indices);

     TEST_EQUALITY(indices.size(),indicator.size());

     for(std::size_t i=0;i<indicator.size();i++) {
       if(indicator[i]==1) {
         out << "Searching (1) " << i << " " << indices[i] << " ";
         // not filtered, should not be in all_filtered array
         TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices[i])==all_filtered.end());
       }
       else if(indicator[i]==0) {
         out << "Searching (0) " << i << " " << indices[i] << " ";
         // filtered, should be in all_filtered array
         TEST_ASSERT(std::find(all_filtered.begin(),all_filtered.end(),indices[i])!=all_filtered.end());
       }
       else {
         // protect that neither one nor zero is set (also report a useful error)
         TEST_ASSERT(indicator[i]==0 || indicator[i]==1);
       }
     }
   }

   out << "testing getFilteredOwnedAndGhostedIndices" << std::endl;
   {
     std::vector<panzer::GlobalOrdinal> indices_f;
     filtered_ugi.getFilteredOwnedAndGhostedIndices(indices_f);

     std::vector<panzer::GlobalOrdinal> indices;
     filtered_ugi.getOwnedAndGhostedIndices(indices);

     // check that the size didn't grow
     TEST_ASSERT(indices_f.size()<indices.size());

     // sort each set of owned indices, then do a set difference and 
     // crop the difference indices
     std::vector<panzer::GlobalOrdinal> diff;
     diff.resize(indices.size()); 
     std::sort(indices.begin(),indices.end());
     std::sort(indices_f.begin(),indices_f.end());
     std::vector<panzer::GlobalOrdinal>::iterator it =
       std::set_difference(indices.begin(),indices.end(),indices_f.begin(),indices_f.end(),diff.begin());
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

}
