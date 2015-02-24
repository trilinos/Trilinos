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
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "Phalanx_KokkosUtilities.hpp"

using panzer_stk_classic::CoordMatcher;

namespace panzer {

  TEUCHOS_UNIT_TEST(periodic_bcs, sorted_permutation)
  {
    PHX::InitializeKokkosDevice();

     std::vector<double> vec(5.0);
     std::vector<std::size_t> permute;
     vec[0] = 0.0; 
     vec[1] = 4.0; 
     vec[2] = 2.0; 
     vec[3] = 3.0; 
     vec[4] = 1.0; 

     panzer_stk_classic::sorted_permutation(vec,permute);

     TEST_EQUALITY(permute.size(),5);
     for(std::size_t i=0;i<permute.size();i++)  
        TEST_EQUALITY(vec[permute[i]],(double) i);

     PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getSideIdsAndCoords)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // run tests (for both nodes and edges)
    /////////////////////////////////////////////
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<std::size_t> & sideIds_edge = *idsAndCoords_edge.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;
    std::vector<Tuple<double,3> > & sideCoords_edge = *idsAndCoords_edge.second;

    TEST_EQUALITY(sideIds.size(),5);
    TEST_EQUALITY(sideIds_edge.size(),4);

    TEST_EQUALITY(sideCoords.size(),5);
    TEST_EQUALITY(sideCoords_edge.size(),4);

    std::vector<std::size_t> permute;
    panzer_stk_classic::sorted_permutation(sideIds,permute);
    for(std::size_t i=0;i<permute.size();i++) {
       std::size_t p = permute[i];

       TEST_EQUALITY(sideIds[p]-1,i*13);

       TEST_FLOATING_EQUALITY(sideCoords[p][0],0.0,1e-14);
       TEST_FLOATING_EQUALITY(sideCoords[p][1],i*1.0/4.0,1e-14);
    }

    // Rather than compute the edge numbers, we can do everything in terms of node numbers
    //   For any edge, find the two nodes associated with it using its coordinates
    //   Use their numbers for comparison, since we know what those should be

    std::vector<std::size_t> permute_edge;
    panzer_stk_classic::sorted_permutation(sideIds_edge,permute_edge);
    for(std::size_t i=0;i<permute_edge.size();i++) {
       std::size_t p = permute_edge[i];

       Tuple<double,3> coord = sideCoords_edge[p]; // coordinate of edge

       int node_l = -1;   // number for node below edge
       int node_u = -1;   // number for node above edge
       int flag   = 0;    // flag should be 2 after this loop

       for(int j=0;j<sideCoords.size();j++) {
         if ((std::abs(sideCoords[j][0] - coord[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord[1]-1.0/8.0)) < 1e-14)){
            node_l = j;
            flag++;
         }
         if ((std::abs(sideCoords[j][0] - coord[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord[1]+1.0/8.0)) < 1e-14)){
            node_u = j;
            flag++;
         }
       }

       TEST_EQUALITY(flag,2);
       TEST_EQUALITY(sideIds[node_l]-1,i*13);
       TEST_EQUALITY(sideIds[node_u]-1,(i+1)*13);

       TEST_FLOATING_EQUALITY(sideCoords_edge[p][0],0.0,1e-14);
       TEST_FLOATING_EQUALITY(sideCoords_edge[p][1],i*1.0/4.0+1.0/8.0,1e-14);
    }

    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getLocalSideIds)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",8);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"left");
    RCP<std::vector<std::size_t> > locallyRequiredIds_edge = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"left","edge");
    if(rank==0) {
       TEST_EQUALITY(locallyRequiredIds->size(),2);
       TEST_EQUALITY(locallyRequiredIds_edge->size(),1);
    }
    else {
       TEST_EQUALITY(locallyRequiredIds->size(),0);
       TEST_EQUALITY(locallyRequiredIds_edge->size(),0);
    }

    PHX::FinalizeKokkosDevice();
  }
 
  TEUCHOS_UNIT_TEST(periodic_bcs, getLocallyMatchedSideIds)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 
    int procCnt = Comm.NumProc(); 

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }
 
    // run tests
    /////////////////////////////////////////////
    
    // We need Ids and Coords on right so we can map edges to nodes
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
    std::vector<std::size_t> & sideIds_edge = *idsAndCoords_edge.first;
    std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;
    std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
    std::vector<Tuple<double,3> > & sideCoords_edge = *idsAndCoords_edge.second;
    std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

    CoordMatcher matcher(1);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedIds 
          = panzer_stk_classic::periodic_helpers::getLocallyMatchedSideIds(sideIds,sideCoords,*mesh,"right",matcher);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedIds_edge 
          = panzer_stk_classic::periodic_helpers::getLocallyMatchedSideIds(sideIds_edge,sideCoords_edge,*mesh,"right",matcher,"edge");

    if(rank==procCnt-1) {
       for(std::size_t i=0;i<matchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*matchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }

       std::vector<std::size_t> permute_edge;
       panzer_stk_classic::sorted_permutation(sideIds_edge,permute_edge);
       for(std::size_t i=0;i<matchedIds_edge->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*matchedIds_edge)[i];

          // Get coordinates for matched edges
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag0 = 0;
          for(std::size_t j=0;j<sideIds.size();j++){
            if(pair.first == sideIds_edge[j]){
              coord_left = sideCoords_edge[j];
              flag0++;
            }
            if(pair.second == sideIds_edge_right[j]){
              coord_right = sideCoords_edge_right[j];
              flag0++;
            }
          }
          TEST_EQUALITY(flag0,2);

          // Get coordinates of associated nodes
          int left_l  = -1;   
          int left_u  = -1; 
          int right_l = -1;   
          int right_u = -1; 
          int flag1   = 0;

          for(std::size_t j=0;j<sideCoords.size();j++) {
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]-1.0/8.0)) < 1e-14)){
               left_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]+1.0/8.0)) < 1e-14)){
               left_u = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/8.0)) < 1e-14)){
               right_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/8.0)) < 1e-14)){
               right_u = j;
               flag1++;
            }
          }
          TEST_EQUALITY(flag1,4);

          // Test equivalence of node numbers
          TEST_EQUALITY(sideIds[left_l],sideIds_right[right_l]-12);
          TEST_EQUALITY(sideIds[left_u],sideIds_right[right_u]-12);
       }
    }
    else {
       TEST_EQUALITY(matchedIds->size(),0)
    }
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getGlobalPairing)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // run tests
    /////////////////////////////////////////////
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > locallyMatchedIds;
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > locallyMatchedIds_edge;
    
       // next line requires global communication
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds = *idsAndCoords.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge = *idsAndCoords_edge.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge = *idsAndCoords_edge.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

       CoordMatcher matcher(1);
       locallyMatchedIds = panzer_stk_classic::periodic_helpers::getLocallyMatchedSideIds(sideIds,sideCoords,*mesh,"right",matcher);
       locallyMatchedIds_edge = panzer_stk_classic::periodic_helpers::getLocallyMatchedSideIds(sideIds_edge,sideCoords_edge,*mesh,"right",matcher,"edge");
    

    Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"left");
    Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds_edge = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"left","edge");

    // next line requires communication
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
          = panzer_stk_classic::periodic_helpers::getGlobalPairing(*locallyRequiredIds,*locallyMatchedIds,*mesh,false);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
          = panzer_stk_classic::periodic_helpers::getGlobalPairing(*locallyRequiredIds_edge,*locallyMatchedIds_edge,*mesh,false);

    if(rank==0) {
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

          // Get coordinates for matched edges
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag0 = 0;
          for(std::size_t j=0;j<sideIds.size();j++){
            if(pair.first == sideIds_edge[j]){
              coord_left = sideCoords_edge[j];
              flag0++;
            }
            if(pair.second == sideIds_edge_right[j]){
              coord_right = sideCoords_edge_right[j];
              flag0++;
            }
          }
          TEST_EQUALITY(flag0,2);

          // Get coordinates of associated nodes
          int left_l  = -1;   
          int left_u  = -1; 
          int right_l = -1;   
          int right_u = -1; 
          int flag1   = 0;

          for(std::size_t j=0;j<sideCoords.size();j++) {
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]-1.0/8.0)) < 1e-14)){
               left_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords[j][1] - (coord_left[1]+1.0/8.0)) < 1e-14)){
               left_u = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/8.0)) < 1e-14)){
               right_l = j;
               flag1++;
            }
            if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/8.0)) < 1e-14)){
               right_u = j;
               flag1++;
            }
          }
          TEST_EQUALITY(flag1,4);

          // Test equivalence of node numbers
          TEST_EQUALITY(sideIds[left_l],sideIds_right[right_l]-12);
          TEST_EQUALITY(sideIds[left_u],sideIds_right[right_u]-12);
       }
    }
    else {
       TEST_EQUALITY(globallyMatchedIds->size(),0);
    }
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, matchPeriodicSides)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    // run tests

    // Nodes
    {
       CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk_classic::periodic_helpers::matchPeriodicSides("left","right",*mesh,matcher);
   
       // match left & right sides
       if(rank==0) {
          for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
             TEST_EQUALITY(pair.first,pair.second-12);
          }
       }
       else
          TEST_EQUALITY(globallyMatchedIds->size(),0);
    }

    // Edges
    {
       CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
             = panzer_stk_classic::periodic_helpers::matchPeriodicSides("left","right",*mesh,matcher,"edge");
   
       // match left & right sides
       if(rank==0) {
          for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_left;
             Tuple<double,3> coord_right;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_left.size();j++){
               if(pair.first == sideIds_edge_left[j]){
                 coord_left = sideCoords_edge_left[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_right[j]){
                 coord_right = sideCoords_edge_right[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int left_l  = -1;   
             int left_u  = -1; 
             int right_l = -1;   
             int right_u = -1; 
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/8.0)) < 1e-14)){
                  left_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/8.0)) < 1e-14)){
                  left_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/8.0)) < 1e-14)){
                  right_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/8.0)) < 1e-14)){
                  right_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-12);
             TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-12);
          }
       }
       else
          TEST_EQUALITY(globallyMatchedIds_edge->size(),0);
    }

    // Nodes
    {
       CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","bottom",*mesh,matcher);

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 
   
       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // Edges
    {
       CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge
             = panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","bottom",*mesh,matcher,"edge");

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"top","edge");

       TEST_EQUALITY(globallyMatchedIds_edge->size(),locallyRequiredIds->size()); 
   
       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_top;
             Tuple<double,3> coord_bottom;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_top.size();j++){
               if(pair.first == sideIds_edge_top[j]){
                 coord_top = sideCoords_edge_top[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_bottom[j]){
                 coord_bottom = sideCoords_edge_bottom[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int top_l  = -1;   
             int top_u  = -1; 
             int bottom_l = -1;   
             int bottom_u = -1; 
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+52);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+52);
       }
    }


    // test a failure case!
    {
       CoordMatcher matcherX(0);
       CoordMatcher matcherY(1);
       
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("left","bottom",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherY),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("bottom","left",*mesh,matcherY),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("left","bottom",*mesh,matcherX,"edge"),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherY,"edge"),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherX,"edge"),std::logic_error);
       TEST_THROW(panzer_stk_classic::periodic_helpers::matchPeriodicSides("bottom","left",*mesh,matcherY,"edge"),std::logic_error);
    }
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    // Nodes
    {
       CoordMatcher matcher(0);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> pMatch 
             = panzer_stk_classic::buildPeriodicBC_Matcher("top","bottom",matcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds = pMatch->getMatchedPair(*mesh);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // Edges
    {
       CoordMatcher matcher(0);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> pMatch 
             = panzer_stk_classic::buildPeriodicBC_Matcher("top","bottom",matcher,"edge");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge = pMatch->getMatchedPair(*mesh);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk_classic::periodic_helpers::getLocalSideIds(*mesh,"top","edge");

       TEST_EQUALITY(globallyMatchedIds_edge->size(),locallyRequiredIds->size()); 

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

             // Get coordinates for matched edges
             Tuple<double,3> coord_top;
             Tuple<double,3> coord_bottom;
             int flag0 = 0;
             for(std::size_t j=0;j<sideIds_top.size();j++){
               if(pair.first == sideIds_edge_top[j]){
                 coord_top = sideCoords_edge_top[j];
                 flag0++;
               }
               if(pair.second == sideIds_edge_bottom[j]){
                 coord_bottom = sideCoords_edge_bottom[j];
                 flag0++;
               }
             }
             TEST_EQUALITY(flag0,2);

             // Get coordinates of associated nodes
             int top_l  = -1;   
             int top_u  = -1; 
             int bottom_l = -1;   
             int bottom_u = -1; 
             int flag1   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag1++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/24.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag1++;
               }
             }
             TEST_EQUALITY(flag1,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+52);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+52);
       }
    }

    // test a failure case!
    {
       CoordMatcher matcherX(0);
       CoordMatcher matcherY(1);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> pMatch;

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("left","bottom",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("top","right",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("top","right",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("bottom","left",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("left","bottom",matcherX,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("top","right",matcherX,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("top","right",matcherY,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk_classic::buildPeriodicBC_Matcher("bottom","left",matcherY,"edge");
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);
    }
    
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher_multi)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    {
       CoordMatcher xmatcher(0);
       CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> tb_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("top","bottom",xmatcher);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> lr_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("left","right",ymatcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       globallyMatchedIds = tb_Match->getMatchedPair(*mesh);
       globallyMatchedIds = lr_Match->getMatchedPair(*mesh,globallyMatchedIds);

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          out << pair.first << " " << pair.second << std::endl;
       }

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          if(pair.first==1 || pair.first==19)
          {   TEST_EQUALITY(pair.second,9); }
          else if(pair.first==10)
          {   TEST_EQUALITY(pair.second,18); }
          else
          {   TEST_EQUALITY(pair.second,pair.first-18); }
       }

    }
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher_multi_edge)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_top = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_bottom = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"bottom","edge");
       std::vector<std::size_t> & sideIds_top = *idsAndCoords_top.first;
       std::vector<std::size_t> & sideIds_bottom = *idsAndCoords_bottom.first;
       std::vector<std::size_t> & sideIds_edge_top = *idsAndCoords_edge_top.first;
       std::vector<std::size_t> & sideIds_edge_bottom = *idsAndCoords_edge_bottom.first;
       std::vector<Tuple<double,3> > & sideCoords_top = *idsAndCoords_top.second;
       std::vector<Tuple<double,3> > & sideCoords_bottom = *idsAndCoords_bottom.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_top = *idsAndCoords_edge_top.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_bottom = *idsAndCoords_edge_bottom.second;

    {
       CoordMatcher xmatcher(0);
       CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> tb_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("top","bottom",xmatcher,"edge");
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> lr_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("left","right",ymatcher,"edge");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds_edge;
       globallyMatchedIds_edge = tb_Match->getMatchedPair(*mesh);
       globallyMatchedIds_edge = lr_Match->getMatchedPair(*mesh,globallyMatchedIds_edge);

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds_edge->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds_edge)[i];

          // Get coordinates for matched edges on top and bottom
          Tuple<double,3> coord_top;
          Tuple<double,3> coord_bottom;
          int flag_tb = 0;
          for(std::size_t j=0;j<sideIds_top.size();j++){
            if(pair.first == sideIds_edge_top[j]){
              coord_top = sideCoords_edge_top[j];
              flag_tb++;
            }
            if(pair.second == sideIds_edge_bottom[j]){
              coord_bottom = sideCoords_edge_bottom[j];
              flag_tb++;
            }
          }

          // Get coordinates for matched edges on left and right
          Tuple<double,3> coord_left;
          Tuple<double,3> coord_right;
          int flag_lr = 0;
          for(std::size_t j=0;j<sideIds_left.size();j++){
            if(pair.first == sideIds_edge_left[j]){
              coord_left = sideCoords_edge_left[j];
              flag_lr++;
            }
            if(pair.second == sideIds_edge_right[j]){
              coord_right = sideCoords_edge_right[j];
              flag_lr++;
            }
          }
          TEST_EQUALITY(flag_tb+flag_lr,2);

          // If node is on top
          if(flag_tb == 2) {
             int top_l  = -1;   
             int top_u  = -1; 
             int bottom_l = -1;   
             int bottom_u = -1; 
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_top.size();j++) {
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_top[j][0] - (coord_top[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_top[j][1] - coord_top[1]) < 1e-14)){
                  top_u = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]-1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_bottom[j][0] - (coord_bottom[0]+1.0/16.0)) < 1e-14) && (std::abs(sideCoords_bottom[j][1] - coord_bottom[1]) < 1e-14)){
                  bottom_u = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_top[top_l],sideIds_bottom[bottom_l]+18);
             TEST_EQUALITY(sideIds_top[top_u],sideIds_bottom[bottom_u]+18);
          }
          // If node is on left
          else if(flag_lr == 2) {
             int left_l  = -1;   
             int left_u  = -1; 
             int right_l = -1;   
             int right_u = -1; 
             int flag   = 0;

             for(std::size_t j=0;j<sideCoords_left.size();j++) {
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14)){
                  left_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14)){
                  left_u = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14)){
                  right_l = j;
                  flag++;
               }
               if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14)){
                  right_u = j;
                  flag++;
               }
             }
             TEST_EQUALITY(flag,4);

             // Test equivalence of node numbers
             TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-8);
             TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-8);
          }

       }

    }
    PHX::FinalizeKokkosDevice();
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher_nodes_and_edges)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    PHX::InitializeKokkosDevice();

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk_classic::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk_classic::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_left = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"left","edge");
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords_edge_right = panzer_stk_classic::periodic_helpers::getSideIdsAndCoords(*mesh,"right","edge");
       std::vector<std::size_t> & sideIds_left = *idsAndCoords_left.first;
       std::vector<std::size_t> & sideIds_right = *idsAndCoords_right.first;
       std::vector<std::size_t> & sideIds_edge_left = *idsAndCoords_edge_left.first;
       std::vector<std::size_t> & sideIds_edge_right = *idsAndCoords_edge_right.first;
       std::vector<Tuple<double,3> > & sideCoords_left = *idsAndCoords_left.second;
       std::vector<Tuple<double,3> > & sideCoords_right = *idsAndCoords_right.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_left = *idsAndCoords_edge_left.second;
       std::vector<Tuple<double,3> > & sideCoords_edge_right = *idsAndCoords_edge_right.second;


    {
       CoordMatcher ymatcher(1);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> node_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("left","right",ymatcher);
       Teuchos::RCP<const panzer_stk_classic::PeriodicBC_MatcherBase> edge_Match 
             = panzer_stk_classic::buildPeriodicBC_Matcher("left","right",ymatcher,"edge");

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds;
       globallyMatchedIds = node_Match->getMatchedPair(*mesh);
       globallyMatchedIds = edge_Match->getMatchedPair(*mesh,globallyMatchedIds);


       // match left & right sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {

          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];

          // Is this a node or edge pairing?
          if(pair.first < 28){ // Node
            TEST_EQUALITY(pair.first,pair.second-8);
          } else {             //Edge

            // Get coordinates for matched edges on left and right
            Tuple<double,3> coord_left;
            Tuple<double,3> coord_right;
            int flag_lr = 0;
            for(std::size_t j=0;j<sideIds_edge_left.size();j++){
              if(pair.first == sideIds_edge_left[j]){
                coord_left = sideCoords_edge_left[j];
                flag_lr++;
              }
              if(pair.second == sideIds_edge_right[j]){
                coord_right = sideCoords_edge_right[j];
                flag_lr++;
              }
            }
            TEST_EQUALITY(flag_lr,2);

            int left_l  = -1;   
            int left_u  = -1; 
            int right_l = -1;   
            int right_u = -1; 
            int flag   = 0;

            for(std::size_t j=0;j<sideCoords_left.size();j++) {
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]-1.0/4.0)) < 1e-14)){
                 left_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_left[j][0] - coord_left[0]) < 1e-14) && (std::abs(sideCoords_left[j][1] - (coord_left[1]+1.0/4.0)) < 1e-14)){
                 left_u = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]-1.0/4.0)) < 1e-14)){
                 right_l = j;
                 flag++;
              }
              if ((std::abs(sideCoords_right[j][0] - coord_right[0]) < 1e-14) && (std::abs(sideCoords_right[j][1] - (coord_right[1]+1.0/4.0)) < 1e-14)){
                 right_u = j;
                 flag++;
              }
            }
            TEST_EQUALITY(flag,4);
            // Test equivalence of node numbers
            TEST_EQUALITY(sideIds_left[left_l],sideIds_right[right_l]-8);
            TEST_EQUALITY(sideIds_left[left_u],sideIds_right[right_u]-8);
          }
       }

    }
    
    PHX::FinalizeKokkosDevice();
  }

}
