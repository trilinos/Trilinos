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

using panzer_stk::CoordMatcher;

namespace panzer {

  TEUCHOS_UNIT_TEST(periodic_bcs, sorted_permutation)
  {
     std::vector<double> vec(5.0);
     std::vector<std::size_t> permute;
     vec[0] = 0.0; 
     vec[1] = 4.0; 
     vec[2] = 2.0; 
     vec[3] = 3.0; 
     vec[4] = 1.0; 

     panzer_stk::sorted_permutation(vec,permute);

     TEST_EQUALITY(permute.size(),5);
     for(std::size_t i=0;i<permute.size();i++)  
        TEST_EQUALITY(vec[permute[i]],(double) i);
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getSideIdsAndCoords)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
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
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

    TEST_EQUALITY(sideIds.size(),5);

    TEST_EQUALITY(sideCoords.size(),5);

    std::vector<std::size_t> permute;
    panzer_stk::sorted_permutation(sideIds,permute);
    for(std::size_t i=0;i<permute.size();i++) {
       std::size_t p = permute[i];

       TEST_EQUALITY(sideIds[p]-1,i*13);

       TEST_FLOATING_EQUALITY(sideCoords[p][0],0.0,1e-14);
       TEST_FLOATING_EQUALITY(sideCoords[p][1],i*1.0/4.0,1e-14);
    }

  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getLocalSideIds)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",8);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"left");
    if(rank==0) {
       TEST_EQUALITY(locallyRequiredIds->size(),2);
    }
    else {
       TEST_EQUALITY(locallyRequiredIds->size(),0);
    }
  }
 
  TEUCHOS_UNIT_TEST(periodic_bcs, getLocallyMatchedSideIds)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 
    int procCnt = Comm.NumProc(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
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
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

    CoordMatcher matcher(1);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedIds 
          = panzer_stk::periodic_helpers::getLocallyMatchedSideIds(sideIds,sideCoords,*mesh,"right",matcher);

    if(rank==procCnt-1) {
       for(std::size_t i=0;i<matchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*matchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
    }
    else {
       TEST_EQUALITY(matchedIds->size(),0)
    }
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getGlobalPairing)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
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
    {
       // next line requires global communication
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"left");
       std::vector<std::size_t> & sideIds = *idsAndCoords.first;
       std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

       CoordMatcher matcher(1);
       locallyMatchedIds = panzer_stk::periodic_helpers::getLocallyMatchedSideIds(sideIds,sideCoords,*mesh,"right",matcher);
    }

    Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"left");

    // next line requires communication
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
          = panzer_stk::periodic_helpers::getGlobalPairing(*locallyRequiredIds,*locallyMatchedIds,*mesh,false);

    if(rank==0) {
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
    }
    else {
       TEST_EQUALITY(globallyMatchedIds->size(),0);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, matchPeriodicSides)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
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
    {
       CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk::periodic_helpers::matchPeriodicSides("left","right",*mesh,matcher);
   
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

    {
       CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = panzer_stk::periodic_helpers::matchPeriodicSides("top","bottom",*mesh,matcher);

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 
   
       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // test a failure case!
    {
       CoordMatcher matcherX(0);
       CoordMatcher matcherY(1);
       
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSides("left","bottom",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherY),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSides("top","right",*mesh,matcherX),std::logic_error);
       TEST_THROW(panzer_stk::periodic_helpers::matchPeriodicSides("bottom","left",*mesh,matcherY),std::logic_error);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    {
       CoordMatcher matcher(0);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch 
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",matcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds = pMatch->getMatchedPair(*mesh);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // test a failure case!
    {
       CoordMatcher matcherX(0);
       CoordMatcher matcherY(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch;

       pMatch = panzer_stk::buildPeriodicBC_Matcher("left","bottom",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherX);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("top","right",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);

       pMatch = panzer_stk::buildPeriodicBC_Matcher("bottom","left",matcherY);
       TEST_THROW(pMatch->getMatchedPair(*mesh),std::logic_error);
    }
    
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, PeriodicBC_Matcher_multi)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
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
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> tb_Match 
             = panzer_stk::buildPeriodicBC_Matcher("top","bottom",xmatcher);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> lr_Match 
             = panzer_stk::buildPeriodicBC_Matcher("left","right",ymatcher);

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

  }

}
