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

// STL includes
#include <iostream>
#include <vector>
#include <set>

// Teuchos includes
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

//Teuchos Testing
#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_Tuple.hpp>

//Tpetra includes
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"

// Intrepid includes
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"

// Panzer includes
#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"

// Panzer_STK includes
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Panzer_SingleBlockDOFManager.cpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#endif

#include <iostream>

typedef int LO;
typedef int GO;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::rcp;
typedef Intrepid::FieldContainer<double> FieldContainer;

namespace {

  TEUCHOS_UNIT_TEST( SingleBlockDOFManager_Test, BasicCreation ){
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",5);
    pl->set("Y Elements",5);
    
    panzer_stk::SquareQuadMeshFactory factory; 
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    RCP<panzer::SingleBlockDOFManager<LO,GO> > my_SingleBlockDOFManager = Teuchos::rcp(new panzer::SingleBlockDOFManager<LO,GO>());
    TEST_EQUALITY(my_SingleBlockDOFManager->getComm(),Teuchos::null);

    RCP<panzer::ConnManager<LO,GO> > conn = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<Intrepid::Basis<double,FieldContainer> > basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer>);
     
    RCP< const panzer::FieldPattern> pattern = Teuchos::rcp(new panzer::IntrepidFieldPattern(basis));

    std::vector<std::string> names;
    names.push_back("Velocity");
    names.push_back("Temperature");
    names.push_back("Radiation Levels");
    my_SingleBlockDOFManager->addField(names[0], pattern);
    my_SingleBlockDOFManager->addField(names[1], pattern);
    my_SingleBlockDOFManager->addField(names[2], pattern);

    my_SingleBlockDOFManager->setConnManager(conn, MPI_COMM_WORLD);

    my_SingleBlockDOFManager->buildGlobalUnknowns();
    //Now that we have created the SingleBlockDOFManager, we can ensure it was created correctly.
    TEST_EQUALITY(my_SingleBlockDOFManager->getConnManager(),conn);


    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Velocity"),0);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Temperature"),1);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Radiation Levels"),2);

    TEST_EQUALITY(my_SingleBlockDOFManager->getNumFields(), 3);

    const std::vector<GO> & vel_offests = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",0); 
    const std::vector<GO> & tem_offests = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",1); 
    const std::vector<GO> & rad_offests = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",2); 

    TEST_EQUALITY(vel_offests.size(),tem_offests.size());
    TEST_EQUALITY(tem_offests.size(),rad_offests.size());
  }

  TEUCHOS_UNIT_TEST( SingleBlockDOFManager_Test, ReorderingFields ){
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",5);
    pl->set("Y Elements",5);
    
    panzer_stk::SquareQuadMeshFactory factory; 
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    RCP<panzer::SingleBlockDOFManager<LO,GO> > my_SingleBlockDOFManager = Teuchos::rcp(new panzer::SingleBlockDOFManager<LO,GO>());
    TEST_EQUALITY(my_SingleBlockDOFManager->getComm(),Teuchos::null);

    RCP<panzer::ConnManager<LO,GO> > conn = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<Intrepid::Basis<double,FieldContainer> > basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer>);
     
    RCP< const panzer::FieldPattern> pattern = Teuchos::rcp(new panzer::IntrepidFieldPattern(basis));

    my_SingleBlockDOFManager->setConnManager(conn, MPI_COMM_WORLD);

    std::vector<std::string> names;
    names.push_back("Velocity");
    names.push_back("Temperature");
    names.push_back("Radiation Levels");
    names.push_back("Frequency");
    names.push_back("Pitch");
    names.push_back("Yaw");
    names.push_back("Brightness");
    names.push_back("Momentum");
    my_SingleBlockDOFManager->addField(names[0], pattern);
    my_SingleBlockDOFManager->addField(names[1], pattern);
    my_SingleBlockDOFManager->addField(names[2], pattern);
    my_SingleBlockDOFManager->addField(names[3], pattern);
    my_SingleBlockDOFManager->addField(names[4], pattern);
    my_SingleBlockDOFManager->addField(names[5], pattern);
    my_SingleBlockDOFManager->addField(names[6], pattern);
    my_SingleBlockDOFManager->addField(names[7], pattern);


    TEST_EQUALITY(my_SingleBlockDOFManager->getNumFields(),8);

    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Velocity"),0);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Temperature"),1);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Radiation Levels"),2);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Frequency"),3);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Pitch"),4);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Yaw"),5);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Brightness"),6);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Momentum"),7);


    std::vector<std::string> invalid_names;
    invalid_names.push_back("Verocity");
    invalid_names.push_back("Temperature");
    invalid_names.push_back("Irradiation Levels");
    invalid_names.push_back("Infrequency");
    invalid_names.push_back("Off-Pitch");
    invalid_names.push_back("Yawn");
    invalid_names.push_back("Brightness");
    invalid_names.push_back("Momentum");

    TEST_ASSERT(my_SingleBlockDOFManager->validFieldOrder(names));
    TEST_ASSERT(!(my_SingleBlockDOFManager->validFieldOrder(invalid_names)));
    TEST_THROW(my_SingleBlockDOFManager->setFieldOrder(invalid_names),std::logic_error);

    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Velocity"),0);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Temperature"),1);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Radiation Levels"),2);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Frequency"),3);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Pitch"),4);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Yaw"),5);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Brightness"),6);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Momentum"),7);

    std::vector<std::string> new_names;
    new_names.push_back("Momentum");
    new_names.push_back("Frequency");
    new_names.push_back("Velocity");
    new_names.push_back("Temperature");
    new_names.push_back("Pitch");
    new_names.push_back("Brightness");
    new_names.push_back("Yaw");
    new_names.push_back("Radiation Levels");

    TEST_ASSERT(my_SingleBlockDOFManager->validFieldOrder(new_names));
    TEST_NOTHROW(my_SingleBlockDOFManager->setFieldOrder(new_names));

    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Velocity"),2);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Temperature"),3);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Radiation Levels"),7);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Frequency"),1);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Pitch"),4);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Yaw"),6);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Brightness"),5);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Momentum"),0);

    TEST_NOTHROW(my_SingleBlockDOFManager->addField("Current",pattern));

    TEST_EQUALITY(my_SingleBlockDOFManager->getNumFields(),9);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Current"),8);

    TEST_NOTHROW(my_SingleBlockDOFManager->buildGlobalUnknowns());

    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Current"),8);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Velocity"),2);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Temperature"),3);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Radiation Levels"),7);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Frequency"),1);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Pitch"),4);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Yaw"),6);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Brightness"),5);
    TEST_EQUALITY(my_SingleBlockDOFManager->getFieldNum("Momentum"),0);

    new_names.push_back("Current");

    //You can't reassign anything.
    TEST_THROW(my_SingleBlockDOFManager->setFieldOrder(new_names),std::logic_error);
    TEST_THROW(my_SingleBlockDOFManager->buildGlobalUnknowns(),std::logic_error);
    TEST_THROW(my_SingleBlockDOFManager->addField("Tide",pattern),std::logic_error);
    TEST_THROW(my_SingleBlockDOFManager->setConnManager(conn,MPI_COMM_WORLD),std::logic_error);

    //All of the fields shoudl be in this block.
    TEST_ASSERT(my_SingleBlockDOFManager->fieldInBlock("Pitch","eblock-0_0"));
    TEST_ASSERT(my_SingleBlockDOFManager->fieldInBlock("Yaw","eblock-0_0"));
    TEST_ASSERT(my_SingleBlockDOFManager->fieldInBlock("Velocity","eblock-0_0"));

    std::vector<int> all_values = my_SingleBlockDOFManager->getBlockFieldNumbers("eblock-0_0");
    std::vector<int> compare_values(9);
    compare_values[0]=0;
    compare_values[1]=1;
    compare_values[2]=2;
    compare_values[3]=3;
    compare_values[4]=4;
    compare_values[5]=5;
    compare_values[6]=6;
    compare_values[7]=7;
    compare_values[8]=8;
    bool all_so_far=true;
    for (size_t i = 0; i < all_values.size(); ++i) {
      bool there=false;
      for (size_t j = 0; j < compare_values.size(); ++j) {
        if(all_values[i]==compare_values[j])
          there=true;
      }
      if(all_so_far && !there)
        all_so_far=false;
    }
    TEST_ASSERT(all_so_far);
  }
  

  TEUCHOS_UNIT_TEST( SingleBlockDOFManager_Test, AccuracyOfLocalElementArrays){
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",3);
    pl->set("Y Elements",3);
    
    //stk::ParallelMachine Comm = MPI_COMM_WORLD;

    //int myRank=stk::parallel_machine_rank(Comm);
    //int numProcs=stk::parallel_machine_size(Comm);

    //TEUCHOS_ASSERT(numProcs==2);

    panzer_stk::SquareQuadMeshFactory factory; 
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    RCP<panzer::SingleBlockDOFManager<LO,GO> > my_SingleBlockDOFManager = Teuchos::rcp(new panzer::SingleBlockDOFManager<LO,GO>());
    TEST_EQUALITY(my_SingleBlockDOFManager->getComm(),Teuchos::null);

    RCP<panzer::ConnManager<LO,GO> > conn = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<Intrepid::Basis<double,FieldContainer> > basis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer>);
     
    RCP< const panzer::FieldPattern> pattern = Teuchos::rcp(new panzer::IntrepidFieldPattern(basis));

    my_SingleBlockDOFManager->setConnManager(conn, MPI_COMM_WORLD);

    std::vector<std::string> names;
    names.push_back("uy");
    names.push_back("p");
    names.push_back("ux");
    my_SingleBlockDOFManager->addField(names[0], pattern);
    my_SingleBlockDOFManager->addField(names[1], pattern);
    my_SingleBlockDOFManager->addField(names[2], pattern);
    int numFields=3;

    TEST_NOTHROW(my_SingleBlockDOFManager->buildGlobalUnknowns());

    std::vector<int> myOwned;
    my_SingleBlockDOFManager->getOwnedIndices(myOwned);

    const std::vector<int> & uy_offsets = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",my_SingleBlockDOFManager->getFieldNum("uy"));
    const std::vector<int> & p_offsets = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",my_SingleBlockDOFManager->getFieldNum("p"));
    const std::vector<int> & ux_offsets = my_SingleBlockDOFManager->getGIDFieldOffsets("eblock-0_0",my_SingleBlockDOFManager->getFieldNum("ux"));

    TEST_EQUALITY(uy_offsets.size(),p_offsets.size());
    TEST_EQUALITY(uy_offsets.size(),ux_offsets.size());

    const std::vector<int> localElements = conn->getElementBlock("eblock-0_0");
    for(std::size_t i=0; i<localElements.size(); ++i){
      const int* elementConn = conn->getConnectivity(localElements[i]);
      //TESTING FOR CORRECTNESS OF DIVISIBILITY.
      //AND Relationship between connectivity ID and GIDs.
      int looper=0;
      std::vector<int> gids;
      my_SingleBlockDOFManager->getElementGIDs(i,gids);
      TEST_EQUALITY(gids.size(),(size_t)4*numFields);
      for(std::size_t j=0; j<gids.size(); ++j){
        TEST_ASSERT((gids[j]-looper)%numFields==0);
        TEST_EQUALITY((elementConn[j/numFields])*numFields+looper,gids[j]);
        if(looper+1==numFields)
          looper=0;
        else
          ++looper;
      }
    }
  }


} /*generic namespace*/ 
