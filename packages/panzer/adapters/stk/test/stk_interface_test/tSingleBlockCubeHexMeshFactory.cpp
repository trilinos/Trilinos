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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace panzer_stk {

void test1(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm);
void test2(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm);
void test4(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm);
void test27(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm);

void entityVecToGIDVec(const std::vector<stk::mesh::Entity*> & eVec,
                             std::vector<stk::mesh::EntityId> & gidVec)
{
   gidVec.resize(eVec.size());
   for(std::size_t i=0;i<eVec.size();i++)
      gidVec[i] = eVec[i]->identifier();

   std::sort(gidVec.begin(),gidVec.end());
}

/*
TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, test_four)
{
   MPI_Comm mpiCommWorld = MPI_COMM_WORLD;
   test4(out,success,mpiCommWorld);
}
*/

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, element_counts)
{
   int colors[] = 
      { 1,
        2,2,
        4,4,4,4,
        27, 27, 27, 27, 27, 27, 27, 27, 27,
        27, 27, 27, 27, 27, 27, 27, 27, 27,
        27, 27, 27, 27, 27, 27, 27, 27, 27 };
   int myrank=0, mycolor=0;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   mycolor = colors[myrank];

   MPI_Comm commUT; // comm under test
   MPI_Comm_split(MPI_COMM_WORLD,colors[myrank],0,&commUT);

   int utSize = 0;
   MPI_Comm_size(commUT, &utSize);

   if(utSize!=mycolor) {
      out << "Processor " << myrank << " quiting because there is nothing to test." << std::endl;
      return;
   }

   switch(mycolor) {
   case 1:
      test1(out,success,commUT);
      break;
   case 2:
      test2(out,success,commUT);
      break;
   case 4:
      test4(out,success,commUT);
      break;
   case 27:
      test27(out,success,commUT);
      break;
   };
}

void test1(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,1);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",1);
   pl->set("Y Procs",1);
   pl->set("Z Procs",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));
}

void test2(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,2);
   int rank; MPI_Comm_rank(comm, &rank); 

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",1);
   pl->set("Y Procs",2);
   pl->set("Z Procs",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   std::vector<stk::mesh::Entity*> myElements;
   std::vector<stk::mesh::EntityId> myGids;
   mesh->getMyElements(myElements);
   entityVecToGIDVec(myElements,myGids);

   if(rank==0) {
      TEST_EQUALITY(myGids.size(),20);
      for(std::size_t i=0;i<5;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+1);
         TEST_EQUALITY(myGids[4*i+1],8*i+2);
         TEST_EQUALITY(myGids[4*i+2],8*i+3);
         TEST_EQUALITY(myGids[4*i+3],8*i+4);
      }
   }
   else if(rank==1) {
      TEST_EQUALITY(myGids.size(),20);
      for(std::size_t i=0;i<5;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+5);
         TEST_EQUALITY(myGids[4*i+1],8*i+6);
         TEST_EQUALITY(myGids[4*i+2],8*i+7);
         TEST_EQUALITY(myGids[4*i+3],8*i+8);
      }
   }
   else TEST_ASSERT(false);
}

void test4(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,4);
   int rank; MPI_Comm_rank(comm, &rank); 

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",1);
   pl->set("Y Procs",2);
   pl->set("Z Procs",2);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   std::vector<stk::mesh::Entity*> myElements;
   std::vector<stk::mesh::EntityId> myGids;
   mesh->getMyElements(myElements);
   entityVecToGIDVec(myElements,myGids);

   if(rank==0) {
      TEST_EQUALITY(myGids.size(),12);
      for(std::size_t i=0;i<3;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+1);
         TEST_EQUALITY(myGids[4*i+1],8*i+2);
         TEST_EQUALITY(myGids[4*i+2],8*i+3);
         TEST_EQUALITY(myGids[4*i+3],8*i+4);
      }
   }
   else if(rank==1) {
      TEST_EQUALITY(myGids.size(),12);
      for(std::size_t i=0;i<3;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+5);
         TEST_EQUALITY(myGids[4*i+1],8*i+6);
         TEST_EQUALITY(myGids[4*i+2],8*i+7);
         TEST_EQUALITY(myGids[4*i+3],8*i+8);
      }
   }
   else if(rank==2) {
      TEST_EQUALITY(myGids.size(),8);
      for(std::size_t i=0;i<2;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+25);
         TEST_EQUALITY(myGids[4*i+1],8*i+26);
         TEST_EQUALITY(myGids[4*i+2],8*i+27);
         TEST_EQUALITY(myGids[4*i+3],8*i+28);
      }
   }
   else if(rank==3) {
      TEST_EQUALITY(myGids.size(),8);
      for(std::size_t i=0;i<2;i++) {
         TEST_EQUALITY(myGids[4*i]  ,8*i+29);
         TEST_EQUALITY(myGids[4*i+1],8*i+30);
         TEST_EQUALITY(myGids[4*i+2],8*i+31);
         TEST_EQUALITY(myGids[4*i+3],8*i+32);
      }
   }
   else TEST_ASSERT(false);
}

void test27(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,27);
   int rank; MPI_Comm_rank(comm, &rank); 

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",3);
   pl->set("Y Procs",3);
   pl->set("Z Procs",3);
   pl->set("X Elements",6);
   pl->set("Y Elements",9);
   pl->set("Z Elements",12);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),6*9*12);

   std::vector<stk::mesh::Entity*> myElements;
   std::vector<stk::mesh::EntityId> myGids;
   mesh->getMyElements(myElements);
   entityVecToGIDVec(myElements,myGids);
}

}
