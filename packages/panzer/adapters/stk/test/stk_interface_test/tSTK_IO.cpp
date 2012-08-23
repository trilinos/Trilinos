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
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include "Intrepid_FieldContainer.hpp"

#ifdef HAVE_IOSS

using Teuchos::RCP;
using Teuchos::rcp;

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer_stk {

RCP<STK_Interface> buildMesh(int xElements,int yElements);
RCP<STK_Interface> buildMesh_cells(int xElements,int yElements);
void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds); 

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double));
void assignBlock(FieldContainer & block,FieldContainer & vertices, double val);

double xval(double x,double y) { return x; }
double yval(double x,double y) { return y; }
double block2(double x,double y) { return (x-0.5)*(x-0.5)+y; }

// triangle tests
TEUCHOS_UNIT_TEST(tSTK_IO, fields)
{
   RCP<STK_Interface> mesh = buildMesh(8,8);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds; 
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock1.resize(localIds["eblock-1_0"]->size(),4);
   out << "assigning" << std::endl;

   assignBlock(ublock0,vert0,xval);
   assignBlock(tblock0,vert0,yval);
   assignBlock(tblock1,vert1,block2);

   mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
   mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
   mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);

   out << "write to exodus" << std::endl;

   mesh->writeToExodus("output.exo");
}

TEUCHOS_UNIT_TEST(tSTK_IO, cell_fields)
{
   RCP<STK_Interface> mesh = buildMesh_cells(8,8);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds; 
   buildLocalIds(*mesh,localIds);

   out << "write to exodus" << std::endl;

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock1.resize(localIds["eblock-1_0"]->size(),4);
   out << "assigning" << std::endl;

   assignBlock(ublock0,vert0,xval);
   assignBlock(tblock0,vert0,yval);
   assignBlock(tblock1,vert1,block2);

   mesh->setCellFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
   mesh->setCellFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
   mesh->setCellFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);

   mesh->writeToExodus("output-cells.exo");
}

TEUCHOS_UNIT_TEST(tSTK_IO, exodus_factory_transient_fields)
{
   STK_ExodusReaderFactory factory("meshes/basic.gen");
   RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
   mesh->addSolutionField("u","block_1");
   mesh->addSolutionField("T","block_1");
   mesh->addSolutionField("T","block_2");
   factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds; 
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["block_1"],vert0);
   mesh->getElementVertices(*localIds["block_2"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0.resize(localIds["block_1"]->size(),4);
   tblock0.resize(localIds["block_1"]->size(),4);
   tblock1.resize(localIds["block_2"]->size(),4);

   mesh->setupTransientExodusFile("transient_exo.exo");

   out << "assigning 4.5" << std::endl;
   {
      assignBlock(ublock0,vert0,6.0);
      assignBlock(tblock0,vert0,7.0);
      assignBlock(tblock1,vert1,8.0);

      mesh->setSolutionFieldData("u","block_1",*localIds["block_1"],ublock0);
      mesh->setSolutionFieldData("T","block_1",*localIds["block_1"],tblock0);
      mesh->setSolutionFieldData("T","block_2",*localIds["block_2"],tblock1);
   }

   out << "write to exodus: 4.5" << std::endl;
   mesh->writeToExodus(4.5);
}

TEUCHOS_UNIT_TEST(tSTK_IO, transient_fields)
{
   RCP<STK_Interface> mesh = buildMesh(20,20);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds; 
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock1.resize(localIds["eblock-1_0"]->size(),4);

   mesh->setupTransientExodusFile("transient.exo");

   out << "assigning 3.0" << std::endl;
   {
      assignBlock(ublock0,vert0,1.0);
      assignBlock(tblock0,vert0,2.0);
      assignBlock(tblock1,vert1,3.0);

      mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
      mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
      mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   }

   out << "write to exodus: 3.0" << std::endl;
   mesh->writeToExodus(3.0);
   TEST_EQUALITY(mesh->getCurrentStateTime(),3.0); 

   out << "assigning 4.5" << std::endl;
   {
      assignBlock(ublock0,vert0,6.0);
      assignBlock(tblock0,vert0,7.0);
      assignBlock(tblock1,vert1,8.0);

      mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
      mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
      mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);
   }

   out << "write to exodus: 4.5" << std::endl;
   mesh->writeToExodus(4.5);
   TEST_EQUALITY(mesh->getCurrentStateTime(),4.5); 

   STK_ExodusReaderFactory factory("transient.exo",2);
   RCP<STK_Interface> mesh_read = factory.buildMesh(MPI_COMM_WORLD);
   TEST_EQUALITY(mesh_read->getInitialStateTime(),4.5);
   TEST_EQUALITY(mesh_read->getCurrentStateTime(),0.0); // writeToExodus has not yet been called
}

RCP<STK_Interface> buildMesh(int xElements,int yElements)
{
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",xElements);  // in each block
    pl->set("Y Elements",yElements);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addSolutionField("u","eblock-0_0");
    mesh->addSolutionField("T","eblock-0_0");
    mesh->addSolutionField("T","eblock-1_0");
 
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    return mesh;
}

RCP<STK_Interface> buildMesh_cells(int xElements,int yElements)
{
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",xElements);  // in each block
    pl->set("Y Elements",yElements);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addCellField("u","eblock-0_0");
    mesh->addCellField("T","eblock-0_0");
    mesh->addCellField("T","eblock-1_0");

    mesh->addSolutionField("T","eblock-1_0");
 
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    return mesh;
}

void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds)
{
   // defines ordering of blocks
   std::vector<std::string> blockIds;
   mesh.getElementBlockNames(blockIds);

   std::vector<std::string>::const_iterator idItr;
   for(idItr=blockIds.begin();idItr!=blockIds.end();++idItr) {
      std::string blockId = *idItr;

      localIds[blockId] = Teuchos::rcp(new std::vector<std::size_t>);
      std::vector<std::size_t> & localBlockIds = *localIds[blockId];

      // grab elements on this block
      std::vector<stk::mesh::Entity*> blockElmts;
      mesh.getMyElements(blockId,blockElmts);

      std::vector<stk::mesh::Entity*>::const_iterator itr;
      for(itr=blockElmts.begin();itr!=blockElmts.end();++itr)
         localBlockIds.push_back(mesh.elementLocalId(*itr));

      std::sort(localBlockIds.begin(),localBlockIds.end());
   }
}

void assignBlock(FieldContainer & block,FieldContainer & vertices, double val)
{
   TEUCHOS_ASSERT(block.dimension(0)==vertices.dimension(0));
   TEUCHOS_ASSERT(block.dimension(1)==vertices.dimension(1));

   std::size_t cellCnt = block.dimension(0); 
   std::size_t nodeCnt = block.dimension(1); 

   for(std::size_t cell=0;cell<cellCnt;cell++) {
      for(std::size_t node=0;node<nodeCnt;node++) {
         block(cell,node) = val;
      }
   }
}

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double))
{
   TEUCHOS_ASSERT(block.dimension(0)==vertices.dimension(0));
   TEUCHOS_ASSERT(block.dimension(1)==vertices.dimension(1));

   std::size_t cellCnt = block.dimension(0); 
   std::size_t nodeCnt = block.dimension(1); 

   for(std::size_t cell=0;cell<cellCnt;cell++) {
      for(std::size_t node=0;node<nodeCnt;node++) {
         double x = vertices(cell,node,0);
         double y = vertices(cell,node,1);
         block(cell,node) = func(x,y);
      }
   }
}

}

#endif
