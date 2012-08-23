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

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

namespace panzer_stk {

typedef shards::Quadrilateral<4> QuadTopo;

Teuchos::RCP<STK_Interface> build2DMesh()
{
   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   Teuchos::RCP<STK_Interface> meshPtr = Teuchos::rcp(new STK_Interface(2));
   STK_Interface & mesh = *meshPtr;

   mesh.addElementBlock("quad_elements",ctd);
   mesh.addSideset("Left",side_ctd);
   mesh.addSideset("Right",side_ctd);
   mesh.addSideset("Top",side_ctd);
   mesh.addSideset("Bottom",side_ctd);

   mesh.initialize(MPI_COMM_WORLD);
   mesh.beginModification();
      std::vector<double> coord(2);
      stk::mesh::Part * block = mesh.getElementBlockPart("quad_elements");
      { 
         // Add four coordinates
         //
         //    4 ---- 3
         //    |      |
         //    |      |
         //    1 ---- 2 
         //
      
         coord[0] = 0.0; coord[1] = 0.0;
         mesh.addNode(1,coord); 

         coord[0] = 1.0; coord[1] = 0.0;
         mesh.addNode(2,coord); 

         coord[0] = 1.0; coord[1] = 1.0;
         mesh.addNode(3,coord); 

         coord[0] = 0.0; coord[1] = 1.0;
         mesh.addNode(4,coord); 

         // add an element
         std::vector<stk::mesh::EntityId> nodes;
         for(std::size_t i=1;i<5;i++)
            nodes.push_back(i);

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(1,nodes);
         mesh.addElement(ed,block);
      }

      { 
         // Add four coordinates
         //
         //    3 ---- 6
         //    |      |
         //    |      |
         //    2 ---- 5 
         //
      
         coord[0] = 2.0; coord[1] = 0.5;
         mesh.addNode(5,coord); 

         coord[0] = 2.1; coord[1] = 1.5;
         mesh.addNode(6,coord); 

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 5;
         nodes[1] = 6;
         nodes[2] = 3;
         nodes[3] = 2;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(2,nodes);
         mesh.addElement(ed,block);
      }

      { 
         // Add four coordinates
         //
         //    8 ---- 7
         //    |      |
         //    |      |
         //    4 ---- 3 
         //
      
         coord[0] = 1.0; coord[1] = 2.5;
         mesh.addNode(7,coord); 

         coord[0] = 0.1; coord[1] = 2.0;
         mesh.addNode(8,coord); 

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 4;
         nodes[1] = 3;
         nodes[2] = 7;
         nodes[3] = 8;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(3,nodes);
         mesh.addElement(ed,block);
      }

   mesh.endModification();
   mesh.buildLocalElementIDs();

   return meshPtr;
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKInterface, interface_test)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif
   RCP<Epetra_Comm> comm = rcpFromRef(Comm);

   STK_Interface mesh(2);

   TEST_EQUALITY(mesh.getDimension(),2);

   mesh.addElementBlock("0",ctd);
   mesh.addSideset("Inflow",side_ctd);
   mesh.addSideset("Outflow",side_ctd);
   mesh.addSideset("Top",side_ctd);
   mesh.addSideset("Bottom",side_ctd);

   TEST_EQUALITY(mesh.getDimension(),2);
   TEST_EQUALITY(mesh.getNumSidesets(),4);
   TEST_EQUALITY(mesh.getNumElementBlocks(),1);

   TEST_ASSERT(not mesh.isModifiable());
   mesh.initialize(MPI_COMM_WORLD);

   TEST_ASSERT(not mesh.isModifiable());
   mesh.beginModification();
      TEST_ASSERT(mesh.isModifiable());

      std::vector<double> coord(2);
      stk::mesh::Part * block = mesh.getElementBlockPart("0");
      { 
         // Add four coordinates
         //
         //    4 ---- 3
         //    |      |
         //    |      |
         //    1 ---- 2 
         //
      
         coord[0] = 0.0; coord[1] = 0.0;
         mesh.addNode(1,coord); 

         coord[0] = 1.0; coord[1] = 0.0;
         mesh.addNode(2,coord); 

         coord[0] = 1.0; coord[1] = 1.0;
         mesh.addNode(3,coord); 

         coord[0] = 0.0; coord[1] = 1.0;
         mesh.addNode(4,coord); 

         // add an element
         std::vector<stk::mesh::EntityId> nodes;
         for(std::size_t i=1;i<5;i++)
            nodes.push_back(i);

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(1,nodes);
         mesh.addElement(ed,block);
      }

      { 
         // Add four coordinates
         //
         //    3 ---- 6
         //    |      |
         //    |      |
         //    2 ---- 5 
         //
      
         coord[0] = 2.0; coord[1] = 0.5;
         mesh.addNode(5,coord); 

         coord[0] = 2.1; coord[1] = 1.5;
         mesh.addNode(6,coord); 

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 5;
         nodes[1] = 6;
         nodes[2] = 3;
         nodes[3] = 2;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(2,nodes);
         mesh.addElement(ed,block);
      }

   mesh.endModification();
   TEST_ASSERT(not mesh.isModifiable());

   stk::mesh::EntityRank nodeRank = mesh.getNodeRank();
   stk::mesh::EntityRank sideRank = mesh.getSideRank();
   stk::mesh::EntityRank elmtRank = mesh.getElementRank();

   TEST_EQUALITY(mesh.getEntityCounts(nodeRank),6);
   TEST_EQUALITY(mesh.getEntityCounts(sideRank),0);
   TEST_EQUALITY(mesh.getEntityCounts(elmtRank),2);

   #ifdef HAVE_IOSS
      TEST_ASSERT(mesh.isWritable());
      TEST_NOTHROW(mesh.writeToExodus("simplemesh.exo"));
   #else
      TEST_ASSERT(not mesh.isWritable());
      TEST_THROW(mesh.writeToExodus("simplemesh.exo"),std::logic_error);
   #endif 

   const double * coords = 0;

   coords = mesh.getNodeCoordinates(5);
   TEST_FLOATING_EQUALITY(coords[0],2.0,1e-14);
   TEST_FLOATING_EQUALITY(coords[1],0.5,1e-14);

   coords = mesh.getNodeCoordinates(2);
   TEST_FLOATING_EQUALITY(coords[0],1.0,1e-14);
   TEST_EQUALITY(coords[1],0.0);

   coords = mesh.getNodeCoordinates(1);
   TEST_EQUALITY(coords[0],0.0);
   TEST_EQUALITY(coords[1],0.0);

   TEST_EQUALITY(mesh.getMaxEntityId(nodeRank),6);
   TEST_EQUALITY(mesh.getMaxEntityId(elmtRank),2);
}

class CompareID {
public:
   CompareID(stk::mesh::EntityId id) : id_(id) {}

   bool operator()(stk::mesh::Entity * e)
   { return e->identifier()==id_; }
   stk::mesh::EntityId id_;
};

TEUCHOS_UNIT_TEST(tSTKInterface, node_sharing_test)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
 
   RCP<STK_Interface> mesh = build2DMesh();
 
   if(mesh->isWritable())
      mesh->writeToExodus("simplemesh.exo");
 
   {
      std::vector<stk::mesh::Entity*> elements;
      mesh->getElementsSharingNode(2,elements);
 
      TEST_EQUALITY(elements.size(),2); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(1))!=elements.end()); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(2))!=elements.end()); 
   }
 
   {
      std::vector<stk::mesh::Entity*> elements;
      mesh->getElementsSharingNode(4,elements);
 
      TEST_EQUALITY(elements.size(),2); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(1))!=elements.end()); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(3))!=elements.end()); 
   }
 
   {
      std::vector<stk::mesh::Entity*> elements;
      mesh->getElementsSharingNode(3,elements);
 
      TEST_EQUALITY(elements.size(),3); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(1))!=elements.end()); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(2))!=elements.end()); 
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(3))!=elements.end()); 
   }

   {
     std::vector<stk::mesh::EntityId> nodes;
     nodes.push_back(3);
     nodes.push_back(4);

     std::vector<stk::mesh::Entity*> elements;
     mesh->getElementsSharingNodes(nodes,elements);

     TEST_EQUALITY(elements.size(),2); 
     TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(1))!=elements.end()); 
     TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(3))!=elements.end()); 
   }
    
   {
      std::vector<stk::mesh::EntityId> nodes;
      nodes.push_back(1);
      nodes.push_back(5);

      std::vector<stk::mesh::Entity*> elements;
      mesh->getElementsSharingNodes(nodes,elements);

      TEST_EQUALITY(elements.size(),0); 
   }
}

TEUCHOS_UNIT_TEST(tSTKInterface, subcellIndices)
{
   using Teuchos::RCP;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   stk::mesh::EntityRank nodeRank = mesh->getNodeRank();
   stk::mesh::EntityRank sideRank = mesh->getSideRank();

   mesh->buildSubcells();

   std::vector<stk::mesh::EntityId> subcells;

   TEST_THROW(mesh->getSubcellIndices(nodeRank,9,subcells),std::logic_error);

   // get nodes
   mesh->getSubcellIndices(nodeRank,3,subcells);
   TEST_EQUALITY(subcells.size(),4);
   TEST_EQUALITY(subcells[0],4);
   TEST_EQUALITY(subcells[1],3);
   TEST_EQUALITY(subcells[2],7);
   TEST_EQUALITY(subcells[3],8);

   // get edges
   mesh->getSubcellIndices(sideRank,3,subcells);
   TEST_EQUALITY(subcells.size(),4);
   //TEST_EQUALITY(subcells[0],20);
   //TEST_EQUALITY(subcells[1],23);
   //TEST_EQUALITY(subcells[2],56);
   //TEST_EQUALITY(subcells[3],32);
}

TEUCHOS_UNIT_TEST(tSTKInterface, local_ids)
{
   using Teuchos::RCP;

   std::vector<stk::mesh::Entity*> elements;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   stk::mesh::EntityRank nodeRank = mesh->getNodeRank();

   mesh->getMyElements(elements);

   // loop over all elements of mesh
   for(std::size_t elmI=0;elmI<elements.size();++elmI) {
      stk::mesh::Entity * elem = elements[elmI];
      std::size_t localId = mesh->elementLocalId(elem);

      stk::mesh::PairIterRelation relations = elem->relations(nodeRank);
      stk::mesh::PairIterRelation::iterator itr;
      stk::mesh::Entity * node = 0;
      for(itr=relations.begin();itr!=relations.end();++itr) {
         if(itr->identifier()==0) {
            node = itr->entity();
            break;
         }
      }

      TEST_ASSERT(node!=0);

      // based on first node check local id of element
      switch(node->identifier()) {
      case 1:
         TEST_EQUALITY(localId,0);
         break;
      case 5:
         TEST_EQUALITY(localId,1);
         break;
      case 4:
         TEST_EQUALITY(localId,2);
         break;
      default:
         TEST_ASSERT(false);
      }
   }
}

TEUCHOS_UNIT_TEST(tSTKInterface, edgeAddTest)
{
   using Teuchos::RCP;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   stk::mesh::EntityRank nodeRank = mesh->getNodeRank();
   stk::mesh::EntityRank sideRank = mesh->getSideRank();
   stk::mesh::EntityRank elmtRank = mesh->getElementRank();

   mesh->buildSubcells();

   if(mesh->isWritable())
      mesh->writeToExodus("simplemesh_wedges.exo");

   TEST_EQUALITY(mesh->getEntityCounts(nodeRank),8);
   TEST_EQUALITY(mesh->getEntityCounts(sideRank),10);
   TEST_EQUALITY(mesh->getEntityCounts(elmtRank),3);
}

}
