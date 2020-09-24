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
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#include "Kokkos_DynRankView.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tExodusEdgeBlock, edge_count)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   
   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   std::size_t xelems=2;
   std::size_t yelems=4;
   std::size_t zelems=5;
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",(int)xelems);
   pl->set("Y Elements",(int)yelems);
   pl->set("Z Elements",(int)zelems);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable())
      mesh->writeToExodus("EdgeBlock1.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),xelems*yelems*zelems);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),xelems*yelems*(zelems+1)+xelems*zelems*(yelems+1)+yelems*zelems*(xelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),xelems*(yelems+1)*(zelems+1)+yelems*(xelems+1)*(zelems+1)+zelems*(xelems+1)*(yelems+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(yelems+1)*(xelems+1)*(zelems+1));
   
   std::vector<stk::mesh::Entity> all_edges;
   mesh->getAllEdges(panzer_stk::STK_Interface::edgeBlockString, all_edges);
   TEST_EQUALITY(all_edges.size(),xelems*(yelems+1)*(zelems+1)
                             +yelems*(xelems+1)*(zelems+1)
                             +zelems*(xelems+1)*(yelems+1));

   std::vector<stk::mesh::Entity> my_edges;
   if (numprocs==1) {
     // all edges belong to rank0, so getMyEdges() is equivalent to getAllEdges()
     mesh->getMyEdges(panzer_stk::STK_Interface::edgeBlockString, my_edges);
     TEST_EQUALITY(my_edges.size(),all_edges.size());
     TEST_EQUALITY(my_edges.size(),xelems*(yelems+1)*(zelems+1)
                                  +yelems*(xelems+1)*(zelems+1)
                                  +zelems*(xelems+1)*(yelems+1));
   }
   else if(numprocs==2 && rank==0) {
     // rank0 owns all edges in it's half of the mesh including the 
     // edges on the face shared with rank1
     std::size_t my_xelems=xelems-1;
     std::size_t my_yelems=yelems;
     std::size_t my_zelems=zelems;
     mesh->getMyEdges(panzer_stk::STK_Interface::edgeBlockString, my_edges);
     TEST_EQUALITY(my_edges.size(),my_xelems*(my_yelems+1)*(my_zelems+1)
                                  +my_yelems*(my_xelems+1)*(my_zelems+1)
                                  +my_zelems*(my_xelems+1)*(my_yelems+1));
   }
   else if(numprocs==2 && rank==1) {
     // rank1 doesn't own the edges on the face shared with rank0
     std::size_t my_xelems=xelems-1;
     std::size_t my_yelems=yelems;
     std::size_t my_zelems=zelems;
     mesh->getMyEdges(panzer_stk::STK_Interface::edgeBlockString, my_edges);
     TEST_EQUALITY(my_edges.size(),my_xelems*(my_yelems+1)*(my_zelems+1)
                                  +my_yelems*(my_xelems+1)*(my_zelems+1)
                                  +my_zelems*(my_xelems+1)*(my_yelems+1)
                                  -my_yelems*(my_zelems+1)               // remove edges owned by rank0
                                  -my_zelems*(my_yelems+1));             // remove edges owned by rank0
   }
   else {
     // fail!
     TEST_ASSERT(false && "This test must run with either 1 or 2 ranks.");
   }
}

TEUCHOS_UNIT_TEST(tExodusEdgeBlock, add_edge_field)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   mesh->addEdgeField("edge_field_1", "eblock-0_0_0");
   mesh->addEdgeField("edge_field_2", "eblock-0_0_0");
   
   stk::mesh::Field<double> * edge_field_1 = mesh->getEdgeField("edge_field_1", "eblock-0_0_0");
   stk::mesh::Field<double> * edge_field_2 = mesh->getEdgeField("edge_field_2", "eblock-0_0_0");

   std::vector<stk::mesh::Entity> edges;
   mesh->getAllEdges(panzer_stk::STK_Interface::edgeBlockString, edges);
   for(auto edge : edges) {
     double* data = stk::mesh::field_data(*edge_field_1, edge);
     // set the edge's field value to edge's entity ID
     *data = mesh->getBulkData()->identifier(edge);
   }
   for(auto edge : edges) {
     double* data = stk::mesh::field_data(*edge_field_2, edge);
     // set the edge's field value to edge's entity ID * 2
     *data = 2*mesh->getBulkData()->identifier(edge);
   }

   if(mesh->isWritable())
      mesh->writeToExodus("EdgeBlock2.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   mesh->getAllEdges(panzer_stk::STK_Interface::edgeBlockString, edges);
   TEST_EQUALITY(edges.size(),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
}

TEUCHOS_UNIT_TEST(tExodusEdgeBlock, set_edge_field_data)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   // need to initialize the list of local edge IDs
   // which is used by edgeLocalId() below
   mesh->buildLocalEdgeIDs();

   mesh->addEdgeField("edge_field_3", "eblock-0_0_0");
   mesh->addEdgeField("edge_field_4", "eblock-0_0_0");
   
   std::vector<stk::mesh::Entity> edges;
   mesh->getMyEdges(panzer_stk::STK_Interface::edgeBlockString, edges);

   Kokkos::DynRankView<double,PHX::Device> edgeValues;
   edgeValues = Kokkos::createDynRankView(edgeValues,"edgeValues",edges.size());

   std::vector<std::size_t> edgeIds;
   for(auto edge : edges) {
     edgeIds.push_back(mesh->edgeLocalId(edge));
   }
   sort(edgeIds.begin(),edgeIds.end());

   for(std::size_t i=0;i<edgeIds.size();i++) {
     edgeValues(i) = 3*mesh->edgeGlobalId(edgeIds[i]);
   }
   mesh->setEdgeFieldData("edge_field_3",
                          "eblock-0_0_0",
                          edgeIds,
                          edgeValues);

   for(std::size_t i=0;i<edgeIds.size();i++) {
     edgeValues(i) = 4*mesh->edgeGlobalId(edgeIds[i]);
   }
   mesh->setEdgeFieldData("edge_field_4",
                          "eblock-0_0_0",
                          edgeIds,
                          edgeValues);

   if(mesh->isWritable())
      mesh->writeToExodus("EdgeBlock3.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));

   mesh->getAllEdges(panzer_stk::STK_Interface::edgeBlockString, edges);
   TEST_EQUALITY(edges.size(),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
}

}
