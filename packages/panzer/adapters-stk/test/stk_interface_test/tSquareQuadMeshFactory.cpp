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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"
#include "Ioss_EdgeBlock.h"

namespace panzer_stk {

inline bool XOR(bool A,bool B)
{ return ! ( (A && B) || ( !A && !B)); }

class LocalIdCompare {
public:
   LocalIdCompare(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   bool operator()(stk::mesh::Entity a,stk::mesh::Entity b) const
   { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b); }

private:
   Teuchos::RCP<const STK_Interface> mesh_;
};

/*
static const double * getNode(const Teuchos::RCP<const STK_Interface> & mesh, stk::mesh::Entity element,int id)
{
   std::vector<stk::mesh::EntityId> nodeIds;
   getNodeIds(mesh->getNodeRank(),element,nodeIds);

   return mesh->getNodeCoordinates(nodeIds[id]);
}
*/

void edge_block_test_helper(Teuchos::FancyOStream &out,
                            bool &success,
                            Teuchos::RCP<Teuchos::ParameterList> pl,
                            std::string exodus_filename,
                            uint32_t expected_edge_block_count)
{
   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus(exodus_filename.c_str());

   {
   Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus",
                                                     exodus_filename.c_str(),
                                                     Ioss::READ_MODEL);
   TEST_ASSERT(db_io);

   Ioss::Region region(db_io);
   TEST_ASSERT(db_io->ok() == true);

   auto all_edge_blocks = region.get_edge_blocks();
   TEST_ASSERT(all_edge_blocks.size() == expected_edge_block_count);

   if (expected_edge_block_count == 1) {
      std::vector<stk::mesh::Entity> edges;
      mesh->getMyEdges(edges);
      TEST_ASSERT(all_edge_blocks[0]->entity_count() == (int64_t)edges.size());
   }
   }
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, periodic_input)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",3);
   Teuchos::ParameterList & pbcs = pl->sublist("Periodic BCs");
   pbcs.set<int>("Count",1);
   pbcs.set("Periodic Condition 1","x-coord left;right");

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),1);
}

// This test was modified to its current lame state when the
// construction of the local element IDs was automated in the
// STK_Interface. (Independent of order of addition in the mesh

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, local_ids)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",3);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);

   std::string strBlock0="eblock-0_0", strBlock1="eblock-1_0";
   std::vector<stk::mesh::Entity> block0, block1;

   mesh->getMyElements(strBlock0,block0);
   mesh->getMyElements(strBlock1,block1);

   std::sort(block0.begin(),block0.end(),LocalIdCompare(mesh));
   std::sort(block1.begin(),block1.end(),LocalIdCompare(mesh));

   if(numprocs==1) {
      TEST_EQUALITY(block0.size(),6);
      TEST_EQUALITY(block1.size(),6);

/*
      {
         const double * coords = getNode(mesh,block0[0],0);
         TEST_FLOATING_EQUALITY(coords[0],0.0,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],0.0,1e-10);
      }

      {
         const double * coords = getNode(mesh,block0[2],1);
         TEST_FLOATING_EQUALITY(coords[0],0.25,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],2.0/3.0,1e-10);
      }
*/
   }
   else if(numprocs==2 && rank==0) {
      TEST_EQUALITY(block0.size(),3);
      TEST_EQUALITY(block1.size(),3);

/*
      {
         const double * coords = getNode(mesh,block0[0],0);
         TEST_FLOATING_EQUALITY(coords[0],0.0,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],0.0,1e-10);
      }

      {
         const double * coords = getNode(mesh,block0[2],1);
         TEST_FLOATING_EQUALITY(coords[0],0.25,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],2.0/3.0,1e-10);
      }
*/
   }
   else if(numprocs==2 && rank==1) {
      TEST_EQUALITY(block0.size(),3);
      TEST_EQUALITY(block1.size(),3);

/*
      {
         const double * coords = getNode(mesh,block0[0],0);
         TEST_FLOATING_EQUALITY(coords[0],0.25,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],0.0,1e-10);
      }

      {
         const double * coords = getNode(mesh,block0[2],1);
         TEST_FLOATING_EQUALITY(coords[0],0.5,1e-10);
         TEST_FLOATING_EQUALITY(coords[1],2.0/3.0,1e-10);
      }
*/
   }
   else {
      // fail!
      TEST_ASSERT(false);
   }
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   SquareQuadMeshFactory factory;
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   if(mesh->isWritable())
      mesh->writeToExodus("SquareQuad.exo");

   // minimal requirements
   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),25);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),60);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),36);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   if(numprocs==1) {
      out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;
      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),3,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),8,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[1]);
         TEST_EQUALITY(elmt0[3],elmt1[0]);

         mesh->getSubcellIndices(mesh->getSideRank(),3,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),8,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),7,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),13,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),14,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),19,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[1]);
         TEST_EQUALITY(elmt0[3],elmt1[0]);

         mesh->getSubcellIndices(mesh->getSideRank(),14,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),19,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),17,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),18,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[1],elmt1[0]);
         TEST_EQUALITY(elmt0[2],elmt1[3]);

         mesh->getSubcellIndices(mesh->getSideRank(),17,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),18,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[1],elmt1[3]);
      }
   }
   else if(numprocs==2 && rank==0) {
      out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;
      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),3,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),8,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[1]);
         TEST_EQUALITY(elmt0[3],elmt1[0]);

         mesh->getSubcellIndices(mesh->getSideRank(),3,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),8,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),17,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),18,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[1],elmt1[0]);
         TEST_EQUALITY(elmt0[2],elmt1[3]);

         mesh->getSubcellIndices(mesh->getSideRank(),17,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),18,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[1],elmt1[3]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),7,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),13,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }
   }
   else if(numprocs==2 && rank==1) {
      out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;
      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),14,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),19,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[1]);
         TEST_EQUALITY(elmt0[3],elmt1[0]);

         mesh->getSubcellIndices(mesh->getSideRank(),14,elmt0);
         mesh->getSubcellIndices(mesh->getSideRank(),19,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }

      {
         std::vector<stk::mesh::EntityId> elmt0, elmt1;

         mesh->getSubcellIndices(mesh->getNodeRank(),3,elmt0);
         mesh->getSubcellIndices(mesh->getNodeRank(),9,elmt1);

         TEST_EQUALITY(elmt0.size(),4);
         TEST_EQUALITY(elmt0.size(),elmt1.size());
         TEST_EQUALITY(elmt0[2],elmt1[0]);
      }
   }
   else {
      // fail!
      TEST_ASSERT(false);
   }

   TEST_EQUALITY(mesh->getMaxEntityId(mesh->getNodeRank()),36);
   TEST_EQUALITY(mesh->getMaxEntityId(mesh->getElementRank()),25);
}

// triangle tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, multi_xblock)
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",3);
   pl->set("X Elements",6);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   if(mesh->isWritable())
      mesh->writeToExodus("SquareQuad_Blocked.exo");

   TEST_EQUALITY(mesh->getNumElementBlocks(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(12+1)*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),12*(12+1)+12*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*12);

   std::vector<stk::mesh::Entity> myElements;
   mesh->getMyElements(myElements);

   TEST_EQUALITY(myElements.size(), (std::size_t) 12*12/numprocs);
}

// triangle tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, side_elmt_access)
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   if(numprocs>2)
      TEUCHOS_ASSERT(false);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",3);
   pl->set("X Elements",6);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   stk::mesh::BulkData & bulkData = *mesh->getBulkData();

   TEST_EQUALITY(mesh->getNumElementBlocks(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(12+1)*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),12*(12+1)+12*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*12);

   std::vector<std::string> blockNames;
   mesh->getElementBlockNames(blockNames);
   TEST_EQUALITY(blockNames.size(),mesh->getNumElementBlocks());
   TEST_EQUALITY(blockNames[0],"eblock-0_0");

   {
      std::vector<stk::mesh::Entity> myElements;
      mesh->getMyElements(blockNames[0],myElements);

      TEST_EQUALITY((int) myElements.size(),24/numprocs);

      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[0]),1+rank*3);
      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[1]),2+rank*3);
      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[2]),3+rank*3);

      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[24/numprocs-3]),40-(numprocs-1)*(1-rank)*3);
      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[24/numprocs-2]),41-(numprocs-1)*(1-rank)*3);
      TEST_EQUALITY((int) mesh->elementGlobalId(myElements[24/numprocs-1]),42-(numprocs-1)*(1-rank)*3);
   }

   {
      {
         std::vector<stk::mesh::Entity> mySides;

         mySides.clear();
         mesh->getMySides("left",mySides);
         TEST_EQUALITY((int) mySides.size(),(rank==0 ? 12 : 0 ));

         mySides.clear();
         mesh->getMySides("right",mySides);
         TEST_EQUALITY((int) mySides.size(),( XOR(rank==0,numprocs==2) ? 12 : 0 ));

         mySides.clear();
         mesh->getMySides("top",mySides);
         TEST_EQUALITY((int) mySides.size(),(numprocs==1 ? 12 : 6));

         mySides.clear();
         mesh->getMySides("bottom",mySides);
         TEST_EQUALITY((int) mySides.size(),(numprocs==1 ? 12 : 6));
      }

      std::vector<std::string> sidesets(4);
      sidesets[0] = "left";
      sidesets[1] = "right";
      sidesets[2] = "top";
      sidesets[3] = "bottom";
      std::vector<std::string>::const_iterator sItr;
      for(sItr=sidesets.begin();sItr!=sidesets.end();++sItr) {
         std::vector<stk::mesh::Entity> mySides;
         mesh->getMySides(*sItr,mySides);

         std::vector<stk::mesh::Entity>::iterator itr;
         for(itr=mySides.begin();itr!=mySides.end();++itr) {
            stk::mesh::Entity side = *itr;

            TEST_EQUALITY(bulkData.entity_rank(side),mesh->getSideRank());
            TEST_EQUALITY(bulkData.num_elements(side),1);
            TEST_EQUALITY(bulkData.num_nodes(side),2);
         }

      }
   }

   {
      {
         std::vector<stk::mesh::Entity> mySides;

         mySides.clear();
         mesh->getMySides("left","eblock-0_0",mySides);
         TEST_EQUALITY((int) mySides.size(),(rank==0 ? 4 : 0 ));

         mySides.clear();
         mesh->getMySides("right","eblock-1_0",mySides);
         TEST_EQUALITY((int) mySides.size(),( XOR(rank==0,numprocs==2) ? 4 : 0 ));

         mySides.clear();
         mesh->getMySides("top","eblock-0_2",mySides);
         TEST_EQUALITY((int) mySides.size(),6/numprocs);

         mySides.clear();
         mesh->getMySides("bottom","eblock-0_0",mySides);
         TEST_EQUALITY((int) mySides.size(),6/numprocs);
      }

      std::vector<std::string> sidesets(4);
      std::vector<std::string> eblocks(4);
      sidesets[0] = "left";   eblocks[0] = "eblock-0_0";
      sidesets[1] = "right";  eblocks[1] = "eblock-1_0";
      sidesets[2] = "top";    eblocks[2] = "eblock-0_2";
      sidesets[3] = "bottom"; eblocks[3] = "eblock-0_0";
      for(std::size_t i=0;i<sidesets.size();++i) {
         std::vector<stk::mesh::Entity> mySides;
         mesh->getMySides(sidesets[i],eblocks[i],mySides);

         std::vector<stk::mesh::Entity>::iterator itr;
         for(itr=mySides.begin();itr!=mySides.end();++itr) {
            stk::mesh::Entity side = *itr;

            TEST_EQUALITY(bulkData.entity_rank(side),mesh->getSideRank());
            TEST_EQUALITY(bulkData.num_elements(side),1);
            TEST_EQUALITY(bulkData.num_nodes(side),2);
         }

      }
   }
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, check_ss)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",1);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
   out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   stk::mesh::Selector ownedAndGhosted = mesh->getMetaData()->locally_owned_part()
                                       | mesh->getMetaData()->globally_shared_part();
   stk::mesh::Selector topPart = *mesh->getMetaData()->get_part("top","Big error!")    & ownedAndGhosted;
   stk::mesh::Selector btmPart = *mesh->getMetaData()->get_part("bottom","Big error!") & ownedAndGhosted;
   stk::mesh::Selector lftPart = *mesh->getMetaData()->get_part("left","Big error!")   & ownedAndGhosted;
   stk::mesh::Selector rhtPart = *mesh->getMetaData()->get_part("right","Big error!")  & ownedAndGhosted;

   const std::vector<stk::mesh::Bucket*> & nodeBuckets = mesh->getBulkData()->buckets(mesh->getNodeRank());

   unsigned lftCnt = stk::mesh::count_selected_entities(lftPart,nodeBuckets);
   unsigned rhtCnt = stk::mesh::count_selected_entities(rhtPart,nodeBuckets);
   unsigned topCnt = stk::mesh::count_selected_entities(topPart,nodeBuckets);
   unsigned btmCnt = stk::mesh::count_selected_entities(btmPart,nodeBuckets);

   if(numprocs==1) {
      TEST_EQUALITY(lftCnt,2);
      TEST_EQUALITY(rhtCnt,2);
      TEST_EQUALITY(topCnt,5);
      TEST_EQUALITY(btmCnt,5);
   }
   else if(rank==0) {
      TEST_EQUALITY(lftCnt,2);
      TEST_EQUALITY(rhtCnt,0);
      TEST_EQUALITY(topCnt,4);
      TEST_EQUALITY(btmCnt,4);
   }
   else if(rank==1) {
      TEST_EQUALITY(lftCnt,0);
      TEST_EQUALITY(rhtCnt,2);
      TEST_EQUALITY(topCnt,4);
      TEST_EQUALITY(btmCnt,4);
   }
   else
      TEUCHOS_ASSERT(false);
}

void test1(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm);
void test2(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm);
void test4(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm);

using Teuchos::RCP;

void entityVecToGIDVec(RCP<STK_Interface> mesh,
                       const std::vector<stk::mesh::Entity> & eVec,
                             std::vector<stk::mesh::EntityId> & gidVec)
{
   gidVec.resize(eVec.size());
   for(std::size_t i=0;i<eVec.size();i++)
      gidVec[i] = mesh->elementGlobalId(eVec[i]);

   std::sort(gidVec.begin(),gidVec.end());
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, sideset_nodeset)
{
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",3);
   pl->set("X Elements",20);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   std::vector<std::string> sidesets, nodesets;

   mesh->getSidesetNames(sidesets);
   mesh->getNodesetNames(nodesets);

   TEST_EQUALITY(sidesets.size(),7);
   TEST_EQUALITY(nodesets.size(),2);

   std::vector<stk::mesh::Entity> nodes;
   std::vector<stk::mesh::Entity> nodes_o;
   mesh->getMyNodes("lower_left","eblock-0_0",nodes);
   mesh->getMyNodes("origin","eblock-0_0",nodes_o);
   if(rank==0) {
      {
         std::vector<std::size_t> localNodeIds;
         std::vector<stk::mesh::Entity> elements;

         TEST_EQUALITY(nodes.size(),1);
         workset_utils::getNodeElements(*mesh,"eblock-0_0",nodes,localNodeIds,elements);

         TEST_EQUALITY(localNodeIds.size(),1);
         TEST_EQUALITY(elements.size(),1);
         TEST_EQUALITY(mesh->elementGlobalId(elements[0]),1);
         TEST_EQUALITY(localNodeIds[0],0);
      }
      {
         std::vector<std::size_t> localNodeIds;
         std::vector<stk::mesh::Entity> elements;

         TEST_EQUALITY(nodes.size(),1);
         workset_utils::getNodeElements(*mesh,"eblock-0_0",nodes_o,localNodeIds,elements);

         TEST_EQUALITY(localNodeIds.size(),1);
         TEST_EQUALITY(elements.size(),1);
         TEST_EQUALITY(mesh->elementGlobalId(elements[0]),1);
         TEST_EQUALITY(localNodeIds[0],0);
      }
   }
   else {
      TEST_EQUALITY(nodes.size(),0);
   }
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, element_counts)
{
   int colors[] =
      { 1,
        2,2,
        4,4,4,4 };
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
   };
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory,rebalance)
{
  // Temporarily disabled until we get rebalance into Trilinos
   // int size = 0;
   // MPI_Comm_size(MPI_COMM_WORLD, &size);

   // RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   // pl->set("X Elements",5);
   // pl->set("Y Elements",5);

   // SquareQuadMeshFactory factory;
   // factory.setParameterList(pl);
   // RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   // TEST_ASSERT(mesh!=Teuchos::null);

   // Teuchos::ParameterList params;

   // mesh->rebalance(params);
   // mesh->buildLocalElementIDs();

   // // check to make sure that user specified parameters work
   // Teuchos::ParameterList zoltan_params;
   // zoltan_params.set("ZOLTAN DEBUG LEVEL","10");
   // mesh->rebalance(zoltan_params);
   // mesh->buildLocalElementIDs();

   // // check the size for the repartitioned mesh
   // if(size==2) {
   //   std::vector<stk::mesh::Entity> elements;
   //   mesh->getMyElements(elements);
   //   TEST_ASSERT(elements.size()==12 || elements.size()==13);
   // }
}

void test1(Teuchos::FancyOStream &out, bool &success, MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,1);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",1);
   pl->set("Y Procs",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),2);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1));
}

void test2(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,2);
   int rank; MPI_Comm_rank(comm, &rank);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",1);
   pl->set("Y Procs",2);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),2);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1));

   std::vector<stk::mesh::Entity> myElements;
   std::vector<stk::mesh::EntityId> myGids;
   mesh->getMyElements(myElements);
   entityVecToGIDVec(mesh,myElements,myGids);

   if(rank==0) {
      TEST_EQUALITY(myGids.size(),4);
      for(std::size_t i=0;i<2;i++) {
         TEST_EQUALITY(myGids[2*i]  ,2*i+1);
         TEST_EQUALITY(myGids[2*i+1],2*i+2);
      }
   }
   else if(rank==1) {
      TEST_EQUALITY(myGids.size(),4);
      for(std::size_t i=0;i<2;i++) {
         TEST_EQUALITY(myGids[2*i]  ,2*i+5);
         TEST_EQUALITY(myGids[2*i+1],2*i+6);
      }
   }
   else TEST_ASSERT(false);
}

void test4(Teuchos::FancyOStream &out, bool &success,MPI_Comm & comm)
{
   int size; MPI_Comm_size(comm, &size); TEST_EQUALITY(size,4);
   int rank; MPI_Comm_rank(comm, &rank);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Procs",2);
   pl->set("Y Procs",2);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);

   SquareQuadMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(comm);
   TEST_ASSERT(mesh!=Teuchos::null);

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),2);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*5+4*3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1));

   // if(mesh->isWritable());
   //    mesh->writeToExodus("test.exo");

   std::vector<stk::mesh::Entity> myElements;
   std::vector<stk::mesh::EntityId> myGids;
   mesh->getMyElements(myElements);
   entityVecToGIDVec(mesh,myElements,myGids);

   if(rank==0) {
      TEST_EQUALITY(myGids.size(),2);
      TEST_EQUALITY(myGids[0],1);
      TEST_EQUALITY(myGids[1],3);
   }
   else if(rank==1) {
      TEST_EQUALITY(myGids.size(),2);
      TEST_EQUALITY(myGids[0],2);
      TEST_EQUALITY(myGids[1],4);
   }
   else if(rank==2) {
      TEST_EQUALITY(myGids.size(),2);
      TEST_EQUALITY(myGids[0],5);
      TEST_EQUALITY(myGids[1],7);
   }
   else if(rank==3) {
      TEST_EQUALITY(myGids.size(),2);
      TEST_EQUALITY(myGids[0],6);
      TEST_EQUALITY(myGids[1],8);
   }
   else TEST_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, default_edge_blocks)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);

   std::size_t expected_edge_block_count = 0;

   edge_block_test_helper(out, success, pl,
                          "SquareQuad_default_edge_blocks.exo",
                          expected_edge_block_count);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, create_edge_blocks_pl)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Create Edge Blocks",true);

   std::size_t expected_edge_block_count = 1;

   edge_block_test_helper(out, success, pl,
                          "SquareQuad_create_edge_blocks_pl.exo",
                          expected_edge_block_count);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshFactory, multiblock_create_edge_blocks_pl)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 2, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Create Edge Blocks",true);

   std::size_t expected_edge_block_count = 1;

   edge_block_test_helper(out, success, pl,
                          "SquareQuad_multiblock_create_edge_blocks_pl.exo",
                          expected_edge_block_count);
}
}
