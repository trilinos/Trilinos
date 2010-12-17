#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

namespace panzer_stk {

inline bool XOR(bool A,bool B)
{ return ! ( (A && B) || ( !A && !B)); }

class LocalIdCompare {
public:
   LocalIdCompare(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   bool operator()(stk::mesh::Entity * a,stk::mesh::Entity * b) const 
   { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b); }

private:
   Teuchos::RCP<const STK_Interface> mesh_;
};

static void getNodeIds(stk::mesh::EntityRank nodeRank,const stk::mesh::Entity * element,std::vector<stk::mesh::EntityId> & nodeIds)
{
   stk::mesh::PairIterRelation nodeRel = element->relations(nodeRank);

   stk::mesh::PairIterRelation::iterator itr;
   for(itr=nodeRel.begin();itr!=nodeRel.end();++itr)
      nodeIds.push_back(itr->entity()->identifier());
}

static const double * getNode(const Teuchos::RCP<const STK_Interface> & mesh, const stk::mesh::Entity * element,int id)
{
   std::vector<stk::mesh::EntityId> nodeIds;
   getNodeIds(mesh->getNodeRank(),element,nodeIds);

   return mesh->getNodeCoordinates(nodeIds[id]); 
}

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

   std::string strBlock0="eblock-0_0", strBlock1="eblock-1_0";
   std::vector<stk::mesh::Entity*> block0, block1;

   mesh->getMyElements(strBlock0,block0);
   mesh->getMyElements(strBlock1,block1);

   std::sort(block0.begin(),block0.end(),LocalIdCompare(mesh));
   std::sort(block1.begin(),block1.end(),LocalIdCompare(mesh));

   if(numprocs==1) {
      TEST_EQUALITY(block0.size(),6);
      TEST_EQUALITY(block1.size(),6);

      for(std::size_t i=0;i<block0.size();i++) {
         TEST_EQUALITY(mesh->elementLocalId(block0[i]),i);
         TEST_EQUALITY(mesh->elementLocalId(block1[i]),i+6);
      }

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
   }
   else if(numprocs==2 && rank==0) {
      TEST_EQUALITY(block0.size(),3);
      TEST_EQUALITY(block1.size(),3);

      for(std::size_t i=0;i<block0.size();i++) {
         TEST_EQUALITY(mesh->elementLocalId(block0[i]),i);
         TEST_EQUALITY(mesh->elementLocalId(block1[i]),i+3);
      }

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
   }
   else if(numprocs==2 && rank==1) {
      TEST_EQUALITY(block0.size(),3);
      TEST_EQUALITY(block1.size(),3);

      for(std::size_t i=0;i<block0.size();i++) {
         TEST_EQUALITY(mesh->elementLocalId(block0[i]),i);
         TEST_EQUALITY(mesh->elementLocalId(block1[i]),i+3);
      }

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
 
   if(mesh->isWritable());
      mesh->writeToExodus("Square.exo");

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
   if(mesh->isWritable());
      mesh->writeToExodus("Square_Blocked.exo");

   TEST_EQUALITY(mesh->getNumElementBlocks(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(12+1)*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),12*(12+1)+12*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*12);

   std::vector<stk::mesh::Entity*> myElements;
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

   TEST_EQUALITY(mesh->getNumElementBlocks(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(12+1)*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),12*(12+1)+12*(12+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*12);

   std::vector<std::string> blockNames;
   mesh->getElementBlockNames(blockNames);
   TEST_EQUALITY(blockNames.size(),mesh->getNumElementBlocks());
   TEST_EQUALITY(blockNames[0],"eblock-0_0");

   {
      std::vector<stk::mesh::Entity*> myElements;
      mesh->getMyElements(blockNames[0],myElements);

      TEST_EQUALITY((int) myElements.size(),24/numprocs);
      
      TEST_EQUALITY((int) myElements[0]->identifier(),1+rank*3); 
      TEST_EQUALITY((int) myElements[1]->identifier(),2+rank*3); 
      TEST_EQUALITY((int) myElements[2]->identifier(),3+rank*3); 

      TEST_EQUALITY((int) myElements[24/numprocs-3]->identifier(),40-(numprocs-1)*(1-rank)*3); 
      TEST_EQUALITY((int) myElements[24/numprocs-2]->identifier(),41-(numprocs-1)*(1-rank)*3); 
      TEST_EQUALITY((int) myElements[24/numprocs-1]->identifier(),42-(numprocs-1)*(1-rank)*3); 
   }

   {
      {
         std::vector<stk::mesh::Entity*> mySides;

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
         std::vector<stk::mesh::Entity*> mySides;
         mesh->getMySides(*sItr,mySides);
   
         std::vector<stk::mesh::Entity*>::iterator itr;
         for(itr=mySides.begin();itr!=mySides.end();++itr) {
            stk::mesh::Entity * side = *itr;
            stk::mesh::PairIterRelation relations = side->relations(mesh->getNodeRank());
            stk::mesh::EntityId n0 = relations[0].entity()->identifier();
            stk::mesh::EntityId n1 = relations[1].entity()->identifier();
   
            TEST_EQUALITY(side->entity_rank(),mesh->getSideRank());
            TEST_EQUALITY((int) side->relations().size(),3);
            TEST_EQUALITY(relations.size(),2);
            // TEST_EQUALITY(side->identifier(),mesh->getEdgeId(n0,n1));
         }
         
      }
   }

   {
      {
         std::vector<stk::mesh::Entity*> mySides;

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
         std::vector<stk::mesh::Entity*> mySides;
         mesh->getMySides(sidesets[i],eblocks[i],mySides);
   
         std::vector<stk::mesh::Entity*>::iterator itr;
         for(itr=mySides.begin();itr!=mySides.end();++itr) {
            stk::mesh::Entity * side = *itr;
            stk::mesh::PairIterRelation relations = side->relations(mesh->getNodeRank());
            stk::mesh::EntityId n0 = relations[0].entity()->identifier();
            stk::mesh::EntityId n1 = relations[1].entity()->identifier();
   
            TEST_EQUALITY(side->entity_rank(),mesh->getSideRank());
            TEST_EQUALITY((int) side->relations().size(),3);
            TEST_EQUALITY(relations.size(),2);
            // TEST_EQUALITY(side->identifier(),mesh->getEdgeId(n0,n1));
         }
         
      }
   }
}

}
