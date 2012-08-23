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

#include <Panzer_STK_MultiBlockMeshFactory.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

MultiBlockMeshFactory::MultiBlockMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
MultiBlockMeshFactory::~MultiBlockMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> MultiBlockMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> MultiBlockMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   RCP<STK_Interface> mesh = rcp(new STK_Interface(2));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   return mesh;
}

void MultiBlockMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges
   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   // now that edges are built, sidets can be added
   addSideSets(mesh);
}

//! From ParameterListAcceptor
void MultiBlockMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> MultiBlockMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<int>("X Elements",2);
      defaultParams->set<int>("Y Elements",2);
   }

   return defaultParams;
}

void MultiBlockMeshFactory::initializeWithDefaults()
{
   nXElems_ = 2;
   nYElems_ = 2;
}

void MultiBlockMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Quadrilateral<4> QuadTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   // build meta data
   //mesh.setDimension(2);
   for(int bx=0;bx<2;bx++) {
      // add this element block
      std::stringstream ebPostfix;
      ebPostfix << "-" << "0_" <<  bx;

      // add element blocks
      mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
   }

   // add sidesets 
   mesh.addSideset("left",side_ctd);
   mesh.addSideset("right",side_ctd);
   mesh.addSideset("top",side_ctd);
   mesh.addSideset("bottom",side_ctd);
}

void MultiBlockMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();

   if(machRank_==0) {
      buildBlock(parallelMach,0,0,mesh);
   }
   else if(machRank_==1) {
      buildBlock(parallelMach,0,1,mesh);
   }
   else TEUCHOS_ASSERT(false);

   mesh.endModification();
}

void MultiBlockMeshFactory::buildBlock(stk::ParallelMachine parallelMach,int xBlock,int yBlock,STK_Interface & mesh) const
{
   int myXElems_start = (machRank_==0 ? 0 : 2);
   int myXElems_end  = (machRank_==0 ? 1 : 3);
   int myYElems_start = 0;
   int myYElems_end  = 1;
   int totalXElems = nXElems_*2;
   int totalYElems = nYElems_*1;

   double x0_ = 0.0;
   double xf_ = 1.0;
   double y0_ = 0.0;
   double yf_ = 1.0;
   double deltaX = (xf_-x0_)/double(totalXElems);
   double deltaY = (yf_-y0_)/double(totalYElems);
 
   std::vector<double> coord(2,0.0);

   // build the nodes
   for(int nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = double(nx)*deltaX+x0_;
      for(int ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = double(ny)*deltaY+y0_;

         mesh.addNode(ny*(totalXElems+1)+nx+1,coord);
      }
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock << "_" << yBlock;
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(int nx=myXElems_start;nx<myXElems_end;++nx) {
      for(int ny=myYElems_start;ny<myYElems_end;++ny) {
         stk::mesh::EntityId gid = totalXElems*ny+nx+1;
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = nx+1+ny*(totalXElems+1);
         nodes[1] = nodes[0]+1;
         nodes[2] = nodes[1]+(totalXElems+1);
         nodes[3] = nodes[2]-1;

         RCP<ElementDescriptor> ed = rcp(new ElementDescriptor(gid,nodes));
         mesh.addElement(ed,block);
      }
   }
}

std::pair<int,int> MultiBlockMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const
{
   unsigned int minElements = nXElems_/size;
   unsigned int extra = nXElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   int nume=0, start=0;
   if(rank<extra) {
      nume  = minElements+1;
      start = rank*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(rank-extra)*minElements;
   }

   return std::make_pair(start+nXElems_*xBlock,nume);
}

std::pair<int,int> MultiBlockMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const
{
   int start = yBlock*nYElems_;

   return std::make_pair(start,nYElems_);
}

const stk::mesh::Relation * MultiBlockMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}

void MultiBlockMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*2;
   std::size_t totalYElems = nYElems_*1;

   // get all part vectors
   stk::mesh::Part * left = mesh.getSideset("left");
   stk::mesh::Part * right = mesh.getSideset("right");
   stk::mesh::Part * top = mesh.getSideset("top");
   stk::mesh::Part * bottom = mesh.getSideset("bottom");

   std::vector<stk::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   // loop over elements adding edges to sidesets
   std::vector<stk::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity * element = (*itr);
      stk::mesh::EntityId gid = element->identifier();      
      stk::mesh::PairIterRelation relations = element->relations(mesh.getEdgeRank());

      std::size_t nx,ny;
      ny = (gid-1) / totalXElems;
      nx = gid-ny*totalXElems-1;

      // vertical boundaries
      ///////////////////////////////////////////

      if(nx+1==totalXElems) { 
         stk::mesh::Entity * edge = getRelationByID(1,relations)->entity();

         // on the right
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,right);
      }

      if(nx==0) {
         stk::mesh::Entity * edge = getRelationByID(3,relations)->entity();

         // on the left
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,left);
      }

      // horizontal boundaries
      ///////////////////////////////////////////

      if(ny==0) {
         stk::mesh::Entity * edge = getRelationByID(0,relations)->entity();

         // on the bottom
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,bottom);
      }

      if(ny+1==totalYElems) {
         stk::mesh::Entity * edge = getRelationByID(2,relations)->entity();

         // on the top
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,top);
      }
   }

   mesh.endModification();
}

} // end panzer_stk
