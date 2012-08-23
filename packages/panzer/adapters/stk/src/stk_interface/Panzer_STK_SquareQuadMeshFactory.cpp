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

#include <Panzer_STK_SquareQuadMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Panzer_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

SquareQuadMeshFactory::SquareQuadMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
SquareQuadMeshFactory::~SquareQuadMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> SquareQuadMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareQuadMeshFactory::buildMesh()");

   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> SquareQuadMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareQuadMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(2));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   if(xProcs_==-1) {
      // default x only decomposition
      xProcs_ = machSize_; 
      yProcs_ = 1;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(int(machSize_)!=xProcs_*yProcs_,std::logic_error,
                      "Cannot build SquareQuadMeshFactory, the product of \"X Procs\" and \"Y Procs\""
                      " must equal the number of processors.");
   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
 
   return mesh;
}

void SquareQuadMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareQuadMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges
   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   // now that edges are built, sidsets can be added
   addSideSets(mesh);

   // add nodesets
   addNodeSets(mesh);
}

//! From ParameterListAcceptor
void SquareQuadMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   x0_ = paramList->get<double>("X0"); 
   y0_ = paramList->get<double>("Y0"); 

   xf_ = paramList->get<double>("Xf"); 
   yf_ = paramList->get<double>("Yf"); 

   xBlocks_ = paramList->get<int>("X Blocks");
   yBlocks_ = paramList->get<int>("Y Blocks");

   nXElems_ = paramList->get<int>("X Elements");
   nYElems_ = paramList->get<int>("Y Elements");

   xProcs_ = paramList->get<int>("X Procs");
   yProcs_ = paramList->get<int>("Y Procs");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> SquareQuadMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<double>("X0",0.0);
      defaultParams->set<double>("Y0",0.0);

      defaultParams->set<double>("Xf",1.0);
      defaultParams->set<double>("Yf",1.0);

      defaultParams->set<int>("X Blocks",1);
      defaultParams->set<int>("Y Blocks",1);

      defaultParams->set<int>("X Procs",-1);
      defaultParams->set<int>("Y Procs",1);

      defaultParams->set<int>("X Elements",5);
      defaultParams->set<int>("Y Elements",5);

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void SquareQuadMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);
}

void SquareQuadMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Quadrilateral<4> QuadTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   // build meta data
   //mesh.setDimension(2);
   for(int bx=0;bx<xBlocks_;bx++) {
      for(int by=0;by<yBlocks_;by++) {

         // add this element block
         {
            std::stringstream ebPostfix;
            ebPostfix << "-" << bx << "_" << by;

            // add element blocks
            mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
         }

      }
   }

   // add sidesets 
   mesh.addSideset("left",side_ctd);
   mesh.addSideset("right",side_ctd);
   mesh.addSideset("top",side_ctd);
   mesh.addSideset("bottom",side_ctd);

   // add nodesets
   mesh.addNodeset("lower_left");
}

void SquareQuadMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();
      // build each block
      for(int xBlock=0;xBlock<xBlocks_;xBlock++) {
         for(int yBlock=0;yBlock<yBlocks_;yBlock++) {
            buildBlock(parallelMach,xBlock,yBlock,mesh);
         }
      }
   mesh.endModification();
}

void SquareQuadMeshFactory::buildBlock(stk::ParallelMachine parallelMach,int xBlock,int yBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<int,int> sizeAndStartX = determineXElemSizeAndStart(xBlock,xProcs_,machRank_);
   std::pair<int,int> sizeAndStartY = determineYElemSizeAndStart(yBlock,yProcs_,machRank_);

   int myXElems_start = sizeAndStartX.first;
   int myXElems_end  = myXElems_start+sizeAndStartX.second;
   int myYElems_start = sizeAndStartY.first;
   int myYElems_end  = myYElems_start+sizeAndStartY.second;
   int totalXElems = nXElems_*xBlocks_;
   int totalYElems = nYElems_*yBlocks_;

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

std::pair<int,int> SquareQuadMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const
{
   std::size_t xProcLoc = procTuple_[0];
   unsigned int minElements = nXElems_/size;
   unsigned int extra = nXElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   int nume=0, start=0;
   if(xProcLoc<extra) {
      nume  = minElements+1;
      start = xProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(xProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nXElems_*xBlock,nume);
}

std::pair<int,int> SquareQuadMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const
{
   std::size_t yProcLoc = procTuple_[1];
   unsigned int minElements = nYElems_/size;
   unsigned int extra = nYElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   int nume=0, start=0;
   if(yProcLoc<extra) {
      nume  = minElements+1;
      start = yProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(yProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nYElems_*yBlock,nume);
}

const stk::mesh::Relation * SquareQuadMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}

void SquareQuadMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;

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

void SquareQuadMeshFactory::addNodeSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // get all part vectors
   stk::mesh::Part * lower_left = mesh.getNodeset("lower_left");

   std::vector<stk::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   if(machRank_==0) 
   {
      // add zero node to lower_left node set
      stk::mesh::Entity * node = bulkData->get_entity(mesh.getNodeRank(),1);
      mesh.addEntityToNodeset(*node,lower_left);
   }

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,2> SquareQuadMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0;

   j = procRank/xProcs_; 
   procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j);
}

} // end panzer_stk
