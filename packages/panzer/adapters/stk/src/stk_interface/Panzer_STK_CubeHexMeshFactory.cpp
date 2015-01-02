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

#include <Panzer_STK_CubeHexMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Panzer_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk_classic {

CubeHexMeshFactory::CubeHexMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
CubeHexMeshFactory::~CubeHexMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> CubeHexMeshFactory::buildMesh(stk_classic::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeHexMeshFactory::buildMesh()");

   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> CubeHexMeshFactory::buildUncommitedMesh(stk_classic::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeHexMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

   machRank_ = stk_classic::parallel_machine_rank(parallelMach);
   machSize_ = stk_classic::parallel_machine_size(parallelMach);

   if(xProcs_==-1) {
      // default x only decomposition
      xProcs_ = machSize_; 
      yProcs_ = 1;
      zProcs_ = 1;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(int(machSize_)!=xProcs_*yProcs_*zProcs_,std::logic_error,
                      "Cannot build CubeHexMeshFactory, the product of \"X Procs\", \"Y Procs\", and \"Z Procs\""
                      " must equal the number of processors.");
   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);
 
   mesh->addPeriodicBCs(periodicBCVec_);

   return mesh;
}

void CubeHexMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk_classic::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeHexMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setOutputToRootOnly(0);
   out.setShowProcRank(true);

   // finish up the edges and faces
   if(buildSubcells_) {
      mesh.buildSubcells();

      out << "CubeHexMesh: Building sub cells" << std::endl;
   }
   else {
      addSides(mesh);

      out << "CubeHexMesh: NOT building sub cells" << std::endl;
   }

   mesh.buildLocalElementIDs();

   // now that edges are built, side and node sets can be added
   addSideSets(mesh);
   addNodeSets(mesh);

   // calls Stk_MeshFactory::rebalance
   this->rebalance(mesh);
}

//! From ParameterListAcceptor
void CubeHexMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   x0_ = paramList->get<double>("X0"); 
   y0_ = paramList->get<double>("Y0"); 
   z0_ = paramList->get<double>("Z0"); 

   xf_ = paramList->get<double>("Xf"); 
   yf_ = paramList->get<double>("Yf"); 
   zf_ = paramList->get<double>("Zf"); 

   xBlocks_ = paramList->get<int>("X Blocks");
   yBlocks_ = paramList->get<int>("Y Blocks");
   zBlocks_ = paramList->get<int>("Z Blocks");

   xProcs_ = paramList->get<int>("X Procs");
   yProcs_ = paramList->get<int>("Y Procs");
   zProcs_ = paramList->get<int>("Z Procs");

   nXElems_ = paramList->get<int>("X Elements");
   nYElems_ = paramList->get<int>("Y Elements");
   nZElems_ = paramList->get<int>("Z Elements");

   buildSubcells_ = paramList->get<bool>("Build Subcells");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> CubeHexMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<double>("X0",0.0);
      defaultParams->set<double>("Y0",0.0);
      defaultParams->set<double>("Z0",0.0);

      defaultParams->set<double>("Xf",1.0);
      defaultParams->set<double>("Yf",1.0);
      defaultParams->set<double>("Zf",1.0);

      defaultParams->set<int>("X Blocks",1);
      defaultParams->set<int>("Y Blocks",1);
      defaultParams->set<int>("Z Blocks",1);

      defaultParams->set<int>("X Procs",-1);
      defaultParams->set<int>("Y Procs",1);
      defaultParams->set<int>("Z Procs",1);

      defaultParams->set<int>("X Elements",5);
      defaultParams->set<int>("Y Elements",5);
      defaultParams->set<int>("Z Elements",5);

      defaultParams->set<bool>("Build Subcells",true);

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void CubeHexMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);
}

void CubeHexMeshFactory::buildMetaData(stk_classic::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Hexahedron<8> HexTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

   // build meta data
   //mesh.setDimension(2);
   for(int bx=0;bx<xBlocks_;bx++) {
      for(int by=0;by<yBlocks_;by++) {
         for(int bz=0;bz<zBlocks_;bz++) {

            std::stringstream ebPostfix;
            ebPostfix << "-" << bx << "_" << by << "_" << bz;

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
   mesh.addSideset("front",side_ctd);
   mesh.addSideset("back",side_ctd);

   mesh.addNodeset("origin");
}

void CubeHexMeshFactory::buildElements(stk_classic::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();
      // build each block
      for(int xBlock=0;xBlock<xBlocks_;xBlock++) {
         for(int yBlock=0;yBlock<yBlocks_;yBlock++) {
            for(int zBlock=0;zBlock<zBlocks_;zBlock++) {
               buildBlock(parallelMach,xBlock,yBlock,zBlock,mesh);
            }
         }
      }
   mesh.endModification();
}

void CubeHexMeshFactory::buildBlock(stk_classic::ParallelMachine parallelMach,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<panzer::Ordinal64,panzer::Ordinal64> sizeAndStartX = determineXElemSizeAndStart(xBlock,xProcs_,machRank_);
   std::pair<panzer::Ordinal64,panzer::Ordinal64> sizeAndStartY = determineYElemSizeAndStart(yBlock,yProcs_,machRank_);
   std::pair<panzer::Ordinal64,panzer::Ordinal64> sizeAndStartZ = determineZElemSizeAndStart(zBlock,zProcs_,machRank_);

   panzer::Ordinal64 myXElems_start = sizeAndStartX.first;
   panzer::Ordinal64 myXElems_end  = myXElems_start+sizeAndStartX.second;
   panzer::Ordinal64 myYElems_start = sizeAndStartY.first;
   panzer::Ordinal64 myYElems_end  = myYElems_start+sizeAndStartY.second;
   panzer::Ordinal64 myZElems_start = sizeAndStartZ.first;
   panzer::Ordinal64 myZElems_end  = myZElems_start+sizeAndStartZ.second;

   panzer::Ordinal64 totalXElems = nXElems_*xBlocks_;
   panzer::Ordinal64 totalYElems = nYElems_*yBlocks_;
   panzer::Ordinal64 totalZElems = nZElems_*zBlocks_;

   double deltaX = (xf_-x0_)/double(totalXElems);
   double deltaY = (yf_-y0_)/double(totalYElems);
   double deltaZ = (zf_-z0_)/double(totalZElems);
 
   std::vector<double> coord(3,0.0);

   // build the nodes
   for(panzer::Ordinal64 nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = double(nx)*deltaX+x0_;
      for(panzer::Ordinal64 ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = double(ny)*deltaY+y0_;
         for(panzer::Ordinal64 nz=myZElems_start;nz<myZElems_end+1;++nz) {
            coord[2] = double(nz)*deltaZ+z0_;

            mesh.addNode(nz*(totalYElems+1)*(totalXElems+1)+ny*(totalXElems+1)+nx+1,coord);
         }
      }
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock << "_" << yBlock << "_" << zBlock;
   stk_classic::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(panzer::Ordinal64 nx=myXElems_start;nx<myXElems_end;++nx) {
      for(panzer::Ordinal64 ny=myYElems_start;ny<myYElems_end;++ny) {
         for(panzer::Ordinal64 nz=myZElems_start;nz<myZElems_end;++nz) {
            stk_classic::mesh::EntityId gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1;
            std::vector<stk_classic::mesh::EntityId> nodes(8);
            nodes[0] = nx+1+ny*(totalXElems+1) +nz*(totalYElems+1)*(totalXElems+1);
            nodes[1] = nodes[0]+1;              
            nodes[2] = nodes[1]+(totalXElems+1);
            nodes[3] = nodes[2]-1;              
            nodes[4] = nodes[0]+(totalYElems+1)*(totalXElems+1);
            nodes[5] = nodes[1]+(totalYElems+1)*(totalXElems+1);
            nodes[6] = nodes[2]+(totalYElems+1)*(totalXElems+1);
            nodes[7] = nodes[3]+(totalYElems+1)*(totalXElems+1);
   
            RCP<ElementDescriptor> ed = rcp(new ElementDescriptor(gid,nodes));
            mesh.addElement(ed,block);
         }
      }
   }
}

std::pair<panzer::Ordinal64,panzer::Ordinal64> CubeHexMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const
{
   std::size_t xProcLoc = procTuple_[0];
   panzer::Ordinal64 minElements = nXElems_/size;
   panzer::Ordinal64 extra = nXElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::Ordinal64 nume=0, start=0;
   if(panzer::Ordinal64(xProcLoc)<extra) {
      nume  = minElements+1;
      start = xProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(xProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nXElems_*xBlock,nume);
}

std::pair<panzer::Ordinal64,panzer::Ordinal64> CubeHexMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const
{
   // int start = yBlock*nYElems_;
   // return std::make_pair(start,nYElems_);

   std::size_t yProcLoc = procTuple_[1];
   panzer::Ordinal64 minElements = nYElems_/size;
   panzer::Ordinal64 extra = nYElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::Ordinal64 nume=0, start=0;
   if(panzer::Ordinal64(yProcLoc)<extra) {
      nume  = minElements+1;
      start = yProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(yProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nYElems_*yBlock,nume);
}

std::pair<panzer::Ordinal64,panzer::Ordinal64> CubeHexMeshFactory::determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int rank) const
{
   // int start = zBlock*nZElems_;
   // return std::make_pair(start,nZElems_);
   std::size_t zProcLoc = procTuple_[2];
   panzer::Ordinal64 minElements = nZElems_/size;
   panzer::Ordinal64 extra = nZElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::Ordinal64 nume=0, start=0;
   if(zProcLoc<Teuchos::as<std::size_t>(extra)) {
      nume  = minElements+1;
      start = zProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(zProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nZElems_*zBlock,nume);
}

// for use with addSideSets
const stk_classic::mesh::Relation * CubeHexMeshFactory::getRelationByID(unsigned ID,stk_classic::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) {
      if(relations[i].identifier()==ID)
         return &relations[i];
   }

   return 0;
}

// this adds side entities only (does not inject them into side sets)
void CubeHexMeshFactory::addSides(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;
   std::size_t totalZElems = nZElems_*zBlocks_;

   std::vector<stk_classic::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   stk_classic::mesh::EntityId offset[6];
   offset[0] = 0;
   offset[1] = offset[0] + totalXElems*totalZElems;
   offset[2] = offset[1] + totalYElems*totalZElems;
   offset[3] = offset[2] + totalXElems*totalZElems;
   offset[4] = offset[3] + totalYElems*totalZElems;
   offset[5] = offset[4] + totalXElems*totalYElems;

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk_classic::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk_classic::mesh::Entity * element = (*itr);
      stk_classic::mesh::EntityId gid = element->identifier();      
      stk_classic::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

      std::size_t nx,ny,nz;
      nz = (gid-1) / (totalXElems*totalYElems);
      gid = (gid-1)-nz*(totalXElems*totalYElems);
      ny = gid / totalXElems;
      nx = gid-ny*totalXElems;

      std::vector<stk_classic::mesh::Part*> parts;

      if(nz==0) {
         // on the back
         stk_classic::mesh::EntityId eid = (1+nx+ny*totalXElems)+offset[4];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,4);
      }
      if(nz+1==totalZElems) {
         // on the front
         stk_classic::mesh::EntityId eid = (1+nx+ny*totalXElems)+offset[5];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,5);
      }

      if(ny==0) {
         // on the bottom 
         stk_classic::mesh::EntityId eid = (1+nx+nz*totalXElems)+offset[0];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,0);
      }
      if(ny+1==totalYElems) {
         // on the top
         stk_classic::mesh::EntityId eid = (1+nx+nz*totalXElems)+offset[2];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,2);
      }

      if(nx==0) {
         // on the left
         stk_classic::mesh::EntityId eid = (1+ny+nz*totalYElems)+offset[3];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,3);
      }
      if(nx+1==totalXElems) {
         // on the right
         stk_classic::mesh::EntityId eid = (1+ny+nz*totalYElems)+offset[1];
         stk_classic::mesh::Entity & side = mesh.getBulkData()->declare_entity(mesh.getSideRank(),eid,parts);
         mesh.getBulkData()->declare_relation(*element,side,1);
      }
   }

   mesh.endModification();
}

void CubeHexMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;
   std::size_t totalZElems = nZElems_*zBlocks_;

   // get all part vectors
   stk_classic::mesh::Part * left = mesh.getSideset("left");
   stk_classic::mesh::Part * right = mesh.getSideset("right");
   stk_classic::mesh::Part * top = mesh.getSideset("top");
   stk_classic::mesh::Part * bottom = mesh.getSideset("bottom");
   stk_classic::mesh::Part * front = mesh.getSideset("front");
   stk_classic::mesh::Part * back = mesh.getSideset("back");

   std::vector<stk_classic::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk_classic::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk_classic::mesh::Entity * element = (*itr);
      stk_classic::mesh::EntityId gid = element->identifier();      
      stk_classic::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

      std::size_t nx,ny,nz;
      nz = (gid-1) / (totalXElems*totalYElems);
      gid = (gid-1)-nz*(totalXElems*totalYElems);
      ny = gid / totalXElems;
      nx = gid-ny*totalXElems;

      if(nz==0) {
         stk_classic::mesh::Entity * side = getRelationByID(4,relations)->entity();

         // on the back
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,back);
      }
      if(nz+1==totalZElems) {
         stk_classic::mesh::Entity * side = getRelationByID(5,relations)->entity();

         // on the front
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,front);
      }

      if(ny==0) {
         stk_classic::mesh::Entity * side = getRelationByID(0,relations)->entity();

         // on the bottom 
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,bottom);
      }
      if(ny+1==totalYElems) {
         stk_classic::mesh::Entity * side = getRelationByID(2,relations)->entity();

         // on the top
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,top);
      }

      if(nx==0) {
         stk_classic::mesh::Entity * side = getRelationByID(3,relations)->entity();

         // on the left
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,left);
      }
      if(nx+1==totalXElems) {
         stk_classic::mesh::Entity * side = getRelationByID(1,relations)->entity();

         // on the right
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,right);
      }
   }

   mesh.endModification();
}

void CubeHexMeshFactory::addNodeSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // get all part vectors
   stk_classic::mesh::Part * origin = mesh.getNodeset("origin");

   Teuchos::RCP<stk_classic::mesh::BulkData> bulkData = mesh.getBulkData();
   if(machRank_==0) 
   {
      // add zero node to origin node set
      stk_classic::mesh::Entity * node = bulkData->get_entity(mesh.getNodeRank(),1);
      mesh.addEntityToNodeset(*node,origin);
   }

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,3> CubeHexMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0,k=0;

   k = procRank/(xProcs_*yProcs_); procRank = procRank % (xProcs_*yProcs_);
   j = procRank/xProcs_;           procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j,k);
}

} // end panzer_stk
