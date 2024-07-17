// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_STK_CubeTetMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <PanzerAdaptersSTK_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

CubeTetMeshFactory::CubeTetMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
CubeTetMeshFactory::~CubeTetMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> CubeTetMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
  PANZER_FUNC_TIME_MONITOR("panzer::CubeTetMeshFactory::buildMesh()");

   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> CubeTetMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeTetMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   if (xProcs_ == -1 && yProcs_ == -1 && zProcs_ == -1) {
     // copied from galeri
     xProcs_ = yProcs_ = zProcs_ = Teuchos::as<int>(pow(Teuchos::as<double>(machSize_), 0.333334));

     if (xProcs_ * yProcs_ * zProcs_ != Teuchos::as<int>(machSize_))  {
       // Simple method to find a set of processor assignments
       xProcs_ = yProcs_ = zProcs_ = 1;

       // This means that this works correctly up to about maxFactor^3
       // processors.
       const int maxFactor = 50;

       int ProcTemp = machSize_;
       int factors[maxFactor];
       for (int jj = 0; jj < maxFactor; jj++) factors[jj] = 0;
       for (int jj = 2; jj < maxFactor; jj++) {
         bool flag = true;
         while (flag) {
           int temp = ProcTemp/jj;
           if (temp*jj == ProcTemp) {
             factors[jj]++;
             ProcTemp = temp;

           } else {
             flag = false;
           }
         }
       }
       xProcs_ = ProcTemp;
       for (int jj = maxFactor-1; jj > 0; jj--) {
         while (factors[jj] != 0) {
           if      ((xProcs_ <= yProcs_) && (xProcs_ <= zProcs_)) xProcs_ = xProcs_*jj;
           else if ((yProcs_ <= xProcs_) && (yProcs_ <= zProcs_)) yProcs_ = yProcs_*jj;
           else                               zProcs_ = zProcs_*jj;
           factors[jj]--;
         }
       }
     }

   } else if(xProcs_==-1) {
      // default x only decomposition
      xProcs_ = machSize_;
      yProcs_ = 1;
      zProcs_ = 1;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(int(machSize_)!=xProcs_*yProcs_*zProcs_,std::logic_error,
                      "Cannot build CubeTetMeshFactory, the product of \"X Procs\", \"Y Procs\", and \"Z Procs\""
                      " must equal the number of processors.");
   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void CubeTetMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeTetMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges and faces
   mesh.buildSubcells();
   mesh.buildLocalElementIDs();
   if(createEdgeBlocks_) {
      mesh.buildLocalEdgeIDs();
   }
   if(createFaceBlocks_) {
      mesh.buildLocalFaceIDs();
   }

   // now that edges are built, sidets can be added
   addSideSets(mesh);
   addNodeSets(mesh);

   mesh.beginModification();
   if(createEdgeBlocks_) {
      addEdgeBlocks(mesh);
   }
   if(createFaceBlocks_) {
      addFaceBlocks(mesh);
   }
   mesh.endModification();

   // calls Stk_MeshFactory::rebalance
   this->rebalance(mesh);
}

//! From ParameterListAcceptor
void CubeTetMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
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

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");
   createFaceBlocks_ = paramList->get<bool>("Create Face Blocks");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> CubeTetMeshFactory::getValidParameters() const
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

      // default to false for backward compatibility
      defaultParams->set<bool>("Create Edge Blocks",false,"Create edge blocks in the mesh");
      defaultParams->set<bool>("Create Face Blocks",false,"Create face blocks in the mesh");

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void CubeTetMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);

   /* This is a tet mesh factory so all elements in all element blocks
    * will be tet4.  This means that all the edges will be line2 and
    * all the faces will be tri3.  The edge and face block names are
    * hard coded to reflect this.
    */
   edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;
   faceBlockName_ = "tri_3_"+panzer_stk::STK_Interface::faceBlockString;
}

void CubeTetMeshFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
   typedef shards::Tetrahedron<4> TetTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<TetTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);
   const CellTopologyData * edge_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);
   const CellTopologyData * face_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

   // build meta data
   //mesh.setDimension(2);
   for(int bx=0;bx<xBlocks_;bx++) {
      for(int by=0;by<yBlocks_;by++) {
         for(int bz=0;bz<zBlocks_;bz++) {

            std::stringstream ebPostfix;
            ebPostfix << "-" << bx << "_" << by << "_" << bz;

            // add element blocks
            mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
            if(createEdgeBlocks_) {
               mesh.addEdgeBlock("eblock"+ebPostfix.str(),
                                 edgeBlockName_,
                                 edge_ctd);
            }
            if(createFaceBlocks_) {
               mesh.addFaceBlock("eblock"+ebPostfix.str(),
                                 faceBlockName_,
                                 face_ctd);
            }
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

void CubeTetMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
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

void CubeTetMeshFactory::buildBlock(stk::ParallelMachine /* parallelMach */,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<int,int> sizeAndStartX = determineXElemSizeAndStart(xBlock,xProcs_,machRank_);
   std::pair<int,int> sizeAndStartY = determineYElemSizeAndStart(yBlock,yProcs_,machRank_);
   std::pair<int,int> sizeAndStartZ = determineZElemSizeAndStart(zBlock,zProcs_,machRank_);

   int myXElems_start = sizeAndStartX.first;
   int myXElems_end  = myXElems_start+sizeAndStartX.second;
   int myYElems_start = sizeAndStartY.first;
   int myYElems_end  = myYElems_start+sizeAndStartY.second;
   int myZElems_start = sizeAndStartZ.first;
   int myZElems_end  = myZElems_start+sizeAndStartZ.second;

   int totalXElems = nXElems_*xBlocks_;
   int totalYElems = nYElems_*yBlocks_;
   int totalZElems = nZElems_*zBlocks_;

   double deltaX = (xf_-x0_)/double(totalXElems);
   double deltaY = (yf_-y0_)/double(totalYElems);
   double deltaZ = (zf_-z0_)/double(totalZElems);

   std::vector<double> coord(3,0.0);

   // build the nodes
   for(int nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = this->getMeshCoord(nx, deltaX, x0_);
      for(int ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = this->getMeshCoord(ny, deltaY, y0_);
         for(int nz=myZElems_start;nz<myZElems_end+1;++nz) {
            coord[2] = this->getMeshCoord(nz, deltaZ, z0_);

            mesh.addNode(nz*(totalYElems+1)*(totalXElems+1)+ny*(totalXElems+1)+nx+1,coord);
         }
      }
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock << "_" << yBlock << "_" << zBlock;
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(int nx=myXElems_start;nx<myXElems_end;++nx) {
      for(int ny=myYElems_start;ny<myYElems_end;++ny) {
         for(int nz=myZElems_start;nz<myZElems_end;++nz) {

            std::vector<stk::mesh::EntityId> nodes(8);
            nodes[0] = nx+1+ny*(totalXElems+1) +nz*(totalYElems+1)*(totalXElems+1);
            nodes[1] = nodes[0]+1;
            nodes[2] = nodes[1]+(totalXElems+1);
            nodes[3] = nodes[2]-1;
            nodes[4] = nodes[0]+(totalYElems+1)*(totalXElems+1);
            nodes[5] = nodes[1]+(totalYElems+1)*(totalXElems+1);
            nodes[6] = nodes[2]+(totalYElems+1)*(totalXElems+1);
            nodes[7] = nodes[3]+(totalYElems+1)*(totalXElems+1);

            buildTetsOnHex(Teuchos::tuple(totalXElems,totalYElems,totalZElems),
                           Teuchos::tuple(nx,ny,nz),
                           block,nodes,mesh);
         }
      }
   }
}

void CubeTetMeshFactory::buildTetsOnHex(const Teuchos::Tuple<int,3> & meshDesc,
                                        const Teuchos::Tuple<int,3> & element,
                                        stk::mesh::Part * block,
                                        const std::vector<stk::mesh::EntityId> & h_nodes,
                                        STK_Interface & mesh) const
{
   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setShowProcRank(true);
   out.setOutputToRootOnly(-1);

   int totalXElems = meshDesc[0]; int totalYElems = meshDesc[1]; int totalZElems = meshDesc[2];
   int nx = element[0]; int ny = element[1]; int nz = element[2];

   stk::mesh::EntityId hex_id = totalXElems*totalYElems*nz+totalXElems*ny+nx+1;
   stk::mesh::EntityId gid_0 = 12*(hex_id-1)+1;
   std::vector<stk::mesh::EntityId> nodes(4);

   // add centroid node
   stk::mesh::EntityId centroid = 0;
   {
      stk::mesh::EntityId largestNode = (totalXElems+1)*(totalYElems+1)*(totalZElems+1);
      centroid = hex_id+largestNode;

      // compute average of coordinates
      std::vector<double> coord(3,0.0);
      for(std::size_t i=0;i<h_nodes.size();i++) {
         const double * node_coord = mesh.getNodeCoordinates(h_nodes[i]);
         coord[0] += node_coord[0];
         coord[1] += node_coord[1];
         coord[2] += node_coord[2];
      }
      coord[0] /= 8.0;
      coord[1] /= 8.0;
      coord[2] /= 8.0;

      mesh.addNode(centroid,coord);
   }

   //
   int idSet[][3] = { { 0, 1, 2}, // back
                      { 0, 2, 3},
                      { 0, 5, 1}, // bottom
                      { 0, 4, 5},
                      { 0, 7, 4}, // left
                      { 0, 3, 7},
                      { 6, 1, 5}, // right
                      { 6, 2, 1},
                      { 6, 3, 2}, // top
                      { 6, 7, 3},
                      { 6, 4, 7}, // front
                      { 6, 5, 4} };

   for(int i=0;i<12;i++) {
      nodes[0] = h_nodes[idSet[i][0]];
      nodes[1] = h_nodes[idSet[i][1]];
      nodes[2] = h_nodes[idSet[i][2]];
      nodes[3] = centroid;

      // add element to mesh
      mesh.addElement(rcp(new ElementDescriptor(gid_0+i,nodes)),block);
   }
}

std::pair<int,int> CubeTetMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int /* rank */) const
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

std::pair<int,int> CubeTetMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int /* rank */) const
{
   // int start = yBlock*nYElems_;
   // return std::make_pair(start,nYElems_);

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

std::pair<int,int> CubeTetMeshFactory::determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int /* rank */) const
{
   // int start = zBlock*nZElems_;
   // return std::make_pair(start,nZElems_);
   std::size_t zProcLoc = procTuple_[2];
   unsigned int minElements = nZElems_/size;
   unsigned int extra = nZElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   int nume=0, start=0;
   if(zProcLoc<extra) {
      nume  = minElements+1;
      start = zProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(zProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nZElems_*zBlock,nume);
}

void CubeTetMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();
   const stk::mesh::EntityRank side_rank = mesh.getSideRank();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;
   std::size_t totalZElems = nZElems_*zBlocks_;

   // get all part vectors
   stk::mesh::Part * left = mesh.getSideset("left");
   stk::mesh::Part * right = mesh.getSideset("right");
   stk::mesh::Part * top = mesh.getSideset("top");
   stk::mesh::Part * bottom = mesh.getSideset("bottom");
   stk::mesh::Part * front = mesh.getSideset("front");
   stk::mesh::Part * back = mesh.getSideset("back");

   std::vector<stk::mesh::Entity> localElmts;
   mesh.getMyElements(localElmts);

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk::mesh::Entity>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity element = (*itr);
      stk::mesh::EntityId gid = mesh.elementGlobalId(element);

      // get hex global id
      stk::mesh::EntityId h_gid = (gid-1)/12+1;
      stk::mesh::EntityId t_offset = gid - (12*(h_gid-1)+1);

      std::size_t nx,ny,nz;
      nz = (h_gid-1) / (totalXElems*totalYElems);
      h_gid = (h_gid-1)-nz*(totalXElems*totalYElems);
      ny = h_gid / totalXElems;
      nx = h_gid-ny*totalXElems;

      if(nz==0 && (t_offset==0 || t_offset==1)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the back
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,back);
      }
      if(nz+1==totalZElems && (t_offset==10 || t_offset==11)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the front
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,front);
      }

      if(ny==0 && (t_offset==2 || t_offset==3)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the bottom
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,bottom);
      }
      if(ny+1==totalYElems && (t_offset==8 || t_offset==9)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the top
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,top);
      }

      if(nx==0 && (t_offset==4 || t_offset==5)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the left
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,left);
      }
      if(nx+1==totalXElems && (t_offset==6 || t_offset==7)) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the right
         if(mesh.entityOwnerRank(side)==machRank_)
            mesh.addEntityToSideset(side,right);
      }
   }

   mesh.endModification();
}

void CubeTetMeshFactory::addNodeSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // get all part vectors
   stk::mesh::Part * origin = mesh.getNodeset("origin");

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   if(machRank_==0)
   {
      // add zero node to origin node set
      stk::mesh::Entity node = bulkData->get_entity(mesh.getNodeRank(),1);
      mesh.addEntityToNodeset(node,origin);
   }

   mesh.endModification();
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void CubeTetMeshFactory::addEdgeBlocks(STK_Interface & mesh) const
{
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   stk::mesh::Part * edge_block = mesh.getEdgeBlock(edgeBlockName_);

   stk::mesh::Selector owned_block = metaData->locally_owned_part();

   std::vector<stk::mesh::Entity> edges;
   bulkData->get_entities(mesh.getEdgeRank(), owned_block, edges);
   mesh.addEntitiesToEdgeBlock(edges, edge_block);
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void CubeTetMeshFactory::addFaceBlocks(STK_Interface & mesh) const
{
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   stk::mesh::Part * face_block = mesh.getFaceBlock(faceBlockName_);

   stk::mesh::Selector owned_block = metaData->locally_owned_part();

   std::vector<stk::mesh::Entity> faces;
   bulkData->get_entities(mesh.getFaceRank(), owned_block, faces);
   mesh.addEntitiesToFaceBlock(faces, face_block);
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,3> CubeTetMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0,k=0;

   k = procRank/(xProcs_*yProcs_); procRank = procRank % (xProcs_*yProcs_);
   j = procRank/xProcs_;           procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j,k);
}

} // end panzer_stk
