// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_STK_CubeHexMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <PanzerAdaptersSTK_config.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

CubeHexMeshFactory::CubeHexMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
CubeHexMeshFactory::~CubeHexMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> CubeHexMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
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

Teuchos::RCP<STK_Interface> CubeHexMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::CubeHexMeshFactory::buildUncomittedMesh()");

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
                      "Cannot build CubeHexMeshFactory, the product of \"X Procs\", \"Y Procs\", and \"Z Procs\""
                      " must equal the number of processors.");
   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void CubeHexMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
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
   if(createEdgeBlocks_) {
      mesh.buildLocalEdgeIDs();
   }
   if(createFaceBlocks_) {
      mesh.buildLocalFaceIDs();
   }

   mesh.beginModification();

   // now that edges are built, side and node sets can be added
   addSideSets(mesh);
   addNodeSets(mesh);
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

   buildInterfaceSidesets_ = paramList->get<bool>("Build Interface Sidesets");

   buildSubcells_ = paramList->get<bool>("Build Subcells");

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");
   createFaceBlocks_ = paramList->get<bool>("Create Face Blocks");
   if (not buildSubcells_ && createEdgeBlocks_) {
      Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
      out.setOutputToRootOnly(0);
      out.setShowProcRank(true);

      out << "CubeHexMesh: NOT creating edge blocks because building sub cells disabled" << std::endl;
      createEdgeBlocks_ = false;
   }
   if (not buildSubcells_ && createFaceBlocks_) {
      Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
      out.setOutputToRootOnly(0);
      out.setShowProcRank(true);

      out << "CubeHexMesh: NOT creating face blocks because building sub cells disabled" << std::endl;
      createFaceBlocks_ = false;
   }

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
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

      defaultParams->set<bool>("Build Interface Sidesets",false);

      defaultParams->set<bool>("Build Subcells",true);

      // default to false for backward compatibility
      defaultParams->set<bool>("Create Edge Blocks",false,"Create edge blocks in the mesh");
      defaultParams->set<bool>("Create Face Blocks",false,"Create face blocks in the mesh");

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

   /* This is a hex mesh factory so all elements in all element blocks
    * will be hex8.  This means that all the edges will be line2 and
    * all the faces will be quad4.  The edge and face block names are
    * hard coded to reflect this.
    */
   edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;
   faceBlockName_ = "quad_4_"+panzer_stk::STK_Interface::faceBlockString;
}

void CubeHexMeshFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
   typedef shards::Hexahedron<8> HexTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();
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

   if(buildInterfaceSidesets_) {
     for(int bx=1;bx<xBlocks_;bx++) {
       std::stringstream ss;
       ss << "vertical_" << bx-1;
       mesh.addSideset(ss.str(),side_ctd);
     }
     for(int by=1;by<yBlocks_;by++) {
       std::stringstream ss;
       ss << "horizontal_" << by-1;
       mesh.addSideset(ss.str(),side_ctd);
     }
     for(int bz=1;bz<zBlocks_;bz++) {
       std::stringstream ss;
       ss << "transverse_" << bz-1;
       mesh.addSideset(ss.str(),side_ctd);
     }
   }

   // add nodesets
   mesh.addNodeset("origin");
}

void CubeHexMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
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

void CubeHexMeshFactory::buildBlock(stk::ParallelMachine /* parallelMach */,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> sizeAndStartX = determineXElemSizeAndStart(xBlock,xProcs_,machRank_);
   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> sizeAndStartY = determineYElemSizeAndStart(yBlock,yProcs_,machRank_);
   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> sizeAndStartZ = determineZElemSizeAndStart(zBlock,zProcs_,machRank_);

   panzer::GlobalOrdinal myXElems_start = sizeAndStartX.first;
   panzer::GlobalOrdinal myXElems_end  = myXElems_start+sizeAndStartX.second;
   panzer::GlobalOrdinal myYElems_start = sizeAndStartY.first;
   panzer::GlobalOrdinal myYElems_end  = myYElems_start+sizeAndStartY.second;
   panzer::GlobalOrdinal myZElems_start = sizeAndStartZ.first;
   panzer::GlobalOrdinal myZElems_end  = myZElems_start+sizeAndStartZ.second;

   panzer::GlobalOrdinal totalXElems = nXElems_*xBlocks_;
   panzer::GlobalOrdinal totalYElems = nYElems_*yBlocks_;
   panzer::GlobalOrdinal totalZElems = nZElems_*zBlocks_;

   double deltaX = (xf_-x0_)/double(totalXElems);
   double deltaY = (yf_-y0_)/double(totalYElems);
   double deltaZ = (zf_-z0_)/double(totalZElems);

   std::vector<double> coord(3,0.0);

   // build the nodes
   for(panzer::GlobalOrdinal nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = this->getMeshCoord(nx, deltaX, x0_);
      for(panzer::GlobalOrdinal ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = this->getMeshCoord(ny, deltaY, y0_);
         for(panzer::GlobalOrdinal nz=myZElems_start;nz<myZElems_end+1;++nz) {
            coord[2] = this->getMeshCoord(nz, deltaZ, z0_);

            mesh.addNode(nz*(totalYElems+1)*(totalXElems+1)+ny*(totalXElems+1)+nx+1,coord);
         }
      }
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock << "_" << yBlock << "_" << zBlock;
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(panzer::GlobalOrdinal nx=myXElems_start;nx<myXElems_end;++nx) {
      for(panzer::GlobalOrdinal ny=myYElems_start;ny<myYElems_end;++ny) {
         for(panzer::GlobalOrdinal nz=myZElems_start;nz<myZElems_end;++nz) {
            stk::mesh::EntityId gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1;
            std::vector<stk::mesh::EntityId> nodes(8);
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

std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> CubeHexMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int /* rank */) const
{
   std::size_t xProcLoc = procTuple_[0];
   panzer::GlobalOrdinal minElements = nXElems_/size;
   panzer::GlobalOrdinal extra = nXElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::GlobalOrdinal nume=0, start=0;
   if(panzer::GlobalOrdinal(xProcLoc)<extra) {
      nume  = minElements+1;
      start = xProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(xProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nXElems_*xBlock,nume);
}

std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> CubeHexMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int /* rank */) const
{
   // int start = yBlock*nYElems_;
   // return std::make_pair(start,nYElems_);

   std::size_t yProcLoc = procTuple_[1];
   panzer::GlobalOrdinal minElements = nYElems_/size;
   panzer::GlobalOrdinal extra = nYElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::GlobalOrdinal nume=0, start=0;
   if(panzer::GlobalOrdinal(yProcLoc)<extra) {
      nume  = minElements+1;
      start = yProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(yProcLoc-extra)*minElements;
   }

   return std::make_pair(start+nYElems_*yBlock,nume);
}

std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> CubeHexMeshFactory::determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int /* rank */) const
{
   // int start = zBlock*nZElems_;
   // return std::make_pair(start,nZElems_);
   std::size_t zProcLoc = procTuple_[2];
   panzer::GlobalOrdinal minElements = nZElems_/size;
   panzer::GlobalOrdinal extra = nZElems_ - minElements*size;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   panzer::GlobalOrdinal nume=0, start=0;
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

// this adds side entities only (does not inject them into side sets)
void CubeHexMeshFactory::addSides(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;
   std::size_t totalZElems = nZElems_*zBlocks_;

   std::vector<stk::mesh::Entity> localElmts;
   mesh.getMyElements(localElmts);

   stk::mesh::EntityId offset[6];
   offset[0] = 0;
   offset[1] = offset[0] + totalXElems*totalZElems;
   offset[2] = offset[1] + totalYElems*totalZElems;
   offset[3] = offset[2] + totalXElems*totalZElems;
   offset[4] = offset[3] + totalYElems*totalZElems;
   offset[5] = offset[4] + totalXElems*totalYElems;

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk::mesh::Entity>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity element = (*itr);
      stk::mesh::EntityId gid = mesh.elementGlobalId(element);

      std::size_t nx,ny,nz;
      nz = (gid-1) / (totalXElems*totalYElems);
      gid = (gid-1)-nz*(totalXElems*totalYElems);
      ny = gid / totalXElems;
      nx = gid-ny*totalXElems;

      std::vector<stk::mesh::Part*> parts;

      if(nz==0) {
         // on the back
         mesh.getBulkData()->declare_element_side(element, 4, parts);
      }
      if(nz+1==totalZElems) {
         // on the front
         mesh.getBulkData()->declare_element_side(element, 5, parts);
      }

      if(ny==0) {
         // on the bottom
         mesh.getBulkData()->declare_element_side(element, 0, parts);
      }
      if(ny+1==totalYElems) {
         // on the top
         mesh.getBulkData()->declare_element_side(element, 2, parts);
      }

      if(nx==0) {
         // on the left
         mesh.getBulkData()->declare_element_side(element, 3, parts);
      }
      if(nx+1==totalXElems) {
         // on the right
         mesh.getBulkData()->declare_element_side(element, 1, parts);
      }
   }

   mesh.endModification();
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void CubeHexMeshFactory::addSideSets(STK_Interface & mesh) const
{
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

   std::vector<stk::mesh::Part*> vertical;
   std::vector<stk::mesh::Part*> horizontal;
   std::vector<stk::mesh::Part*> transverse;

   if(buildInterfaceSidesets_) {
     for(int bx=1;bx<xBlocks_;bx++) {
       std::stringstream ss;
       ss << "vertical_" << bx-1;
       vertical.push_back(mesh.getSideset(ss.str()));
     }
     for(int by=1;by<yBlocks_;by++) {
       std::stringstream ss;
       ss << "horizontal_" << by-1;
       horizontal.push_back(mesh.getSideset(ss.str()));
     }
     for(int bz=1;bz<zBlocks_;bz++) {
       std::stringstream ss;
       ss << "transverse_" << bz-1;
       transverse.push_back(mesh.getSideset(ss.str()));
     }
   }

   std::vector<stk::mesh::Entity> localElmts;
   mesh.getMyElements(localElmts);

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk::mesh::Entity>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity element = (*itr);
      stk::mesh::EntityId gid = mesh.elementGlobalId(element);

      std::size_t nx,ny,nz;
      nz = (gid-1) / (totalXElems*totalYElems);
      gid = (gid-1)-nz*(totalXElems*totalYElems);
      ny = gid / totalXElems;
      nx = gid-ny*totalXElems;

      if(nz % nZElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 4);

         // on the back
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(nz==0) {
	     mesh.addEntityToSideset(side,back);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = nz/nZElems_-1;
	       mesh.addEntityToSideset(side,transverse[index]);
	     }
	   }
	 }
      }
      if((nz+1) % nZElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 5);

         // on the front
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(nz+1==totalZElems) {
	     mesh.addEntityToSideset(side,front);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = (nz+1)/nZElems_-1;
	       mesh.addEntityToSideset(side,transverse[index]);
	     }
	   }
	 }
      }

      if(ny % nYElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 0);

         // on the bottom
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(ny==0) {
	     mesh.addEntityToSideset(side,bottom);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = ny/nYElems_-1;
	       mesh.addEntityToSideset(side,horizontal[index]);
	     }
	   }
	 }
      }
      if((ny+1) % nYElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 2);

         // on the top
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(ny+1==totalYElems) {
	     mesh.addEntityToSideset(side,top);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = (ny+1)/nYElems_-1;
	       mesh.addEntityToSideset(side,horizontal[index]);
	     }
	   }
	 }
      }

      if(nx % nXElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 3);

         // on the left
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(nx==0) {
	     mesh.addEntityToSideset(side,left);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = nx/nXElems_-1;
	       mesh.addEntityToSideset(side,vertical[index]);
	     }
	   }
	 }
      }
      if((nx+1) % nXElems_==0) {
         stk::mesh::Entity side = mesh.findConnectivityById(element, side_rank, 1);

         // on the right
         if(mesh.entityOwnerRank(side)==machRank_) {
	   if(nx+1==totalXElems) {
	     mesh.addEntityToSideset(side,right);
	   } else {
	     if(buildInterfaceSidesets_) {
	       int index = (nx+1)/nXElems_-1;
	       mesh.addEntityToSideset(side,vertical[index]);
	     }
	   }
	 }
      }
   }
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void CubeHexMeshFactory::addNodeSets(STK_Interface & mesh) const
{
   // get all part vectors
   stk::mesh::Part * origin = mesh.getNodeset("origin");

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   if(machRank_==0)
   {
      // add zero node to origin node set
      stk::mesh::Entity node = bulkData->get_entity(mesh.getNodeRank(),1);
      mesh.addEntityToNodeset(node,origin);
   }
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void CubeHexMeshFactory::addEdgeBlocks(STK_Interface & mesh) const
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
void CubeHexMeshFactory::addFaceBlocks(STK_Interface & mesh) const
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
Teuchos::Tuple<std::size_t,3> CubeHexMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0,k=0;

   k = procRank/(xProcs_*yProcs_); procRank = procRank % (xProcs_*yProcs_);
   j = procRank/xProcs_;           procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j,k);
}

} // end panzer_stk
