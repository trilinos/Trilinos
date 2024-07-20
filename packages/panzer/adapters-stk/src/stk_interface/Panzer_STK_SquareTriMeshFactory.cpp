// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Panzer_STK_SquareTriMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <PanzerAdaptersSTK_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

SquareTriMeshFactory::SquareTriMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
SquareTriMeshFactory::~SquareTriMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> SquareTriMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareTriMeshFactory::buildMesh()");

   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> SquareTriMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareTriMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(2));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   if (xProcs_ == -1 && yProcs_ == -1) {
     // copied from galeri
     xProcs_ = yProcs_ = Teuchos::as<int>(pow(Teuchos::as<double>(machSize_), 0.5));

     if (xProcs_ * yProcs_ != Teuchos::as<int>(machSize_))  {
       // Simple method to find a set of processor assignments
       xProcs_ = yProcs_ = 1;

       // This means that this works correctly up to about maxFactor^2
       // processors.
       const int maxFactor = 100;

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
           if      (xProcs_ <= yProcs_) xProcs_ = xProcs_*jj;
           else                         yProcs_ = yProcs_*jj;
           factors[jj]--;
         }
       }
     }

   } else if(xProcs_==-1) {
      // default x only decomposition
      xProcs_ = machSize_;
      yProcs_ = 1;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(int(machSize_)!=xProcs_*yProcs_,std::logic_error,
                      "Cannot build SquareTriMeshFactory, the product of \"X Procs\" and \"Y Procs\""
                      " must equal the number of processors.");
   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void SquareTriMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::SquareTriMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges
   mesh.buildSubcells();
   mesh.buildLocalElementIDs();
   if(createEdgeBlocks_) {
      mesh.buildLocalEdgeIDs();
   }

   // now that edges are built, sidets can be added
   addSideSets(mesh);

   // add nodesets
   addNodeSets(mesh);

   if(createEdgeBlocks_) {
      addEdgeBlocks(mesh);
   }

   // calls Stk_MeshFactory::rebalance
   this->rebalance(mesh);
}

//! From ParameterListAcceptor
void SquareTriMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
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

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> SquareTriMeshFactory::getValidParameters() const
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

      // default to false for backward compatibility
      defaultParams->set<bool>("Create Edge Blocks",false,"Create edge blocks in the mesh");

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void SquareTriMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);

   /* This is a tri mesh factory so all elements in all element blocks
    * will be tri3.  This means that all the edges will be line2.
    * The edge block name is hard coded to reflect this.
    */
   edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;
}

void SquareTriMeshFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
   typedef shards::Triangle<> TriTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<TriTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);
   const CellTopologyData * edge_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

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
            if(createEdgeBlocks_) {
               mesh.addEdgeBlock("eblock"+ebPostfix.str(),
                                 edgeBlockName_,
                                 edge_ctd);
            }
         }

      }
   }

   // add sidesets
   mesh.addSideset("left",side_ctd);
   mesh.addSideset("right",side_ctd);
   mesh.addSideset("top",side_ctd);
   mesh.addSideset("bottom",side_ctd);

   // add nodesets
   mesh.addNodeset("origin");
}

void SquareTriMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
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

void SquareTriMeshFactory::buildBlock(stk::ParallelMachine /* parallelMach */, int xBlock, int yBlock, STK_Interface& mesh) const
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
      coord[0] = this->getMeshCoord(nx, deltaX, x0_);
      for(int ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = this->getMeshCoord(ny, deltaY, y0_);

         mesh.addNode(ny*(totalXElems+1)+nx+1,coord);
      }
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock << "_" << yBlock;
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(int nx=myXElems_start;nx<myXElems_end;++nx) {
      for(int ny=myYElems_start;ny<myYElems_end;++ny) {
         stk::mesh::EntityId gid_a = 2*(totalXElems*ny+nx+1)-1;
         stk::mesh::EntityId gid_b = gid_a+1;
         std::vector<stk::mesh::EntityId> nodes(3);
         stk::mesh::EntityId sw,se,ne,nw;
         sw = nx+1+ny*(totalXElems+1);
         se = sw+1;
         ne = se+(totalXElems+1);
         nw = ne-1;

         nodes[0] = sw;
         nodes[1] = se;
         nodes[2] = ne;
         mesh.addElement(rcp(new ElementDescriptor(gid_a,nodes)),block);

         nodes[0] = sw;
         nodes[1] = ne;
         nodes[2] = nw;
         mesh.addElement(rcp(new ElementDescriptor(gid_b,nodes)),block);
      }
   }
}

std::pair<int,int> SquareTriMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int /* rank */) const
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

std::pair<int,int> SquareTriMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int /* rank */) const
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

void SquareTriMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;
   std::size_t totalYElems = nYElems_*yBlocks_;

   // get all part vectors
   stk::mesh::Part * left = mesh.getSideset("left");
   stk::mesh::Part * right = mesh.getSideset("right");
   stk::mesh::Part * top = mesh.getSideset("top");
   stk::mesh::Part * bottom = mesh.getSideset("bottom");

   std::vector<stk::mesh::Entity> localElmts;
   mesh.getMyElements(localElmts);

   // loop over elements adding edges to sidesets
   std::vector<stk::mesh::Entity>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity element = (*itr);
      stk::mesh::EntityId gid = mesh.elementGlobalId(element);

      bool lower = (gid%2 != 0);
      std::size_t block = lower ? (gid+1)/2 : gid/2;
      std::size_t nx,ny;
      ny = (block-1) / totalXElems;
      nx = block-ny*totalXElems-1;

      // vertical boundaries
      ///////////////////////////////////////////

      if(nx+1==totalXElems && lower) {
         stk::mesh::Entity edge = mesh.findConnectivityById(element, stk::topology::EDGE_RANK, 1);

         // on the right
         if(mesh.entityOwnerRank(edge)==machRank_)
            mesh.addEntityToSideset(edge,right);
      }

      if(nx==0 && !lower) {
         stk::mesh::Entity edge = mesh.findConnectivityById(element, stk::topology::EDGE_RANK, 2);

         // on the left
         if(mesh.entityOwnerRank(edge)==machRank_)
            mesh.addEntityToSideset(edge,left);
      }

      // horizontal boundaries
      ///////////////////////////////////////////

      if(ny==0 && lower) {
         stk::mesh::Entity edge = mesh.findConnectivityById(element, stk::topology::EDGE_RANK, 0);

         // on the bottom
         if(mesh.entityOwnerRank(edge)==machRank_)
            mesh.addEntityToSideset(edge,bottom);
      }

      if(ny+1==totalYElems && !lower) {
         stk::mesh::Entity edge = mesh.findConnectivityById(element, stk::topology::EDGE_RANK, 1);

         // on the top
         if(mesh.entityOwnerRank(edge)==machRank_)
            mesh.addEntityToSideset(edge,top);
      }
   }

   mesh.endModification();
}

void SquareTriMeshFactory::addNodeSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // get all part vectors
   stk::mesh::Part * origin = mesh.getNodeset("origin");

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   if(machRank_==0)
   {
      stk::mesh::Entity node = bulkData->get_entity(mesh.getNodeRank(),1);

      // add zero node to origin node set
      mesh.addEntityToNodeset(node,origin);
   }

   mesh.endModification();
}

void SquareTriMeshFactory::addEdgeBlocks(STK_Interface & mesh) const
{
   mesh.beginModification();

   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   stk::mesh::Part * edge_block = mesh.getEdgeBlock(edgeBlockName_);

   stk::mesh::Selector owned_block = metaData->locally_owned_part();

   std::vector<stk::mesh::Entity> edges;
   bulkData->get_entities(mesh.getEdgeRank(), owned_block, edges);
   mesh.addEntitiesToEdgeBlock(edges, edge_block);

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,2> SquareTriMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0;

   j = procRank/xProcs_;
   procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j);
}

} // end panzer_stk
