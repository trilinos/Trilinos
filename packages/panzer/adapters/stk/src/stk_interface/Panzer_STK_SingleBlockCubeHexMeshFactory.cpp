#include <Panzer_STK_SingleBlockCubeHexMeshFactory.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

SingleBlockCubeHexMeshFactory::SingleBlockCubeHexMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
SingleBlockCubeHexMeshFactory::~SingleBlockCubeHexMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> SingleBlockCubeHexMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> SingleBlockCubeHexMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   TEST_FOR_EXCEPTION(machSize_!=xProcs_*yProcs_*zProcs_,std::logic_error,
                      "Cannot build SingleBlockCubeHexMeshFactory, the product of \"X Procs\", \"Y Procs\", and \"Z Procs\""
                      " must equal the number of processors.");

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);
 
   mesh->addPeriodicBCs(periodicBCVec_);

   return mesh;
}

void SingleBlockCubeHexMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // compute tuple to construct offsets
   procTuple_ = procRankToProcTuple(machRank_);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges and faces
   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   // now that edges are built, sidets can be added
   addSideSets(mesh);
}

//! From ParameterListAcceptor
void SingleBlockCubeHexMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   x0_ = paramList->get<double>("X0"); 
   y0_ = paramList->get<double>("Y0"); 
   z0_ = paramList->get<double>("Z0"); 

   xf_ = paramList->get<double>("Xf"); 
   yf_ = paramList->get<double>("Yf"); 
   zf_ = paramList->get<double>("Zf"); 

   xProcs_ = paramList->get<int>("X Procs");
   yProcs_ = paramList->get<int>("Y Procs");
   zProcs_ = paramList->get<int>("Z Procs");

   nXElems_ = paramList->get<int>("X Elements");
   nYElems_ = paramList->get<int>("Y Elements");
   nZElems_ = paramList->get<int>("Z Elements");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> SingleBlockCubeHexMeshFactory::getValidParameters() const
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

      defaultParams->set<int>("X Procs",1);
      defaultParams->set<int>("Y Procs",1);
      defaultParams->set<int>("Z Procs",1);

      defaultParams->set<int>("X Elements",5);
      defaultParams->set<int>("Y Elements",5);
      defaultParams->set<int>("Z Elements",5);

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void SingleBlockCubeHexMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);
}

void SingleBlockCubeHexMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Hexahedron<8> HexTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

   // add element blocks
   mesh.addElementBlock(std::string("eblock-0_0"),ctd);

   // add sidesets 
   mesh.addSideset("left",side_ctd);
   mesh.addSideset("right",side_ctd);
   mesh.addSideset("top",side_ctd);
   mesh.addSideset("bottom",side_ctd);
   mesh.addSideset("front",side_ctd);
   mesh.addSideset("back",side_ctd);
}

void SingleBlockCubeHexMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();
      // build each block
      buildBlock(parallelMach,mesh);
   mesh.endModification();
}

void SingleBlockCubeHexMeshFactory::buildBlock(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<int,int> sizeAndStartX = determineXElemSizeAndStart();
   std::pair<int,int> sizeAndStartY = determineYElemSizeAndStart();
   std::pair<int,int> sizeAndStartZ = determineZElemSizeAndStart();

   int myXElems_start = sizeAndStartX.first;
   int myXElems_end  = myXElems_start+sizeAndStartX.second;
   int myYElems_start = sizeAndStartY.first;
   int myYElems_end  = myYElems_start+sizeAndStartY.second;
   int myZElems_start = sizeAndStartZ.first;
   int myZElems_end  = myZElems_start+sizeAndStartZ.second;

   int totalXElems = nXElems_;
   int totalYElems = nYElems_;
   int totalZElems = nZElems_;

   double deltaX = (xf_-x0_)/double(totalXElems);
   double deltaY = (yf_-y0_)/double(totalYElems);
   double deltaZ = (zf_-z0_)/double(totalZElems);
 
   std::vector<double> coord(3,0.0);

   // build the nodes
   for(int nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = double(nx)*deltaX+x0_;
      for(int ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = double(ny)*deltaY+y0_;
         for(int nz=myZElems_start;nz<myZElems_end+1;++nz) {
            coord[2] = double(nz)*deltaZ+z0_;

            mesh.addNode(nz*(totalYElems+1)*(totalXElems+1)+ny*(totalXElems+1)+nx+1,coord);
         }
      }
   }

   stk::mesh::Part * block = mesh.getElementBlockPart("eblock-0_0");

   // build the elements
   for(int nx=myXElems_start;nx<myXElems_end;++nx) {
      for(int ny=myYElems_start;ny<myYElems_end;++ny) {
         for(int nz=myZElems_start;nz<myZElems_end;++nz) {
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

std::pair<int,int> SingleBlockCubeHexMeshFactory::determineXElemSizeAndStart() const
{
   std::size_t xProcLoc = procTuple_[0];
   std::size_t minElements = nXElems_/xProcs_;
   std::size_t extra = nXElems_ - minElements*xProcs_;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting X index and number of elements
   std::size_t nume=0, start=0;
   if(xProcLoc<extra) {
      nume  = minElements+1;
      start = xProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(xProcLoc-extra)*minElements;
   }

   return std::make_pair(start,nume);
}

std::pair<int,int> SingleBlockCubeHexMeshFactory::determineYElemSizeAndStart() const
{
   std::size_t yProcLoc = procTuple_[1];
   std::size_t minElements = nYElems_/yProcs_;
   std::size_t extra = nYElems_ - minElements*yProcs_;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting index and number of elements
   std::size_t nume=0, start=0;
   if(yProcLoc<extra) {
      nume  = minElements+1;
      start = yProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(yProcLoc-extra)*minElements;
   }

   return std::make_pair(start,nume);
}

std::pair<int,int> SingleBlockCubeHexMeshFactory::determineZElemSizeAndStart() const
{
   std::size_t zProcLoc = procTuple_[2];
   std::size_t minElements = nZElems_/zProcs_;
   std::size_t extra = nZElems_ - minElements*zProcs_;

   TEUCHOS_ASSERT(minElements>0);

   // first "extra" elements get an extra column of elements
   // this determines the starting index and number of elements
   std::size_t nume=0, start=0;
   if(zProcLoc<extra) {
      nume  = minElements+1;
      start = zProcLoc*(minElements+1);
   }
   else {
      nume  = minElements;
      start = extra*(minElements+1)+(zProcLoc-extra)*minElements;
   }

   return std::make_pair(start,nume);
}

// for use with addSideSets
const stk::mesh::Relation * SingleBlockCubeHexMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}

void SingleBlockCubeHexMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_;
   std::size_t totalYElems = nYElems_;
   std::size_t totalZElems = nZElems_;

   // get all part vectors
   stk::mesh::Part * left = mesh.getSideset("left");
   stk::mesh::Part * right = mesh.getSideset("right");
   stk::mesh::Part * top = mesh.getSideset("top");
   stk::mesh::Part * bottom = mesh.getSideset("bottom");
   stk::mesh::Part * front = mesh.getSideset("front");
   stk::mesh::Part * back = mesh.getSideset("back");

   std::vector<stk::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   // gid = totalXElems*totalYElems*nz+totalXElems*ny+nx+1

   // loop over elements adding sides to sidesets
   std::vector<stk::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity * element = (*itr);
      stk::mesh::EntityId gid = element->identifier();      
      stk::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

      std::size_t nx,ny,nz;
      nz = (gid-1) / (totalXElems*totalYElems);
      gid = (gid-1)-nz*(totalXElems*totalYElems);
      ny = gid / totalXElems;
      nx = gid-ny*totalXElems;

      if(nz==0) {
         stk::mesh::Entity * side = getRelationByID(4,relations)->entity();

         // on the back
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,back);
      }
      if(nz+1==totalZElems) {
         stk::mesh::Entity * side = getRelationByID(5,relations)->entity();

         // on the front
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,front);
      }

      if(ny==0) {
         stk::mesh::Entity * side = getRelationByID(0,relations)->entity();

         // on the bottom 
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,bottom);
      }
      if(ny+1==totalYElems) {
         stk::mesh::Entity * side = getRelationByID(2,relations)->entity();

         // on the top
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,top);
      }

      if(nx==0) {
         stk::mesh::Entity * side = getRelationByID(3,relations)->entity();

         // on the left
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,left);
      }
      if(nx+1==totalXElems) {
         stk::mesh::Entity * side = getRelationByID(1,relations)->entity();

         // on the right
         if(side->owner_rank()==machRank_)
            mesh.addEntityToSideset(*side,right);
      }
   }

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,3> SingleBlockCubeHexMeshFactory::procRankToProcTuple(std::size_t procRank) const
{

   std::size_t i=0,j=0,k=0;

   k = procRank/(xProcs_*yProcs_); procRank = procRank % (xProcs_*yProcs_);
   j = procRank/xProcs_;           procRank = procRank % xProcs_;
   i = procRank;

   return Teuchos::tuple(i,j,k);
}

} // end panzer_stk
