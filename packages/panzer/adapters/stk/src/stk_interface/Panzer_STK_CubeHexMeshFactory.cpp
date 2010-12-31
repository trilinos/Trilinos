#include <Panzer_STK_CubeHexMeshFactory.hpp>

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
   RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);
 
  return mesh;
}

void CubeHexMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   buildElements(parallelMach,mesh);

   // finish up the edges
   mesh.buildSubcells();

   // now that edges are built, sidets can be added
   // addSideSets(mesh);
}

//! From ParameterListAcceptor
void CubeHexMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters());

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

   nXElems_ = paramList->get<int>("X Elements");
   nYElems_ = paramList->get<int>("Y Elements");
   nZElems_ = paramList->get<int>("Z Elements");
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

      defaultParams->set<int>("X Elements",5);
      defaultParams->set<int>("Y Elements",5);
      defaultParams->set<int>("Z Elements",5);
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

void CubeHexMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Hexahedron<8> HexTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();

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
   mesh.addSideset("left");
   mesh.addSideset("right");
   mesh.addSideset("top");
   mesh.addSideset("bottom");
   mesh.addSideset("front");
   mesh.addSideset("back");
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

void CubeHexMeshFactory::buildBlock(stk::ParallelMachine parallelMach,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<int,int> sizeAndStartX = determineXElemSizeAndStart(xBlock,machSize_,machRank_);
   std::pair<int,int> sizeAndStartY = determineYElemSizeAndStart(yBlock,machSize_,machRank_);
   std::pair<int,int> sizeAndStartZ = determineZElemSizeAndStart(zBlock,machSize_,machRank_);

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
      coord[0] = double(nx)*deltaX+x0_;
      for(int ny=myYElems_start;ny<myYElems_end+1;++ny) {
         coord[1] = double(ny)*deltaY+y0_;
         for(int nz=myZElems_start;nz<myZElems_end+1;++nz) {
            coord[2] = double(nz)*deltaZ+z0_;

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
         for(int nz=myZElems_start;nz<myYElems_end;++nz) {
            stk::mesh::EntityId gid = totalXElems*totalZElems*nz+totalXElems*ny+nx+1;
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

std::pair<int,int> CubeHexMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const
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

std::pair<int,int> CubeHexMeshFactory::determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const
{
   int start = yBlock*nYElems_;

   return std::make_pair(start,nYElems_);
}

std::pair<int,int> CubeHexMeshFactory::determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int rank) const
{
   int start = zBlock*nZElems_;

   return std::make_pair(start,nZElems_);
}

// for use with addSideSets
const stk::mesh::Relation * CubeHexMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}

void CubeHexMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   int totalXElems = nXElems_*xBlocks_;
   int totalYElems = nYElems_*yBlocks_;

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
      stk::mesh::PairIterRelation relations = element->relations(stk::mesh::Edge);
      
      // vertical boundaries
      if(gid % totalXElems==0) { 
         stk::mesh::Entity * edge = getRelationByID(1,relations)->entity();

         // on the right
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,right);
      }
      else if((gid-1) % totalXElems ==0) {
         stk::mesh::Entity * edge = getRelationByID(3,relations)->entity();

         // on the left
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,left);
      }

      // horizontal boundaries
      if(gid <= (std::size_t) totalXElems) {
         stk::mesh::Entity * edge = getRelationByID(0,relations)->entity();

         // on the bottom
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,bottom);
      }
      else if(gid > (std::size_t) totalXElems*totalYElems-totalXElems ) {
         stk::mesh::Entity * edge = getRelationByID(2,relations)->entity();

         // on the top
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,top);
      }
   }

   mesh.endModification();
}

} // end panzer_stk
