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

#include <Panzer_STK_LineMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Panzer_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

LineMeshFactory::LineMeshFactory()
{
   initializeWithDefaults();
}

//! Destructor
LineMeshFactory::~LineMeshFactory()
{
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> LineMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::LineMeshFactory::buildMesh()");

   // build all meta data
   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

Teuchos::RCP<STK_Interface> LineMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::LineMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(1));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   procTuple_ = procRankToProcTuple(machRank_);

   // build meta information: blocks and side set setups
   buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
 
   return mesh;
}

void LineMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::LineMeshFactory::completeMeshConstruction()");

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
void LineMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   x0_ = paramList->get<double>("X0"); 
   xf_ = paramList->get<double>("Xf"); 
   xBlocks_ = paramList->get<int>("X Blocks");
   nXElems_ = paramList->get<int>("X Elements");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> LineMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<double>("X0",0.0);
      defaultParams->set<double>("Xf",1.0);
      defaultParams->set<int>("X Blocks",1);
      defaultParams->set<int>("X Elements",5);

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void LineMeshFactory::initializeWithDefaults()
{
   // get valid parameters
   RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

   // set that parameter list
   setParameterList(validParams);
}

void LineMeshFactory::buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const
{
   typedef shards::Line<2> LineTopo;
   const CellTopologyData * ctd = shards::getCellTopologyData<LineTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(0,0);

   // build meta data
   //mesh.setDimension(2);
   for(int bx=0;bx<xBlocks_;bx++) {
      // add this element block
      {
         std::stringstream ebPostfix;
         ebPostfix << "-" << bx;

         // add element blocks
         mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
      }
   }

   // add sidesets 
   mesh.addSideset("left",side_ctd);
   mesh.addSideset("right",side_ctd);
}

void LineMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();
      // build each block
      for(int xBlock=0;xBlock<xBlocks_;xBlock++) {
         buildBlock(parallelMach,xBlock,mesh);
      }
   mesh.endModification();
}

void LineMeshFactory::buildBlock(stk::ParallelMachine parallelMach,int xBlock,STK_Interface & mesh) const
{
   // grab this processors rank and machine size
   std::pair<int,int> sizeAndStartX = determineXElemSizeAndStart(xBlock,machSize_,machRank_);

   int myXElems_start = sizeAndStartX.first;
   int myXElems_end  = myXElems_start+sizeAndStartX.second;
   int totalXElems = nXElems_*xBlocks_;

   double deltaX = (xf_-x0_)/double(totalXElems);

   // build the nodes
   std::vector<double> coord(1,0.0);
   for(int nx=myXElems_start;nx<myXElems_end+1;++nx) {
      coord[0] = double(nx)*deltaX+x0_;
      mesh.addNode(nx+1,coord);
   }

   std::stringstream blockName;
   blockName << "eblock-" << xBlock;
   stk::mesh::Part * block = mesh.getElementBlockPart(blockName.str());

   // build the elements
   for(int nx=myXElems_start;nx<myXElems_end;++nx) {
      stk::mesh::EntityId gid = nx+1;
      std::vector<stk::mesh::EntityId> nodes(2);
      nodes[0] = nx+1;
      nodes[1] = nodes[0]+1;

      RCP<ElementDescriptor> ed = rcp(new ElementDescriptor(gid,nodes));
      mesh.addElement(ed,block);
   }
}

std::pair<int,int> LineMeshFactory::determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const
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

const stk::mesh::Relation * LineMeshFactory::getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const
{
   for(std::size_t i=0;i<relations.size();i++) 
      if(relations[i].identifier()==ID)
         return &relations[i];

   return 0;
}

void LineMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   std::size_t totalXElems = nXElems_*xBlocks_;

   // get all part vectors
   stk::mesh::Part * left = mesh.getSideset("left");
   stk::mesh::Part * right = mesh.getSideset("right");

   std::vector<stk::mesh::Entity*> localElmts;
   mesh.getMyElements(localElmts);

   // loop over elements adding edges to sidesets
   std::vector<stk::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
      stk::mesh::Entity * element = (*itr);
      stk::mesh::EntityId gid = element->identifier();      
      stk::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

      std::size_t nx = gid-1;

      // vertical boundaries
      ///////////////////////////////////////////

      if(nx+1==totalXElems) { 
         stk::mesh::Entity * edge = getRelationByID(1,relations)->entity();

         // on the right
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,right);
      }

      if(nx==0) {
         stk::mesh::Entity * edge = getRelationByID(0,relations)->entity();

         // on the left
         if(edge->owner_rank()==machRank_)
            mesh.addEntityToSideset(*edge,left);
      }
   }

   mesh.endModification();
}

//! Convert processor rank to a tuple
Teuchos::Tuple<std::size_t,2> LineMeshFactory::procRankToProcTuple(std::size_t procRank) const
{
   std::size_t i=0,j=0;

   j = procRank/machSize_; 
   procRank = procRank % machSize_;
   i = procRank;

   return Teuchos::tuple(i,j);
}

} // end panzer_stk
