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

#include <Panzer_STK_CustomMeshFactory.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Panzer_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk_classic {

  CustomMeshFactory::CustomMeshFactory()
  {
    initializeWithDefaults();
  }

  //! Destructor
  CustomMeshFactory::~CustomMeshFactory()
  {
  }

  //! Build the mesh object
  Teuchos::RCP<STK_Interface> 
  CustomMeshFactory::buildMesh(stk_classic::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::buildMesh()");

    // build all meta data
    RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

    // commit meta data
    mesh->initialize(parallelMach);

    // build bulk data
    completeMeshConstruction(*mesh,parallelMach);

    return mesh;
  }

  Teuchos::RCP<STK_Interface> 
  CustomMeshFactory::buildUncommitedMesh(stk_classic::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::buildUncomittedMesh()");

    RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

    machRank_ = stk_classic::parallel_machine_rank(parallelMach);
    machSize_ = stk_classic::parallel_machine_size(parallelMach);

    // add blocks and side sets (global geometry setup)
    buildMetaData(parallelMach,*mesh);
 
    mesh->addPeriodicBCs(periodicBCVec_);

    return mesh;
  }

  void 
  CustomMeshFactory::completeMeshConstruction(STK_Interface &mesh,
                                              stk_classic::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::completeMeshConstruction()");

    if(not mesh.isInitialized())
      mesh.initialize(parallelMach);
   
    // add node and element information
    buildElements(parallelMach,mesh);

    // build edges and faces; fyi: addSides(mesh) builds only edges
    mesh.buildSubcells();  
    mesh.buildLocalElementIDs();

    // now that edges are built, side and node sets can be added
    addSideSets(mesh);

    // calls Stk_MeshFactory::rebalance
    this->rebalance(mesh);
  }

  //! From ParameterListAcceptor
  void 
  CustomMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> &paramList)
  {
    paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

    setMyParamList(paramList);

    Dimension_ = paramList->get<int>("Dimension");

    NumBlocks_ = paramList->get<int>("NumBlocks");

    NumNodesPerProc_ = paramList->get<int>("NumNodesPerProc");
    Nodes_  = paramList->get<int*>("Nodes");

    Coords_ = paramList->get<double*>("Coords");

    NumElementsPerProc_ = paramList->get<int>("NumElementsPerProc");
    Elements_ = paramList->get<int*>("Elements");

    BlockIDs_ = paramList->get<int*>("BlockIDs");   
    Element2Nodes_ = paramList->get<int*>("Element2Nodes");   

    // read in periodic boundary conditions
    parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_);
  }

  //! From ParameterListAcceptor
  Teuchos::RCP<const Teuchos::ParameterList> 
  CustomMeshFactory::getValidParameters() const
  {
    static RCP<Teuchos::ParameterList> defaultParams;

    // fill with default values
    if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<int>("Dimension",3);

      defaultParams->set<int>("NumBlocks",0);
      
      defaultParams->set<int>("NumNodesPerProc",0);
      defaultParams->set<int*>("Nodes",NULL);

      defaultParams->set<double*>("Coords",NULL);

      defaultParams->set<int>("NumElementsPerProc",0);
      defaultParams->set<int*>("Elements",NULL);

      defaultParams->set<int*>("BlockIDs",NULL);
      defaultParams->set<int*>("Element2Nodes",NULL);
      
      Teuchos::ParameterList &bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
    }

    return defaultParams;
  }

  void 
  CustomMeshFactory::initializeWithDefaults()
  {
    // get valid parameters
    RCP<Teuchos::ParameterList> validParams = rcp(new Teuchos::ParameterList(*getValidParameters()));

    // set that parameter list
    setParameterList(validParams);
  }

  void 
  CustomMeshFactory::buildMetaData(stk_classic::ParallelMachine parallelMach, 
                                   STK_Interface &mesh) const
  {
    typedef shards::Hexahedron<8> HexTopo;
    const CellTopologyData * ctd = shards::getCellTopologyData<HexTopo>();
    const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

    for (int blk=0;blk<NumBlocks_;++blk) {
      std::stringstream ebPostfix;
      ebPostfix << "-" <<  blk ;

      // add element blocks
      mesh.addElementBlock("eblock"+ebPostfix.str(),ctd);
    }

    // add sidesets 
    mesh.addSideset("left",side_ctd);
    mesh.addSideset("right",side_ctd);
    mesh.addSideset("top",side_ctd);
    mesh.addSideset("bottom",side_ctd);
    mesh.addSideset("front",side_ctd);
    mesh.addSideset("back",side_ctd);
  }

  void 
  CustomMeshFactory::buildElements(stk_classic::ParallelMachine parallelMach,
                                   STK_Interface &mesh) const
  {
    mesh.beginModification();

    const int dim = mesh.getDimension();

    // build the nodes
    std::vector<double> coord(dim,0.0);
    for (int i=0;i<NumNodesPerProc_;++i) {
      for (int k=0;k<dim;++k)
        coord[k] = Coords_[i*dim+k];
      mesh.addNode(Nodes_[i], coord);
    }
   
    // build the elements
    for (int i=0;i<NumElementsPerProc_;++i) { 

      // get block by its name
      std::stringstream block_name;
      block_name << "eblock-" << BlockIDs_[i];

      stk_classic::mesh::Part *block = mesh.getElementBlockPart(block_name.str());

      // construct element and its nodal connectivity
      stk_classic::mesh::EntityId elt = Elements_[i];
      std::vector<stk_classic::mesh::EntityId> elt2nodes(8);
     
      for (int k=0;k<8;++k)
        elt2nodes[k] = Element2Nodes_[i*8+k];
     
      RCP<ElementDescriptor> ed = rcp(new ElementDescriptor(elt,elt2nodes));
      mesh.addElement(ed,block);
    }

    mesh.endModification();
  }

  void 
  CustomMeshFactory::addSideSets(STK_Interface &mesh) const
  {
    mesh.beginModification();

    // get all part vectors
    // stk_classic::mesh::Part *left = mesh.getSideset("left");
    // stk_classic::mesh::Part *right = mesh.getSideset("right");
    // stk_classic::mesh::Part *top = mesh.getSideset("top");
    // stk_classic::mesh::Part *bottom = mesh.getSideset("bottom");
    // stk_classic::mesh::Part *front = mesh.getSideset("front");
    // stk_classic::mesh::Part *back = mesh.getSideset("back");

    std::vector<stk_classic::mesh::Entity*> localElmts;
    mesh.getMyElements(localElmts);

    // loop over elements adding sides to sidesets
    // std::vector<stk_classic::mesh::Entity*>::const_iterator itr;
    // for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
    //   stk_classic::mesh::Entity * element = (*itr);
    //   stk_classic::mesh::EntityId gid = element->identifier();      
    //   stk_classic::mesh::PairIterRelation relations = element->relations(mesh.getSideRank());

    //   std::size_t nx,ny,nz;
    //   nz = (gid-1) / (totalXElems*totalYElems);
    //   gid = (gid-1)-nz*(totalXElems*totalYElems);
    //   ny = gid / totalXElems;
    //   nx = gid-ny*totalXElems;

    //   if(nz==0) {
    //     stk_classic::mesh::Entity * side = getRelationByID(4,relations)->entity();

    //     // on the back
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,back);
    //   }
    //   if(nz+1==totalZElems) {
    //     stk_classic::mesh::Entity * side = getRelationByID(5,relations)->entity();

    //     // on the front
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,front);
    //   }

    //   if(ny==0) {
    //     stk_classic::mesh::Entity * side = getRelationByID(0,relations)->entity();

    //     // on the bottom 
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,bottom);
    //   }
    //   if(ny+1==totalYElems) {
    //     stk_classic::mesh::Entity * side = getRelationByID(2,relations)->entity();

    //     // on the top
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,top);
    //   }

    //   if(nx==0) {
    //     stk_classic::mesh::Entity * side = getRelationByID(3,relations)->entity();

    //     // on the left
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,left);
    //   }
    //   if(nx+1==totalXElems) {
    //     stk_classic::mesh::Entity * side = getRelationByID(1,relations)->entity();

    //     // on the right
    //     if(side->owner_rank()==machRank_)
    //       mesh.addEntityToSideset(*side,right);
    //   }
    // }

    mesh.endModification();
  }

} // end panzer_stk
