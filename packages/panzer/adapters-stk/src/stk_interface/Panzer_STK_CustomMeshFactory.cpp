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
#include <PanzerAdaptersSTK_config.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

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
  CustomMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::buildMesh()");

    // build all meta data
    RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

    // commit meta data; allocate field variables
    mesh->initialize(parallelMach);

    // build bulk data
    completeMeshConstruction(*mesh,parallelMach);

    return mesh;
  }

  Teuchos::RCP<STK_Interface>
  CustomMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::buildUncomittedMesh()");

    RCP<STK_Interface> mesh = rcp(new STK_Interface(3));

    machRank_ = stk::parallel_machine_rank(parallelMach);
    machSize_ = stk::parallel_machine_size(parallelMach);

    // add blocks and side sets (global geometry setup)
    buildMetaData(*mesh);

    mesh->addPeriodicBCs(periodicBCVec_);

    return mesh;
  }

  void
  CustomMeshFactory::completeMeshConstruction(STK_Interface &mesh,
                                              stk::ParallelMachine parallelMach) const
  {
    PANZER_FUNC_TIME_MONITOR("panzer::CustomMeshFactory::completeMeshConstruction()");

    if (not mesh.isInitialized())
      mesh.initialize(parallelMach);

    // add node and element information
    buildElements(mesh);

    // build edges and faces; fyi: addSides(mesh) builds only edges
    mesh.buildSubcells();
    mesh.buildLocalElementIDs();



    // now that edges are built, side and node sets can be added
    addSideSets(mesh);

    // set solution fields
    fillSolutionFieldData(mesh);

    // calls Stk_MeshFactory::rebalance
    //this->rebalance(mesh);
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

    BlockIDs_ = paramList->get<int*>("BlockIDs");
    Element2Nodes_ = paramList->get<int*>("Element2Nodes");

    OffsetToGlobalElementIDs_ = paramList->get<int>("OffsetToGlobalElementIDs");

    // solution fields from tramonto
    ChargeDensity_ = paramList->get<double*>("ChargeDensity");
    ElectricPotential_ = paramList->get<double*>("ElectricPotential");

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

      defaultParams->set<int*>("BlockIDs",NULL);
      defaultParams->set<int*>("Element2Nodes",NULL);

      defaultParams->set<int>("OffsetToGlobalElementIDs", 0);

      defaultParams->set<double*>("ChargeDensity",NULL);
      defaultParams->set<double*>("ElectricPotential",NULL);

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
  CustomMeshFactory::buildMetaData(STK_Interface &mesh) const
  {
    typedef shards::Hexahedron<8> HexTopo;
    const CellTopologyData *ctd = shards::getCellTopologyData<HexTopo>();
    const CellTopologyData *side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(2,0);

    for (int blk=0;blk<NumBlocks_;++blk) {
      std::stringstream block_id;
      block_id << "eblock-" << blk;

      // add element blocks
      mesh.addElementBlock(block_id.str(),ctd);

      mesh.addSolutionField("CHARGE_DENSITY",block_id.str());
      mesh.addSolutionField("ELECTRIC_POTENTIAL",block_id.str());
    }

    mesh.addSideset("left",  side_ctd);
    mesh.addSideset("right", side_ctd);
    mesh.addSideset("top",   side_ctd);
    mesh.addSideset("bottom",side_ctd);
    mesh.addSideset("front", side_ctd);
    mesh.addSideset("back",  side_ctd);

    mesh.addSideset("wall",  side_ctd);
  }

  void
  CustomMeshFactory::buildElements(STK_Interface &mesh) const
  {
    mesh.beginModification();

    const int dim = mesh.getDimension();

    // build the nodes
    std::vector<double> coords(dim,0.0);
    for (int i=0;i<NumNodesPerProc_;++i) {
      for (int k=0;k<dim;++k)
        coords[k] = Coords_[i*dim+k];
      mesh.addNode(Nodes_[i], coords);
    }

    // build the elements
    std::vector<std::string> block_ids;
    mesh.getElementBlockNames(block_ids);  

    for (int i=0;i<NumElementsPerProc_;++i) {

      // get block by its name
      std::stringstream block_id;
      block_id << "eblock-" << BlockIDs_[i];

      stk::mesh::Part *block = mesh.getElementBlockPart(block_id.str());

      // construct element and its nodal connectivity
      stk::mesh::EntityId elt = i + OffsetToGlobalElementIDs_;
      std::vector<stk::mesh::EntityId> elt2nodes(8);

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
    stk::mesh::BulkData& bulkData = *mesh.getBulkData();
    const stk::mesh::EntityRank sideRank = mesh.getSideRank();

    // get all part vectors
    stk::mesh::Part *box[6];
    box[0] = mesh.getSideset("front");
    box[1] = mesh.getSideset("right");
    box[2] = mesh.getSideset("back");
    box[3] = mesh.getSideset("left");
    box[4] = mesh.getSideset("bottom");
    box[5] = mesh.getSideset("top");

    stk::mesh::Part *wall = mesh.getSideset("wall");

    std::vector<stk::mesh::Entity> elements;
    mesh.getMyElements(elements);

    // loop over elements adding sides to sidesets
    for (std::vector<stk::mesh::Entity>::const_iterator
           itr=elements.begin();itr!=elements.end();++itr) {
      stk::mesh::Entity element = (*itr);
      const size_t numSides = bulkData.num_connectivity(element, sideRank);
      stk::mesh::Entity const* relations = bulkData.begin(element, sideRank);

      // loop over side id checking element neighbors
      for (std::size_t i=0;i<numSides;++i) {
        stk::mesh::Entity side = relations[i];
        const size_t numNeighbors = bulkData.num_elements(side);
        stk::mesh::Entity const* neighbors = bulkData.begin_elements(side);

        if (numNeighbors == 1) {
          if (mesh.entityOwnerRank(side) == machRank_)
            mesh.addEntityToSideset(side, box[i]);
        }
        else if (numNeighbors == 2) {
          std::string neig_block_id_0 = mesh.containingBlockId(neighbors[0]);
          std::string neig_block_id_1 = mesh.containingBlockId(neighbors[1]);
          if (neig_block_id_0 != neig_block_id_1 && mesh.entityOwnerRank(side) == machRank_)
            mesh.addEntityToSideset(side, wall);
        }
        else {
          // runtime exception
        }
      }
    }

    mesh.endModification();
  }

  void
  CustomMeshFactory::fillSolutionFieldData(STK_Interface &mesh) const
  {
    for (int blk=0;blk<NumBlocks_;++blk) {

      std::stringstream block_id;
      block_id << "eblock-" << blk;
      
      // elements in this processor for this block
      std::vector<stk::mesh::Entity> elements;    
      mesh.getMyElements(block_id.str(), elements);

      // size of elements in the current block
      std::size_t n_elements = elements.size();
      
      // build local element index
      std::vector<std::size_t> local_ids;
      for (std::vector<stk::mesh::Entity>::const_iterator
             itr=elements.begin();itr!=elements.end();++itr) 
        local_ids.push_back(mesh.elementLocalId(*itr));

      // re-index solution fields in the same order of local_ids
      std::vector<double> charge_density_by_local_ids, electric_potential_by_local_ids;
      for (std::vector<stk::mesh::Entity>::const_iterator
             itr=elements.begin();itr!=elements.end();++itr) {
        int q = mesh.elementGlobalId(*itr) - OffsetToGlobalElementIDs_;
        for (int k=0;k<8;++k) {
          int loc = q*8 + k;
          charge_density_by_local_ids.push_back(ChargeDensity_[loc]);
          electric_potential_by_local_ids.push_back(ElectricPotential_[loc]);
        }
      }

      // wrap the buffer with a proper container
      FieldContainer charge_density(n_elements, 8, &charge_density_by_local_ids[0]),
        electric_potential(n_elements, 8, &electric_potential_by_local_ids[0]);

      // write out to stk mesh
      mesh.setSolutionFieldData("CHARGE_DENSITY",
                                block_id.str(),
                                local_ids,
                                charge_density, 1.0);
      
      mesh.setSolutionFieldData("ELECTRIC_POTENTIAL",
                                block_id.str(),
                                local_ids,
                                electric_potential, 1.0);
    }
  }
} // end panzer_stk
