// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_QuadraticToLinearMeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp" // for plist validation
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/DumpMeshInfo.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

QuadraticToLinearMeshFactory::QuadraticToLinearMeshFactory(const std::string& quadMeshFileName,
                                                 stk::ParallelMachine mpi_comm,
                                                 const bool print_debug)
  : createEdgeBlocks_(false),
    print_debug_(print_debug)
{
  panzer_stk::STK_ExodusReaderFactory factory;
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("File Name",quadMeshFileName);
  factory.setParameterList(pl);
  quadMesh_ = factory.buildMesh(mpi_comm);

  // TODO for currently supported input topologies, this should always be true
  // but may need to be generalized in the future
  edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;

   // check the conversion is supported
   // and get output topology
   this->getOutputTopology(); 
}

QuadraticToLinearMeshFactory::QuadraticToLinearMeshFactory(const Teuchos::RCP<panzer_stk::STK_Interface>& quadMesh,
                                                 const bool print_debug)
  : quadMesh_(quadMesh),
    createEdgeBlocks_(false),
    print_debug_(print_debug)
{
  // TODO for currently supported input topologies, this should always be true
  // but may need to be generalized in the future
  edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;

   // check the conversion is supported
   // and get output topology
   this->getOutputTopology(); 
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> QuadraticToLinearMeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::QuadraticToLinearMeshFactory::buildMesh()");

   // Make sure the Comms match between the meshes
   {
     int result = MPI_UNEQUAL;
     MPI_Comm_compare(parallelMach, quadMesh_->getBulkData()->parallel(), &result);
     TEUCHOS_ASSERT(result != MPI_UNEQUAL);
   }

   // build all meta data
   RCP<STK_Interface> mesh = this->buildUncommitedMesh(parallelMach);

   // commit meta data
   mesh->initialize(parallelMach);

   // build bulk data
   this->completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

void QuadraticToLinearMeshFactory::getOutputTopology()
{
  bool errFlag = false;

  std::vector<std::string> eblock_names;
  quadMesh_->getElementBlockNames(eblock_names);

  // check that we have a supported topology
  auto inputTopo = quadMesh_->getCellTopology(eblock_names[0]);
  if (std::find(supportedInputTopos_.begin(),
                supportedInputTopos_.end(),*inputTopo) == supportedInputTopos_.end()) errFlag = true;
  TEUCHOS_TEST_FOR_EXCEPTION(errFlag,std::logic_error,
    "ERROR :: Input topology " << inputTopo->getName() << " currently unsupported by QuadraticToLinearMeshFactory!");

  // check that the topology is the same over blocks
  // not sure this is 100% foolproof
  for (auto & eblock : eblock_names) {
    auto cellTopo = quadMesh_->getCellTopology(eblock);
    if (*cellTopo != *inputTopo) errFlag = true;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(errFlag, std::logic_error, 
    "ERROR :: The mesh has different topologies on different blocks!");

  outputTopoData_ = outputTopoMap_[inputTopo->getName()];

  nDim_ = outputTopoData_->dimension;
  nNodes_ = outputTopoData_->node_count;

  return;
}

Teuchos::RCP<STK_Interface> QuadraticToLinearMeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::QuadraticToLinearMeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(nDim_));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   // build meta information: blocks and side set setups
   this->buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void QuadraticToLinearMeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::QuadraticToLinearMeshFactory::completeMeshConstruction()");

   if(not mesh.isInitialized())
      mesh.initialize(parallelMach);

   // add node and element information
   this->buildElements(parallelMach,mesh);

   // finish up the edges
#ifndef ENABLE_UNIFORM
   mesh.buildSubcells();
#endif
   mesh.buildLocalElementIDs();
   if(createEdgeBlocks_) {
      mesh.buildLocalEdgeIDs();
   }

   // now that edges are built, sidesets can be added
#ifndef ENABLE_UNIFORM
   this->addSideSets(mesh);
#endif

   // add nodesets
   this->addNodeSets(mesh);

   // TODO this functionality may be untested
   if(createEdgeBlocks_) {
      this->addEdgeBlocks(mesh);
   }

   // Copy field data
   this->copyCellFieldData(mesh);

   // calls Stk_MeshFactory::rebalance
   // this->rebalance(mesh);
}

//! From ParameterListAcceptor
void QuadraticToLinearMeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   // offsetGIDs_ = (paramList->get<std::string>("Offset mesh GIDs above 32-bit int limit") == "ON") ? true : false;

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> QuadraticToLinearMeshFactory::getValidParameters() const
{
   static RCP<Teuchos::ParameterList> defaultParams;

   // fill with default values
   if(defaultParams == Teuchos::null) {
      defaultParams = rcp(new Teuchos::ParameterList);

      defaultParams->set<std::string>("Offset mesh GIDs above 32-bit int limit", "OFF",
        "If 64-bit GIDs are supported, the mesh element and node global indices will start at a value greater than 32-bit limit.",
        rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("OFF", "ON"))));

      // default to false for backward compatibility
      defaultParams->set<bool>("Create Edge Blocks",false,"Create edge blocks in the mesh");

      Teuchos::ParameterList & bcs = defaultParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
   }

   return defaultParams;
}

void QuadraticToLinearMeshFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
  if (print_debug_) {
    std::cout << "\n\n**** DEBUG: begin printing source quad mesh exodus file metadata ****\n";
    stk::mesh::impl::dump_all_meta_info(*(quadMesh_->getMetaData()), std::cout);
    std::cout << "\n\n**** DEBUG: end printing source quad mesh exodus file metadata ****\n";
  }

  // Required topologies
  auto ctd = outputTopoData_;
  // will only work for 2D and 3D meshes
  const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(nDim_-1,0);

  // Add in element blocks
  std::vector<std::string> element_block_names;
  quadMesh_->getElementBlockNames(element_block_names);
  {
    for (const auto& n : element_block_names)
       mesh.addElementBlock(n,ctd);
  }

  // Add in sidesets
  {
    std::vector<std::string> sideset_names;
    quadMesh_->getSidesetNames(sideset_names);
    for (const auto& n : sideset_names)
      mesh.addSideset(n,side_ctd);
  }

  // Add in nodesets
  {
    std::vector<std::string> nodeset_names;
    quadMesh_->getNodesetNames(nodeset_names);
    for (const auto& n : nodeset_names)
      mesh.addNodeset(n);
  }

  if(createEdgeBlocks_) {
    const CellTopologyData * edge_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);
    std::vector<std::string> element_block_names2;
    quadMesh_->getElementBlockNames(element_block_names2);
    for (const auto& block_name : element_block_names2)
      mesh.addEdgeBlock(block_name,edgeBlockName_,edge_ctd);
  }

  // Add in nodal fields
  {
    const auto& fields = quadMesh_->getMetaData()->get_fields(mesh.getNodeRank());
    for (const auto& f : fields) {
      if (print_debug_)
        std::cout << "Field=" << f->name() << ", rank=" << f->entity_rank() << std::endl;

      // Cull out the coordinate fields. That is a default field in
      // stk that is automatically created by stk.
      if (f->name() != "coordinates") {
        for (const auto& n : element_block_names)
          mesh.addSolutionField(f->name(),n);
      }
    }
  }

  // Add in element fields
  {
    const auto& fields = quadMesh_->getMetaData()->get_fields(mesh.getElementRank());
    for (const auto& f : fields) {
      if (print_debug_)
        std::cout << "Add Cell Field=" << f->name() << ", rank=" << f->entity_rank() << std::endl;

      for (const auto& n : element_block_names)
        mesh.addCellField(f->name(),n);
    }
  }

  // NOTE: skipping edge and face fields for now. Can't error out since sidesets count as edge fields.
  // TEUCHOS_TEST_FOR_EXCEPTION(quadMesh_->getMetaData()->get_fields(mesh.getEdgeRank()).size() != 0,std::runtime_error,
  //                            "ERROR: the Quad8 mesh contains edge fields\""
  //                            << quadMesh_->getMetaData()->get_fields(mesh.getEdgeRank())[0]->name()
  //                            << "\". Edge fields are not supported in Quad8ToQuad4!");

  if (print_debug_) {
    std::cout << "\n\n**** DEBUG: begin printing source linear mesh exodus file metadata ****\n";
    stk::mesh::impl::dump_all_meta_info(*(mesh.getMetaData()), std::cout);
    std::cout << "\n\n**** DEBUG: end printing source linear mesh exodus file metadata ****\n";
  }
}

void QuadraticToLinearMeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();

   auto metadata = mesh.getMetaData();
   auto bulkdata = mesh.getBulkData();

   // Loop over element blocks
   std::vector<std::string> block_names;
   quadMesh_->getElementBlockNames(block_names);
   for (const auto& block_name : block_names) {

     // Get the elements on this process
     std::vector<stk::mesh::Entity> elements;
     quadMesh_->getMyElements(block_name,elements);

     if (print_debug_) {
       std::cout << "*************************************************" << std::endl;
       std::cout << "block_name=" << block_name << ", num_my_elements=" << elements.size() << std::endl;
       std::cout << "*************************************************" << std::endl;
     }

     for (const auto& element : elements) {

       const auto element_gid = quadMesh_->getBulkData()->identifier(element);

       if (print_debug_) {
         std::cout << "rank=" << machRank_
                   << ", block=" << block_name
                   << ", element entity_id=" << element
                   << ", gid=" << element_gid << std::endl;
       }

       // Register nodes with the mesh
       std::vector<stk::mesh::EntityId> nodes(nNodes_);
       for (size_t i=0; i < nNodes_; ++i) {
         // NOTE: this assumes that the linear topology is nested in
         // the quadratic topology as an extended topology - i.e. the first n
         // nodes of the quadratic topology are the corresponding linear
         // nodes. Shards topologies enforce this through the concept
         // of extended topologies.
         stk::mesh::Entity node_entity = quadMesh_->findConnectivityById(element,mesh.getNodeRank(),i);
         TEUCHOS_ASSERT(node_entity.is_local_offset_valid());
         const auto node_gid = quadMesh_->getBulkData()->identifier(node_entity);
         const double* node_coords = quadMesh_->getNodeCoordinates(node_entity);
         std::vector<double> coords_vec(nDim_);
         for (size_t j=0; j < nDim_; ++j) coords_vec[j] = node_coords[j]; 
         mesh.addNode(node_gid,coords_vec);
         nodes[i]=node_gid;
         if (print_debug_) {
           if (nDim_==2) {
             std::cout << "elem gid=" << quadMesh_->getBulkData()->identifier(element)
                       << ", node_gid=" << node_gid << ", (" 
                       << coords_vec[0] << "," << coords_vec[1] << ")\n";
           } else {
             std::cout << "elem gid=" << quadMesh_->getBulkData()->identifier(element)
                       << ", node_gid=" << node_gid << ", (" 
                       << coords_vec[0] << "," << coords_vec[1] << "," << coords_vec[2] << ")\n";
           }
         }
       }

       // Register element with the element block
       auto element_descriptor = panzer_stk::buildElementDescriptor(element_gid,nodes);
       auto element_block_part = mesh.getMetaData()->get_part(block_name);
       mesh.addElement(element_descriptor,element_block_part);
     }
   }
   mesh.endModification();
}

void QuadraticToLinearMeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // Loop over sidesets
   std::vector<std::string> sideset_names;
   quadMesh_->getSidesetNames(sideset_names);
   for (const auto& sideset_name : sideset_names) {

     stk::mesh::Part* sideset_part = mesh.getSideset(sideset_name);

     std::vector<stk::mesh::Entity> q_sides;
     quadMesh_->getMySides(sideset_name,q_sides);

     // Loop over edges
     for (const auto q_ent : q_sides) {
       // The edge numbering scheme uses the element/node gids, so it
       // should be consistent between the quadratic and linear meshes
       // since we used the same gids. We use this fact to populate
       // the quadratic sidesets.
       stk::mesh::EntityId ent_gid = quadMesh_->getBulkData()->identifier(q_ent);
       stk::mesh::Entity lin_ent = mesh.getBulkData()->get_entity(mesh.getSideRank(),ent_gid);
       mesh.addEntityToSideset(lin_ent,sideset_part);
     }
   }

   mesh.endModification();
}

void QuadraticToLinearMeshFactory::addNodeSets(STK_Interface & mesh) const
{}

void QuadraticToLinearMeshFactory::addEdgeBlocks(STK_Interface & mesh) const
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

void QuadraticToLinearMeshFactory::copyCellFieldData(STK_Interface & lin_mesh) const
{
  // Vector of pointers to field data
  const auto& fields = lin_mesh.getMetaData()->get_fields(lin_mesh.getElementRank());

  // loop over fields and add the data to the new mesh.
  for (const auto& field : fields) {

    if (print_debug_)
      std::cout << "Adding field values for \"" << *field << "\" to the linear mesh!\n";

    auto field_name = field->name();

    // Divide into scalar and vector fields, ignoring any other types
    // for now.
    if (field->type_is<double>() &&
        field_name != "coordinates" &&
        field_name != "PROC_ID" &&
        field_name != "LOAD_BAL") {

      // Loop over element blocks
      std::vector<std::string> block_names;
      quadMesh_->getElementBlockNames(block_names);
      for (const auto& block : block_names) {

        auto* lin_field = lin_mesh.getCellField(field_name,block);
        TEUCHOS_ASSERT(lin_field != nullptr);
        // The q mesh doesn't have the field names set up, so a query
        // fails. Go to stk directly in this case.
        auto* q_field = quadMesh_->getMetaData()->get_field(quadMesh_->getElementRank(),field_name);
#ifdef PANZER_DEBUG
        TEUCHOS_ASSERT(q_field != nullptr);
#endif

        // Get the elements for this block.
        std::vector<stk::mesh::Entity> lin_elements;
        lin_mesh.getMyElements(block,lin_elements);
        std::vector<stk::mesh::Entity> q_elements;
        quadMesh_->getMyElements(block,q_elements);
        TEUCHOS_ASSERT(lin_elements.size() == q_elements.size());

        for (size_t i=0; i < lin_elements.size(); ++i) {
#ifdef PANZER_DEBUG
          TEUCHOS_ASSERT(lin_mesh.getBulkData()->identifier(lin_elements[i]) ==
                         quadMesh_->getBulkData()->identifier(q_elements[i]));
#endif

          double* const lin_val = static_cast<double*>(stk::mesh::field_data(*lin_field,lin_elements[i]));
          const double* const q_val = static_cast<double*>(stk::mesh::field_data(*q_field,q_elements[i]));
          *lin_val = *q_val;

          if (print_debug_) {
            std::cout << "field=" << field_name << ", block=" << block
                      << ", lin_e(" << lin_mesh.getBulkData()->identifier(lin_elements[i])  << ") = " << *lin_val
                      << ", q_e(" << quadMesh_->getBulkData()->identifier(q_elements[i])  << ") = " << *q_val
                      << std::endl;
          }

        }
      }
    }
  }
}

} // end panzer_stk
