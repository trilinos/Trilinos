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
#include "Panzer_STK_Quad8ToQuad4MeshFactory.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp" // for plist validation
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/DumpMeshInfo.hpp>

// #define ENABLE_UNIFORM

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

Quad8ToQuad4MeshFactory::Quad8ToQuad4MeshFactory(const std::string& quad8MeshFileName,
                                                 stk::ParallelMachine mpi_comm,
                                                 const bool print_debug)
  : createEdgeBlocks_(false),
    print_debug_(print_debug)
{
  panzer_stk::STK_ExodusReaderFactory factory;
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("File Name",quad8MeshFileName);
  factory.setParameterList(pl);
  quad8Mesh_ = factory.buildMesh(mpi_comm);

  edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;
}

Quad8ToQuad4MeshFactory::Quad8ToQuad4MeshFactory(const Teuchos::RCP<panzer_stk::STK_Interface>& quad8Mesh,
                                                 const bool print_debug)
  : quad8Mesh_(quad8Mesh),
    createEdgeBlocks_(false),
    print_debug_(print_debug)
{
  edgeBlockName_ = "line_2_"+panzer_stk::STK_Interface::edgeBlockString;
}

//! Build the mesh object
Teuchos::RCP<STK_Interface> Quad8ToQuad4MeshFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::Quad8ToQuad4MeshFactory::buildMesh()");

   // Make sure the Quad8 and Quad4 Comms match
   {
     int result = MPI_UNEQUAL;
     MPI_Comm_compare(parallelMach, quad8Mesh_->getBulkData()->parallel(), &result);
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

Teuchos::RCP<STK_Interface> Quad8ToQuad4MeshFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::Quad8ToQuad4MeshFactory::buildUncomittedMesh()");

   RCP<STK_Interface> mesh = rcp(new STK_Interface(2));

   machRank_ = stk::parallel_machine_rank(parallelMach);
   machSize_ = stk::parallel_machine_size(parallelMach);

   // build meta information: blocks and side set setups
   this->buildMetaData(parallelMach,*mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void Quad8ToQuad4MeshFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::Quad8ToQuad4MeshFactory::completeMeshConstruction()");

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

   // now that edges are built, sidsets can be added
#ifndef ENABLE_UNIFORM
   this->addSideSets(mesh);
#endif

   // add nodesets
   this->addNodeSets(mesh);

   if(createEdgeBlocks_) {
      this->addEdgeBlocks(mesh);
   }

   // Copy field data
   this->copyCellFieldData(mesh);

   // calls Stk_MeshFactory::rebalance
   // this->rebalance(mesh);
}

//! From ParameterListAcceptor
void Quad8ToQuad4MeshFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   paramList->validateParametersAndSetDefaults(*getValidParameters(),0);

   setMyParamList(paramList);

   // offsetGIDs_ = (paramList->get<std::string>("Offset mesh GIDs above 32-bit int limit") == "ON") ? true : false;

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> Quad8ToQuad4MeshFactory::getValidParameters() const
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

void Quad8ToQuad4MeshFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
  if (print_debug_) {
    std::cout << "\n\n**** DEBUG: begin printing source Quad8 exodus file metadata ****\n";
    stk::mesh::impl::dump_all_meta_info(*(quad8Mesh_->getMetaData()), std::cout);
    std::cout << "\n\n**** DEBUG: end printing source Quad8 exodus file metadata ****\n";
  }

  // Required topologies
  using QuadTopo = shards::Quadrilateral<4>;
  const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
  const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

  // Add in element blocks
  std::vector<std::string> element_block_names;
  quad8Mesh_->getElementBlockNames(element_block_names);
  {
    for (const auto& n : element_block_names)
       mesh.addElementBlock(n,ctd);
  }

  // Add in sidesets
  {
    std::vector<std::string> sideset_names;
    quad8Mesh_->getSidesetNames(sideset_names);
    for (const auto& n : sideset_names)
      mesh.addSideset(n,side_ctd);
  }

  // Add in nodesets
  {
    std::vector<std::string> nodeset_names;
    quad8Mesh_->getNodesetNames(nodeset_names);
    for (const auto& n : nodeset_names)
      mesh.addNodeset(n);
  }

  if(createEdgeBlocks_) {
    const CellTopologyData * edge_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);
    std::vector<std::string> element_block_names2;
    quad8Mesh_->getElementBlockNames(element_block_names2);
    for (const auto& block_name : element_block_names2)
      mesh.addEdgeBlock(block_name,edgeBlockName_,edge_ctd);
  }

  // Add in nodal fields
  {
    const auto& fields = quad8Mesh_->getMetaData()->get_fields(mesh.getNodeRank());
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
    const auto& fields = quad8Mesh_->getMetaData()->get_fields(mesh.getElementRank());
    for (const auto& f : fields) {
      if (print_debug_)
        std::cout << "Add Cell Field=" << f->name() << ", rank=" << f->entity_rank() << std::endl;

      for (const auto& n : element_block_names)
        mesh.addCellField(f->name(),n);
    }
  }

  // NOTE: skipping edge and face fields for now. Can't error out since sidesets count as edge fields.
  // TEUCHOS_TEST_FOR_EXCEPTION(quad8Mesh_->getMetaData()->get_fields(mesh.getEdgeRank()).size() != 0,std::runtime_error,
  //                            "ERROR: the Quad8 mesh contains edge fields\""
  //                            << quad8Mesh_->getMetaData()->get_fields(mesh.getEdgeRank())[0]->name()
  //                            << "\". Edge fields are not supported in Quad8ToQuad4!");

  if (print_debug_) {
    std::cout << "\n\n**** DEBUG: begin printing source Quad4 exodus file metadata ****\n";
    stk::mesh::impl::dump_all_meta_info(*(mesh.getMetaData()), std::cout);
    std::cout << "\n\n**** DEBUG: end printing source Quad4 exodus file metadata ****\n";
  }
}

void Quad8ToQuad4MeshFactory::buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const
{
   mesh.beginModification();

   auto metadata = mesh.getMetaData();
   auto bulkdata = mesh.getBulkData();

   // Loop over element blocks
   std::vector<std::string> block_names;
   quad8Mesh_->getElementBlockNames(block_names);
   for (const auto& block_name : block_names) {

     // Get the elements on this process
     std::vector<stk::mesh::Entity> elements;
     quad8Mesh_->getMyElements(block_name,elements);

     if (print_debug_) {
       std::cout << "*************************************************" << std::endl;
       std::cout << "block_name=" << block_name << ", num_my_elements=" << elements.size() << std::endl;
       std::cout << "*************************************************" << std::endl;
     }

     for (const auto& element : elements) {

       const auto element_gid = quad8Mesh_->getBulkData()->identifier(element);

       if (print_debug_) {
         std::cout << "rank=" << machRank_
                   << ", block=" << block_name
                   << ", element entity_id=" << element
                   << ", gid=" << element_gid << std::endl;
       }

       // Register nodes with the mesh
       std::vector<stk::mesh::EntityId> nodes(4);
       for (int i=0; i < 4; ++i) {
         // NOTE: this assumes that the Quad4 topology is nested in
         // the Quad8 as an extended topology - i.e. the first four
         // nodes of the quad8 topology are the corresponding quad4
         // nodes. Shards topologies enforce this throught the concept
         // of extended topologies.
         stk::mesh::Entity node_entity = quad8Mesh_->findConnectivityById(element,mesh.getNodeRank(),i);
         TEUCHOS_ASSERT(node_entity.is_local_offset_valid());
         const auto node_gid = quad8Mesh_->getBulkData()->identifier(node_entity);
         const double* node_coords = quad8Mesh_->getNodeCoordinates(node_entity);
         std::vector<double> coords_vec(2);
         coords_vec[0] = node_coords[0];
         coords_vec[1] = node_coords[1];
         mesh.addNode(node_gid,coords_vec);
         nodes[i]=node_gid;
         if (print_debug_) {
           std::cout << "elem gid=" << quad8Mesh_->getBulkData()->identifier(element)
                     << ", node_gid=" << node_gid << ", (" << coords_vec[0] << "," << coords_vec[1] << ")\n";
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

void Quad8ToQuad4MeshFactory::addSideSets(STK_Interface & mesh) const
{
   mesh.beginModification();

   // Loop over sidesets
   std::vector<std::string> sideset_names;
   quad8Mesh_->getSidesetNames(sideset_names);
   for (const auto& sideset_name : sideset_names) {

     stk::mesh::Part* sideset_part = mesh.getSideset(sideset_name);

     std::vector<stk::mesh::Entity> q8_sides;
     quad8Mesh_->getMySides(sideset_name,q8_sides);

     // Loop over edges
     for (const auto q8_edge : q8_sides) {
       // The edge numbering scheme uses the element/node gids, so it
       // should be consistent between the quad8 and quad4 meshes
       // since we used the same gids. We use this fact to populate
       // the quad4 sidesets.
       stk::mesh::EntityId edge_gid = quad8Mesh_->getBulkData()->identifier(q8_edge);
       stk::mesh::Entity q4_edge = mesh.getBulkData()->get_entity(mesh.getEdgeRank(),edge_gid);
       mesh.addEntityToSideset(q4_edge,sideset_part);
     }
   }

   mesh.endModification();
}

void Quad8ToQuad4MeshFactory::addNodeSets(STK_Interface & mesh) const
{}

void Quad8ToQuad4MeshFactory::addEdgeBlocks(STK_Interface & mesh) const
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

void Quad8ToQuad4MeshFactory::copyCellFieldData(STK_Interface & q4_mesh) const
{
  // Vector of pointers to field data
  const auto& fields = q4_mesh.getMetaData()->get_fields(q4_mesh.getElementRank());

  // loop over fields and add the data to the new mesh.
  for (const auto& field : fields) {

    if (print_debug_)
      std::cout << "Adding field values for \"" << *field << "\" to the Quad4 mesh!\n";

    auto field_name = field->name();

    // Divide into scalar and vector fields, ignoring any other types
    // for now.
    if (field->type_is<double>() &&
        field_name != "coordinates" &&
        field_name != "PROC_ID" &&
        field_name != "LOAD_BAL") {

      // Loop over element blocks
      std::vector<std::string> block_names;
      quad8Mesh_->getElementBlockNames(block_names);
      for (const auto& block : block_names) {

        auto* q4_field = q4_mesh.getCellField(field_name,block);
        TEUCHOS_ASSERT(q4_field != nullptr);
        // The q8 mesh doesn't have the field names set up, so a query
        // fails. Go to stk directly in this case.
        auto* q8_field = quad8Mesh_->getMetaData()->get_field(quad8Mesh_->getElementRank(),field_name);
#ifdef PANZER_DEBUG
        TEUCHOS_ASSERT(q8_field != nullptr);
#endif

        // Get the elements for this block.
        std::vector<stk::mesh::Entity> q4_elements;
        q4_mesh.getMyElements(block,q4_elements);
        std::vector<stk::mesh::Entity> q8_elements;
        quad8Mesh_->getMyElements(block,q8_elements);
        TEUCHOS_ASSERT(q4_elements.size() == q8_elements.size());

        for (size_t i=0; i < q4_elements.size(); ++i) {
#ifdef PANZER_DEBUG
          TEUCHOS_ASSERT(q4_mesh.getBulkData()->identifier(q4_elements[i]) ==
                         quad8Mesh_->getBulkData()->identifier(q8_elements[i]));
#endif

          double* const q4_val = static_cast<double*>(stk::mesh::field_data(*q4_field,q4_elements[i]));
          const double* const q8_val = static_cast<double*>(stk::mesh::field_data(*q8_field,q8_elements[i]));
          *q4_val = *q8_val;

          if (print_debug_) {
            std::cout << "field=" << field_name << ", block=" << block
                      << ", q4e(" << q4_mesh.getBulkData()->identifier(q4_elements[i])  << ") = " << *q4_val
                      << ", q8e(" << quad8Mesh_->getBulkData()->identifier(q8_elements[i])  << ") = " << *q8_val
                      << std::endl;
          }

        }
      }
    }
  }
}

} // end panzer_stk
