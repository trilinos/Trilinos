// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <PanzerAdaptersSTK_config.hpp>

#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Interface.hpp"

#ifdef PANZER_HAVE_IOSS

#include <Ionit_Initializer.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_EdgeBlock.h>
#include <Ioss_FaceBlock.h>
#include <Ioss_Region.h>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#ifdef PANZER_HAVE_UMR
#include <Ioumr_DatabaseIO.hpp>
#endif

#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace panzer_stk {

int getMeshDimension(const std::string & meshStr,
                     stk::ParallelMachine parallelMach,
                     const std::string & typeStr)
{
  stk::io::StkMeshIoBroker meshData(parallelMach);
  meshData.use_simple_fields();
  meshData.property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", false));
  meshData.add_mesh_database(meshStr, fileTypeToIOSSType(typeStr), stk::io::READ_MESH);
  meshData.create_input_mesh();
  return Teuchos::as<int>(meshData.meta_data_ptr()->spatial_dimension());
}

std::string fileTypeToIOSSType(const std::string & fileType)
{
  std::string IOSSType;
  if      (fileType=="Exodus")
    IOSSType = "exodusii";
#ifdef PANZER_HAVE_UMR
  else if (fileType=="Exodus Refinement")
    IOSSType = "Refinement";
#endif
  else if (fileType=="Pamgen")
    IOSSType = "pamgen";
  
  return IOSSType;
}

STK_ExodusReaderFactory::STK_ExodusReaderFactory()
  : fileName_(""), fileType_(""), restartIndex_(0), userMeshScaling_(false), keepPerceptData_(false),
    keepPerceptParentElements_(false), rebalancing_("default"),
meshScaleFactor_(0.0), levelsOfRefinement_(0),
    createEdgeBlocks_(false), createFaceBlocks_(false), geometryName_("")
{ }

STK_ExodusReaderFactory::STK_ExodusReaderFactory(const std::string & fileName,
                                                 const int restartIndex)
  : fileName_(fileName), fileType_("Exodus"), restartIndex_(restartIndex), userMeshScaling_(false),
    keepPerceptData_(false), keepPerceptParentElements_(false), rebalancing_("default"),
    meshScaleFactor_(0.0), levelsOfRefinement_(0), createEdgeBlocks_(false), createFaceBlocks_(false), geometryName_("")
{ }

Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::buildMesh()");

   using Teuchos::RCP;
   using Teuchos::rcp;

   RCP<STK_Interface> mesh = buildUncommitedMesh(parallelMach);

   // in here you would add your fields...but it is better to use
   // the two step construction

   // this calls commit on meta data
   mesh->initialize(parallelMach,false,doPerceptRefinement());

   completeMeshConstruction(*mesh,parallelMach);

   return mesh;
}

/** This builds all the meta data of the mesh. Does not call metaData->commit.
  * Allows user to add solution fields and other pieces. The mesh can be "completed"
  * by calling <code>completeMeshConstruction</code>.
  */
Teuchos::RCP<STK_Interface> STK_ExodusReaderFactory::buildUncommitedMesh(stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::buildUncomittedMesh()");

   using Teuchos::RCP;
   using Teuchos::rcp;

   // read in meta data
   stk::io::StkMeshIoBroker* meshData = new stk::io::StkMeshIoBroker(parallelMach);
   meshData->use_simple_fields();
   meshData->property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", false));

   // add in "FAMILY_TREE" entity for doing refinement
   std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
   entity_rank_names.push_back("FAMILY_TREE");
   meshData->set_rank_name_vector(entity_rank_names);

#ifdef PANZER_HAVE_UMR
   // this line registers Ioumr with Ioss
   Ioumr::IOFactory::factory();

   meshData->property_add(Ioss::Property("GEOMETRY_FILE", geometryName_));
   meshData->property_add(Ioss::Property("NUMBER_REFINEMENTS", levelsOfRefinement_));
#endif

   meshData->add_mesh_database(fileName_, fileTypeToIOSSType(fileType_), stk::io::READ_MESH);

   meshData->create_input_mesh();
   RCP<stk::mesh::MetaData> metaData = Teuchos::rcp(meshData->meta_data_ptr());

   RCP<STK_Interface> mesh = rcp(new STK_Interface(metaData));
   mesh->initializeFromMetaData();
   mesh->instantiateBulkData(parallelMach);
   meshData->set_bulk_data(Teuchos::get_shared_ptr(mesh->getBulkData()));

   // read in other transient fields, these will be useful later when
   // trying to read other fields for use in solve
   meshData->add_all_mesh_fields_as_input_fields();

   // store mesh data pointer for later use in initializing
   // bulk data
   mesh->getMetaData()->declare_attribute_with_delete(meshData);

   // build element blocks
   registerElementBlocks(*mesh,*meshData);
   registerSidesets(*mesh);
   registerNodesets(*mesh);

   if (createEdgeBlocks_) {
      registerEdgeBlocks(*mesh,*meshData);
   }
   if (createFaceBlocks_ && mesh->getMetaData()->spatial_dimension() > 2) {
      registerFaceBlocks(*mesh,*meshData);
   }

   buildMetaData(parallelMach, *mesh);

   mesh->addPeriodicBCs(periodicBCVec_);
   mesh->setBoundingBoxSearchFlag(useBBoxSearch_);

   return mesh;
}

void STK_ExodusReaderFactory::completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::STK_ExodusReaderFactory::completeMeshConstruction()");

   using Teuchos::RCP;
   using Teuchos::rcp;


   if(not mesh.isInitialized()) {
     mesh.initialize(parallelMach,true,doPerceptRefinement());
   }

   // grab mesh data pointer to build the bulk data
   stk::mesh::MetaData & metaData = *mesh.getMetaData();
   stk::mesh::BulkData & bulkData = *mesh.getBulkData();
   stk::io::StkMeshIoBroker * meshData =
     const_cast<stk::io::StkMeshIoBroker *>(metaData.get_attribute<stk::io::StkMeshIoBroker>());
         // if const_cast is wrong ... why does it feel so right?
         // I believe this is safe since we are basically hiding this object under the covers
         // until the mesh construction can be completed...below I cleanup the object myself.
   TEUCHOS_ASSERT(metaData.remove_attribute(meshData));
      // remove the MeshData attribute

   // build mesh bulk data
   meshData->populate_bulk_data();

   if (doPerceptRefinement()) {
     const bool deleteParentElements = !keepPerceptParentElements_;
     mesh.refineMesh(levelsOfRefinement_,deleteParentElements);
   }

   // The following section of code is applicable if mesh scaling is
   // turned on from the input file.
   if (userMeshScaling_)
   {
     stk::mesh::Field<double>* coord_field = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");

     stk::mesh::Selector select_all_local = metaData.locally_owned_part() | metaData.globally_shared_part();
     stk::mesh::BucketVector const& my_node_buckets = bulkData.get_buckets(stk::topology::NODE_RANK, select_all_local);

     int mesh_dim = mesh.getDimension();

     // Scale the mesh
     const double inv_msf = 1.0/meshScaleFactor_;
     for (size_t i=0; i < my_node_buckets.size(); ++i)
     {
       stk::mesh::Bucket& b = *(my_node_buckets[i]);
       double* coordinate_data = field_data( *coord_field, b );

       for (size_t j=0; j < b.size(); ++j) {
         for (int k=0; k < mesh_dim; ++k) {
           coordinate_data[mesh_dim*j + k] *= inv_msf;
         }
       }
     }
   }

   // put in a negative index and (like python) the restart will be from the back
   // (-1 is the last time step)
   int restartIndex = restartIndex_;
   if(restartIndex<0) {
     std::pair<int,double> lastTimeStep = meshData->get_input_ioss_region()->get_max_time();
     restartIndex = 1+restartIndex+lastTimeStep.first;
   }

   // populate mesh fields with specific index
   meshData->read_defined_input_fields(restartIndex);

   mesh.buildSubcells();
   mesh.buildLocalElementIDs();

   mesh.beginModification();
   if (createEdgeBlocks_) {
      mesh.buildLocalEdgeIDs();
      addEdgeBlocks(mesh);
   }
   if (createFaceBlocks_ && mesh.getMetaData()->spatial_dimension() > 2) {
      mesh.buildLocalFaceIDs();
      addFaceBlocks(mesh);
   }
   mesh.endModification();

   if (userMeshScaling_) {
     stk::mesh::Field<double>* coord_field = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");
     std::vector< const stk::mesh::FieldBase *> fields;
     fields.push_back(coord_field);

     stk::mesh::communicate_field_data(bulkData, fields);
   }

   if(restartIndex>0) // process_input_request is a no-op if restartIndex<=0 ... thus there would be no inital time
      mesh.setInitialStateTime(meshData->get_input_ioss_region()->get_state_time(restartIndex));
   else
      mesh.setInitialStateTime(0.0); // no initial time to speak, might as well use 0.0

   // clean up mesh data object
   delete meshData;

   if(rebalancing_ == "default")
     // calls Stk_MeshFactory::rebalance
     this->rebalance(mesh);
   else if(rebalancing_ != "none")
   {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
				"ERROR: Rebalancing was not set to a valid choice");
   }
}

//! From ParameterListAcceptor
void STK_ExodusReaderFactory::setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList)
{
   TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(!paramList->isParameter("File Name"),
        Teuchos::Exceptions::InvalidParameterName,
        "Error, the parameter {name=\"File Name\","
        "type=\"string\""
        "\nis required in parameter (sub)list \""<< paramList->name() <<"\"."
        "\n\nThe parsed parameter parameter list is: \n" << paramList->currentParametersString()
   );

   // Set default values here. Not all the params should be set so this
   // has to be done manually as opposed to using
   // validateParametersAndSetDefaults().
   if(!paramList->isParameter("Restart Index"))
     paramList->set<int>("Restart Index", -1);

   if(!paramList->isParameter("File Type"))
     paramList->set("File Type", "Exodus");

   if(!paramList->isSublist("Periodic BCs"))
     paramList->sublist("Periodic BCs");

   Teuchos::ParameterList& p_bcs = paramList->sublist("Periodic BCs");
   if (!p_bcs.isParameter("Count"))
     p_bcs.set<int>("Count", 0);

   if(!paramList->isParameter("Levels of Uniform Refinement"))
     paramList->set<int>("Levels of Uniform Refinement", 0);

   if(!paramList->isParameter("Keep Percept Data"))
     paramList->set<bool>("Keep Percept Data", false);

   if(!paramList->isParameter("Keep Percept Parent Elements"))
     paramList->set<bool>("Keep Percept Parent Elements", false);

   if(!paramList->isParameter("Rebalancing"))
     paramList->set<std::string>("Rebalancing", "default");

   if(!paramList->isParameter("Create Edge Blocks"))
     // default to false to prevent massive exodiff test failures
     paramList->set<bool>("Create Edge Blocks", false);

   if(!paramList->isParameter("Create Face Blocks"))
     // default to false to prevent massive exodiff test failures
     paramList->set<bool>("Create Face Blocks", false);

   if(!paramList->isParameter("Geometry File Name"))
     paramList->set("Geometry File Name", "");

   paramList->validateParameters(*getValidParameters(),0);

   setMyParamList(paramList);

   fileName_ = paramList->get<std::string>("File Name");

   geometryName_ = paramList->get<std::string>("Geometry File Name");

   restartIndex_ = paramList->get<int>("Restart Index");

   fileType_ = paramList->get<std::string>("File Type");

   // get any mesh scale factor
   if (paramList->isParameter("Scale Factor"))
   {
     meshScaleFactor_ = paramList->get<double>("Scale Factor");
     userMeshScaling_ = true;
   }

   // read in periodic boundary conditions
   parsePeriodicBCList(Teuchos::rcpFromRef(paramList->sublist("Periodic BCs")),periodicBCVec_,useBBoxSearch_);

   keepPerceptData_ = paramList->get<bool>("Keep Percept Data");

   keepPerceptParentElements_ = paramList->get<bool>("Keep Percept Parent Elements");

   rebalancing_ = paramList->get<std::string>("Rebalancing");

   levelsOfRefinement_ = paramList->get<int>("Levels of Uniform Refinement");

   createEdgeBlocks_ = paramList->get<bool>("Create Edge Blocks");
   createFaceBlocks_ = paramList->get<bool>("Create Face Blocks");
}

//! From ParameterListAcceptor
Teuchos::RCP<const Teuchos::ParameterList> STK_ExodusReaderFactory::getValidParameters() const
{
   static Teuchos::RCP<Teuchos::ParameterList> validParams;

   if(validParams==Teuchos::null) {
      validParams = Teuchos::rcp(new Teuchos::ParameterList);
      validParams->set<std::string>("File Name","<file name not set>","Name of exodus file to be read",
                                    Teuchos::rcp(new Teuchos::FileNameValidator));
      validParams->set<std::string>("Geometry File Name","<file name not set>","Name of geometry file for refinement", 
                                    Teuchos::rcp(new Teuchos::FileNameValidator));

      validParams->set<int>("Restart Index",-1,"Index of solution to read in",
			    Teuchos::rcp(new Teuchos::AnyNumberParameterEntryValidator(Teuchos::AnyNumberParameterEntryValidator::PREFER_INT,Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes(true))));

      validParams->set<std::string>("File Type", "Exodus",
        "Choose input file type - either \"Exodus\", \"Exodus Refinement\" or \"Pamgen\"",
        rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("Exodus", "Pamgen"
#ifdef PANZER_HAVE_UMR
                                                                     ,"Exodus Refinement"
#endif        
        ))));

      validParams->set<double>("Scale Factor", 1.0, "Scale factor to apply to mesh after read",
                               Teuchos::rcp(new Teuchos::AnyNumberParameterEntryValidator(Teuchos::AnyNumberParameterEntryValidator::PREFER_DOUBLE,Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes(true))));

      Teuchos::ParameterList & bcs = validParams->sublist("Periodic BCs");
      bcs.set<int>("Count",0); // no default periodic boundary conditions
      
      validParams->set("Levels of Uniform Refinement",0,"Number of levels of inline uniform mesh refinement");

      validParams->set("Keep Percept Data",false,"Keep the Percept mesh after uniform refinement is applied");

      validParams->set("Keep Percept Parent Elements",false,"Keep the parent element information in the Percept data");

      validParams->set("Rebalancing","default","The type of rebalancing to be performed on the mesh after creation (default, none)");

      // default to false for backward compatibility
      validParams->set("Create Edge Blocks",false,"Create or copy edge blocks in the mesh");
      validParams->set("Create Face Blocks",false,"Create or copy face blocks in the mesh");
   }

   return validParams.getConst();
}

void STK_ExodusReaderFactory::registerElementBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> femMetaData = mesh.getMetaData();

   // here we use the Ioss interface because they don't add
   // "bonus" element blocks and its easier to determine
   // "real" element blocks versus STK-only blocks
   const Ioss::ElementBlockContainer & elem_blocks = meshData.get_input_ioss_region()->get_element_blocks();
   for(Ioss::ElementBlockContainer::const_iterator itr=elem_blocks.begin();itr!=elem_blocks.end();++itr) {
      Ioss::GroupingEntity * entity = *itr;
      const std::string & name = entity->name();

      const stk::mesh::Part * part = femMetaData->get_part(name);
      shards::CellTopology cellTopo = stk::mesh::get_cell_topology(femMetaData->get_topology(*part));
      const CellTopologyData * ct = cellTopo.getCellTopologyData();

      TEUCHOS_ASSERT(ct!=0);
      mesh.addElementBlock(part->name(),ct);

      if (createEdgeBlocks_) {
         createUniqueEdgeTopologyMap(mesh, part);
      }
      if (createFaceBlocks_ && mesh.getMetaData()->spatial_dimension() > 2) {
         createUniqueFaceTopologyMap(mesh, part);
      }
   }
}

template <typename SetType>
void buildSetNames(const SetType & setData,std::vector<std::string> & names)
{
   // pull out all names for this set
   for(typename SetType::const_iterator itr=setData.begin();itr!=setData.end();++itr) {
      Ioss::GroupingEntity * entity = *itr;
      names.push_back(entity->name());
   }
}

void STK_ExodusReaderFactory::registerSidesets(STK_Interface & mesh) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      const stk::mesh::PartVector & subsets = part->subsets();
      shards::CellTopology cellTopo = stk::mesh::get_cell_topology(metaData->get_topology(*part));
      const CellTopologyData * ct = cellTopo.getCellTopologyData();

      // if a side part ==> this is a sideset: now storage is recursive
      // on part contains all sub parts with consistent topology
      if(part->primary_entity_rank()==mesh.getSideRank() && ct==0 && subsets.size()>0) {
         TEUCHOS_TEST_FOR_EXCEPTION(subsets.size()!=1,std::runtime_error,
                            "STK_ExodusReaderFactory::registerSidesets error - part \"" << part->name() <<
                            "\" has more than one subset");

         // grab cell topology and name of subset part
         const stk::mesh::Part * ss_part = subsets[0];
         shards::CellTopology ss_cellTopo = stk::mesh::get_cell_topology(metaData->get_topology(*ss_part));
         const CellTopologyData * ss_ct = ss_cellTopo.getCellTopologyData();

         // only add subset parts that have no topology
         if(ss_ct!=0)
            mesh.addSideset(part->name(),ss_ct);
      }
   }
}

void STK_ExodusReaderFactory::registerNodesets(STK_Interface & mesh) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::mesh::PartVector & parts = metaData->get_parts();

   stk::mesh::PartVector::const_iterator partItr;
   for(partItr=parts.begin();partItr!=parts.end();++partItr) {
      const stk::mesh::Part * part = *partItr;
      shards::CellTopology cellTopo = stk::mesh::get_cell_topology(metaData->get_topology(*part));
      const CellTopologyData * ct = cellTopo.getCellTopologyData();

      // if a side part ==> this is a sideset: now storage is recursive
      // on part contains all sub parts with consistent topology
      if(part->primary_entity_rank()==mesh.getNodeRank() && ct==0) {

         // only add subset parts that have no topology
         if(part->name()!=STK_Interface::nodesString)
            mesh.addNodeset(part->name());
      }
   }
}

void STK_ExodusReaderFactory::registerEdgeBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   /* For each edge block in the exodus file, check it's topology
    * against the list of edge topologies for each element block.
    * If it matches, add the edge block for that element block.
    * This will add the edge block as a subset of the element
    * block in the STK mesh.
    */
   const Ioss::EdgeBlockContainer & edge_blocks = meshData.get_input_ioss_region()->get_edge_blocks();
   for(Ioss::EdgeBlockContainer::const_iterator ebc_iter=edge_blocks.begin();ebc_iter!=edge_blocks.end();++ebc_iter) {
      Ioss::GroupingEntity * entity = *ebc_iter;
      const stk::mesh::Part * edgeBlockPart = metaData->get_part(entity->name());
      const stk::topology edgeBlockTopo = metaData->get_topology(*edgeBlockPart);

      for (auto ebuet_iter : elemBlockUniqueEdgeTopologies_) {
         std::string elemBlockName = ebuet_iter.first;
         std::vector<stk::topology> uniqueEdgeTopologies = ebuet_iter.second;

         auto find_result = std::find(uniqueEdgeTopologies.begin(), uniqueEdgeTopologies.end(), edgeBlockTopo);
         if (find_result != uniqueEdgeTopologies.end()) {
            mesh.addEdgeBlock(elemBlockName, edgeBlockPart->name(), edgeBlockTopo);
         }
      }
   }
}

void STK_ExodusReaderFactory::registerFaceBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const
{
   using Teuchos::RCP;

   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   /* For each face block in the exodus file, check it's topology
    * against the list of face topologies for each element block.
    * If it matches, add the face block for that element block.
    * This will add the face block as a subset of the element
    * block in the STK mesh.
    */
   const Ioss::FaceBlockContainer & face_blocks = meshData.get_input_ioss_region()->get_face_blocks();
   for(Ioss::FaceBlockContainer::const_iterator fbc_itr=face_blocks.begin();fbc_itr!=face_blocks.end();++fbc_itr) {
      Ioss::GroupingEntity * entity = *fbc_itr;
      const stk::mesh::Part * faceBlockPart = metaData->get_part(entity->name());
      const stk::topology faceBlockTopo = metaData->get_topology(*faceBlockPart);

      for (auto ebuft_iter : elemBlockUniqueFaceTopologies_) {
         std::string elemBlockName = ebuft_iter.first;
         std::vector<stk::topology> uniqueFaceTopologies = ebuft_iter.second;

         auto find_result = std::find(uniqueFaceTopologies.begin(), uniqueFaceTopologies.end(), faceBlockTopo);
         if (find_result != uniqueFaceTopologies.end()) {
            mesh.addFaceBlock(elemBlockName, faceBlockPart->name(), faceBlockTopo);
         }
      }
   }
}

bool topo_less (stk::topology &i,stk::topology &j) { return (i.value() < j.value()); }
bool topo_equal (stk::topology &i,stk::topology &j) { return (i.value() == j.value()); }

void STK_ExodusReaderFactory::createUniqueEdgeTopologyMap(STK_Interface & mesh, const stk::mesh::Part *elemBlockPart) const
{
   using Teuchos::RCP;

   /* For a given element block, add it's edge topologies to a vector,
    * sort it, dedupe it and save it to the "unique topo" map.
    */
   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::topology elemBlockTopo = metaData->get_topology(*elemBlockPart);

   std::vector<stk::topology> edge_topologies;
   for (unsigned i=0;i<elemBlockTopo.num_edges();i++) {
       edge_topologies.push_back(elemBlockTopo.edge_topology(i));
   }
   std::sort(edge_topologies.begin(), edge_topologies.end(), topo_less);
   std::vector<stk::topology>::iterator new_end;
   new_end = std::unique(edge_topologies.begin(), edge_topologies.end(), topo_equal);
   edge_topologies.resize( std::distance(edge_topologies.begin(),new_end) );

   elemBlockUniqueEdgeTopologies_[elemBlockPart->name()] = edge_topologies;
}

void STK_ExodusReaderFactory::createUniqueFaceTopologyMap(STK_Interface & mesh, const stk::mesh::Part *elemBlockPart) const
{
   using Teuchos::RCP;

   /* For a given element block, add it's face topologies to a vector,
    * sort it, dedupe it and save it to the "unique topo" map.
    */
   RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   const stk::topology elemBlockTopo = metaData->get_topology(*elemBlockPart);

   std::vector<stk::topology> face_topologies;
   for (unsigned i=0;i<elemBlockTopo.num_faces();i++) {
      face_topologies.push_back(elemBlockTopo.face_topology(i));
   }
   std::sort(face_topologies.begin(), face_topologies.end(), topo_less);
   std::vector<stk::topology>::iterator new_end;
   new_end = std::unique(face_topologies.begin(), face_topologies.end(), topo_equal);
   face_topologies.resize( std::distance(face_topologies.begin(),new_end) );

   elemBlockUniqueFaceTopologies_[elemBlockPart->name()] = face_topologies;
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void STK_ExodusReaderFactory::addEdgeBlocks(STK_Interface & mesh) const
{
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   /* For each element block, iterate over it's edge topologies.
    * For each edge topology, get the matching edge block and
    * add all edges of that topology to the edge block.
    */
   for (auto iter : elemBlockUniqueEdgeTopologies_) {
      std::string elemBlockName = iter.first;
      std::vector<stk::topology> uniqueEdgeTopologies = iter.second;

      for (auto topo : uniqueEdgeTopologies ) {
         const stk::mesh::Part * elemBlockPart = metaData->get_part(elemBlockName);
         const stk::mesh::Part & edgeTopoPart  = metaData->get_topology_root_part(topo);

         stk::mesh::Selector owned_block;
         owned_block  = *elemBlockPart;
         owned_block &= edgeTopoPart;
         owned_block &= metaData->locally_owned_part();

         std::string edge_block_name = mkBlockName(panzer_stk::STK_Interface::edgeBlockString, topo.name());
         stk::mesh::Part * edge_block = mesh.getEdgeBlock(edge_block_name);

         std::vector<stk::mesh::Entity> all_edges_for_topo;
         bulkData->get_entities(mesh.getEdgeRank(),owned_block,all_edges_for_topo);
         mesh.addEntitiesToEdgeBlock(all_edges_for_topo, edge_block);
      }
   }
}

// Pre-Condition: call beginModification() before entry
// Post-Condition: call endModification() after exit
void STK_ExodusReaderFactory::addFaceBlocks(STK_Interface & mesh) const
{
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();

   /* For each element block, iterate over it's face topologies.
    * For each face topology, get the matching face block and
    * add all faces of that topology to the face block.
    */
   for (auto iter : elemBlockUniqueFaceTopologies_) {
      std::string elemBlockName = iter.first;
      std::vector<stk::topology> uniqueFaceTopologies = iter.second;

      for (auto topo : uniqueFaceTopologies ) {
         const stk::mesh::Part * elemBlockPart = metaData->get_part(elemBlockName);
         const stk::mesh::Part & faceTopoPart  = metaData->get_topology_root_part(topo);

         stk::mesh::Selector owned_block;
         owned_block  = *elemBlockPart;
         owned_block &= faceTopoPart;
         owned_block &= metaData->locally_owned_part();

         std::string face_block_name = mkBlockName(panzer_stk::STK_Interface::faceBlockString, topo.name());
         stk::mesh::Part * face_block = mesh.getFaceBlock(face_block_name);

         std::vector<stk::mesh::Entity> all_faces_for_topo;
         bulkData->get_entities(mesh.getFaceRank(),owned_block,all_faces_for_topo);
         mesh.addEntitiesToFaceBlock(all_faces_for_topo, face_block);
      }
   }
}

void STK_ExodusReaderFactory::buildMetaData(stk::ParallelMachine /* parallelMach */, STK_Interface & mesh) const
{
   if (createEdgeBlocks_) {
      /* For each element block, iterate over it's edge topologies.
       * For each edge topology, create an edge block for that topology.
       */
      for (auto iter : elemBlockUniqueEdgeTopologies_) {
         std::string elemBlockName = iter.first;
         std::vector<stk::topology> uniqueEdgeTopologies = iter.second;

         for (auto topo : uniqueEdgeTopologies ) {
            std::string edge_block_name = mkBlockName(panzer_stk::STK_Interface::edgeBlockString, topo.name());
            mesh.addEdgeBlock(elemBlockName, edge_block_name, topo);
         }
      }
   }
   if (createFaceBlocks_ && mesh.getMetaData()->spatial_dimension() > 2) {
      /* For each element block, iterate over it's face topologies.
       * For each face topology, create a face block for that topology.
       */
      for (auto iter : elemBlockUniqueFaceTopologies_) {
         std::string elemBlockName = iter.first;
         std::vector<stk::topology> uniqueFaceTopologies = iter.second;

         for (auto topo : uniqueFaceTopologies ) {
            std::string face_block_name = mkBlockName(panzer_stk::STK_Interface::faceBlockString, topo.name());
            mesh.addFaceBlock(elemBlockName, face_block_name, topo);
         }
      }
   }
}

bool STK_ExodusReaderFactory::doPerceptRefinement() const
{
  return (fileType_!="Exodus Refinement") && (levelsOfRefinement_ > 0);
}

std::string STK_ExodusReaderFactory::mkBlockName(std::string base, std::string topo_name) const
{
   std::string name;
   name = topo_name+"_"+base;
   std::transform(name.begin(), name.end(), name.begin(),
                  [](const char c)
                  { return char(std::tolower(c)); });
   return name;
}

}

#endif
