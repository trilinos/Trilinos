// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <PanzerAdaptersSTK_config.hpp>
#include <Panzer_STK_Interface.hpp>

#include <Teuchos_as.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <limits>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>

// #include <stk_rebalance/Rebalance.hpp>
// #include <stk_rebalance/Partition.hpp>
// #include <stk_rebalance/ZoltanPartition.hpp>
// #include <stk_rebalance_utils/RebalanceUtils.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#ifdef PANZER_HAVE_IOSS
#include <Ionit_Initializer.h>
#include <stk_io/IossBridge.hpp>
#include <stk_io/WriteMesh.hpp>
#endif

#ifdef PANZER_HAVE_PERCEPT
#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#endif

#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include <set>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

ElementDescriptor::ElementDescriptor() {}
ElementDescriptor::ElementDescriptor(stk::mesh::EntityId gid,const std::vector<stk::mesh::EntityId> & nodes)
   : gid_(gid), nodes_(nodes) {}
ElementDescriptor::~ElementDescriptor() {}

/** Constructor function for building the element descriptors.
  */
Teuchos::RCP<ElementDescriptor>
buildElementDescriptor(stk::mesh::EntityId elmtId,std::vector<stk::mesh::EntityId> & nodes)
{
   return Teuchos::rcp(new ElementDescriptor(elmtId,nodes));
}

const std::string STK_Interface::coordsString = "coordinates";
const std::string STK_Interface::nodesString = "nodes";
const std::string STK_Interface::edgesString = "edges";
const std::string STK_Interface::facesString = "faces";
const std::string STK_Interface::edgeBlockString = "edge_block";
const std::string STK_Interface::faceBlockString = "face_block";

STK_Interface::STK_Interface()
   : dimension_(0), initialized_(false), currentLocalId_(0), initialStateTime_(0.0), currentStateTime_(0.0), useFieldCoordinates_(false)
{
  metaData_ = rcp(new stk::mesh::MetaData());
  metaData_->use_simple_fields();
}

STK_Interface::STK_Interface(Teuchos::RCP<stk::mesh::MetaData> metaData)
  : dimension_(0), initialized_(false), currentLocalId_(0), initialStateTime_(0.0), currentStateTime_(0.0), useFieldCoordinates_(false)
{
  metaData_ = metaData;
  metaData_->use_simple_fields();
}

STK_Interface::STK_Interface(unsigned dim)
   : dimension_(dim), initialized_(false), currentLocalId_(0), useFieldCoordinates_(false)
{
   std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
   entity_rank_names.push_back("FAMILY_TREE");

   metaData_ = rcp(new stk::mesh::MetaData(dimension_,entity_rank_names));
   metaData_->use_simple_fields();

   initializeFromMetaData();
}

void STK_Interface::addSideset(const std::string & name,const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk::mesh::Part * sideset = metaData_->get_part(name);
   if(sideset==nullptr)
      sideset = &metaData_->declare_part_with_topology(name,
         stk::mesh::get_topology(shards::CellTopology(ctData), dimension_));
   sidesets_.insert(std::make_pair(name,sideset));
}

void STK_Interface::addNodeset(const std::string & name)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk::mesh::Part * nodeset = metaData_->get_part(name);
   if(nodeset==nullptr) {
      const CellTopologyData * ctData = shards::getCellTopologyData<shards::Node>();
      nodeset = &metaData_->declare_part_with_topology(name,
         stk::mesh::get_topology(shards::CellTopology(ctData), dimension_));
   }
   nodesets_.insert(std::make_pair(name,nodeset));
}

void STK_Interface::addSolutionField(const std::string & fieldName,const std::string & blockId)
{
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToSolution_.find(key)==fieldNameToSolution_.end()) {
      SolutionFieldType * field = metaData_->get_field<double>(stk::topology::NODE_RANK, fieldName);
      if(field==0)
         field = &metaData_->declare_field<double>(stk::topology::NODE_RANK, fieldName);
      if ( initialized_ )  {
        metaData_->enable_late_fields();
        stk::mesh::put_field_on_mesh(*field, metaData_->universal_part(), nullptr);
      }
      fieldNameToSolution_[key] = field;
   }
}

void STK_Interface::addCellField(const std::string & fieldName,const std::string & blockId)
{
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToCellField_.find(key)==fieldNameToCellField_.end()) {
      SolutionFieldType * field = metaData_->get_field<double>(stk::topology::ELEMENT_RANK, fieldName);
      if(field==0)
         field = &metaData_->declare_field<double>(stk::topology::ELEMENT_RANK, fieldName);

      if ( initialized_ )  {
        metaData_->enable_late_fields();
        stk::mesh::put_field_on_mesh(*field, metaData_->universal_part(), nullptr);
      }
      fieldNameToCellField_[key] = field;
   }
}

void STK_Interface::addEdgeField(const std::string & fieldName,const std::string & blockId)
{
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToEdgeField_.find(key)==fieldNameToEdgeField_.end()) {
      SolutionFieldType * field = metaData_->get_field<double>(stk::topology::EDGE_RANK, fieldName);
      if(field==0) {
         field = &metaData_->declare_field<double>(stk::topology::EDGE_RANK, fieldName);
      }

      if ( initialized_ )  {
        metaData_->enable_late_fields();
        stk::mesh::put_field_on_mesh(*field, metaData_->universal_part(), nullptr);
      }
      fieldNameToEdgeField_[key] = field;
   }
}

void STK_Interface::addFaceField(const std::string & fieldName,const std::string & blockId)
{
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToFaceField_.find(key)==fieldNameToFaceField_.end()) {
      SolutionFieldType * field = metaData_->get_field<double>(stk::topology::FACE_RANK, fieldName);
      if(field==0) {
         field = &metaData_->declare_field<double>(stk::topology::FACE_RANK, fieldName);
      }

      if ( initialized_ )  {
        metaData_->enable_late_fields();
        stk::mesh::put_field_on_mesh(*field, metaData_->universal_part(), nullptr);
      }
      fieldNameToFaceField_[key] = field;
   }
}

void STK_Interface::addMeshCoordFields(const std::string & blockId,
                                       const std::vector<std::string> & coordNames,
                                       const std::string & dispPrefix)
{
   TEUCHOS_ASSERT(dimension_!=0);
   TEUCHOS_ASSERT(dimension_==coordNames.size());
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");

   // we only allow one alternative coordinate field
   TEUCHOS_TEST_FOR_EXCEPTION(meshCoordFields_.find(blockId)!=meshCoordFields_.end(),std::invalid_argument,
                              "STK_Interface::addMeshCoordFields: Can't set more than one set of coordinate "
                              "fields for element block \""+blockId+"\".");

   // Note that there is a distinction between the key which is used for lookups
   // and the field that lives on the mesh, which is used for printing the displacement.

   // just copy the coordinate names
   meshCoordFields_[blockId] = coordNames;

   // must fill in the displacement fields
   std::vector<std::string> & dispFields = meshDispFields_[blockId];
   dispFields.resize(dimension_);

   for(unsigned i=0;i<dimension_;i++) {
      std::pair<std::string,std::string> key = std::make_pair(coordNames[i],blockId);
      std::string dispName = dispPrefix+coordNames[i];

      dispFields[i] = dispName; // record this field as a
                                // displacement field

      // add & declare field if not already added...currently assuming linears
      if(fieldNameToSolution_.find(key)==fieldNameToSolution_.end()) {

         SolutionFieldType * field = metaData_->get_field<double>(stk::topology::NODE_RANK, dispName);
         if(field==0) {
            field = &metaData_->declare_field<double>(stk::topology::NODE_RANK, dispName);
         }
         fieldNameToSolution_[key] = field;
      }
   }
}

void STK_Interface::addInformationRecords(const std::vector<std::string> & info_records)
{
   informationRecords_.insert(info_records.begin(), info_records.end());
}

void STK_Interface::initialize(stk::ParallelMachine parallelMach,bool setupIO,
                               const bool buildRefinementSupport)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0); // no zero dimensional meshes!

   if(mpiComm_==Teuchos::null)
      mpiComm_ = getSafeCommunicator(parallelMach);

   procRank_ = stk::parallel_machine_rank(*mpiComm_->getRawMpiComm());

   // associating the field with a part: universal part!
   stk::mesh::put_field_on_mesh( *coordinatesField_ , metaData_->universal_part(), getDimension(), nullptr);
   stk::mesh::put_field_on_mesh( *edgesField_ , metaData_->universal_part(), getDimension(), nullptr);
   if (dimension_ > 2)
     stk::mesh::put_field_on_mesh( *facesField_ , metaData_->universal_part(), getDimension(), nullptr);
   stk::mesh::put_field_on_mesh( *processorIdField_ , metaData_->universal_part(), nullptr);
   stk::mesh::put_field_on_mesh( *loadBalField_ , metaData_->universal_part(), nullptr);

   initializeFieldsInSTK(fieldNameToSolution_, setupIO);
   initializeFieldsInSTK(fieldNameToCellField_, setupIO);
   initializeFieldsInSTK(fieldNameToEdgeField_, setupIO);
   initializeFieldsInSTK(fieldNameToFaceField_, setupIO);

#ifdef PANZER_HAVE_IOSS
   if(setupIO) {
      // setup Exodus file IO
      /////////////////////////////////////////

      // add element blocks
      {
         std::map<std::string, stk::mesh::Part*>::iterator itr;
         for(itr=elementBlocks_.begin();itr!=elementBlocks_.end();++itr)
            if(!stk::io::is_part_io_part(*itr->second))
               stk::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add edge blocks
      {
         std::map<std::string, stk::mesh::Part*>::iterator itr;
         for(itr=edgeBlocks_.begin();itr!=edgeBlocks_.end();++itr)
            if(!stk::io::is_part_edge_block_io_part(*itr->second))
               stk::io::put_edge_block_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add face blocks
      {
         std::map<std::string, stk::mesh::Part*>::iterator itr;
         for(itr=faceBlocks_.begin();itr!=faceBlocks_.end();++itr)
            if(!stk::io::is_part_face_block_io_part(*itr->second))
               stk::io::put_face_block_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add side sets
      {
         std::map<std::string, stk::mesh::Part*>::iterator itr;
         for(itr=sidesets_.begin();itr!=sidesets_.end();++itr)
            if(!stk::io::is_part_io_part(*itr->second))
               stk::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add node sets
      {
         std::map<std::string, stk::mesh::Part*>::iterator itr;
         for(itr=nodesets_.begin();itr!=nodesets_.end();++itr)
            if(!stk::io::is_part_io_part(*itr->second))
               stk::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add nodes
      if(!stk::io::is_part_io_part(*nodesPart_))
         stk::io::put_io_part_attribute(*nodesPart_);

      stk::io::set_field_role(*coordinatesField_, Ioss::Field::MESH);
      stk::io::set_field_output_type(*coordinatesField_, stk::io::FieldOutputType::VECTOR_3D);
      stk::io::set_field_role(*edgesField_, Ioss::Field::MESH);
      stk::io::set_field_output_type(*edgesField_, stk::io::FieldOutputType::VECTOR_3D);
      if (dimension_ > 2) {
         stk::io::set_field_role(*facesField_, Ioss::Field::MESH);
         stk::io::set_field_output_type(*facesField_, stk::io::FieldOutputType::VECTOR_3D);
      }
      stk::io::set_field_role(*processorIdField_, Ioss::Field::TRANSIENT);
      // stk::io::set_field_role(*loadBalField_, Ioss::Field::TRANSIENT);

   }
#endif

   if (buildRefinementSupport) {
#ifdef PANZER_HAVE_PERCEPT
     refinedMesh_ = Teuchos::rcp(new percept::PerceptMesh(this->getMetaData().get(),
                                                          this->getBulkData().get(),
                                                          true));

     breakPattern_ = Teuchos::rcp(new percept::URP_Heterogeneous_3D(*refinedMesh_));
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "ERROR: Uniform refinement requested. This requires the Percept package to be enabled in Trilinos!");
#endif
   }

   if(bulkData_==Teuchos::null)
      instantiateBulkData(*mpiComm_->getRawMpiComm());

   metaData_->commit();

   initialized_ = true;
}

void STK_Interface::initializeFieldsInSTK(const std::map<std::pair<std::string,std::string>,SolutionFieldType*> & nameToField,
                                          bool setupIO)
{
   std::set<SolutionFieldType*> uniqueFields;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*>::const_iterator fieldIter;
   for(fieldIter=nameToField.begin();fieldIter!=nameToField.end();++fieldIter)
      uniqueFields.insert(fieldIter->second); // this makes setting up IO easier!

   {
      std::set<SolutionFieldType*>::const_iterator uniqueFieldIter;
      for(uniqueFieldIter=uniqueFields.begin();uniqueFieldIter!=uniqueFields.end();++uniqueFieldIter)
        stk::mesh::put_field_on_mesh(*(*uniqueFieldIter),metaData_->universal_part(),nullptr);
   }

#ifdef PANZER_HAVE_IOSS
   if(setupIO) {
      // add solution fields
      std::set<SolutionFieldType*>::const_iterator uniqueFieldIter;
      for(uniqueFieldIter=uniqueFields.begin();uniqueFieldIter!=uniqueFields.end();++uniqueFieldIter)
         stk::io::set_field_role(*(*uniqueFieldIter), Ioss::Field::TRANSIENT);
   }
#endif
}

void STK_Interface::instantiateBulkData(stk::ParallelMachine parallelMach)
{
   TEUCHOS_ASSERT(bulkData_==Teuchos::null);
   if(mpiComm_==Teuchos::null)
      mpiComm_ = getSafeCommunicator(parallelMach);

   std::unique_ptr<stk::mesh::BulkData> bulkUPtr = stk::mesh::MeshBuilder(*mpiComm_->getRawMpiComm()).create(Teuchos::get_shared_ptr(metaData_));
   bulkData_ = rcp(bulkUPtr.release());
}

void STK_Interface::beginModification()
{
   TEUCHOS_TEST_FOR_EXCEPTION(bulkData_==Teuchos::null,std::logic_error,
                      "STK_Interface: Must call \"initialized\" or \"instantiateBulkData\" before \"beginModification\"");

   bulkData_->modification_begin();
}

void STK_Interface::endModification()
{
   TEUCHOS_TEST_FOR_EXCEPTION(bulkData_==Teuchos::null,std::logic_error,
                      "STK_Interface: Must call \"initialized\" or \"instantiateBulkData\" before \"endModification\"");

   // TODO: Resolving sharing this way comes at a cost in performance. The STK
   // team has decided that users need to declare their own sharing. We should
   // find where shared entities are being created in Panzer and declare it.
   // Once this is done, the extra code below can be deleted.

    stk::CommSparse comm(bulkData_->parallel());

    for (int phase=0;phase<2;++phase) {
      for (int i=0;i<bulkData_->parallel_size();++i) {
        if ( i != bulkData_->parallel_rank() ) {
          const stk::mesh::BucketVector& buckets = bulkData_->buckets(stk::topology::NODE_RANK);
          for (size_t j=0;j<buckets.size();++j) {
            const stk::mesh::Bucket& bucket = *buckets[j];
            if ( bucket.owned() ) {
              for (size_t k=0;k<bucket.size();++k) {
                stk::mesh::EntityKey key = bulkData_->entity_key(bucket[k]);
                comm.send_buffer(i).pack<stk::mesh::EntityKey>(key);
              }
            }
          }
        }
      }

      if (phase == 0 ) {
        comm.allocate_buffers();
      }
      else {
        comm.communicate();
      }
    }

    for (int i=0;i<bulkData_->parallel_size();++i) {
      if ( i != bulkData_->parallel_rank() ) {
        while(comm.recv_buffer(i).remaining()) {
          stk::mesh::EntityKey key;
          comm.recv_buffer(i).unpack<stk::mesh::EntityKey>(key);
          stk::mesh::Entity node = bulkData_->get_entity(key);
          if ( bulkData_->is_valid(node) ) {
            bulkData_->add_node_sharing(node, i);
          }
        }
      }
    }


    bulkData_->modification_end();

    buildEntityCounts();
    buildMaxEntityIds();
}

void STK_Interface::addNode(stk::mesh::EntityId gid, const std::vector<double> & coord)
{
   TEUCHOS_TEST_FOR_EXCEPTION(not isModifiable(),std::logic_error,
                      "STK_Interface::addNode: STK_Interface must be modifiable to add a node");
   TEUCHOS_TEST_FOR_EXCEPTION(coord.size()!=getDimension(),std::logic_error,
                      "STK_Interface::addNode: number of coordinates in vector must mation dimension");
   TEUCHOS_TEST_FOR_EXCEPTION(gid==0,std::logic_error,
                      "STK_Interface::addNode: STK has STUPID restriction of no zero GIDs, pick something else");
   stk::mesh::EntityRank nodeRank = getNodeRank();

   stk::mesh::Entity node = bulkData_->declare_entity(nodeRank,gid,nodesPartVec_);

   // set coordinate vector
   double * fieldCoords = stk::mesh::field_data(*coordinatesField_,node);
   for(std::size_t i=0;i<coord.size();++i)
      fieldCoords[i] = coord[i];
}

void STK_Interface::addEntityToSideset(stk::mesh::Entity entity,stk::mesh::Part * sideset)
{
   std::vector<stk::mesh::Part*> sidesetV;
   sidesetV.push_back(sideset);

   bulkData_->change_entity_parts(entity,sidesetV);
}

void STK_Interface::addEntityToNodeset(stk::mesh::Entity entity,stk::mesh::Part * nodeset)
{
   std::vector<stk::mesh::Part*> nodesetV;
   nodesetV.push_back(nodeset);

   bulkData_->change_entity_parts(entity,nodesetV);
}

void STK_Interface::addEntityToEdgeBlock(stk::mesh::Entity entity,stk::mesh::Part * edgeblock)
{
   std::vector<stk::mesh::Part*> edgeblockV;
   edgeblockV.push_back(edgeblock);

   bulkData_->change_entity_parts(entity,edgeblockV);
}
void STK_Interface::addEntitiesToEdgeBlock(std::vector<stk::mesh::Entity> entities,stk::mesh::Part * edgeblock)
{
   if (entities.size() > 0) {
      std::vector<stk::mesh::Part*> edgeblockV;
      edgeblockV.push_back(edgeblock);

      bulkData_->change_entity_parts(entities,edgeblockV);
   }
}

void STK_Interface::addEntityToFaceBlock(stk::mesh::Entity entity,stk::mesh::Part * faceblock)
{
   std::vector<stk::mesh::Part*> faceblockV;
   faceblockV.push_back(faceblock);

   bulkData_->change_entity_parts(entity,faceblockV);
}
void STK_Interface::addEntitiesToFaceBlock(std::vector<stk::mesh::Entity> entities,stk::mesh::Part * faceblock)
{
   if (entities.size() > 0) {
      std::vector<stk::mesh::Part*> faceblockV;
      faceblockV.push_back(faceblock);

      bulkData_->change_entity_parts(entities,faceblockV);
   }
}

void STK_Interface::addElement(const Teuchos::RCP<ElementDescriptor> & ed,stk::mesh::Part * block)
{
   std::vector<stk::mesh::Part*> blockVec;
   blockVec.push_back(block);

   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::EntityRank nodeRank = getNodeRank();
   stk::mesh::Entity element = bulkData_->declare_entity(elementRank,ed->getGID(),blockVec);

   // build relations that give the mesh structure
   const std::vector<stk::mesh::EntityId> & nodes = ed->getNodes();
   for(std::size_t i=0;i<nodes.size();++i) {
      // add element->node relation
      stk::mesh::Entity node = bulkData_->get_entity(nodeRank,nodes[i]);
      TEUCHOS_ASSERT(bulkData_->is_valid(node));
      bulkData_->declare_relation(element,node,i);
   }

   ProcIdData * procId = stk::mesh::field_data(*processorIdField_,element);
   procId[0] = Teuchos::as<ProcIdData>(procRank_);
}


void STK_Interface::addEdges()
{
   // loop over elements
   stk::mesh::EntityRank edgeRank = getEdgeRank();
   std::vector<stk::mesh::Entity> localElmts;
   getMyElements(localElmts);
   std::vector<stk::mesh::Entity>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
     stk::mesh::Entity element = (*itr);
     stk::mesh::EntityId gid = bulkData_->identifier(element);
     std::vector<stk::mesh::EntityId> subcellIds;
     getSubcellIndices(edgeRank,gid,subcellIds);

     for(std::size_t i=0;i<subcellIds.size();++i) {
       stk::mesh::Entity edge = bulkData_->get_entity(edgeRank,subcellIds[i]);
       stk::mesh::Entity const* relations = bulkData_->begin_nodes(edge);

       double * node_coord_1 = stk::mesh::field_data(*coordinatesField_,relations[0]);
       double * node_coord_2 = stk::mesh::field_data(*coordinatesField_,relations[1]);

       // set coordinate vector
       double * edgeCoords = stk::mesh::field_data(*edgesField_,edge);
       for(std::size_t j=0;j<getDimension();++j)
          edgeCoords[j] = (node_coord_1[j]+node_coord_2[j])/2.0;
     }
   }
}

void STK_Interface::addFaces()
{
  // loop over elements
  stk::mesh::EntityRank faceRank = getFaceRank();
  std::vector<stk::mesh::Entity> localElmts;
  getMyElements(localElmts);
  std::vector<stk::mesh::Entity>::const_iterator itr;
  for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
    stk::mesh::Entity element = (*itr);
    stk::mesh::EntityId gid = elementGlobalId(element);
    std::vector<stk::mesh::EntityId> subcellIds;
    getSubcellIndices(faceRank,gid,subcellIds);

    for(std::size_t i=0;i<subcellIds.size();++i) {
      stk::mesh::Entity face = bulkData_->get_entity(faceRank,subcellIds[i]);
      stk::mesh::Entity const* relations = bulkData_->begin_nodes(face);
      const size_t num_relations = bulkData_->num_nodes(face);

      // set coordinate vector
      double * faceCoords = stk::mesh::field_data(*facesField_,face);
      for(std::size_t j=0;j<getDimension();++j){
        faceCoords[j] = 0.0;
        for(std::size_t k=0;k<num_relations;++k)
          faceCoords[j] += stk::mesh::field_data(*coordinatesField_,relations[k])[j];
        faceCoords[j] /= double(num_relations);
      }
    }
  }
}

stk::mesh::Entity STK_Interface::findConnectivityById(stk::mesh::Entity src, stk::mesh::EntityRank tgt_rank, unsigned rel_id) const
{
  const size_t num_rels = bulkData_->num_connectivity(src, tgt_rank);
  stk::mesh::Entity const* relations = bulkData_->begin(src, tgt_rank);
  stk::mesh::ConnectivityOrdinal const* ordinals = bulkData_->begin_ordinals(src, tgt_rank);
  for (size_t i = 0; i < num_rels; ++i) {
    if (ordinals[i] == static_cast<stk::mesh::ConnectivityOrdinal>(rel_id)) {
      return relations[i];
    }
  }

  return stk::mesh::Entity();
}

///////////////////////////////////////////////////////////////////////////////
//
//  writeToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
writeToExodus(const std::string& filename,
              const bool append)
{
  setupExodusFile(filename,append);
  writeToExodus(0.0);
} // end of writeToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  setupExodusFile()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
setupExodusFile(const std::string& filename,
                const bool append,
                const bool append_after_restart_time,
                const double restart_time)
{
  using std::runtime_error;
  using stk::io::StkMeshIoBroker;
  using stk::mesh::FieldVector;
  using stk::ParallelMachine;
  using Teuchos::rcp;
  PANZER_FUNC_TIME_MONITOR("STK_Interface::setupExodusFile(filename)");
#ifdef PANZER_HAVE_IOSS
  TEUCHOS_ASSERT(not mpiComm_.is_null())

  ParallelMachine comm = *mpiComm_->getRawMpiComm();
  meshData_ = rcp(new StkMeshIoBroker(comm));
  meshData_->set_bulk_data(Teuchos::get_shared_ptr(bulkData_));
  Ioss::PropertyManager props;
  props.add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", "FALSE"));
  if (append) {
    if (append_after_restart_time) {
      meshIndex_ = meshData_->create_output_mesh(filename, stk::io::APPEND_RESULTS,
                                                 props, restart_time);
    }
    else // Append results to the end of the file
      meshIndex_ = meshData_->create_output_mesh(filename, stk::io::APPEND_RESULTS, props);
  }
  else
    meshIndex_ = meshData_->create_output_mesh(filename, stk::io::WRITE_RESULTS, props);
  const FieldVector& fields = metaData_->get_fields();
  for (size_t i(0); i < fields.size(); ++i) {
    // Do NOT add MESH type stk fields to exodus io, but do add everything
    // else. This allows us to avoid having to catch a throw for
    // re-registering coordinates, sidesets, etc... Note that some
    // fields like LOAD_BAL don't always have a role assigned, so for
    // roles that point to null, we need to register them as well.
    auto role = stk::io::get_field_role(*fields[i]);
    if (role != nullptr) {
      if (*role != Ioss::Field::MESH)
        meshData_->add_field(meshIndex_, *fields[i]);
    } else {
      meshData_->add_field(meshIndex_, *fields[i]);
    }
  }

  // convert the set to a vector
  std::vector<std::string> deduped_info_records(informationRecords_.begin(),informationRecords_.end());
  meshData_->add_info_records(meshIndex_, deduped_info_records);
#else
  TEUCHOS_ASSERT(false)
#endif
} // end of setupExodusFile()

///////////////////////////////////////////////////////////////////////////////
//
//  writeToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
writeToExodus(
  double timestep)
{
  using std::cout;
  using std::endl;
  using Teuchos::FancyOStream;
  using Teuchos::rcpFromRef;
  PANZER_FUNC_TIME_MONITOR("STK_Interface::writeToExodus(timestep)");
#ifdef PANZER_HAVE_IOSS
  if (not meshData_.is_null())
  {
    currentStateTime_ = timestep;
    globalToExodus(GlobalVariable::ADD);
    meshData_->begin_output_step(meshIndex_, currentStateTime_);
    meshData_->write_defined_output_fields(meshIndex_);
    globalToExodus(GlobalVariable::WRITE);
    meshData_->end_output_step(meshIndex_);
  }
  else // if (meshData_.is_null())
  {
    FancyOStream out(rcpFromRef(cout));
    out.setOutputToRootOnly(0);
    out << "WARNING:  Exodus I/O has been disabled or not setup properly; "
        << "not writing to Exodus." << endl;
  } // end if (meshData_.is_null()) or not
#else
  TEUCHOS_ASSERT(false);
#endif
} // end of writeToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  globalToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
globalToExodus(
  const GlobalVariable& flag)
{
  using std::invalid_argument;
  using std::string;
  using Teuchos::Array;

  // Loop over all the global variables to be added to the Exodus output file.
  // For each global variable, we determine the data type, and then add or
  // write it accordingly, depending on the value of flag.
  for (auto i = globalData_.begin(); i != globalData_.end(); ++i)
  {
    const string& name = globalData_.name(i);

    // Integers.
    if (globalData_.isType<int>(name))
    {
      const auto& value = globalData_.get<int>(name);
      if (flag == GlobalVariable::ADD)
      {
        try
        {
          meshData_->add_global(meshIndex_, name, value,
            stk::util::ParameterType::INTEGER);
        }
        catch (...)
        {
          return;
        }
      }
      else // if (flag == GlobalVariable::WRITE)
        meshData_->write_global(meshIndex_, name, value,
          stk::util::ParameterType::INTEGER);
    }

    // Doubles.
    else if (globalData_.isType<double>(name))
    {
      const auto& value = globalData_.get<double>(name);
      if (flag == GlobalVariable::ADD)
      {
        try
        {
          meshData_->add_global(meshIndex_, name, value,
            stk::util::ParameterType::DOUBLE);
        }
        catch (...)
        {
          return;
        }
      }
      else // if (flag == GlobalVariable::WRITE)
        meshData_->write_global(meshIndex_, name, value,
          stk::util::ParameterType::DOUBLE);
    }

    // Vectors of integers.
    else if (globalData_.isType<Array<int>>(name))
    {
      const auto& value = globalData_.get<Array<int>>(name).toVector();
      if (flag == GlobalVariable::ADD)
      {
        try
        {
          meshData_->add_global(meshIndex_, name, value,
            stk::util::ParameterType::INTEGERVECTOR);
        }
        catch (...)
        {
          return;
        }
      }
      else // if (flag == GlobalVariable::WRITE)
        meshData_->write_global(meshIndex_, name, value,
          stk::util::ParameterType::INTEGERVECTOR);
    }

    // Vectors of doubles.
    else if (globalData_.isType<Array<double>>(name))
    {
      const auto& value = globalData_.get<Array<double>>(name).toVector();
      if (flag == GlobalVariable::ADD)
      {
        try
        {
          meshData_->add_global(meshIndex_, name, value,
            stk::util::ParameterType::DOUBLEVECTOR);
        }
        catch (...)
        {
          return;
        }
      }
      else // if (flag == GlobalVariable::WRITE)
        meshData_->write_global(meshIndex_, name, value,
          stk::util::ParameterType::DOUBLEVECTOR);
    }

    // If the data type is something else, throw an exception.
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, invalid_argument,
        "STK_Interface::globalToExodus():  The global variable to be added "  \
        "to the Exodus output file is of an invalid type.  Valid types are "  \
        "int and double, along with std::vectors of those types.")
  } // end loop over globalData_
} // end of globalToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  addGlobalToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
addGlobalToExodus(
  const std::string& key,
  const int&         value)
{
  globalData_.set(key, value);
} // end of addGlobalToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  addGlobalToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
addGlobalToExodus(
  const std::string& key,
  const double&      value)
{
  globalData_.set(key, value);
} // end of addGlobalToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  addGlobalToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
addGlobalToExodus(
  const std::string&      key,
  const std::vector<int>& value)
{
  using Teuchos::Array;
  globalData_.set(key, Array<int>(value));
} // end of addGlobalToExodus()

///////////////////////////////////////////////////////////////////////////////
//
//  addGlobalToExodus()
//
///////////////////////////////////////////////////////////////////////////////
void
STK_Interface::
addGlobalToExodus(
  const std::string&         key,
  const std::vector<double>& value)
{
  using Teuchos::Array;
  globalData_.set(key, Array<double>(value));
} // end of addGlobalToExodus()

bool STK_Interface::isWritable() const
{
   #ifdef PANZER_HAVE_IOSS
      return true;
   #else
      return false;
   #endif
}

void STK_Interface::getElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity> & elements) const
{
   stk::mesh::EntityRank nodeRank = getNodeRank();

   // get all relations for node
   stk::mesh::Entity node = bulkData_->get_entity(nodeRank,nodeId);
   const size_t numElements = bulkData_->num_elements(node);
   stk::mesh::Entity const* relations = bulkData_->begin_elements(node);

   // extract elements sharing nodes
   elements.insert(elements.end(), relations, relations + numElements);
}

void STK_Interface::getOwnedElementsSharingNode(stk::mesh::Entity node,std::vector<stk::mesh::Entity> & elements,
                                                std::vector<int> & relIds) const
{
   // get all relations for node
   const size_t numElements = bulkData_->num_elements(node);
   stk::mesh::Entity const* relations = bulkData_->begin_elements(node);
   stk::mesh::ConnectivityOrdinal const* rel_ids = bulkData_->begin_element_ordinals(node);

   // extract elements sharing nodes
   for (size_t i = 0; i < numElements; ++i) {
      stk::mesh::Entity element = relations[i];

     // if owned by this processor
      if(bulkData_->parallel_owner_rank(element) == static_cast<int>(procRank_)) {
         elements.push_back(element);
         relIds.push_back(rel_ids[i]);
      }
   }
}

void STK_Interface::getOwnedElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity> & elements,
                                                                           std::vector<int> & relIds, unsigned int matchType) const
{
   stk::mesh::EntityRank rank;
   if(matchType == 0)
     rank = getNodeRank();
   else if(matchType == 1)
     rank = getEdgeRank();
   else if(matchType == 2)
     rank = getFaceRank();
   else
     TEUCHOS_ASSERT(false);

   stk::mesh::Entity node = bulkData_->get_entity(rank,nodeId);

   getOwnedElementsSharingNode(node,elements,relIds);
}

void STK_Interface::getElementsSharingNodes(const std::vector<stk::mesh::EntityId> nodeIds,std::vector<stk::mesh::Entity> & elements) const
{
   std::vector<stk::mesh::Entity> current;

   getElementsSharingNode(nodeIds[0],current); // fill it with elements touching first node
   std::sort(current.begin(),current.end());   // sort for intersection on the pointer

   // find intersection with remaining nodes
   for(std::size_t n=1;n<nodeIds.size();++n) {
      // get elements associated with next node
      std::vector<stk::mesh::Entity> nextNode;
      getElementsSharingNode(nodeIds[n],nextNode); // fill it with elements touching first node
      std::sort(nextNode.begin(),nextNode.end());   // sort for intersection on the pointer ID

      // intersect next node elements with current element list
      std::vector<stk::mesh::Entity> intersection(std::min(nextNode.size(),current.size()));
      std::vector<stk::mesh::Entity>::const_iterator endItr
            = std::set_intersection(current.begin(),current.end(),
                                    nextNode.begin(),nextNode.end(),
                                    intersection.begin());
      std::size_t newLength = endItr-intersection.begin();
      intersection.resize(newLength);

      // store intersection
      current.clear();
      current = intersection;
   }

   // return the elements computed
   elements = current;
}

void STK_Interface::getNodeIdsForElement(stk::mesh::Entity element,std::vector<stk::mesh::EntityId> & nodeIds) const
{
  stk::mesh::Entity const* nodeRel = getBulkData()->begin_nodes(element);
  const size_t numNodes = getBulkData()->num_nodes(element);

  nodeIds.reserve(numNodes);
  for(size_t i = 0; i < numNodes; ++i) {
    nodeIds.push_back(elementGlobalId(nodeRel[i]));
  }
}

void STK_Interface::buildEntityCounts()
{
   entityCounts_.clear();
   stk::mesh::comm_mesh_counts(*bulkData_,entityCounts_);
}

void STK_Interface::buildMaxEntityIds()
{
   // developed to mirror "comm_mesh_counts" in stk_mesh/base/Comm.cpp

   const auto entityRankCount =  metaData_->entity_rank_count();
   const size_t   commCount        = 10; // entityRankCount

   TEUCHOS_ASSERT(entityRankCount<10);

   // stk::ParallelMachine mach = bulkData_->parallel();
   stk::ParallelMachine mach = *mpiComm_->getRawMpiComm();

   std::vector<stk::mesh::EntityId> local(commCount,0);

   // determine maximum ID for this processor for each entity type
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();
   for(stk::mesh::EntityRank i=stk::topology::NODE_RANK;
       i < static_cast<stk::mesh::EntityRank>(entityRankCount); ++i) {
      std::vector<stk::mesh::Entity> entities;

      stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(i),entities);

      // determine maximum ID for this processor
      std::vector<stk::mesh::Entity>::const_iterator itr;
      for(itr=entities.begin();itr!=entities.end();++itr) {
         stk::mesh::EntityId id = bulkData_->identifier(*itr);
         if(id>local[i])
            local[i] = id;
      }
   }

   // get largest IDs across processors
   stk::all_reduce(mach,stk::ReduceMax<10>(&local[0]));
   maxEntityId_.assign(local.begin(),local.begin()+entityRankCount+1);
}

std::size_t STK_Interface::getEntityCounts(unsigned entityRank) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(entityRank>=entityCounts_.size(),std::logic_error,
                      "STK_Interface::getEntityCounts: Entity counts do not include rank: " << entityRank);

   return entityCounts_[entityRank];
}

stk::mesh::EntityId STK_Interface::getMaxEntityId(unsigned entityRank) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(entityRank>=maxEntityId_.size(),std::logic_error,
                      "STK_Interface::getMaxEntityId: Max entity ids do not include rank: " << entityRank);

   return maxEntityId_[entityRank];
}

void STK_Interface::buildSubcells()
{
   stk::mesh::PartVector emptyPartVector;
   stk::mesh::create_adjacent_entities(*bulkData_,emptyPartVector);

   buildEntityCounts();
   buildMaxEntityIds();

   addEdges();
   if (dimension_ > 2)
     addFaces();
}

const double * STK_Interface::getNodeCoordinates(stk::mesh::EntityId nodeId) const
{
   stk::mesh::Entity node = bulkData_->get_entity(getNodeRank(),nodeId);
   return stk::mesh::field_data(*coordinatesField_,node);
}

const double * STK_Interface::getNodeCoordinates(stk::mesh::Entity node) const
{
   return stk::mesh::field_data(*coordinatesField_,node);
}

void STK_Interface::getSubcellIndices(unsigned entityRank,stk::mesh::EntityId elementId,
                                      std::vector<stk::mesh::EntityId> & subcellIds) const
{
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::Entity cell = bulkData_->get_entity(elementRank,elementId);

   TEUCHOS_TEST_FOR_EXCEPTION(!bulkData_->is_valid(cell),std::logic_error,
                      "STK_Interface::getSubcellIndices: could not find element requested (GID = " << elementId << ")");

   const size_t numSubcells = bulkData_->num_connectivity(cell, static_cast<stk::mesh::EntityRank>(entityRank));
   stk::mesh::Entity const* subcells = bulkData_->begin(cell, static_cast<stk::mesh::EntityRank>(entityRank));
   subcellIds.clear();
   subcellIds.resize(numSubcells,0);

   // loop over relations and fill subcell vector
   for(size_t i = 0; i < numSubcells; ++i) {
      stk::mesh::Entity subcell = subcells[i];
      subcellIds[i] = bulkData_->identifier(subcell);
   }
}

void STK_Interface::getMyElements(std::vector<stk::mesh::Entity> & elements) const
{
   // setup local ownership
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();

   // grab elements
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getMyElements(const std::string & blockID,std::vector<stk::mesh::Entity> & elements) const
{
   stk::mesh::Part * elementBlock = getElementBlockPart(blockID);

   TEUCHOS_TEST_FOR_EXCEPTION(elementBlock==0,std::logic_error,"Could not find element block \"" << blockID << "\"");

   // setup local ownership
   // stk::mesh::Selector block = *elementBlock;
   stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & (*elementBlock);

   // grab elements
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getNeighborElements(std::vector<stk::mesh::Entity> & elements) const
{
   // setup local ownership
   stk::mesh::Selector neighborBlock = (!metaData_->locally_owned_part());

   // grab elements
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::get_selected_entities(neighborBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getNeighborElements(const std::string & blockID,std::vector<stk::mesh::Entity> & elements) const
{
   stk::mesh::Part * elementBlock = getElementBlockPart(blockID);

   TEUCHOS_TEST_FOR_EXCEPTION(elementBlock==0,std::logic_error,"Could not find element block \"" << blockID << "\"");

   // setup local ownership
   stk::mesh::Selector neighborBlock = (!metaData_->locally_owned_part()) & (*elementBlock);

   // grab elements
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::get_selected_entities(neighborBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getMyEdges(std::vector<stk::mesh::Entity> & edges) const
{
   // setup local ownership
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();

   // grab elements
   stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(getEdgeRank()),edges);
}

void STK_Interface::getMyEdges(const std::string & edgeBlockName,std::vector<stk::mesh::Entity> & edges) const
{
   stk::mesh::Part * edgeBlockPart = getEdgeBlock(edgeBlockName);
   TEUCHOS_TEST_FOR_EXCEPTION(edgeBlockPart==0,std::logic_error,
                      "Unknown edge block \"" << edgeBlockName << "\"");

   stk::mesh::Selector edge_block = *edgeBlockPart;
   stk::mesh::Selector owned_block = metaData_->locally_owned_part() & edge_block;

   // grab elements
   stk::mesh::get_selected_entities(owned_block,bulkData_->buckets(getEdgeRank()),edges);
}

void STK_Interface::getMyEdges(const std::string & edgeBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & edges) const
{
   stk::mesh::Part * edgeBlockPart = getEdgeBlock(edgeBlockName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(edgeBlockPart==0,EdgeBlockException,
                      "Unknown edge block \"" << edgeBlockName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector edge_block = *edgeBlockPart;
   stk::mesh::Selector element_block = *elmtPart;
   stk::mesh::Selector owned_block = metaData_->locally_owned_part() & element_block & edge_block;

   // grab elements
   stk::mesh::get_selected_entities(owned_block,bulkData_->buckets(getEdgeRank()),edges);
}

void STK_Interface::getAllEdges(const std::string & edgeBlockName,std::vector<stk::mesh::Entity> & edges) const
{
   stk::mesh::Part * edgeBlockPart = getEdgeBlock(edgeBlockName);
   TEUCHOS_TEST_FOR_EXCEPTION(edgeBlockPart==0,std::logic_error,
                      "Unknown edge block \"" << edgeBlockName << "\"");

   stk::mesh::Selector edge_block = *edgeBlockPart;

   // grab elements
   stk::mesh::get_selected_entities(edge_block,bulkData_->buckets(getEdgeRank()),edges);
}

void STK_Interface::getAllEdges(const std::string & edgeBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & edges) const
{
   stk::mesh::Part * edgeBlockPart = getEdgeBlock(edgeBlockName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(edgeBlockPart==0,EdgeBlockException,
                      "Unknown edge block \"" << edgeBlockName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector edge_block = *edgeBlockPart;
   stk::mesh::Selector element_block = *elmtPart;
   stk::mesh::Selector element_edge_block = element_block & edge_block;

   // grab elements
   stk::mesh::get_selected_entities(element_edge_block,bulkData_->buckets(getEdgeRank()),edges);
}

void STK_Interface::getMyFaces(std::vector<stk::mesh::Entity> & faces) const
{
   // setup local ownership
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();

   // grab elements
   stk::mesh::EntityRank faceRank = getFaceRank();
   stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(faceRank),faces);
}

void STK_Interface::getMyFaces(const std::string & faceBlockName,std::vector<stk::mesh::Entity> & faces) const
{
   stk::mesh::Part * faceBlockPart = getFaceBlock(faceBlockName);
   TEUCHOS_TEST_FOR_EXCEPTION(faceBlockPart==0,std::logic_error,
                      "Unknown face block \"" << faceBlockName << "\"");

   stk::mesh::Selector face_block = *faceBlockPart;
   stk::mesh::Selector owned_block = metaData_->locally_owned_part() & face_block;

   // grab elements
   stk::mesh::get_selected_entities(owned_block,bulkData_->buckets(getFaceRank()),faces);
}

void STK_Interface::getMyFaces(const std::string & faceBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & faces) const
{
   stk::mesh::Part * faceBlockPart = getFaceBlock(faceBlockName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(faceBlockPart==0,FaceBlockException,
                      "Unknown face block \"" << faceBlockName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector face_block = *faceBlockPart;
   stk::mesh::Selector element_block = *elmtPart;
   stk::mesh::Selector owned_block = metaData_->locally_owned_part() & element_block & face_block;

   // grab elements
   stk::mesh::get_selected_entities(owned_block,bulkData_->buckets(getFaceRank()),faces);
}

void STK_Interface::getAllFaces(const std::string & faceBlockName,std::vector<stk::mesh::Entity> & faces) const
{
   stk::mesh::Part * faceBlockPart = getFaceBlock(faceBlockName);
   TEUCHOS_TEST_FOR_EXCEPTION(faceBlockPart==0,std::logic_error,
                      "Unknown face block \"" << faceBlockName << "\"");

   stk::mesh::Selector face_block = *faceBlockPart;

   // grab elements
   stk::mesh::get_selected_entities(face_block,bulkData_->buckets(getFaceRank()),faces);
}

void STK_Interface::getAllFaces(const std::string & faceBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & faces) const
{
   stk::mesh::Part * faceBlockPart = getFaceBlock(faceBlockName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(faceBlockPart==0,FaceBlockException,
                      "Unknown face block \"" << faceBlockName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector face_block = *faceBlockPart;
   stk::mesh::Selector element_block = *elmtPart;
   stk::mesh::Selector element_face_block = element_block & face_block;

   // grab elements
   stk::mesh::get_selected_entities(element_face_block,bulkData_->buckets(getFaceRank()),faces);
}

void STK_Interface::getMySides(const std::string & sideName,std::vector<stk::mesh::Entity> & sides) const
{
   stk::mesh::Part * sidePart = getSideset(sideName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                      "Unknown side set \"" << sideName << "\"");

   stk::mesh::Selector side = *sidePart;
   stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & side;

   // grab elements
   stk::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity> & sides) const
{
   stk::mesh::Part * sidePart = getSideset(sideName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,SidesetException,
                      "Unknown side set \"" << sideName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector side = *sidePart;
   stk::mesh::Selector block = *elmtPart;
   stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & block & side;

   // grab elements
   stk::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getAllSides(const std::string & sideName,std::vector<stk::mesh::Entity> & sides) const
{
   stk::mesh::Part * sidePart = getSideset(sideName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                      "Unknown side set \"" << sideName << "\"");

   stk::mesh::Selector side = *sidePart;

   // grab elements
   stk::mesh::get_selected_entities(side,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getAllSides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity> & sides) const
{
   stk::mesh::Part * sidePart = getSideset(sideName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,SidesetException,
                      "Unknown side set \"" << sideName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector side = *sidePart;
   stk::mesh::Selector block = *elmtPart;
   stk::mesh::Selector sideBlock = block & side;

   // grab elements
   stk::mesh::get_selected_entities(sideBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getMyNodes(const std::string & nodesetName,const std::string & blockName,std::vector<stk::mesh::Entity> & nodes) const
{
   stk::mesh::Part * nodePart = getNodeset(nodesetName);
   stk::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(nodePart==0,SidesetException,
                      "Unknown node set \"" << nodesetName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk::mesh::Selector nodeset = *nodePart;
   stk::mesh::Selector block = *elmtPart;
   stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & block & nodeset;

   // grab elements
   stk::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getNodeRank()),nodes);
}

void STK_Interface::getElementBlockNames(std::vector<std::string> & names) const
{
   // TEUCHOS_ASSERT(initialized_); // all blocks must have been added

   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator blkItr;   // Element blocks
   for(blkItr=elementBlocks_.begin();blkItr!=elementBlocks_.end();++blkItr)
      names.push_back(blkItr->first);
}

void STK_Interface::getSidesetNames(std::vector<std::string> & names) const
{
   // TEUCHOS_ASSERT(initialized_); // all blocks must have been added

   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator sideItr;   // Side sets
   for(sideItr=sidesets_.begin();sideItr!=sidesets_.end();++sideItr)
      names.push_back(sideItr->first);
}

void STK_Interface::getNodesetNames(std::vector<std::string> & names) const
{
   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator nodeItr;   // Node sets
   for(nodeItr=nodesets_.begin();nodeItr!=nodesets_.end();++nodeItr)
      names.push_back(nodeItr->first);
}

void STK_Interface::getEdgeBlockNames(std::vector<std::string> & names) const
{
   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator edgeBlockItr;   // Edge blocks
   for(edgeBlockItr=edgeBlocks_.begin();edgeBlockItr!=edgeBlocks_.end();++edgeBlockItr)
      names.push_back(edgeBlockItr->first);
}

void STK_Interface::getFaceBlockNames(std::vector<std::string> & names) const
{
   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator faceBlockItr;   // Face blocks
   for(faceBlockItr=faceBlocks_.begin();faceBlockItr!=faceBlocks_.end();++faceBlockItr)
      names.push_back(faceBlockItr->first);
}

std::size_t STK_Interface::elementLocalId(stk::mesh::Entity elmt) const
{
   return elementLocalId(bulkData_->identifier(elmt));
   // const std::size_t * fieldCoords = stk::mesh::field_data(*localIdField_,*elmt);
   // return fieldCoords[0];
}

std::size_t STK_Interface::elementLocalId(stk::mesh::EntityId gid) const
{
   // stk::mesh::EntityRank elementRank = getElementRank();
   // stk::mesh::Entity elmt = bulkData_->get_entity(elementRank,gid);
   // TEUCHOS_ASSERT(elmt->owner_rank()==procRank_);
   // return elementLocalId(elmt);
   std::unordered_map<stk::mesh::EntityId,std::size_t>::const_iterator itr = localIDHash_.find(gid);
   TEUCHOS_ASSERT(itr!=localIDHash_.end());
   return itr->second;
}

bool STK_Interface::isEdgeLocal(stk::mesh::Entity edge) const
{
   return isEdgeLocal(bulkData_->identifier(edge));
}

bool STK_Interface::isEdgeLocal(stk::mesh::EntityId gid) const
{
   std::unordered_map<stk::mesh::EntityId,std::size_t>::const_iterator itr = localEdgeIDHash_.find(gid);
   if (itr==localEdgeIDHash_.end()) {
     return false;
   }
   return true;
}

std::size_t STK_Interface::edgeLocalId(stk::mesh::Entity edge) const
{
   return edgeLocalId(bulkData_->identifier(edge));
}

std::size_t STK_Interface::edgeLocalId(stk::mesh::EntityId gid) const
{
   std::unordered_map<stk::mesh::EntityId,std::size_t>::const_iterator itr = localEdgeIDHash_.find(gid);
   TEUCHOS_ASSERT(itr!=localEdgeIDHash_.end());
   return itr->second;
}

bool STK_Interface::isFaceLocal(stk::mesh::Entity face) const
{
   return isFaceLocal(bulkData_->identifier(face));
}

bool STK_Interface::isFaceLocal(stk::mesh::EntityId gid) const
{
   std::unordered_map<stk::mesh::EntityId,std::size_t>::const_iterator itr = localFaceIDHash_.find(gid);
   if (itr==localFaceIDHash_.end()) {
     return false;
   }
   return true;
}

std::size_t STK_Interface::faceLocalId(stk::mesh::Entity face) const
{
   return faceLocalId(bulkData_->identifier(face));
}

std::size_t STK_Interface::faceLocalId(stk::mesh::EntityId gid) const
{
   std::unordered_map<stk::mesh::EntityId,std::size_t>::const_iterator itr = localFaceIDHash_.find(gid);
   TEUCHOS_ASSERT(itr!=localFaceIDHash_.end());
   return itr->second;
}


std::string STK_Interface::containingBlockId(stk::mesh::Entity elmt) const
{
   for(const auto & eb_pair : elementBlocks_)
      if(bulkData_->bucket(elmt).member(*(eb_pair.second)))
         return eb_pair.first;
   return "";
}

stk::mesh::Field<double> * STK_Interface::getSolutionField(const std::string & fieldName,
                                                           const std::string & blockId) const
{
   // look up field in map
   std::map<std::pair<std::string,std::string>, SolutionFieldType*>::const_iterator
         iter = fieldNameToSolution_.find(std::make_pair(fieldName,blockId));

   // check to make sure field was actually found
   TEUCHOS_TEST_FOR_EXCEPTION(iter==fieldNameToSolution_.end(),std::runtime_error,
                      "Solution field name \"" << fieldName << "\" in block ID \"" << blockId << "\" was not found");

   return iter->second;
}

stk::mesh::Field<double> * STK_Interface::getCellField(const std::string & fieldName,
                                                       const std::string & blockId) const
{
   // look up field in map
   std::map<std::pair<std::string,std::string>, SolutionFieldType*>::const_iterator
         iter = fieldNameToCellField_.find(std::make_pair(fieldName,blockId));

   // check to make sure field was actually found
   TEUCHOS_TEST_FOR_EXCEPTION(iter==fieldNameToCellField_.end(),std::runtime_error,
                      "Cell field named \"" << fieldName << "\" in block ID \"" << blockId << "\" was not found");

   return iter->second;
}

stk::mesh::Field<double> * STK_Interface::getEdgeField(const std::string & fieldName,
                                                       const std::string & blockId) const
{
   // look up field in map
   std::map<std::pair<std::string,std::string>, SolutionFieldType*>::const_iterator
         iter = fieldNameToEdgeField_.find(std::make_pair(fieldName,blockId));

   // check to make sure field was actually found
   TEUCHOS_TEST_FOR_EXCEPTION(iter==fieldNameToEdgeField_.end(),std::runtime_error,
                      "Edge field named \"" << fieldName << "\" in block ID \"" << blockId << "\" was not found");

   return iter->second;
}

stk::mesh::Field<double> * STK_Interface::getFaceField(const std::string & fieldName,
                                                       const std::string & blockId) const
{
   // look up field in map
   std::map<std::pair<std::string,std::string>, SolutionFieldType*>::const_iterator
         iter = fieldNameToFaceField_.find(std::make_pair(fieldName,blockId));

   // check to make sure field was actually found
   TEUCHOS_TEST_FOR_EXCEPTION(iter==fieldNameToFaceField_.end(),std::runtime_error,
                      "Face field named \"" << fieldName << "\" in block ID \"" << blockId << "\" was not found");

   return iter->second;
}

Teuchos::RCP<const std::vector<stk::mesh::Entity> > STK_Interface::getElementsOrderedByLID() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   if(orderedElementVector_==Teuchos::null) {
      // safe because essentially this is a call to modify a mutable object
      const_cast<STK_Interface*>(this)->buildLocalElementIDs();
   }

   return orderedElementVector_.getConst();
}

void STK_Interface::addElementBlock(const std::string & name,const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);

   stk::mesh::Part * block = metaData_->get_part(name);
   if(block==0) {
     block = &metaData_->declare_part_with_topology(name, stk::mesh::get_topology(shards::CellTopology(ctData), dimension_));
   }

   // construct cell topology object for this block
   Teuchos::RCP<shards::CellTopology> ct
         = Teuchos::rcp(new shards::CellTopology(ctData));

   // add element block part and cell topology
   elementBlocks_.insert(std::make_pair(name,block));
   elementBlockCT_.insert(std::make_pair(name,ct));
}

Teuchos::RCP<const std::vector<stk::mesh::Entity> > STK_Interface::getEdgesOrderedByLID() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   if(orderedEdgeVector_==Teuchos::null) {
      // safe because essentially this is a call to modify a mutable object
      const_cast<STK_Interface*>(this)->buildLocalEdgeIDs();
   }

   return orderedEdgeVector_.getConst();
}

void STK_Interface::addEdgeBlock(const std::string & elemBlockName,
                                 const std::string & edgeBlockName,
                                 const stk::topology & topology)
{
   TEUCHOS_ASSERT(not initialized_);

   stk::mesh::Part * edge_block = metaData_->get_part(edgeBlockName);
   if(edge_block==0) {
      edge_block = &metaData_->declare_part_with_topology(edgeBlockName, topology);
   }

   /* There is only one edge block for each edge topology, so declare it
    * as a subset of the element block even if it wasn't just created.
    */
   stk::mesh::Part * elem_block = metaData_->get_part(elemBlockName);
   metaData_->declare_part_subset(*elem_block, *edge_block);

   // add edge block part
   edgeBlocks_.insert(std::make_pair(edgeBlockName,edge_block));
}

void STK_Interface::addEdgeBlock(const std::string & elemBlockName,
                                 const std::string & edgeBlockName,
                                 const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);

   addEdgeBlock(elemBlockName,
                edgeBlockName,
                stk::mesh::get_topology(shards::CellTopology(ctData), dimension_));
}

Teuchos::RCP<const std::vector<stk::mesh::Entity> > STK_Interface::getFacesOrderedByLID() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   if(orderedFaceVector_==Teuchos::null) {
      // safe because essentially this is a call to modify a mutable object
      const_cast<STK_Interface*>(this)->buildLocalFaceIDs();
   }

   return orderedFaceVector_.getConst();
}

void STK_Interface::addFaceBlock(const std::string & elemBlockName,
                                 const std::string & faceBlockName,
                                 const stk::topology & topology)
{
   TEUCHOS_ASSERT(not initialized_);

   stk::mesh::Part * face_block = metaData_->get_part(faceBlockName);
   if(face_block==0) {
      face_block = &metaData_->declare_part_with_topology(faceBlockName, topology);
   }

   /* There is only one face block for each edge topology, so declare it
    * as a subset of the element block even if it wasn't just created.
    */
   stk::mesh::Part * elem_block = metaData_->get_part(elemBlockName);
   metaData_->declare_part_subset(*elem_block, *face_block);

   // add face block part
   faceBlocks_.insert(std::make_pair(faceBlockName,face_block));
}

void STK_Interface::addFaceBlock(const std::string & elemBlockName,
                                 const std::string & faceBlockName,
                                 const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);

   addFaceBlock(elemBlockName,
                faceBlockName,
                stk::mesh::get_topology(shards::CellTopology(ctData), dimension_));
}

void STK_Interface::initializeFromMetaData()
{
   dimension_ = metaData_->spatial_dimension();

   // declare coordinates and node parts
   coordinatesField_ = &metaData_->declare_field<double>(stk::topology::NODE_RANK, coordsString);
   edgesField_       = &metaData_->declare_field<double>(stk::topology::EDGE_RANK, edgesString);
   if (dimension_ > 2)
     facesField_       = &metaData_->declare_field<double>(stk::topology::FACE_RANK, facesString);
   processorIdField_ = &metaData_->declare_field<ProcIdData>(stk::topology::ELEMENT_RANK, "PROC_ID");
   loadBalField_     = &metaData_->declare_field<double>(stk::topology::ELEMENT_RANK, "LOAD_BAL");

   // stk::mesh::put_field( *coordinatesField_ , getNodeRank() , metaData_->universal_part() );

   nodesPart_        = &metaData_->declare_part(nodesString,getNodeRank());
   nodesPartVec_.push_back(nodesPart_);
   edgesPart_        = &metaData_->declare_part(edgesString,getEdgeRank());
   edgesPartVec_.push_back(edgesPart_);
   if (dimension_ > 2) {
     facesPart_        = &metaData_->declare_part(facesString,getFaceRank());
     facesPartVec_.push_back(facesPart_);
   }
}

void STK_Interface::buildLocalElementIDs()
{
   currentLocalId_ = 0;

   orderedElementVector_ = Teuchos::null; // forces rebuild of ordered lists

   // might be better (faster) to do this by buckets
   std::vector<stk::mesh::Entity> elements;
   getMyElements(elements);

   for(std::size_t index=0;index<elements.size();++index) {
      stk::mesh::Entity element = elements[index];

      // set processor rank
      ProcIdData * procId = stk::mesh::field_data(*processorIdField_,element);
      procId[0] = Teuchos::as<ProcIdData>(procRank_);

      localIDHash_[bulkData_->identifier(element)] = currentLocalId_;

      currentLocalId_++;
   }

   // copy elements into the ordered element vector
   orderedElementVector_ = Teuchos::rcp(new std::vector<stk::mesh::Entity>(elements));

   elements.clear();
   getNeighborElements(elements);

   for(std::size_t index=0;index<elements.size();++index) {
      stk::mesh::Entity element = elements[index];

      // set processor rank
      ProcIdData * procId = stk::mesh::field_data(*processorIdField_,element);
      procId[0] = Teuchos::as<ProcIdData>(procRank_);

      localIDHash_[bulkData_->identifier(element)] = currentLocalId_;

      currentLocalId_++;
   }

   orderedElementVector_->insert(orderedElementVector_->end(),elements.begin(),elements.end());
}

void STK_Interface::applyElementLoadBalanceWeights()
{
  std::vector<std::string> names;
  getElementBlockNames(names);

  for(std::size_t b=0;b<names.size();b++) {
    // find user specified block weight, otherwise use 1.0
    std::map<std::string,double>::const_iterator bw_itr = blockWeights_.find(names[b]);
    double blockWeight = (bw_itr!=blockWeights_.end()) ? bw_itr->second : 1.0;

    std::vector<stk::mesh::Entity> elements;
    getMyElements(names[b],elements);

    for(std::size_t index=0;index<elements.size();++index) {
      // set local element ID
      double * loadBal = stk::mesh::field_data(*loadBalField_,elements[index]);
      loadBal[0] = blockWeight;
    }
  }
}

void STK_Interface::buildLocalEdgeIDs()
{
   currentLocalId_ = 0;

   orderedEdgeVector_ = Teuchos::null; // forces rebuild of ordered lists

   // might be better (faster) to do this by buckets
   std::vector<stk::mesh::Entity> edges;
   getMyEdges(edges);

   for(std::size_t index=0;index<edges.size();++index) {
      stk::mesh::Entity edge = edges[index];
      localEdgeIDHash_[bulkData_->identifier(edge)] = currentLocalId_;
      currentLocalId_++;
   }

   // copy edges into the ordered edge vector
   orderedEdgeVector_ = Teuchos::rcp(new std::vector<stk::mesh::Entity>(edges));
}

void STK_Interface::buildLocalFaceIDs()
{
   currentLocalId_ = 0;

   orderedFaceVector_ = Teuchos::null; // forces rebuild of ordered lists

   // might be better (faster) to do this by buckets
   std::vector<stk::mesh::Entity> faces;
   getMyFaces(faces);

   for(std::size_t index=0;index<faces.size();++index) {
      stk::mesh::Entity face = faces[index];
      localFaceIDHash_[bulkData_->identifier(face)] = currentLocalId_;
      currentLocalId_++;
   }

   // copy faces into the ordered face vector
   orderedFaceVector_ = Teuchos::rcp(new std::vector<stk::mesh::Entity>(faces));
}

bool
STK_Interface::isMeshCoordField(const std::string & eBlock,
                                const std::string & fieldName,
                                int & axis) const
{
  std::map<std::string,std::vector<std::string> >::const_iterator blkItr = meshCoordFields_.find(eBlock);
  if(blkItr==meshCoordFields_.end()) {
    return false;
  }

  axis = 0;
  for(axis=0;axis<Teuchos::as<int>(blkItr->second.size());axis++) {
    if(blkItr->second[axis]==fieldName)
      break; // found axis, break
  }

  if(axis>=Teuchos::as<int>(blkItr->second.size()))
    return false;

  return true;
}

std::pair<Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >, Teuchos::RCP<std::vector<unsigned int> > >
STK_Interface::getPeriodicNodePairing() const
{
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > vec;
   Teuchos::RCP<std::vector<unsigned int > > type_vec = rcp(new std::vector<unsigned int>);
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & matchers = getPeriodicBCVector();
   const bool & useBBoxSearch = useBoundingBoxSearch();
   std::vector<std::vector<std::string> > matchedSides(3); // (coord,edge,face)

   // build up the vectors by looping over the matched pair
   for(std::size_t m=0;m<matchers.size();m++){
      unsigned int type;
      if(matchers[m]->getType() == "coord")
        type = 0;
      else if(matchers[m]->getType() == "edge")
        type = 1;
      else if(matchers[m]->getType() == "face")
        type = 2;
      else
        TEUCHOS_ASSERT(false);
#ifdef PANZER_HAVE_STKSEARCH

      if (useBBoxSearch) {
         vec = matchers[m]->getMatchedPair(*this,matchedSides[type],vec);
      } else {
         vec = matchers[m]->getMatchedPair(*this,vec);
      }
#else 
      TEUCHOS_TEST_FOR_EXCEPTION(useBBoxSearch,std::logic_error,
          "panzer::STK_Interface::getPeriodicNodePairing(): Requested bounding box search, but "
          "did not compile with STK_SEARCH enabled.");
      vec = matchers[m]->getMatchedPair(*this,vec);
#endif
      type_vec->insert(type_vec->begin(),vec->size()-type_vec->size(),type);
      matchedSides[type].push_back(matchers[m]->getLeftSidesetName());
   }

   return std::make_pair(vec,type_vec);
}

bool STK_Interface::validBlockId(const std::string & blockId) const
{
   std::map<std::string, stk::mesh::Part*>::const_iterator blkItr = elementBlocks_.find(blockId);

   return blkItr!=elementBlocks_.end();
}

void STK_Interface::print(std::ostream & os) const
{
   std::vector<std::string> blockNames, sidesetNames, nodesetNames;

   getElementBlockNames(blockNames);
   getSidesetNames(sidesetNames);
   getNodesetNames(nodesetNames);

   os << "STK Mesh data:\n";
   os << "   Spatial dim = " << getDimension() << "\n";
   if(getDimension()==2)
      os << "   Entity counts (Nodes, Edges, Cells) = ( "
         << getEntityCounts(getNodeRank()) << ", "
         << getEntityCounts(getEdgeRank()) << ", "
         << getEntityCounts(getElementRank()) << " )\n";
   else if(getDimension()==3)
      os << "   Entity counts (Nodes, Edges, Faces, Cells) = ( "
         << getEntityCounts(getNodeRank()) << ", "
         << getEntityCounts(getEdgeRank()) << ", "
         << getEntityCounts(getSideRank()) << ", "
         << getEntityCounts(getElementRank()) << " )\n";
   else
      os << "   Entity counts (Nodes, Cells) = ( "
         << getEntityCounts(getNodeRank()) << ", "
         << getEntityCounts(getElementRank()) << " )\n";

   os << "   Element blocks = ";
   for(std::size_t i=0;i<blockNames.size();i++)
      os << "\"" << blockNames[i] << "\" ";
   os << "\n";
   os << "   Sidesets = ";
   for(std::size_t i=0;i<sidesetNames.size();i++)
      os << "\"" << sidesetNames[i] << "\" ";
   os << "\n";
   os << "   Nodesets = ";
   for(std::size_t i=0;i<nodesetNames.size();i++)
      os << "\"" << nodesetNames[i] << "\" ";
   os << std::endl;

   // print out periodic boundary conditions
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & bcVector
         = getPeriodicBCVector();
   if(bcVector.size()!=0) {
      os << "   Periodic BCs:\n";
      for(std::size_t i=0;i<bcVector.size();i++)
         os << "      " << bcVector[i]->getString() << "\n";
      os << std::endl;
   }
}

void STK_Interface::printMetaData(std::ostream & os) const
{
   std::vector<std::string> blockNames, sidesetNames, nodesetNames;

   getElementBlockNames(blockNames);
   getSidesetNames(sidesetNames);
   getNodesetNames(nodesetNames);

   os << "STK Meta data:\n";
   os << "   Element blocks = ";
   for(std::size_t i=0;i<blockNames.size();i++)
      os << "\"" << blockNames[i] << "\" ";
   os << "\n";
   os << "   Sidesets = ";
   for(std::size_t i=0;i<sidesetNames.size();i++)
      os << "\"" << sidesetNames[i] << "\" ";
   os << "\n";
   os << "   Nodesets = ";
   for(std::size_t i=0;i<nodesetNames.size();i++)
      os << "\"" << nodesetNames[i] << "\" ";
   os << std::endl;

   // print out periodic boundary conditions
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & bcVector
         = getPeriodicBCVector();
   if(bcVector.size()!=0) {
      os << "   Periodic BCs:\n";
      for(std::size_t i=0;i<bcVector.size();i++)
         os << "      " << bcVector[i]->getString() << "\n";
      os << std::endl;
   }

   // print all available fields in meta data
   os << "   Fields = ";
   const stk::mesh::FieldVector & fv = metaData_->get_fields();
   for(std::size_t i=0;i<fv.size();i++)
      os << "\"" << fv[i]->name() << "\" ";
   os << std::endl;
}

Teuchos::RCP<const shards::CellTopology> STK_Interface::getCellTopology(const std::string & eBlock) const
{
   std::map<std::string, Teuchos::RCP<shards::CellTopology> >::const_iterator itr;
   itr = elementBlockCT_.find(eBlock);

   if(itr==elementBlockCT_.end()) {
      std::stringstream ss;
      printMetaData(ss);
      TEUCHOS_TEST_FOR_EXCEPTION(itr==elementBlockCT_.end(),std::logic_error,
                                 "STK_Interface::getCellTopology: No such element block \"" +eBlock +"\" available.\n\n"
                              << "STK Meta Data follows: \n" << ss.str());
   }

   return itr->second;
}

Teuchos::RCP<Teuchos::MpiComm<int> > STK_Interface::getSafeCommunicator(stk::ParallelMachine parallelMach) const
{
   MPI_Comm newComm;
   const int err = MPI_Comm_dup (parallelMach, &newComm);
   TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
     "panzer::STK_Interface: MPI_Comm_dup failed with error \""
     << Teuchos::mpiErrorCodeToString (err) << "\".");

   return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper (newComm,MPI_Comm_free)));
}

void STK_Interface::rebalance(const Teuchos::ParameterList & /* params */)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Rebalance not currently supported");
#if 0
  // make sure weights have been set (a local operation)
  applyElementLoadBalanceWeights();

  stk::mesh::Selector selector(getMetaData()->universal_part());
  stk::mesh::Selector owned_selector(getMetaData()->locally_owned_part());

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out << "Load balance before: " << stk::rebalance::check_balance(*getBulkData(), loadBalField_, getElementRank(), &selector) << std::endl;

  // perform reblance
  Teuchos::ParameterList graph;
  if(params.begin()!=params.end())
    graph.sublist(stk::rebalance::Zoltan::default_parameters_name()) = params;
  stk::rebalance::Zoltan zoltan_partition(*mpiComm_->getRawMpiComm(), getDimension(), graph);
  stk::rebalance::rebalance(*getBulkData(), owned_selector, &getCoordinatesField(), loadBalField_, zoltan_partition);

  out << "Load balance after: " << stk::rebalance::check_balance(*getBulkData(), loadBalField_, getElementRank(), &selector) << std::endl;

  currentLocalId_ = 0;
  orderedElementVector_ = Teuchos::null; // forces rebuild of ordered lists
#endif
}

Teuchos::RCP<const Teuchos::Comm<int> >
STK_Interface::getComm() const
{
  TEUCHOS_ASSERT(this->isInitialized());
  return mpiComm_;
}

void STK_Interface::refineMesh(const int numberOfLevels, const bool deleteParentElements) {
#ifdef PANZER_HAVE_PERCEPT
  TEUCHOS_TEST_FOR_EXCEPTION(numberOfLevels < 1,std::runtime_error,
                             "ERROR: Number of levels for uniform refinement must be greater than 0");
  TEUCHOS_ASSERT(nonnull(refinedMesh_));
  TEUCHOS_ASSERT(nonnull(breakPattern_));

  refinedMesh_->setCoordinatesField();

  percept::UniformRefiner breaker(*refinedMesh_,*breakPattern_);
  breaker.setRemoveOldElements(deleteParentElements);

  for (int i=0; i < numberOfLevels; ++i)
    breaker.doBreak();

#else
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "ERROR: Uniform refinement requested. This requires the Percept package to be enabled in Trilinos!");
#endif
}


} // end namespace panzer_stk
