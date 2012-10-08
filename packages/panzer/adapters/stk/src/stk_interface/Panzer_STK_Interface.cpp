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

#include <Panzer_STK_Interface.hpp>

#include <Teuchos_as.hpp>

#include <limits>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#ifdef HAVE_IOSS
#include <Ionit_Initializer.h>
#include <stk_io/IossBridge.hpp>
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

STK_Interface::STK_Interface()
   : dimension_(0), initialized_(false), currentLocalId_(0), initialStateTime_(0.0), currentStateTime_(0.0)
{
   metaData_ = rcp(new stk::mesh::fem::FEMMetaData());
}

STK_Interface::STK_Interface(unsigned dim)
   : dimension_(dim), initialized_(false), currentLocalId_(0)
{
   std::vector<std::string> entity_rank_names = stk::mesh::fem::entity_rank_names(dimension_);
   entity_rank_names.push_back("FAMILY_TREE");

   metaData_ = rcp(new stk::mesh::fem::FEMMetaData());
   metaData_->FEM_initialize(dimension_,entity_rank_names);

   initializeFromMetaData();
}

void STK_Interface::addSideset(const std::string & name,const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk::mesh::Part * sideset = metaData_->get_part(name);
   if(sideset==NULL)
      sideset = &metaData_->declare_part(name,stk::mesh::fem::CellTopology(ctData)); 
   sidesets_.insert(std::make_pair(name,sideset));
}

void STK_Interface::addNodeset(const std::string & name)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk::mesh::Part * nodeset = metaData_->get_part(name);
   if(nodeset==NULL) {
      const CellTopologyData * ctData = shards::getCellTopologyData<shards::Node>();
      nodeset = &metaData_->declare_part(name,stk::mesh::fem::CellTopology(ctData)); 
   }
   nodesets_.insert(std::make_pair(name,nodeset));
}

void STK_Interface::addSolutionField(const std::string & fieldName,const std::string & blockId) 
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToSolution_.find(key)==fieldNameToSolution_.end()) {
      SolutionFieldType * field = metaData_->get_field<SolutionFieldType>(fieldName);
      if(field==0)
         field = &metaData_->declare_field<SolutionFieldType>(fieldName);     
      fieldNameToSolution_[key] = field;
   }
}

void STK_Interface::addCellField(const std::string & fieldName,const std::string & blockId) 
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_TEST_FOR_EXCEPTION(!validBlockId(blockId),ElementBlockException,
                      "Unknown element block \"" << blockId << "\"");
   std::pair<std::string,std::string> key = std::make_pair(fieldName,blockId);

   // add & declare field if not already added...currently assuming linears
   if(fieldNameToCellField_.find(key)==fieldNameToCellField_.end()) {
      SolutionFieldType * field = metaData_->get_field<SolutionFieldType>(fieldName);
      if(field==0)
         field = &metaData_->declare_field<SolutionFieldType>(fieldName);     
      fieldNameToCellField_[key] = field;
   }
}

void STK_Interface::initialize(stk::ParallelMachine parallelMach,bool setupIO) 
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0); // no zero dimensional meshes!

   if(mpiComm_==Teuchos::null)
      mpiComm_ = getSafeCommunicator(parallelMach);

   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::EntityRank nodeRank = getNodeRank();

   procRank_ = stk::parallel_machine_rank(*mpiComm_->getRawMpiComm());

   // associating the field with a part: universal part!
   stk::mesh::put_field( *coordinatesField_ , nodeRank, metaData_->universal_part(), getDimension());
   stk::mesh::put_field( *processorIdField_ , elementRank, metaData_->universal_part());
   stk::mesh::put_field( *localIdField_ , elementRank, metaData_->universal_part());

   initializeFieldsInSTK(fieldNameToSolution_,nodeRank,setupIO);
   initializeFieldsInSTK(fieldNameToCellField_,elementRank,setupIO);

#ifdef HAVE_IOSS
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
      stk::io::set_field_role(*processorIdField_, Ioss::Field::TRANSIENT);
   }
#endif

   metaData_->commit();
   if(bulkData_==Teuchos::null)
      instantiateBulkData(*mpiComm_->getRawMpiComm());

   initialized_ = true;
}

void STK_Interface::initializeFieldsInSTK(const std::map<std::pair<std::string,std::string>,SolutionFieldType*> & nameToField,
                                          stk::mesh::EntityRank rank,bool setupIO)
{
   std::set<SolutionFieldType*> uniqueFields;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*>::const_iterator fieldIter;
   for(fieldIter=nameToField.begin();fieldIter!=nameToField.end();++fieldIter)
      uniqueFields.insert(fieldIter->second); // this makes setting up IO easier!

   {
      std::set<SolutionFieldType*>::const_iterator uniqueFieldIter;
      for(uniqueFieldIter=uniqueFields.begin();uniqueFieldIter!=uniqueFields.end();++uniqueFieldIter)
         stk::mesh::put_field(*(*uniqueFieldIter),rank,metaData_->universal_part());
   }

#ifdef HAVE_IOSS
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

   bulkData_ = rcp(new stk::mesh::BulkData(stk::mesh::fem::FEMMetaData::get_meta_data(*metaData_),*mpiComm_->getRawMpiComm()));
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

   bulkData_->modification_end();

   buildEntityCounts();
   buildMaxEntityIds();
}

void STK_Interface::addNode(stk::mesh::EntityId gid, const std::vector<double> & coord)
{
   TEUCHOS_TEST_FOR_EXCEPTION(not isModifiable(),std::logic_error,
                      "STK_Interface::addNode: STK_Interface must be modifiable to add a node");
   TEUCHOS_TEST_FOR_EXCEPTION(not coord.size()==getDimension(),std::logic_error,
                      "STK_Interface::addNode: number of coordinates in vector must mation dimension");
   TEUCHOS_TEST_FOR_EXCEPTION(gid==0,std::logic_error,
                      "STK_Interface::addNode: STK has STUPID restriction of no zero GIDs, pick something else");
   stk::mesh::EntityRank nodeRank = getNodeRank();

   stk::mesh::Entity & node = bulkData_->declare_entity(nodeRank,gid,nodesPartVec_);

   // set coordinate vector
   double * fieldCoords = stk::mesh::field_data(*coordinatesField_,node);
   for(std::size_t i=0;i<coord.size();++i)
      fieldCoords[i] = coord[i];
}

void STK_Interface::addEntityToSideset(stk::mesh::Entity & entity,stk::mesh::Part * sideset)
{
   std::vector<stk::mesh::Part*> sidesetV;
   sidesetV.push_back(sideset);

   bulkData_->change_entity_parts(entity,sidesetV);
}

void STK_Interface::addEntityToNodeset(stk::mesh::Entity & entity,stk::mesh::Part * nodeset)
{
   std::vector<stk::mesh::Part*> nodesetV;
   nodesetV.push_back(nodeset);

   bulkData_->change_entity_parts(entity,nodesetV);
}

void STK_Interface::addElement(const Teuchos::RCP<ElementDescriptor> & ed,stk::mesh::Part * block)
{
   std::vector<stk::mesh::Part*> blockVec;
   blockVec.push_back(block);

   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::EntityRank nodeRank = getNodeRank();
   stk::mesh::Entity & element = bulkData_->declare_entity(elementRank,ed->getGID(),blockVec);

   // build relations that give the mesh structure
   const std::vector<stk::mesh::EntityId> & nodes = ed->getNodes();
   for(std::size_t i=0;i<nodes.size();++i) {
      // add element->node relation
      stk::mesh::Entity * node = bulkData_->get_entity(nodeRank,nodes[i]);
      TEUCHOS_ASSERT(node!=0);
      bulkData_->declare_relation(element,*node,i);
   }

   ProcIdData * procId = stk::mesh::field_data(*processorIdField_,element);
   procId[0] = Teuchos::as<ProcIdData>(procRank_);
}

void STK_Interface::writeToExodus(const std::string & filename)
{
   #ifdef HAVE_IOSS
      TEUCHOS_ASSERT(mpiComm_!=Teuchos::null);
      stk::ParallelMachine comm = *mpiComm_->getRawMpiComm();

      Ioss::Init::Initializer io;
      stk::io::MeshData meshData;
      stk::io::create_output_mesh(filename, comm, *bulkData_, meshData);
      stk::io::define_output_fields(meshData,*metaData_);

      stk::io::process_output_request(meshData, *bulkData_, 0.0);
   #else 
      TEUCHOS_ASSERT(false);
   #endif
}

void STK_Interface::setupTransientExodusFile(const std::string & filename)
{
   #ifdef HAVE_IOSS
      TEUCHOS_ASSERT(mpiComm_!=Teuchos::null);
      stk::ParallelMachine comm = *mpiComm_->getRawMpiComm();

      Ioss::Init::Initializer io;
      meshData_ = Teuchos::rcp(new stk::io::MeshData);
      stk::io::create_output_mesh(filename, comm, *bulkData_, *meshData_);
      stk::io::define_output_fields(*meshData_,*metaData_);
   #else 
      TEUCHOS_ASSERT(false);
   #endif
}

void STK_Interface::writeToExodus(double timestep)
{
   #ifdef HAVE_IOSS
      currentStateTime_ = timestep;
      stk::io::process_output_request(*meshData_, *bulkData_, timestep);
   #else 
      TEUCHOS_ASSERT(false);
   #endif
}

bool STK_Interface::isWritable() const
{
   #ifdef HAVE_IOSS
      return true;
   #else
      return false;
   #endif
}

void STK_Interface::getElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity *> & elements) const
{
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::EntityRank nodeRank = getNodeRank();

   // get all relations for node
   stk::mesh::Entity * node = bulkData_->get_entity(nodeRank,nodeId);
   stk::mesh::PairIterRelation relations = node->relations(elementRank);

   // extract elements sharing nodes
   stk::mesh::PairIterRelation::iterator itr;
   for(itr=relations.begin();itr!=relations.end();++itr) {
      elements.push_back(itr->entity());
   }
}

void STK_Interface::getOwnedElementsSharingNode(stk::mesh::Entity * node,std::vector<stk::mesh::Entity *> & elements,
                                                                         std::vector<int> & localNodeId) const
{
   stk::mesh::EntityRank elementRank = getElementRank();

   // get all relations for node
   stk::mesh::PairIterRelation relations = node->relations(elementRank);

   // extract elements sharing nodes
   stk::mesh::PairIterRelation::iterator itr;
   for(itr=relations.begin();itr!=relations.end();++itr) {
      stk::mesh::Entity * element = itr->entity();
      
      // if owned by this processor 
      if(element->owner_rank() == procRank_) {
         elements.push_back(element);
         localNodeId.push_back(itr->identifier());
      }

   }
}

void STK_Interface::getOwnedElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity *> & elements,
                                                                           std::vector<int> & localNodeId) const
{
   stk::mesh::Entity * node = bulkData_->get_entity(getNodeRank(),nodeId);

   getOwnedElementsSharingNode(node,elements,localNodeId);
}

void STK_Interface::getElementsSharingNodes(const std::vector<stk::mesh::EntityId> nodeIds,std::vector<stk::mesh::Entity *> & elements) const
{
   std::vector<stk::mesh::Entity*> current;

   getElementsSharingNode(nodeIds[0],current); // fill it with elements touching first node
   std::sort(current.begin(),current.end());   // sort for intersection on the pointer 

   // find intersection with remaining nodes
   for(std::size_t n=1;n<nodeIds.size();++n) {
      // get elements associated with next node
      std::vector<stk::mesh::Entity*> nextNode;
      getElementsSharingNode(nodeIds[n],nextNode); // fill it with elements touching first node
      std::sort(nextNode.begin(),nextNode.end());   // sort for intersection on the pointer ID

      // intersect next node elements with current element list
      std::vector<stk::mesh::Entity*> intersection(std::min(nextNode.size(),current.size())); 
      std::vector<stk::mesh::Entity*>::const_iterator endItr
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

void STK_Interface::buildEntityCounts()
{
   entityCounts_.clear();
   stk::mesh::comm_mesh_counts(*bulkData_,entityCounts_);
}

void STK_Interface::buildMaxEntityIds()
{
   // developed to mirror "comm_mesh_counts" in stk_mesh/base/Comm.cpp

   const unsigned entityRankCount =  metaData_->entity_rank_count();
   const size_t   commCount        = 10; // entityRankCount

   TEUCHOS_ASSERT(entityRankCount<10);

   // stk::ParallelMachine mach = bulkData_->parallel();
   stk::ParallelMachine mach = *mpiComm_->getRawMpiComm();

   std::vector<stk::mesh::EntityId> local(commCount,0);

   // determine maximum ID for this processor for each entity type
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();
   for(unsigned i=0;i<entityRankCount;i++) {
      std::vector<stk::mesh::Entity*> entities;

      stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(i),entities);

      // determine maximum ID for this processor
      std::vector<stk::mesh::Entity*>::const_iterator itr;  
      for(itr=entities.begin();itr!=entities.end();++itr) {
         stk::mesh::EntityId id = (*itr)->identifier();
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
}

const double * STK_Interface::getNodeCoordinates(stk::mesh::EntityId nodeId) const
{
   stk::mesh::Entity * node = bulkData_->get_entity(getNodeRank(),nodeId);
   return stk::mesh::field_data(*coordinatesField_,*node);
}

const double * STK_Interface::getNodeCoordinates(stk::mesh::Entity * node) const
{
   return stk::mesh::field_data(*coordinatesField_,*node);
}

void STK_Interface::getSubcellIndices(unsigned entityRank,stk::mesh::EntityId elementId,
                                      std::vector<stk::mesh::EntityId> & subcellIds) const                       
{
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::Entity * cell = bulkData_->get_entity(elementRank,elementId);
   
   TEUCHOS_TEST_FOR_EXCEPTION(cell==0,std::logic_error,
                      "STK_Interface::getSubcellIndices: could not find element requested (GID = " << elementId << ")");

   stk::mesh::PairIterRelation subcells = cell->relations(entityRank);
   subcellIds.clear();
   subcellIds.resize(subcells.size(),0);

   // loop over relations and fill subcell vector
   stk::mesh::PairIterRelation::iterator iter;
   for(iter=subcells.begin();iter!=subcells.end();++iter) {
      TEUCHOS_ASSERT(iter->identifier()<subcellIds.size());
      subcellIds[iter->identifier()] = iter->entity()->identifier();
   }
}

void STK_Interface::getMyElements(std::vector<stk::mesh::Entity*> & elements) const
{
   // setup local ownership
   stk::mesh::Selector ownedPart = metaData_->locally_owned_part();

   // grab elements
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::get_selected_entities(ownedPart,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getMyElements(const std::string & blockID,std::vector<stk::mesh::Entity*> & elements) const
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

void STK_Interface::getMySides(const std::string & sideName,std::vector<stk::mesh::Entity*> & sides) const
{
   stk::mesh::Part * sidePart = getSideset(sideName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                      "Unknown side set \"" << sideName << "\"");

   stk::mesh::Selector side = *sidePart;
   stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & side;

   // grab elements
   stk::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity*> & sides) const
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

void STK_Interface::getMyNodes(const std::string & nodesetName,const std::string & blockName,std::vector<stk::mesh::Entity*> & nodes) const
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
   std::map<std::string, stk::mesh::Part*>::const_iterator sideItr;   // Element blocks
   for(sideItr=sidesets_.begin();sideItr!=sidesets_.end();++sideItr) 
      names.push_back(sideItr->first);
}

void STK_Interface::getNodesetNames(std::vector<std::string> & names) const
{
   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk::mesh::Part*>::const_iterator nodeItr;   // Element blocks
   for(nodeItr=nodesets_.begin();nodeItr!=nodesets_.end();++nodeItr) 
      names.push_back(nodeItr->first);
}

std::size_t STK_Interface::elementLocalId(stk::mesh::Entity * elmt) const
{
   const std::size_t * fieldCoords = stk::mesh::field_data(*localIdField_,*elmt);
   return fieldCoords[0];
}

std::size_t STK_Interface::elementLocalId(stk::mesh::EntityId gid) const
{
   stk::mesh::EntityRank elementRank = getElementRank();
   stk::mesh::Entity * elmt = bulkData_->get_entity(elementRank,gid);
   TEUCHOS_ASSERT(elmt->owner_rank()==procRank_);
   return elementLocalId(elmt);
}


std::string STK_Interface::containingBlockId(stk::mesh::Entity * elmt)
{
   std::map<std::string,stk::mesh::Part*>::const_iterator itr;
   for(itr=elementBlocks_.begin();itr!=elementBlocks_.end();++itr)
      if(elmt->bucket().member(*itr->second))
         return itr->first;
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

Teuchos::RCP<const std::vector<stk::mesh::Entity*> > STK_Interface::getElementsOrderedByLID() const
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
      block = &metaData_->declare_part(name,stk::mesh::fem::CellTopology(ctData));
   }

   // construct cell topology object for this block
   Teuchos::RCP<shards::CellTopology> ct
         = Teuchos::rcp(new shards::CellTopology(ctData));

   // add element block part and cell topology
   elementBlocks_.insert(std::make_pair(name,block));
   elementBlockCT_.insert(std::make_pair(name,ct));
}

void STK_Interface::initializeFromMetaData()
{
   dimension_ = metaData_->spatial_dimension();   

   // declare coordinates and node parts
   coordinatesField_ = &metaData_->declare_field<VectorFieldType>(coordsString);
   processorIdField_ = &metaData_->declare_field<ProcIdFieldType>("PROC_ID");
   localIdField_     = &metaData_->declare_field<LocalIdFieldType>("LOCAL_ID");

   // stk::mesh::put_field( *coordinatesField_ , getNodeRank() , metaData_->universal_part() );

   nodesPart_        = &metaData_->declare_part(nodesString,getNodeRank());
   nodesPartVec_.push_back(nodesPart_);
}

void STK_Interface::buildLocalElementIDs()
{
   currentLocalId_ = 0;
   
   orderedElementVector_ = Teuchos::null; // forces rebuild of ordered lists

   // might be better (faster) to do this by buckets
   Teuchos::RCP<std::vector<stk::mesh::Entity*> > elements 
      = Teuchos::rcp(new std::vector<stk::mesh::Entity*>);
   getMyElements(*elements);
 
   for(std::size_t index=0;index<elements->size();++index) {
      stk::mesh::Entity & element = *((*elements)[index]);

      // set processor rank
      ProcIdData * procId = stk::mesh::field_data(*processorIdField_,element);
      procId[0] = Teuchos::as<ProcIdData>(procRank_);

      // set local element ID
      std::size_t * localId = stk::mesh::field_data(*localIdField_,element);
      localId[0] = currentLocalId_;
      currentLocalId_++;
   }

   orderedElementVector_ = elements; 
}

Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > 
STK_Interface::getPeriodicNodePairing() const
{
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > vec;
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & matchers = getPeriodicBCVector();

   // build up the vectors by looping over the matched pair
   for(std::size_t m=0;m<matchers.size();m++)
      vec = matchers[m]->getMatchedPair(*this,vec);

   return vec;
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

} // end namespace panzer_stk
