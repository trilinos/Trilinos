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

#include <Panzer_config.hpp>
#include <Panzer_STK_Interface.hpp>

#include <Teuchos_as.hpp>

#include <limits>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>
#include <stk_rebalance_utils/RebalanceUtils.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#ifdef HAVE_IOSS
#include <Ionit_Initializer.h>
#include <stk_io/IossBridge.hpp>
#endif

#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include <set>

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk_classic {

ElementDescriptor::ElementDescriptor() {}
ElementDescriptor::ElementDescriptor(stk_classic::mesh::EntityId gid,const std::vector<stk_classic::mesh::EntityId> & nodes)
   : gid_(gid), nodes_(nodes) {}
ElementDescriptor::~ElementDescriptor() {}

/** Constructor function for building the element descriptors.
  */ 
Teuchos::RCP<ElementDescriptor> 
buildElementDescriptor(stk_classic::mesh::EntityId elmtId,std::vector<stk_classic::mesh::EntityId> & nodes)
{
   return Teuchos::rcp(new ElementDescriptor(elmtId,nodes));
}

const std::string STK_Interface::coordsString = "coordinates";
const std::string STK_Interface::nodesString = "nodes";
const std::string STK_Interface::edgesString = "edges";

STK_Interface::STK_Interface()
   : dimension_(0), initialized_(false), currentLocalId_(0), initialStateTime_(0.0), currentStateTime_(0.0), useFieldCoordinates_(false), useLowerCase_(false)
{
   metaData_ = rcp(new stk_classic::mesh::fem::FEMMetaData());
}

STK_Interface::STK_Interface(unsigned dim)
   : dimension_(dim), initialized_(false), currentLocalId_(0), useFieldCoordinates_(false), useLowerCase_(false)
{
   std::vector<std::string> entity_rank_names = stk_classic::mesh::fem::entity_rank_names(dimension_);
   entity_rank_names.push_back("FAMILY_TREE");

   metaData_ = rcp(new stk_classic::mesh::fem::FEMMetaData());
   metaData_->FEM_initialize(dimension_,entity_rank_names);

   initializeFromMetaData();
}

void STK_Interface::addSideset(const std::string & name,const CellTopologyData * ctData)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk_classic::mesh::Part * sideset = metaData_->get_part(name);
   if(sideset==NULL)
      sideset = &metaData_->declare_part(name,stk_classic::mesh::fem::CellTopology(ctData)); 
   sidesets_.insert(std::make_pair(name,sideset));
}

void STK_Interface::addNodeset(const std::string & name)
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0);

   stk_classic::mesh::Part * nodeset = metaData_->get_part(name);
   if(nodeset==NULL) {
      const CellTopologyData * ctData = shards::getCellTopologyData<shards::Node>();
      nodeset = &metaData_->declare_part(name,stk_classic::mesh::fem::CellTopology(ctData)); 
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

         SolutionFieldType * field = metaData_->get_field<SolutionFieldType>(dispName);
         if(field==0) {
            field = &metaData_->declare_field<SolutionFieldType>(dispName);     
         }
         fieldNameToSolution_[key] = field;
      }
   }
}

void STK_Interface::initialize(stk_classic::ParallelMachine parallelMach,bool setupIO) 
{
   TEUCHOS_ASSERT(not initialized_);
   TEUCHOS_ASSERT(dimension_!=0); // no zero dimensional meshes!

   if(mpiComm_==Teuchos::null)
      mpiComm_ = getSafeCommunicator(parallelMach);

   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::EntityRank nodeRank = getNodeRank();
   stk_classic::mesh::EntityRank edgeRank = getEdgeRank();

   procRank_ = stk_classic::parallel_machine_rank(*mpiComm_->getRawMpiComm());

   // associating the field with a part: universal part!
   stk_classic::mesh::put_field( *coordinatesField_ , nodeRank, metaData_->universal_part(), getDimension());
   stk_classic::mesh::put_field( *edgesField_ , edgeRank, metaData_->universal_part(), getDimension());
   stk_classic::mesh::put_field( *processorIdField_ , elementRank, metaData_->universal_part());
   stk_classic::mesh::put_field( *loadBalField_ , elementRank, metaData_->universal_part());

   initializeFieldsInSTK(fieldNameToSolution_,nodeRank,setupIO);
   initializeFieldsInSTK(fieldNameToCellField_,elementRank,setupIO);

#ifdef HAVE_IOSS
   if(setupIO) {
      // setup Exodus file IO
      /////////////////////////////////////////

      // add element blocks
      {
         std::map<std::string, stk_classic::mesh::Part*>::iterator itr;
         for(itr=elementBlocks_.begin();itr!=elementBlocks_.end();++itr) 
            if(!stk_classic::io::is_part_io_part(*itr->second))
               stk_classic::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add side sets
      {
         std::map<std::string, stk_classic::mesh::Part*>::iterator itr;
         for(itr=sidesets_.begin();itr!=sidesets_.end();++itr) 
            if(!stk_classic::io::is_part_io_part(*itr->second))
               stk_classic::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }

      // add node sets
      {
         std::map<std::string, stk_classic::mesh::Part*>::iterator itr;
         for(itr=nodesets_.begin();itr!=nodesets_.end();++itr) 
            if(!stk_classic::io::is_part_io_part(*itr->second))
               stk_classic::io::put_io_part_attribute(*itr->second); // this can only be called once per part
      }
   
      // add nodes 
      if(!stk_classic::io::is_part_io_part(*nodesPart_))
	stk_classic::io::put_io_part_attribute(*nodesPart_);

      stk_classic::io::set_field_role(*coordinatesField_, Ioss::Field::MESH);
      stk_classic::io::set_field_role(*edgesField_, Ioss::Field::MESH);
      stk_classic::io::set_field_role(*processorIdField_, Ioss::Field::TRANSIENT);
      // stk_classic::io::set_field_role(*loadBalField_, Ioss::Field::TRANSIENT);
   }
#endif

   metaData_->commit();
   if(bulkData_==Teuchos::null)
      instantiateBulkData(*mpiComm_->getRawMpiComm());

   initialized_ = true;
}

void STK_Interface::initializeFieldsInSTK(const std::map<std::pair<std::string,std::string>,SolutionFieldType*> & nameToField,
                                          stk_classic::mesh::EntityRank rank,bool setupIO)
{
   std::set<SolutionFieldType*> uniqueFields;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*>::const_iterator fieldIter;
   for(fieldIter=nameToField.begin();fieldIter!=nameToField.end();++fieldIter)
      uniqueFields.insert(fieldIter->second); // this makes setting up IO easier!

   {
      std::set<SolutionFieldType*>::const_iterator uniqueFieldIter;
      for(uniqueFieldIter=uniqueFields.begin();uniqueFieldIter!=uniqueFields.end();++uniqueFieldIter)
         stk_classic::mesh::put_field(*(*uniqueFieldIter),rank,metaData_->universal_part());
   }

#ifdef HAVE_IOSS
   if(setupIO) {
      // add solution fields
      std::set<SolutionFieldType*>::const_iterator uniqueFieldIter;
      for(uniqueFieldIter=uniqueFields.begin();uniqueFieldIter!=uniqueFields.end();++uniqueFieldIter)
         stk_classic::io::set_field_role(*(*uniqueFieldIter), Ioss::Field::TRANSIENT);
   }
#endif
}

void STK_Interface::instantiateBulkData(stk_classic::ParallelMachine parallelMach)
{
   TEUCHOS_ASSERT(bulkData_==Teuchos::null);
   if(mpiComm_==Teuchos::null)
      mpiComm_ = getSafeCommunicator(parallelMach);

   bulkData_ = rcp(new stk_classic::mesh::BulkData(stk_classic::mesh::fem::FEMMetaData::get_meta_data(*metaData_),*mpiComm_->getRawMpiComm()));
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

void STK_Interface::addNode(stk_classic::mesh::EntityId gid, const std::vector<double> & coord)
{
   TEUCHOS_TEST_FOR_EXCEPTION(not isModifiable(),std::logic_error,
                      "STK_Interface::addNode: STK_Interface must be modifiable to add a node");
   TEUCHOS_TEST_FOR_EXCEPTION(not coord.size()==getDimension(),std::logic_error,
                      "STK_Interface::addNode: number of coordinates in vector must mation dimension");
   TEUCHOS_TEST_FOR_EXCEPTION(gid==0,std::logic_error,
                      "STK_Interface::addNode: STK has STUPID restriction of no zero GIDs, pick something else");
   stk_classic::mesh::EntityRank nodeRank = getNodeRank();

   stk_classic::mesh::Entity & node = bulkData_->declare_entity(nodeRank,gid,nodesPartVec_);

   // set coordinate vector
   double * fieldCoords = stk_classic::mesh::field_data(*coordinatesField_,node);
   for(std::size_t i=0;i<coord.size();++i)
      fieldCoords[i] = coord[i];
}

void STK_Interface::addEntityToSideset(stk_classic::mesh::Entity & entity,stk_classic::mesh::Part * sideset)
{
   std::vector<stk_classic::mesh::Part*> sidesetV;
   sidesetV.push_back(sideset);

   bulkData_->change_entity_parts(entity,sidesetV);
}

void STK_Interface::addEntityToNodeset(stk_classic::mesh::Entity & entity,stk_classic::mesh::Part * nodeset)
{
   std::vector<stk_classic::mesh::Part*> nodesetV;
   nodesetV.push_back(nodeset);

   bulkData_->change_entity_parts(entity,nodesetV);
}

void STK_Interface::addElement(const Teuchos::RCP<ElementDescriptor> & ed,stk_classic::mesh::Part * block)
{
   std::vector<stk_classic::mesh::Part*> blockVec;
   blockVec.push_back(block);

   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::EntityRank nodeRank = getNodeRank();
   stk_classic::mesh::Entity & element = bulkData_->declare_entity(elementRank,ed->getGID(),blockVec);

   // build relations that give the mesh structure
   const std::vector<stk_classic::mesh::EntityId> & nodes = ed->getNodes();
   for(std::size_t i=0;i<nodes.size();++i) {
      // add element->node relation
      stk_classic::mesh::Entity * node = bulkData_->get_entity(nodeRank,nodes[i]);
      TEUCHOS_ASSERT(node!=0);
      bulkData_->declare_relation(element,*node,i);
   }

   ProcIdData * procId = stk_classic::mesh::field_data(*processorIdField_,element);
   procId[0] = Teuchos::as<ProcIdData>(procRank_);
}


void STK_Interface::addEdges()
{
   // loop over elements
   stk_classic::mesh::EntityRank edgeRank = getEdgeRank();
   stk_classic::mesh::EntityRank nodeRank = getNodeRank();
   std::vector<stk_classic::mesh::Entity*> localElmts;
   getMyElements(localElmts);
   std::vector<stk_classic::mesh::Entity*>::const_iterator itr;
   for(itr=localElmts.begin();itr!=localElmts.end();++itr) {
     stk_classic::mesh::Entity * element = (*itr);
     stk_classic::mesh::EntityId gid = element->identifier();
     std::vector<stk_classic::mesh::EntityId> subcellIds;
     getSubcellIndices(edgeRank,gid,subcellIds);

     for(std::size_t i=0;i<subcellIds.size();++i) {
       stk_classic::mesh::Entity * edge = bulkData_->get_entity(edgeRank,subcellIds[i]);
       stk_classic::mesh::PairIterRelation relations = edge->relations(nodeRank);

       double * node_coord_1 = stk_classic::mesh::field_data(*coordinatesField_,*(relations[0].entity()));
       double * node_coord_2 = stk_classic::mesh::field_data(*coordinatesField_,*(relations[1].entity()));

       // set coordinate vector
       double * edgeCoords = stk_classic::mesh::field_data(*edgesField_,*edge);
       for(std::size_t i=0;i<getDimension();++i)
          edgeCoords[i] = (node_coord_1[i]+node_coord_2[i])/2.0;
     }
   }
}


void STK_Interface::writeToExodus(const std::string & filename)
{
   PANZER_FUNC_TIME_MONITOR("STK_Interface::writeToExodus(filename)");

   #ifdef HAVE_IOSS
      TEUCHOS_ASSERT(mpiComm_!=Teuchos::null);
      stk_classic::ParallelMachine comm = *mpiComm_->getRawMpiComm();

      Ioss::Init::Initializer io;
      stk_classic::io::MeshData meshData;
      stk_classic::io::create_output_mesh(filename, comm, *bulkData_, meshData,getUseLowerCaseForIO());
      stk_classic::io::define_output_fields(meshData,*metaData_);

      stk_classic::io::process_output_request(meshData, *bulkData_, 0.0);
   #else 
      TEUCHOS_ASSERT(false);
   #endif
}

void STK_Interface::setupTransientExodusFile(const std::string & filename)
{
   PANZER_FUNC_TIME_MONITOR("STK_Interface::setupTransientExodusFile(filename)");

   #ifdef HAVE_IOSS
      TEUCHOS_ASSERT(mpiComm_!=Teuchos::null);
      stk_classic::ParallelMachine comm = *mpiComm_->getRawMpiComm();

      Ioss::Init::Initializer io;
      meshData_ = Teuchos::rcp(new stk_classic::io::MeshData);
      stk_classic::io::create_output_mesh(filename, comm, *bulkData_, *meshData_,getUseLowerCaseForIO());
      stk_classic::io::define_output_fields(*meshData_,*metaData_);
   #else 
      TEUCHOS_ASSERT(false);
   #endif
}

void STK_Interface::writeToExodus(double timestep)
{
   PANZER_FUNC_TIME_MONITOR("STK_Interface::writeToExodus(timestep)");

   #ifdef HAVE_IOSS
      if(meshData_!=Teuchos::null) {
        currentStateTime_ = timestep;
        stk_classic::io::process_output_request(*meshData_, *bulkData_, timestep);
      }
      else {
        Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
        out.setOutputToRootOnly(0);
        out << "WARNING: Exodus I/O has been disabled or not setup properly, not writting to Exodus" << std::endl;
      }
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

void STK_Interface::getElementsSharingNode(stk_classic::mesh::EntityId nodeId,std::vector<stk_classic::mesh::Entity *> & elements) const
{
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::EntityRank nodeRank = getNodeRank();

   // get all relations for node
   stk_classic::mesh::Entity * node = bulkData_->get_entity(nodeRank,nodeId);
   stk_classic::mesh::PairIterRelation relations = node->relations(elementRank);

   // extract elements sharing nodes
   stk_classic::mesh::PairIterRelation::iterator itr;
   for(itr=relations.begin();itr!=relations.end();++itr) {
      elements.push_back(itr->entity());
   }
}

void STK_Interface::getOwnedElementsSharingNode(stk_classic::mesh::Entity * node,std::vector<stk_classic::mesh::Entity *> & elements,
                                                                         std::vector<int> & localNodeId) const
{
   stk_classic::mesh::EntityRank elementRank = getElementRank();

   // get all relations for node
   stk_classic::mesh::PairIterRelation relations = node->relations(elementRank);

   // extract elements sharing nodes
   stk_classic::mesh::PairIterRelation::iterator itr;
   for(itr=relations.begin();itr!=relations.end();++itr) {
      stk_classic::mesh::Entity * element = itr->entity();
      
      // if owned by this processor 
      if(element->owner_rank() == procRank_) {
         elements.push_back(element);
         localNodeId.push_back(itr->identifier());
      }

   }
}

void STK_Interface::getOwnedElementsSharingNode(stk_classic::mesh::EntityId nodeId,std::vector<stk_classic::mesh::Entity *> & elements,
                                                                           std::vector<int> & localNodeId, unsigned int matchType) const
{
   stk_classic::mesh::EntityRank rank;
   if(matchType == 0)
     rank = getNodeRank();
   else if(matchType == 1)
     rank = getEdgeRank();
   else
     TEUCHOS_ASSERT(false);

   stk_classic::mesh::Entity * node = bulkData_->get_entity(rank,nodeId);

   getOwnedElementsSharingNode(node,elements,localNodeId);
}

void STK_Interface::getElementsSharingNodes(const std::vector<stk_classic::mesh::EntityId> nodeIds,std::vector<stk_classic::mesh::Entity *> & elements) const
{
   std::vector<stk_classic::mesh::Entity*> current;

   getElementsSharingNode(nodeIds[0],current); // fill it with elements touching first node
   std::sort(current.begin(),current.end());   // sort for intersection on the pointer 

   // find intersection with remaining nodes
   for(std::size_t n=1;n<nodeIds.size();++n) {
      // get elements associated with next node
      std::vector<stk_classic::mesh::Entity*> nextNode;
      getElementsSharingNode(nodeIds[n],nextNode); // fill it with elements touching first node
      std::sort(nextNode.begin(),nextNode.end());   // sort for intersection on the pointer ID

      // intersect next node elements with current element list
      std::vector<stk_classic::mesh::Entity*> intersection(std::min(nextNode.size(),current.size())); 
      std::vector<stk_classic::mesh::Entity*>::const_iterator endItr
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
   stk_classic::mesh::comm_mesh_counts(*bulkData_,entityCounts_);
}

void STK_Interface::buildMaxEntityIds()
{
   // developed to mirror "comm_mesh_counts" in stk_mesh/base/Comm.cpp

   const unsigned entityRankCount =  metaData_->entity_rank_count();
   const size_t   commCount        = 10; // entityRankCount

   TEUCHOS_ASSERT(entityRankCount<10);

   // stk_classic::ParallelMachine mach = bulkData_->parallel();
   stk_classic::ParallelMachine mach = *mpiComm_->getRawMpiComm();

   std::vector<stk_classic::mesh::EntityId> local(commCount,0);

   // determine maximum ID for this processor for each entity type
   stk_classic::mesh::Selector ownedPart = metaData_->locally_owned_part();
   for(unsigned i=0;i<entityRankCount;i++) {
      std::vector<stk_classic::mesh::Entity*> entities;

      stk_classic::mesh::get_selected_entities(ownedPart,bulkData_->buckets(i),entities);

      // determine maximum ID for this processor
      std::vector<stk_classic::mesh::Entity*>::const_iterator itr;  
      for(itr=entities.begin();itr!=entities.end();++itr) {
         stk_classic::mesh::EntityId id = (*itr)->identifier();
         if(id>local[i])
            local[i] = id;
      }
   }

   // get largest IDs across processors
   stk_classic::all_reduce(mach,stk_classic::ReduceMax<10>(&local[0]));
   maxEntityId_.assign(local.begin(),local.begin()+entityRankCount+1); 
}

std::size_t STK_Interface::getEntityCounts(unsigned entityRank) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(entityRank>=entityCounts_.size(),std::logic_error,
                      "STK_Interface::getEntityCounts: Entity counts do not include rank: " << entityRank);
                      
   return entityCounts_[entityRank];
}

stk_classic::mesh::EntityId STK_Interface::getMaxEntityId(unsigned entityRank) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(entityRank>=maxEntityId_.size(),std::logic_error,
                      "STK_Interface::getMaxEntityId: Max entity ids do not include rank: " << entityRank);
                      
   return maxEntityId_[entityRank];
}

void STK_Interface::buildSubcells()
{
   stk_classic::mesh::PartVector emptyPartVector;
   stk_classic::mesh::create_adjacent_entities(*bulkData_,emptyPartVector);

   buildEntityCounts();
   buildMaxEntityIds();

   addEdges();
}

const double * STK_Interface::getNodeCoordinates(stk_classic::mesh::EntityId nodeId) const
{
   stk_classic::mesh::Entity * node = bulkData_->get_entity(getNodeRank(),nodeId);
   return stk_classic::mesh::field_data(*coordinatesField_,*node);
}

const double * STK_Interface::getNodeCoordinates(stk_classic::mesh::Entity * node) const
{
   return stk_classic::mesh::field_data(*coordinatesField_,*node);
}

void STK_Interface::getSubcellIndices(unsigned entityRank,stk_classic::mesh::EntityId elementId,
                                      std::vector<stk_classic::mesh::EntityId> & subcellIds) const                       
{
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::Entity * cell = bulkData_->get_entity(elementRank,elementId);
   
   TEUCHOS_TEST_FOR_EXCEPTION(cell==0,std::logic_error,
                      "STK_Interface::getSubcellIndices: could not find element requested (GID = " << elementId << ")");

   stk_classic::mesh::PairIterRelation subcells = cell->relations(entityRank);
   subcellIds.clear();
   subcellIds.resize(subcells.size(),0);

   // loop over relations and fill subcell vector
   stk_classic::mesh::PairIterRelation::iterator iter;
   for(iter=subcells.begin();iter!=subcells.end();++iter) {
      TEUCHOS_ASSERT(iter->identifier()<subcellIds.size());
      subcellIds[iter->identifier()] = iter->entity()->identifier();
   }
}

void STK_Interface::getMyElements(std::vector<stk_classic::mesh::Entity*> & elements) const
{
   // setup local ownership
   stk_classic::mesh::Selector ownedPart = metaData_->locally_owned_part();

   // grab elements
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::get_selected_entities(ownedPart,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getMyElements(const std::string & blockID,std::vector<stk_classic::mesh::Entity*> & elements) const
{
   stk_classic::mesh::Part * elementBlock = getElementBlockPart(blockID);

   TEUCHOS_TEST_FOR_EXCEPTION(elementBlock==0,std::logic_error,"Could not find element block \"" << blockID << "\"");

   // setup local ownership
   // stk_classic::mesh::Selector block = *elementBlock;
   stk_classic::mesh::Selector ownedBlock = metaData_->locally_owned_part() & (*elementBlock);

   // grab elements
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getNeighborElements(std::vector<stk_classic::mesh::Entity*> & elements) const
{
   // setup local ownership
   stk_classic::mesh::Selector neighborBlock = (!metaData_->locally_owned_part());

   // grab elements
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::get_selected_entities(neighborBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getNeighborElements(const std::string & blockID,std::vector<stk_classic::mesh::Entity*> & elements) const
{
   stk_classic::mesh::Part * elementBlock = getElementBlockPart(blockID);

   TEUCHOS_TEST_FOR_EXCEPTION(elementBlock==0,std::logic_error,"Could not find element block \"" << blockID << "\"");

   // setup local ownership
   stk_classic::mesh::Selector neighborBlock = (!metaData_->locally_owned_part()) & (*elementBlock);

   // grab elements
   stk_classic::mesh::EntityRank elementRank = getElementRank();
   stk_classic::mesh::get_selected_entities(neighborBlock,bulkData_->buckets(elementRank),elements);
}

void STK_Interface::getMySides(const std::string & sideName,std::vector<stk_classic::mesh::Entity*> & sides) const
{
   stk_classic::mesh::Part * sidePart = getSideset(sideName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                      "Unknown side set \"" << sideName << "\"");

   stk_classic::mesh::Selector side = *sidePart;
   stk_classic::mesh::Selector ownedBlock = metaData_->locally_owned_part() & side;

   // grab elements
   stk_classic::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk_classic::mesh::Entity*> & sides) const
{
   stk_classic::mesh::Part * sidePart = getSideset(sideName);
   stk_classic::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,SidesetException,
                      "Unknown side set \"" << sideName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk_classic::mesh::Selector side = *sidePart;
   stk_classic::mesh::Selector block = *elmtPart;
   stk_classic::mesh::Selector ownedBlock = metaData_->locally_owned_part() & block & side;

   // grab elements
   stk_classic::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getSideRank()),sides);
}

void STK_Interface::getMyNodes(const std::string & nodesetName,const std::string & blockName,std::vector<stk_classic::mesh::Entity*> & nodes) const
{
   stk_classic::mesh::Part * nodePart = getNodeset(nodesetName);
   stk_classic::mesh::Part * elmtPart = getElementBlockPart(blockName);
   TEUCHOS_TEST_FOR_EXCEPTION(nodePart==0,SidesetException,
                      "Unknown node set \"" << nodesetName << "\"");
   TEUCHOS_TEST_FOR_EXCEPTION(elmtPart==0,ElementBlockException,
                      "Unknown element block \"" << blockName << "\"");

   stk_classic::mesh::Selector nodeset = *nodePart;
   stk_classic::mesh::Selector block = *elmtPart;
   stk_classic::mesh::Selector ownedBlock = metaData_->locally_owned_part() & block & nodeset;

   // grab elements
   stk_classic::mesh::get_selected_entities(ownedBlock,bulkData_->buckets(getNodeRank()),nodes);
}

void STK_Interface::getElementBlockNames(std::vector<std::string> & names) const
{
   // TEUCHOS_ASSERT(initialized_); // all blocks must have been added

   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk_classic::mesh::Part*>::const_iterator blkItr;   // Element blocks
   for(blkItr=elementBlocks_.begin();blkItr!=elementBlocks_.end();++blkItr) 
      names.push_back(blkItr->first);
}

void STK_Interface::getSidesetNames(std::vector<std::string> & names) const
{
   // TEUCHOS_ASSERT(initialized_); // all blocks must have been added

   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk_classic::mesh::Part*>::const_iterator sideItr;   // Element blocks
   for(sideItr=sidesets_.begin();sideItr!=sidesets_.end();++sideItr) 
      names.push_back(sideItr->first);
}

void STK_Interface::getNodesetNames(std::vector<std::string> & names) const
{
   names.clear();

   // fill vector with automagically ordered string values
   std::map<std::string, stk_classic::mesh::Part*>::const_iterator nodeItr;   // Element blocks
   for(nodeItr=nodesets_.begin();nodeItr!=nodesets_.end();++nodeItr) 
      names.push_back(nodeItr->first);
}

std::size_t STK_Interface::elementLocalId(stk_classic::mesh::Entity * elmt) const
{
   return elementLocalId(elmt->identifier());
   // const std::size_t * fieldCoords = stk_classic::mesh::field_data(*localIdField_,*elmt);
   // return fieldCoords[0];
}

std::size_t STK_Interface::elementLocalId(stk_classic::mesh::EntityId gid) const
{
   // stk_classic::mesh::EntityRank elementRank = getElementRank();
   // stk_classic::mesh::Entity * elmt = bulkData_->get_entity(elementRank,gid);
   // TEUCHOS_ASSERT(elmt->owner_rank()==procRank_);
   // return elementLocalId(elmt);
   boost::unordered_map<stk_classic::mesh::EntityId,std::size_t>::const_iterator itr = localIDHash_.find(gid);
   TEUCHOS_ASSERT(itr!=localIDHash_.end());
   return itr->second;
}


std::string STK_Interface::containingBlockId(stk_classic::mesh::Entity * elmt)
{
   std::map<std::string,stk_classic::mesh::Part*>::const_iterator itr;
   for(itr=elementBlocks_.begin();itr!=elementBlocks_.end();++itr)
      if(elmt->bucket().member(*itr->second))
         return itr->first;
   return "";
}

stk_classic::mesh::Field<double> * STK_Interface::getSolutionField(const std::string & fieldName,
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

stk_classic::mesh::Field<double> * STK_Interface::getCellField(const std::string & fieldName,
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

Teuchos::RCP<const std::vector<stk_classic::mesh::Entity*> > STK_Interface::getElementsOrderedByLID() const
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

   stk_classic::mesh::Part * block = metaData_->get_part(name);
   if(block==0) {
      block = &metaData_->declare_part(name,stk_classic::mesh::fem::CellTopology(ctData));
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
   edgesField_       = &metaData_->declare_field<VectorFieldType>(edgesString);
   processorIdField_ = &metaData_->declare_field<ProcIdFieldType>("PROC_ID");
   loadBalField_     = &metaData_->declare_field<SolutionFieldType>("LOAD_BAL");

   // stk_classic::mesh::put_field( *coordinatesField_ , getNodeRank() , metaData_->universal_part() );

   nodesPart_        = &metaData_->declare_part(nodesString,getNodeRank());
   nodesPartVec_.push_back(nodesPart_);
   edgesPart_        = &metaData_->declare_part(edgesString,getEdgeRank());
   edgesPartVec_.push_back(edgesPart_);
}

void STK_Interface::buildLocalElementIDs()
{
   currentLocalId_ = 0;
   
   orderedElementVector_ = Teuchos::null; // forces rebuild of ordered lists

   // might be better (faster) to do this by buckets
   std::vector<stk_classic::mesh::Entity*> elements;
   getMyElements(elements);
 
   for(std::size_t index=0;index<elements.size();++index) {
      stk_classic::mesh::Entity & element = *(elements[index]);

      // set processor rank
      ProcIdData * procId = stk_classic::mesh::field_data(*processorIdField_,element);
      procId[0] = Teuchos::as<ProcIdData>(procRank_);

      localIDHash_[element.identifier()] = currentLocalId_;

      currentLocalId_++;
   }

   // copy elements into the ordered element vector
   orderedElementVector_ = Teuchos::rcp(new std::vector<stk_classic::mesh::Entity*>(elements));

   elements.clear();
   getNeighborElements(elements);

   for(std::size_t index=0;index<elements.size();++index) {
      stk_classic::mesh::Entity & element = *(elements[index]);

      // set processor rank
      ProcIdData * procId = stk_classic::mesh::field_data(*processorIdField_,element);
      procId[0] = Teuchos::as<ProcIdData>(procRank_);

      localIDHash_[element.identifier()] = currentLocalId_;

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

    std::vector<stk_classic::mesh::Entity*> elements;
    getMyElements(names[b],elements);

    for(std::size_t index=0;index<elements.size();++index) {
      // set local element ID
      double * loadBal = stk_classic::mesh::field_data(*loadBalField_,*elements[index]);
      loadBal[0] = blockWeight;
    }
  }
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
  for(axis=0;axis<blkItr->second.size();axis++) {
    if(blkItr->second[axis]==fieldName) 
      break; // found axis, break
  }
    
  if(axis>=blkItr->second.size())
    return false;
 
  return true;
}

std::pair<Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >, Teuchos::RCP<std::vector<unsigned int> > >
STK_Interface::getPeriodicNodePairing() const
{
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > vec;
   Teuchos::RCP<std::vector<unsigned int > > type_vec = rcp(new std::vector<unsigned int>);
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & matchers = getPeriodicBCVector();

   // build up the vectors by looping over the matched pair
   for(std::size_t m=0;m<matchers.size();m++){
      vec = matchers[m]->getMatchedPair(*this,vec);
      unsigned int type;
      if(matchers[m]->getType() == "coord")
        type = 0;
      else if(matchers[m]->getType() == "edge")
        type = 1;
      else
        TEUCHOS_ASSERT(false);
      type_vec->insert(type_vec->begin(),vec->size()-type_vec->size(),type);
   }

   return std::make_pair(vec,type_vec);
}

bool STK_Interface::validBlockId(const std::string & blockId) const
{
   std::map<std::string, stk_classic::mesh::Part*>::const_iterator blkItr = elementBlocks_.find(blockId);

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
   const stk_classic::mesh::FieldVector & fv = metaData_->get_fields(); 
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

Teuchos::RCP<Teuchos::MpiComm<int> > STK_Interface::getSafeCommunicator(stk_classic::ParallelMachine parallelMach) const
{
   MPI_Comm newComm;
   const int err = MPI_Comm_dup (parallelMach, &newComm);
   TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
     "panzer::STK_Interface: MPI_Comm_dup failed with error \""
     << Teuchos::mpiErrorCodeToString (err) << "\".");

   return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper (newComm,MPI_Comm_free)));
}

void STK_Interface::rebalance(const Teuchos::ParameterList & params)
{
  // make sure weights have been set (a local operation)
  applyElementLoadBalanceWeights();

  stk_classic::mesh::Selector selector(getMetaData()->universal_part());
  stk_classic::mesh::Selector owned_selector(getMetaData()->locally_owned_part());

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  out << "Load balance before: " << stk_classic::rebalance::check_balance(*getBulkData(), loadBalField_, getElementRank(), &selector) << std::endl;

  // perform reblance
  Teuchos::ParameterList graph;
  if(params.begin()!=params.end())
    graph.sublist(stk_classic::rebalance::Zoltan::default_parameters_name()) = params;
  stk_classic::rebalance::Zoltan zoltan_partition(*mpiComm_->getRawMpiComm(), getDimension(), graph);
  stk_classic::rebalance::rebalance(*getBulkData(), owned_selector, &getCoordinatesField(), loadBalField_, zoltan_partition);

  out << "Load balance after: " << stk_classic::rebalance::check_balance(*getBulkData(), loadBalField_, getElementRank(), &selector) << std::endl;

  currentLocalId_ = 0;
  orderedElementVector_ = Teuchos::null; // forces rebuild of ordered lists
}

} // end namespace panzer_stk_classic
