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

#include "Panzer_STK_WorksetFactory.hpp"

// Panzer includes
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_WorksetDescriptor.hpp"

// Panzer STK includes
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_WorksetUtilities.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"

namespace panzer_stk {

namespace {


void getSubcellElements(const panzer_stk::STK_Interface & mesh,
            const std::string & blockId,
            const std::vector<stk::mesh::Entity> & entities,
            std::vector<std::size_t> & localEntityIds,
            std::vector<stk::mesh::Entity> & elements)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;

    const size_t num_rels = bulkData.num_elements(entity);
    stk::mesh::Entity const* relations = bulkData.begin_elements(entity);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(entity);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = relations[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock = bucket.member(*blockPart);
      bool onProc = bucket.member(*ownedPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}

void getUniversalSubcellElements(const panzer_stk::STK_Interface & mesh,
         const std::string & blockId,
         const std::vector<stk::mesh::Entity> & entities,
         std::vector<std::size_t> & localEntityIds,
         std::vector<stk::mesh::Entity> & elements)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;

    const size_t num_rels = bulkData.num_elements(entity);
    stk::mesh::Entity const* element_rels = bulkData.begin_elements(entity);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(entity);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = element_rels[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock = bucket.member(*blockPart);
      bool onProc = bucket.member(*universalPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}


void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId,
                     const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds,
                     std::vector<stk::mesh::Entity> & elements)
{
   getSubcellElements(mesh,blockId,sides,localSideIds,elements);
}

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId_a,
                     const std::string & blockId_b,
                     const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds_a,
                     std::vector<stk::mesh::Entity> & elements_a,
                     std::vector<std::size_t> & localSideIds_b,
                     std::vector<stk::mesh::Entity> & elements_b)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart_a = mesh.getElementBlockPart(blockId_a);
  stk::mesh::Part * blockPart_b = mesh.getElementBlockPart(blockId_b);
  stk::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator sidesItr;
  for(sidesItr=sides.begin();sidesItr!=sides.end();++sidesItr) {
    stk::mesh::Entity side = *sidesItr;

     // these are used below the loop to insert into the appropriate vectors
    stk::mesh::Entity element_a = stk::mesh::Entity(), element_b = stk::mesh::Entity();
    std::size_t entityId_a=0, entityId_b=0;

    const size_t num_rels = bulkData.num_elements(side);
    stk::mesh::Entity const* element_rels = bulkData.begin_elements(side);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(side);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = element_rels[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock_a = bucket.member(*blockPart_a);
      bool inBlock_b = bucket.member(*blockPart_b);
      bool onProc = bucket.member(*ownedPart);
      bool unProc = bucket.member(*universalPart);

      if(inBlock_a && onProc) {
        TEUCHOS_ASSERT(element_a==stk::mesh::Entity()); // sanity check
        element_a = element;
        entityId_a = entityId;
      }
      if(inBlock_b && unProc) {
        TEUCHOS_ASSERT(element_b==stk::mesh::Entity()); // sanity check
        element_b = element;
        entityId_b = entityId;
      }
    }

    if(element_a!=stk::mesh::Entity() && element_b!=stk::mesh::Entity()) {      // add element and Side ID to output vectors
      elements_a.push_back(element_a);
      localSideIds_a.push_back(entityId_a);

      // add element and Side ID to output vectors
      elements_b.push_back(element_b);
      localSideIds_b.push_back(entityId_b);
    }
  }
}

void getNodeElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId,
                     const std::vector<stk::mesh::Entity> & nodes,
                     std::vector<std::size_t> & localNodeIds,
                     std::vector<stk::mesh::Entity> & elements)
{
   getSubcellElements(mesh,blockId,nodes,localNodeIds,elements);
}

void getSubcellEntities(const panzer_stk::STK_Interface & mesh,
            const std::vector<stk::mesh::Entity> & entities,
            std::vector<std::vector<stk::mesh::Entity> > & subcells)
{
  // exit if there is no work to do
  if(entities.size()==0) {
    subcells.clear();
    return;
  }

  stk::mesh::BulkData& bulkData = *mesh.getBulkData();
  stk::mesh::EntityRank master_rank = bulkData.entity_rank(entities[0]);

  std::vector<std::set<stk::mesh::Entity> > subcells_set(master_rank);

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;

    // sanity check, enforcing that there is only one rank
    TEUCHOS_ASSERT(bulkData.entity_rank(entity)==master_rank);

    for(int i=0; i<static_cast<int>(master_rank); i++) {
      stk::mesh::EntityRank const to_rank = static_cast<stk::mesh::EntityRank>(i);
      const size_t num_rels = bulkData.num_connectivity(entity, to_rank);
      stk::mesh::Entity const* relations = bulkData.begin(entity, to_rank);

      // for each relation insert the appropriate entity (into the set
      // which gurantees uniqueness
      for(std::size_t e=0; e<num_rels; ++e) {
        stk::mesh::Entity subcell = relations[e];

        subcells_set[i].insert(subcell);
      }
    }
  }

  // copy unique entities into vector
  subcells.clear();
  subcells.resize(subcells_set.size());
  for(std::size_t i=0;i<subcells_set.size();i++)
    subcells[i].assign(subcells_set[i].begin(),subcells_set[i].end());
}

void getSideElementCascade(const panzer_stk::STK_Interface & mesh,
                           const std::string & blockId,
                           const std::vector<stk::mesh::Entity> & sides,
                           std::vector<std::size_t> & localSubcellDim,
                           std::vector<std::size_t> & localSubcellIds,
                           std::vector<stk::mesh::Entity> & elements)
{
  // This is the alogrithm, for computing the side element
  // cascade. The requirements are that for a particular set of sides
  // we compute all elements and subcells where they touch the side. Note
  // that elements can be and will be repeated within this list.

  std::vector<std::vector<stk::mesh::Entity> > subcells;
  getSubcellEntities(mesh,sides,subcells);
  subcells.push_back(sides);

  // subcells now contains a unique list of faces, edges and nodes that
  // intersect with the sides

  for(std::size_t d=0;d<subcells.size();d++) {
    std::vector<std::size_t> subcellIds;
    std::vector<stk::mesh::Entity> subcellElements;

    // find elements connected to the subcells and their local subcell information
    getSubcellElements(mesh,blockId,subcells[d],subcellIds,subcellElements);

    // sanity check
    TEUCHOS_ASSERT(subcellIds.size()==subcellElements.size());

    // concatenate with found elements
    localSubcellDim.insert(localSubcellDim.end(),subcellElements.size(),d);
    localSubcellIds.insert(localSubcellIds.end(),subcellIds.begin(),subcellIds.end());
    elements.insert(elements.end(),subcellElements.begin(),subcellElements.end());
  }
}

Teuchos::RCP<std::vector<panzer::Workset> >
buildCascadeWorksets(const panzer_stk::STK_Interface & mesh,
                     const panzer::LocalMeshInfo & local_mesh,
                     const panzer::WorksetDescriptor & description,
                     const panzer::WorksetOptions & options)
{
  using Teuchos::RCP;

  // If this sideset doesn't exist, then return nothing
  if(local_mesh.sidesets.find(description.getElementBlock()) == local_mesh.sidesets.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>);
  const auto & sideset_map = local_mesh.sidesets.at(description.getElementBlock());
  if(sideset_map.find(description.getSideset()) == sideset_map.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>);

  const auto & sideset_info = sideset_map.at(description.getSideset());

  std::vector<stk::mesh::Entity> sideEntities;
  mesh.getAllSides(description.getSideset(),description.getElementBlock(),sideEntities);


  std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> > cell_ids;
  {
    std::unordered_map<int,size_t> local_id_map;
    for(size_t i=0; i<sideset_info.local_cells.size(); ++i)
      local_id_map[sideset_info.local_cells(i)] = i;

    std::vector<stk::mesh::Entity> elements;
    std::vector<std::size_t> local_subcell_ids, subcell_dim;
    getSideElementCascade(mesh, description.getElementBlock(),
                          sideEntities,subcell_dim,local_subcell_ids,elements);

    // build local cell_ids, mapped by local side id
    for(std::size_t elm=0;elm<elements.size();++elm)
      cell_ids[std::make_pair(subcell_dim[elm],local_subcell_ids[elm])].push_back(local_id_map.at(mesh.elementLocalId(elements[elm])));
  }

  TEUCHOS_ASSERT(description.getWorksetSize() != panzer::NO_ELEMENTS);

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor

  // worksets to be returned
  Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>);

  // loop over each combination of subcell dimension/index
  for(auto itr=cell_ids.begin(); itr!=cell_ids.end();++itr) {

    // Sells associated with this combination
    const auto & cells = itr->second;
    const size_t num_cells = cells.size();

    // If no cells were found, skip
    if(num_cells == 0)
      continue;

    // Define the maximum workset size
    const size_t max_partition_size = (description.getWorksetSize() < 0) ? num_cells : description.getWorksetSize();

    // Storage for cells owned by sideset
    std::vector<size_t> partition_cells;
    partition_cells.reserve(max_partition_size);

    size_t cell_count = 0;
    while(cell_count < num_cells){

      panzer::LocalMeshPartition partition;

      size_t this_partition_size = max_partition_size;
      if(cell_count + this_partition_size > num_cells)
        this_partition_size = num_cells - cell_count;

      // Not sure if clear wipes the reserve...
      partition_cells.clear();
      for(size_t i=0; i<this_partition_size;++i)
        partition_cells.push_back(cells[cell_count+i]);

      // Setup the partitions for the workset
      panzer::partitioning_utilities::setupSubLocalMeshInfo(sideset_info, partition_cells, partition);

      // Now add some additional info to the partition
      partition.element_block_name = description.getElementBlock();
      partition.sideset_name = description.getSideset();
      partition.cell_topology = sideset_info.cell_topology;
      partition.subcell_dimension = itr->first.first;
      partition.subcell_index = itr->first.second;

      worksets->push_back(panzer::Workset());
      worksets->back().setup(partition,options);

      cell_count += this_partition_size;
    }
  }
  return worksets;
}

}

/** Set mesh
  */
void
WorksetFactory::
setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh)
{ mesh_ = mesh; }


Teuchos::RCP<std::vector<panzer::Workset> >
WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & description) const
{

  // There is a special case of 'cascade' that we need to consider
  if(description.buildCascade()){
    TEUCHOS_ASSERT(not description.sideAssembly());
    TEUCHOS_ASSERT(not description.connectsElementBlocks());

    panzer::WorksetOptions options;
    options.align_side_points_ = false;
    options.side_assembly_ = false;
    options.orientations_ = this->getOrientationsInterface();

    // TODO: To get this running with the non-stk path, we need to first get the full connectivity into the LocalMeshInfo object
    // This is something we need to do anyway to better represent connectivity in the worksets
    return buildCascadeWorksets(*mesh_, getMeshInfo(), description, options);
  }

  // This supports many kinds of worksets - except for cascade...
  return panzer::buildWorksets(getMeshInfo(), description, this->getOrientationsInterface());
}

const panzer::LocalMeshInfo &
WorksetFactory::
getMeshInfo() const
{
  if(mesh_info_.is_null()){
    TEUCHOS_TEST_FOR_EXCEPT(mesh_.is_null());
    mesh_info_ = panzer_stk::generateLocalMeshInfo(*mesh_);
  }
  return *mesh_info_;
}

}
