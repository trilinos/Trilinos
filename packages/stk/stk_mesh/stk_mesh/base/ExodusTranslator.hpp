// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef EXODUSTRANSLATOR_HPP_
#define EXODUSTRANSLATOR_HPP_

#include <stk_topology/topology.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_mesh/base/FindRestriction.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

inline bool is_elem_block(const stk::mesh::Part& part)
{
  return (part.primary_entity_rank() == stk::topology::ELEMENT_RANK && part.id() > 0);
}

inline bool is_elem_block_ptr(const stk::mesh::Part* part) {
  return part!=nullptr && is_elem_block(*part);
}

inline bool has_subset_elem_block(const stk::mesh::Part &part)
{
  const PartVector& subsets = part.subsets();
  return std::any_of(subsets.cbegin(), subsets.cend(), is_elem_block_ptr);
}

}//namespace impl

inline bool is_element_block(const stk::mesh::Part &part)
{
    return (impl::is_elem_block(part) && !impl::has_subset_elem_block(part));
}

inline bool is_node_set(const stk::mesh::Part &part)
{
    return (part.primary_entity_rank() == stk::topology::NODE_RANK && part.id() > 0 );
}

inline bool has_super_set_face_part(const stk::mesh::Part &part)
{
  const PartVector& supersets = part.supersets();
  for(const Part* superset : supersets) {
    if(!stk::mesh::is_auto_declared_part(*superset) &&
       (superset->primary_entity_rank() == stk::topology::FACE_RANK ||
        superset->primary_entity_rank() == stk::topology::EDGE_RANK))
    {
      return true;
    }
  }
  return false;
}

inline bool is_side_set(const stk::mesh::Part &part)
{
  if (part.id() > 0) {
    bool isFacePart = part.primary_entity_rank() == stk::topology::FACE_RANK
            || part.primary_entity_rank() == stk::topology::EDGE_RANK;
    return isFacePart && !has_super_set_face_part(part);
  }
  return false;
}

inline void fill_element_block_parts(const MetaData& meta, stk::topology elemTopo,
                                     PartVector& elemBlockParts,
                                     bool sortById=true)
{
  elemBlockParts.clear();
  const PartVector &parts = meta.get_mesh_parts();
  for (Part* part : parts) {
    if (elemTopo == stk::topology::INVALID_TOPOLOGY) {
      if (is_element_block(*part)) {
        elemBlockParts.push_back(part);
      }   
    }   
    else {
      stk::topology partTopo = part->topology();
      if (is_element_block(*part) && partTopo == elemTopo) {
        elemBlockParts.push_back(part);
      }   
    }   
  }
  if (sortById) {
    stk::util::sort_and_unique(elemBlockParts, stk::mesh::PartLessById());
  }
}

inline stk::mesh::Part* get_element_block_part(const stk::mesh::BulkData& bulkData,
                                               const std::vector<stk::mesh::PartOrdinal>& partOrdinals,
                                               stk::mesh::EntityId elemId)
{
  stk::mesh::Part* elementBlockPart = nullptr;;
  const stk::mesh::PartVector& allParts = bulkData.mesh_meta_data().get_parts();
  unsigned blockCounter = 0;

  for(stk::mesh::PartOrdinal partOrdinal : partOrdinals)
  {
    stk::mesh::Part* part = allParts[partOrdinal];
    if(stk::mesh::is_element_block(*part))
    {
      STK_ThrowRequireMsg(blockCounter==0, "element global-id "<<elemId<<" associated with 2 element-blocks: '"<<elementBlockPart->name()<<"' (id="<<elementBlockPart->id()<<"), and '"<<part->name()<<"' (id="<<part->id()<<"). Each element should only be in 1 element-block.");
      elementBlockPart = part;
      blockCounter++;
    }
  }

  STK_ThrowRequireMsg(elementBlockPart != nullptr, "Failed to find element-block.");
  return elementBlockPart;
}

inline stk::mesh::Part* get_element_block_part(const stk::mesh::BulkData& bulkData,
                                               stk::mesh::Entity element)
{
  const stk::mesh::Bucket& bucket = bulkData.bucket(element);
  STK_ThrowAssertMsg(bucket.entity_rank() == stk::topology::ELEM_RANK, "get_element_block_part must be called with entity of rank stk::topology::ELEM_RANK, not "<<bucket.entity_rank());

  const stk::mesh::PartVector& parts = bucket.supersets();
  stk::mesh::Part* elementBlockPart = nullptr;;
  unsigned blockCounter = 0;
  for(stk::mesh::Part *part : parts)
  {
    if(stk::mesh::is_element_block(*part))
    {
      STK_ThrowRequireMsg(blockCounter==0, "element global-id "<<bulkData.identifier(element)<<" associated with 2 element-blocks: '"<<elementBlockPart->name()<<"' (id="<<elementBlockPart->id()<<"), and '"<<part->name()<<"' (id="<<part->id()<<"). Each element should only be in 1 element-block.");
      elementBlockPart = part;
      blockCounter++;
    }
  }

  STK_ThrowRequireMsg(elementBlockPart != nullptr, "Failed to find element-block for element global-id "<<bulkData.identifier(element));
  return elementBlockPart;
}

class ExodusTranslator
{
public:
    typedef int64_t IdType;

    ExodusTranslator(const stk::mesh::BulkData& b) : bulkData(b) {}

    size_t get_number_global_entities(stk::mesh::EntityRank rank) const
    {
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(bulkData, counts);
        return counts[rank];
    }

    size_t get_number_global_elements() const
    {
        return get_number_global_entities(stk::topology::ELEMENT_RANK);
    }

    size_t get_number_global_nodes() const
    {
        return get_number_global_entities(stk::topology::NODE_RANK);
    }

    size_t get_number_element_blocks() const
    {
        return get_element_block_parts().size();
    }

    size_t get_number_node_sets() const
    {
        return get_node_set_parts().size();
    }

    size_t get_number_side_sets() const
    {
        return get_side_set_parts().size();
    }

    void fill_node_set_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(get_node_set_parts(), ids);
    }

    void fill_node_set_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(get_node_set_parts(), names);
    }

    void fill_side_set_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(get_side_set_parts(), ids);
    }

    void fill_side_set_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(get_side_set_parts(), names);
    }

    void fill_element_block_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(get_element_block_parts(), ids);
    }

    void fill_element_block_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(get_element_block_parts(), names);
    }

    size_t get_local_num_entities_for_id(int setId, stk::mesh::EntityRank rank)
    {
        stk::mesh::Selector localAndShared = bulkData.mesh_meta_data().locally_owned_part() | bulkData.mesh_meta_data().globally_shared_part();
        return get_local_num_entitites_for_id_and_selector(setId, rank, localAndShared);
    }

    size_t get_global_num_entities_for_id(int setId, stk::mesh::EntityRank rank)
    {
        size_t numLocal = get_local_num_entitites_for_id_and_selector(setId, rank, bulkData.mesh_meta_data().locally_owned_part());
        return stk::get_global_sum(bulkData.parallel(), numLocal);
    }

    size_t get_global_num_distribution_factors_in_side_set(int setId) const
    {
        const stk::mesh::Part *part = get_exodus_part_of_rank(setId, bulkData.mesh_meta_data().side_rank());
        stk::mesh::Selector localAndShared = bulkData.mesh_meta_data().locally_owned_part() | bulkData.mesh_meta_data().globally_shared_part();

        std::vector<stk::mesh::Entity> entities;
        stk::mesh::get_selected_entities(*part & localAndShared, bulkData.buckets(bulkData.mesh_meta_data().side_rank()), entities);

        size_t numDF=0;
        for(size_t i = 0; i < entities.size(); i++)
            numDF += bulkData.num_nodes(entities[i]);

        return stk::get_global_sum(bulkData.parallel(), numDF);
    }

    stk::mesh::PartVector get_element_block_parts() const
    {
        stk::mesh::PartVector elementBlockParts;
        fill_element_block_parts(bulkData.mesh_meta_data(), stk::topology::INVALID_TOPOLOGY, elementBlockParts, false);
        return elementBlockParts;
    }

    stk::mesh::PartVector get_side_set_parts() const
    {
        stk::mesh::PartVector sideParts;
        const stk::mesh::PartVector &parts = bulkData.mesh_meta_data().get_mesh_parts();
        for(size_t i = 0; i < parts.size(); i++)
            if(is_side_set(*parts[i]))
                sideParts.push_back(parts[i]);
        return sideParts;
    }

    stk::mesh::PartVector get_node_set_parts() const
    {
        stk::mesh::PartVector nodeParts;
        const stk::mesh::PartVector &parts = bulkData.mesh_meta_data().get_mesh_parts();
        for(size_t i = 0; i < parts.size(); i++)
            if(is_node_set(*parts[i]))
                nodeParts.push_back(parts[i]);
        return nodeParts;
    }

    const stk::mesh::Part* get_exodus_part_of_rank(int id, stk::mesh::EntityRank rank) const
    {
        if(rank == stk::topology::NODE_RANK)
            return get_exodus_part_from_map(get_node_set_parts(), id);
        if(rank == bulkData.mesh_meta_data().side_rank())
            return get_exodus_part_from_map(get_side_set_parts(), id);
        if(rank == stk::topology::ELEM_RANK)
            return get_exodus_part_from_map(get_element_block_parts(), id);
        return nullptr;
    }

protected:
    size_t get_local_num_entitites_for_id_and_selector(int setId, stk::mesh::EntityRank rank, stk::mesh::Selector sel)
    {
        const stk::mesh::Part *part = get_exodus_part_of_rank(setId, rank);
        return stk::mesh::count_selected_entities(*part & sel, bulkData.buckets(rank));
    }

    template<typename ExodusPartMap, typename ExodusIdVector>
    void fill_exodus_ids(const ExodusPartMap& exodusParts, ExodusIdVector &exodusIds) const
    {
        exodusIds.clear();
        exodusIds.reserve(exodusParts.size());

        typename ExodusPartMap::const_iterator iter = exodusParts.begin();
        for(; iter != exodusParts.end(); iter++)
        {
            exodusIds.push_back((*iter)->id());
        }
    }

    template<typename ExodusPartMap>
    void fill_exodus_names(const ExodusPartMap& exodusParts, std::vector<std::string> &exodusNames) const
    {
        exodusNames.clear();
        exodusNames.reserve(exodusParts.size());

        typename ExodusPartMap::const_iterator iter = exodusParts.begin();
        for(; iter != exodusParts.end(); iter++)
        {
            exodusNames.push_back((*iter)->name());
        }
    }

    template<typename ExodusPartMap, typename ExodusId>
    const stk::mesh::Part* get_exodus_part_from_map(const ExodusPartMap& exodusParts, const ExodusId id) const
    {
        stk::mesh::Part *part = NULL;
        typename ExodusPartMap::const_iterator iter = exodusParts.begin();
        for(; iter != exodusParts.end(); iter++)
        {
            if((*iter)->id() == id)
            {
                part = *iter;
                break;
            }
        }
        return part;
    }

    const stk::mesh::BulkData& bulkData;
};

}
}

#endif /* EXODUSTRANSLATOR_HPP_ */
