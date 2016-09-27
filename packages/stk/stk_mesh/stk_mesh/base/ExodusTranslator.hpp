#ifndef EXODUSTRANSLATOR_HPP_
#define EXODUSTRANSLATOR_HPP_

#include <stk_topology/topology.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FindRestriction.hpp>

namespace stk
{
namespace mesh
{

inline bool is_element_block(const stk::mesh::Part &part)
{
    return (part.primary_entity_rank() == stk::topology::ELEMENT_RANK && part.id() > 0 );
}
inline bool is_node_set(const stk::mesh::Part &part)
{
    return (part.primary_entity_rank() == stk::topology::NODE_RANK && part.id() > 0 );
}

inline bool has_super_set_face_part(const stk::mesh::Part &part)
{
    size_t numSuperSets = part.supersets().size();
    for(size_t i = 0; i < numSuperSets; i++)
    {
        if(     !stk::mesh::is_auto_declared_part(*part.supersets()[i]) &&
                (part.supersets()[i]->primary_entity_rank() == stk::topology::FACE_RANK ||
                 part.supersets()[i]->primary_entity_rank() == stk::topology::EDGE_RANK))
        {
            return true;
        }
    }
    return false;
}
inline bool is_side_set(const stk::mesh::Part &part)
{
    bool isFacePart = part.primary_entity_rank() == stk::topology::FACE_RANK
            || part.primary_entity_rank() == stk::topology::EDGE_RANK;
    return isFacePart && !has_super_set_face_part(part) && part.id() > 0;
}

class ExodusTranslator
{
public:
    typedef int64_t IdType;

    ExodusTranslator(const stk::mesh::BulkData& b) : mBulkData(b)
    {
        cache_mesh_constructs();
    }

    size_t get_number_global_entities(stk::mesh::EntityRank rank) const
    {
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(mBulkData, counts);
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
        return mElementBlockParts.size();
    }

    size_t get_number_node_sets() const
    {
        return mNodeSetParts.size();
    }

    size_t get_number_side_sets() const
    {
        return mSideSetParts.size();
    }

    void fill_node_set_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(mNodeSetParts, ids);
    }

    void fill_node_set_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(mNodeSetParts, names);
    }

    void fill_side_set_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(mSideSetParts, ids);
    }

    void fill_side_set_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(mSideSetParts, names);
    }

    void fill_element_block_ids(std::vector<IdType> &ids) const
    {
        fill_exodus_ids(mElementBlockParts, ids);
    }

    void fill_element_block_names(std::vector<std::string> &names) const
    {
        fill_exodus_names(mElementBlockParts, names);
    }

    size_t get_local_num_entities_for_id(int setId, stk::mesh::EntityRank rank)
    {
        stk::mesh::Selector localAndShared = mBulkData.mesh_meta_data().locally_owned_part() | mBulkData.mesh_meta_data().globally_shared_part();
        return get_local_num_entitites_for_id_and_selector(setId, rank, localAndShared);
    }

    size_t get_global_num_entities_for_id(int setId, stk::mesh::EntityRank rank)
    {
        size_t numLocal = get_local_num_entitites_for_id_and_selector(setId, rank, mBulkData.mesh_meta_data().locally_owned_part());
        return stk::get_global_sum(mBulkData.parallel(), numLocal);
    }

    size_t get_global_num_distribution_factors_in_side_set(int setId) const
    {
        const stk::mesh::Part *part = get_exodus_part_of_rank(setId, mBulkData.mesh_meta_data().side_rank());
        stk::mesh::Selector localAndShared = mBulkData.mesh_meta_data().locally_owned_part() | mBulkData.mesh_meta_data().globally_shared_part();

        std::vector<stk::mesh::Entity> entities;
        stk::mesh::get_selected_entities(*part & localAndShared, mBulkData.buckets(mBulkData.mesh_meta_data().side_rank()), entities);

        size_t numDF=0;
        for(size_t i = 0; i < entities.size(); i++)
            numDF += mBulkData.num_nodes(entities[i]);

        return stk::get_global_sum(mBulkData.parallel(), numDF);
    }

    const stk::mesh::PartVector & get_node_set_parts() { return mNodeSetParts; }
    const stk::mesh::PartVector & get_side_set_parts() { return mSideSetParts; }
    const stk::mesh::PartVector & get_element_block_parts() { return mElementBlockParts; }

    const stk::mesh::Part* get_exodus_part_of_rank(int id, stk::mesh::EntityRank rank) const
    {
        if(rank == stk::topology::NODE_RANK)
            return get_exodus_part_from_map(mNodeSetParts, id);
        if(rank == mBulkData.mesh_meta_data().side_rank())
            return get_exodus_part_from_map(mSideSetParts, id);
        if(rank == stk::topology::ELEM_RANK)
            return get_exodus_part_from_map(mElementBlockParts, id);
        return nullptr;
    }

private:
    size_t get_local_num_entitites_for_id_and_selector(int setId, stk::mesh::EntityRank rank, stk::mesh::Selector sel)
    {
        const stk::mesh::Part *part = get_exodus_part_of_rank(setId, rank);
        return stk::mesh::count_selected_entities(*part & sel, mBulkData.buckets(rank));
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

    void cache_mesh_constructs()
    {
        cache_element_blocks();
        cache_side_sets();
        cache_node_sets();
    }

    void cache_element_blocks()
    {
        const stk::mesh::PartVector &parts = mBulkData.mesh_meta_data().get_mesh_parts();
        for(size_t i = 0; i < parts.size(); i++)
            if(is_element_block(*parts[i]))
                mElementBlockParts.push_back(parts[i]);
    }

    void cache_side_sets()
    {
        const stk::mesh::PartVector &parts = mBulkData.mesh_meta_data().get_mesh_parts();
        for(size_t i = 0; i < parts.size(); i++)
        {
            if(is_side_set(*parts[i]))
            {
                mSideSetParts.push_back(parts[i]);
            }
        }
    }

    void cache_node_sets()
    {
        const stk::mesh::PartVector &parts = mBulkData.mesh_meta_data().get_mesh_parts();
        for(size_t i = 0; i < parts.size(); i++)
        {
            if(is_node_set(*parts[i]))
            {
                mNodeSetParts.push_back(parts[i]);
            }
        }
    }

private:
    const stk::mesh::BulkData& mBulkData;
    stk::mesh::PartVector mElementBlockParts;
    stk::mesh::PartVector mSideSetParts;
    stk::mesh::PartVector mNodeSetParts;
};

}
}

#endif /* EXODUSTRANSLATOR_HPP_ */
