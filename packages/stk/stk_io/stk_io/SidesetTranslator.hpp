#ifndef SIDESETTRANSLATOR_HPP_
#define SIDESETTRANSLATOR_HPP_

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "StkIoUtils.hpp"

namespace stk {
namespace io {

inline size_t get_number_sides_in_sideset(const stk::mesh::BulkData& bulk,
                                          int sideset_id,
                                          stk::mesh::Selector selector,
                                          stk::topology stk_element_topology,
                                          const stk::mesh::BucketVector& buckets)
{
    if (bulk.has_sideset_data())
    {
        size_t num_sides = 0;

        const stk::mesh::SideSet& sset = bulk.get_sideset_data(sideset_id);

        for(const stk::mesh::ElemIdSide& elem_and_side : sset)
        {
            stk::mesh::EntityId elem_id = elem_and_side.elem_id;
            int zero_based_side_ord = elem_and_side.side_ordinal;
            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk, elem_id, zero_based_side_ord, bulk.mesh_meta_data().side_rank());
            if(bulk.is_valid(side))
            {
                if(selector(bulk.bucket(side)))
                {
                    stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elem_id);
                    if(stk_element_topology == stk::topology::INVALID_TOPOLOGY ||
                       stk_element_topology == bulk.bucket(element).topology())
                    {
                        ++num_sides;
                    }
                }
            }
        }

        return num_sides;
    }
    else
    {
        return count_selected_entities(selector, buckets);
    }
}

template<typename INT>
void fill_element_and_side_ids(Ioss::GroupingEntity & io,
                               stk::mesh::Part * const part,
                               const stk::mesh::BulkData & bulk_data,
                               stk::topology stk_element_topology,
                               const stk::mesh::Selector *subset_selector,
                               stk::mesh::EntityVector &sides,
                               std::vector<INT>& elem_side_ids)
{
    if (bulk_data.has_sideset_data())
    {
        const stk::mesh::SideSet& sset = bulk_data.get_sideset_data(part->id());
        size_t num_sides = sset.size();
        elem_side_ids.reserve(num_sides*2);

        stk::mesh::Selector selector = *part & bulk_data.mesh_meta_data().locally_owned_part();
        if(subset_selector)
            selector &= *subset_selector;

        for(size_t i=0;i<sset.size();++i)
        {
            stk::mesh::EntityId elem_id = sset[i].elem_id;
            int zero_based_side_ord = sset[i].side_ordinal;
            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_id_side_pair_of_rank(bulk_data, elem_id, zero_based_side_ord, bulk_data.mesh_meta_data().side_rank());
            if(bulk_data.is_valid(side))
            {
                if(selector(bulk_data.bucket(side)))
                {
                    stk::mesh::Entity element = bulk_data.get_entity(stk::topology::ELEM_RANK, elem_id);
                    if(bulk_data.bucket(element).topology() == stk_element_topology)
                    {
                        elem_side_ids.push_back(elem_id);
                        elem_side_ids.push_back(zero_based_side_ord+1);
                        sides.push_back(side);
                    }
                }
            }
        }
    }
    else
    {
        const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(*part);

        stk::mesh::EntityRank type = part_primary_entity_rank(*part);
        size_t num_sides = get_entities(*part, type, bulk_data, sides, false, subset_selector);
        elem_side_ids.reserve(num_sides * 2);

        stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;

        for(size_t i = 0; i < num_sides; ++i)
        {
            std::vector<stk::mesh::Entity> side;
            side.push_back(sides[i]);
            std::vector<stk::mesh::Entity> side_elements;
            std::vector<stk::mesh::Entity> side_nodes(bulk_data.begin_nodes(sides[i]), bulk_data.end_nodes(sides[i]));

            get_entities_through_relations(bulk_data, side_nodes, elem_rank, side_elements);
            const size_t num_side_elem = side_elements.size();

            std::sort(side_elements.begin(), side_elements.end(), stk::mesh::EntityLess(bulk_data));

            stk::mesh::Entity suitable_elem = stk::mesh::Entity();
            stk::mesh::ConnectivityOrdinal suitable_ordinal = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;
            for(size_t j = 0; j < num_side_elem; ++j)
            {
                const stk::mesh::Entity elem = side_elements[j];
                const stk::mesh::Bucket &elemBucket = bulk_data.bucket(elem);
                const bool isSelectingEverything = subset_selector == NULL;
                const bool isElementBeingOutput = (isSelectingEverything || (*subset_selector)(elemBucket))
                                                  && elemBucket.member(meta_data.locally_owned_part());
                if(isElementBeingOutput)
                {
                    const stk::mesh::Entity * elem_sides = bulk_data.begin(elem, type);
                    stk::mesh::ConnectivityOrdinal const * side_ordinal = bulk_data.begin_ordinals(elem, type);
                    const size_t num_elem_sides = bulk_data.num_connectivity(elem, type);
                    for(size_t k = 0; k < num_elem_sides; ++k)
                    {
                        if(elem_sides[k] == side[0])
                        {
                            suitable_elem = elem;
                            suitable_ordinal = side_ordinal[k];
                            break;
                        }
                    }
                }
            }

            if(!bulk_data.is_valid(suitable_elem))
            {
                std::ostringstream oss;
                oss << "ERROR, no suitable element found";
                throw std::runtime_error(oss.str());
            }

            elem_side_ids.push_back(bulk_data.identifier(suitable_elem));
            elem_side_ids.push_back(suitable_ordinal + 1); // Ioss is 1-based, mesh is 0-based.
        }
    }
}

}
}

#endif /* SIDESETTRANSLATOR_HPP_ */
