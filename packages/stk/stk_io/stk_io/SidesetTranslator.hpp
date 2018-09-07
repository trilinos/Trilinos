// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef SIDESETTRANSLATOR_HPP_
#define SIDESETTRANSLATOR_HPP_

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/OutputParams.hpp"

namespace stk {
namespace io {

inline size_t get_number_sides_in_sideset(OutputParams &params,
                                          const stk::mesh::Part &ssPart,
                                          stk::mesh::Selector selector,
                                          stk::topology stk_element_topology,
                                          const stk::mesh::BucketVector& buckets,
                                          const stk::mesh::Part *parentElementBlock = nullptr)
{
    const mesh::BulkData &bulk = params.bulk_data();
    const stk::mesh::Selector *subset_selector = params.get_subset_selector();
    const stk::mesh::Selector *output_selector = params.get_output_selector();

    const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(ssPart);

    if (bulk.does_sideset_exist(parentPart))
    {
        selector &= ( bulk.mesh_meta_data().locally_owned_part() | bulk.mesh_meta_data().globally_shared_part());
        stk::mesh::Selector parentElementSelector =  (parentElementBlock == nullptr) ? stk::mesh::Selector() : *parentElementBlock;

        size_t num_sides = 0;

        const stk::mesh::SideSet& sset = bulk.get_sideset(parentPart);

        for(const stk::mesh::SideSetEntry& elem_and_side : sset)
        {
            stk::mesh::Entity element = elem_and_side.element;
            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulk, element, elem_and_side.side);
            if(bulk.is_valid(side))
            {
                if(selector(bulk.bucket(side)))
                {
                    stk::mesh::Bucket &elementBucket = bulk.bucket(element);

                    if(stk_element_topology == stk::topology::INVALID_TOPOLOGY ||
                       stk_element_topology == elementBucket.topology())
                    {
                        bool selectedByParent = (parentElementBlock == nullptr) ? true : parentElementSelector(elementBucket);
                        bool selectedByBucket = (   subset_selector == nullptr) ? true :    (*subset_selector)(elementBucket);
                        bool selectedByOutput = (   output_selector == nullptr) ? true :    (*output_selector)(elementBucket);

                        if(selectedByBucket && selectedByParent && selectedByOutput)
                        {
                            ++num_sides;
                        }
                    }
                }
            }
        }

        return num_sides;
    }
    else
    {
        selector &= bulk.mesh_meta_data().locally_owned_part();
        return count_selected_entities(selector, buckets);
    }
}

template<typename INT>
void fill_element_and_side_ids(stk::io::OutputParams &params,
                               Ioss::GroupingEntity & io,
                               stk::mesh::Part * const part,
                               const stk::mesh::Part *parentElementBlock,
                               stk::topology stk_element_topology,
                               stk::mesh::EntityVector &sides,
                               std::vector<INT>& elem_side_ids)
{
    const mesh::BulkData &bulk_data = params.bulk_data();
    const stk::mesh::Selector *subset_selector = params.get_subset_selector();
    const stk::mesh::Selector *output_selector = params.get_output_selector();

    const stk::mesh::Part &parentPart = stk::io::get_sideset_parent(*part);

    if (bulk_data.does_sideset_exist(parentPart))
    {
        const stk::mesh::SideSet& sset = bulk_data.get_sideset(parentPart);
        size_t num_sides = sset.size();
        elem_side_ids.reserve(num_sides*2);

        stk::mesh::Selector selector = *part & ( bulk_data.mesh_meta_data().locally_owned_part() | bulk_data.mesh_meta_data().globally_shared_part() );
        if(subset_selector)
            selector &= *subset_selector;
        if(output_selector)
            selector &= *output_selector;

        stk::mesh::Selector parentElementSelector =  (parentElementBlock == nullptr) ? stk::mesh::Selector() : *parentElementBlock;

        for(size_t i=0;i<sset.size();++i)
        {
            stk::mesh::Entity element = sset[i].element;
            stk::mesh::EntityId elemId = bulk_data.identifier(element);
            int zero_based_side_ord = sset[i].side;
            stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair(bulk_data, element, zero_based_side_ord);
            if(bulk_data.is_valid(side))
            {
                if(selector(bulk_data.bucket(side)))
                {
                    stk::mesh::Bucket &elementBucket = bulk_data.bucket(element);
                    if(stk_element_topology == stk::topology::INVALID_TOPOLOGY ||
                       stk_element_topology == elementBucket.topology())
                    {
                        bool selectedByParent = (parentElementBlock == nullptr) ? true : parentElementSelector(elementBucket);
                        bool selectedByBucket = (   subset_selector == nullptr) ? true :    (*subset_selector)(elementBucket);
                        bool selectedByOutput = (   output_selector == nullptr) ? true :    (*output_selector)(elementBucket);

                        if(selectedByBucket && selectedByParent && selectedByOutput)
                        {
                            elem_side_ids.push_back(elemId);
                            elem_side_ids.push_back(zero_based_side_ord+1);
                            sides.push_back(side);
                        }
                    }
                }
            }
        }
    }
    else
    {
        const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(*part);

        stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(*part);
        size_t num_sides = stk::io::get_entities(params, *part, type, sides, false);
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

            stk::mesh::Entity suitable_elem_with_correct_polarity = stk::mesh::Entity();
            stk::mesh::ConnectivityOrdinal suitable_ordinal_with_correct_polarity = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

            bool foundElementWithCorrectPolarity = false;

            for(size_t j = 0; j < num_side_elem; ++j)
            {
                const stk::mesh::Entity elem = side_elements[j];
                const stk::mesh::Bucket &elemBucket = bulk_data.bucket(elem);
                const bool isSelectingEverything = subset_selector == nullptr;
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

                    for(size_t k = 0; k < num_elem_sides; ++k)
                    {
                        if(elem_sides[k] == side[0])
                        {
                            stk::mesh::EquivAndPositive result = stk::mesh::is_side_equivalent_and_positive(bulk_data, elem, side_ordinal[k], side_nodes);

                            bool hasCorrectPolarity = result.is_equiv && result.is_positive;
                            if (hasCorrectPolarity) {
                                foundElementWithCorrectPolarity = true;
                                suitable_elem_with_correct_polarity = elem;
                                suitable_ordinal_with_correct_polarity = side_ordinal[k];
                                break;
                            }
                        }
                    }
                }
            }

            if(foundElementWithCorrectPolarity) {
                suitable_elem = suitable_elem_with_correct_polarity;
                suitable_ordinal = suitable_ordinal_with_correct_polarity;
            }

            ThrowRequireMsg( bulk_data.is_valid(suitable_elem), __FILE__ << ", " << __FUNCTION__ << ", ERROR, no suitable element found");

            elem_side_ids.push_back(bulk_data.identifier(suitable_elem));
            elem_side_ids.push_back(suitable_ordinal + 1); // Ioss is 1-based, mesh is 0-based.
        }
    }
}

}
}

#endif /* SIDESETTRANSLATOR_HPP_ */
