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

#ifndef SIDESETTRANSLATOR_HPP_
#define SIDESETTRANSLATOR_HPP_

#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/EquivalentEntityBlocks.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/SideSetUtil.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/EntityLess.hpp"
#include "stk_mesh/base/Relation.hpp"
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/OutputParams.hpp"

namespace stk {
namespace io {

template<typename INT>
void fill_element_and_side_ids_from_sideset(const stk::mesh::SideSet& sset,
                                            stk::io::OutputParams &params,
                                            const stk::mesh::Part * part,
                                            const stk::mesh::Part *parentElementBlock,
                                            stk::topology stk_element_topology,
                                            stk::mesh::EntityVector &sides,
                                            std::vector<INT>& elem_side_ids,
                                            INT sideOrdOffset = 0)
{
  const mesh::BulkData &bulk_data = params.bulk_data();
  const stk::mesh::Selector *elemSubsetSelector = params.get_subset_selector();
  const stk::mesh::Selector *elemOutputSelector = params.get_output_selector(stk::topology::ELEM_RANK);
  const stk::mesh::Selector *faceOutputSelector = params.get_output_selector(bulk_data.mesh_meta_data().side_rank());

  size_t num_sides = sset.size();
  elem_side_ids.reserve(num_sides*2);

  stk::mesh::Selector selector = *part & construct_sideset_selector(params);
  if (nullptr != faceOutputSelector) {
    selector &= *faceOutputSelector;
  }
  stk::mesh::Selector parentElementSelector =  (parentElementBlock == nullptr) ? stk::mesh::Selector() : *parentElementBlock;
  const stk::mesh::EntityRank sideRank = part->primary_entity_rank();

  unsigned previousBucketId = stk::mesh::INVALID_BUCKET_ID;

  for(size_t i=0;i<sset.size();++i)
  {
    stk::mesh::Entity element = sset[i].element;
    stk::mesh::EntityId elemId = bulk_data.identifier(element);
    int zero_based_side_ord = sset[i].side - sideOrdOffset;
    stk::mesh::Entity side = stk::mesh::get_side_entity_for_elem_side_pair_of_rank(bulk_data, element, zero_based_side_ord, sideRank);
    if(bulk_data.is_valid(side))
    {
      stk::mesh::Bucket &sideBucket = bulk_data.bucket(side);
      bool sideIsSelected = false;

      if (sideBucket.bucket_id() == previousBucketId) {
        sideIsSelected = true;
      } else {
        sideIsSelected = selector(sideBucket);

        if (sideIsSelected) {
          previousBucketId = sideBucket.bucket_id();
        }
      }

      if (sideIsSelected) {
        stk::mesh::Bucket &elementBucket = bulk_data.bucket(element);
        if (elementBucket.owned() && (stk_element_topology == stk::topology::INVALID_TOPOLOGY ||
                                         stk_element_topology == elementBucket.topology())) {

          bool selectedByParent = (parentElementBlock == nullptr) ? true : parentElementSelector(elementBucket);
          bool selectedByBucket = (elemSubsetSelector == nullptr) ? true : (*elemSubsetSelector)(elementBucket);
          bool selectedByOutput = (elemOutputSelector == nullptr) ? true : (*elemOutputSelector)(elementBucket);

          if (selectedByBucket && selectedByParent && selectedByOutput) {
            elem_side_ids.push_back(elemId);
            elem_side_ids.push_back(zero_based_side_ord + 1);
            sides.push_back(side);
          }
        }
      }
    }
  }
}

inline void fill_side_elements_and_nodes(const stk::mesh::BulkData &bulk_data,
                                         stk::mesh::Entity side,
                                         std::vector<stk::mesh::Entity> &side_elements,
                                         std::vector<stk::mesh::Entity> &side_nodes)
{
    side_nodes.assign(bulk_data.begin_nodes(side), bulk_data.end_nodes(side));

    stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
    get_entities_through_relations(bulk_data, side_nodes, elem_rank, side_elements);

    std::sort(side_elements.begin(), side_elements.end(), stk::mesh::EntityLess(bulk_data));
}

inline bool is_correct_sub_topology(const stk::mesh::BulkData &mesh,
                                    stk::mesh::Entity element,
                                    unsigned sideOrdinal,
                                    stk::mesh::Entity side)
{
    stk::topology elemTopology = mesh.bucket(element).topology();
    stk::topology sideTopology = mesh.bucket(side).topology();

    stk::topology subTopology = elemTopology.sub_topology(mesh.mesh_meta_data().side_rank(), sideOrdinal);
    return (sideTopology == subTopology);
}

template<typename INT>
void fill_element_and_side_ids_from_connectivity(stk::io::OutputParams &params,
                                                 const stk::mesh::Part * part,
                                                 const stk::mesh::Part *parentElementBlock,
                                                 stk::topology stk_element_topology,
                                                 stk::mesh::EntityVector &sides,
                                                 std::vector<INT>& elem_side_ids,
                                                 INT sideOrdOffset = 0)
{
    const stk::mesh::BulkData &bulk_data = params.bulk_data();
    const stk::mesh::Selector *subset_selector = params.get_subset_selector();
    const stk::mesh::Selector *output_selector = params.get_output_selector(stk::topology::ELEM_RANK);

    const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(*part);

    stk::mesh::EntityRank type = stk::io::part_primary_entity_rank(*part);
    stk::mesh::EntityVector allSides;
    size_t num_sides = stk::io::get_entities(params, *part, type, allSides, false);
    elem_side_ids.reserve(num_sides * 2);

    std::vector<stk::mesh::Entity> side_elements;
    std::vector<stk::mesh::Entity> side_nodes;
    for(size_t i = 0; i < num_sides; ++i)
    {
        stk::mesh::Entity side = allSides[i];

        fill_side_elements_and_nodes(bulk_data, side, side_elements, side_nodes);

        stk::mesh::Entity suitable_elem = stk::mesh::Entity();
        stk::mesh::ConnectivityOrdinal suitable_ordinal = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

        stk::mesh::Entity suitable_elem_with_correct_polarity = stk::mesh::Entity();
        stk::mesh::ConnectivityOrdinal suitable_ordinal_with_correct_polarity = stk::mesh::INVALID_CONNECTIVITY_ORDINAL;

        bool foundElementWithCorrectPolarity = false;

        for(const stk::mesh::Entity& elem : side_elements)
        {
            const stk::mesh::Bucket &elemBucket = bulk_data.bucket(elem);
            bool selectedByBucket = (   subset_selector == nullptr) ? true :    (*subset_selector)(elemBucket);
            bool selectedByOutput = (   output_selector == nullptr) ? true :    (*output_selector)(elemBucket);
            const bool isElementBeingOutput = (selectedByBucket && selectedByOutput) && elemBucket.member(meta_data.locally_owned_part());

            if(isElementBeingOutput)
            {
                const stk::mesh::Entity * elem_sides = bulk_data.begin(elem, type);
                stk::mesh::ConnectivityOrdinal const * side_ordinal = bulk_data.begin_ordinals(elem, type);
                const size_t num_elem_sides = bulk_data.num_connectivity(elem, type);

                for(size_t k = 0; k < num_elem_sides; ++k)
                {
                    if(elem_sides[k] == side)
                    {
                        suitable_elem = elem;
                        suitable_ordinal = side_ordinal[k];
                        break;
                    }
                }

                if (bulk_data.is_valid(suitable_elem))
                {
                    for(size_t k = 0; k < num_elem_sides; ++k)
                    {
                        if(elem_sides[k] == side)
                        {
                            bool hasCorrectPolarity = false;

                            if(is_correct_sub_topology(bulk_data, elem, side_ordinal[k], side)) {
                                stk::mesh::EquivAndPositive result = stk::mesh::is_side_equivalent_and_positive(bulk_data, elem, side_ordinal[k], side_nodes);
                                hasCorrectPolarity = result.is_equiv && result.is_positive;
                            }

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
        }

        if(foundElementWithCorrectPolarity) {
            suitable_elem = suitable_elem_with_correct_polarity;
            suitable_ordinal = suitable_ordinal_with_correct_polarity;
        }

        if (bulk_data.is_valid(suitable_elem))
        {
            int oneBasedOrdinal = suitable_ordinal - sideOrdOffset + 1;
            elem_side_ids.push_back(bulk_data.identifier(suitable_elem));
            elem_side_ids.push_back(oneBasedOrdinal);
            sides.push_back(side);
        }
    }
}

template<typename INT>
void fill_element_and_side_ids(stk::io::OutputParams &params,
                               const stk::mesh::Part * part,
                               const stk::mesh::Part *parentElementBlock,
                               stk::topology stk_element_topology,
                               stk::mesh::EntityVector &sides,
                               std::vector<INT>& elem_side_ids,
                               INT sideOrdOffset = 0)
{
    const mesh::BulkData &bulk_data = params.bulk_data();
    const stk::mesh::Part &parentPart = stk::mesh::get_sideset_parent(*part);

    if (bulk_data.does_sideset_exist(parentPart))
    {
        const stk::mesh::SideSet& sset = bulk_data.get_sideset(parentPart);
        fill_element_and_side_ids_from_sideset(sset, params, part, parentElementBlock,
                                               stk_element_topology, sides, elem_side_ids, sideOrdOffset);
    }
    else
    {
        fill_element_and_side_ids_from_connectivity(params, part, parentElementBlock,
                                                    stk_element_topology, sides, elem_side_ids, sideOrdOffset);
    }
}


inline size_t get_number_sides_in_sideset(OutputParams &params,
                                          const stk::mesh::Part &ssPart,
                                          stk::topology stk_element_topology,
                                          const stk::mesh::Part *parentElementBlock = nullptr,
                                          int sideOrdOffset = 0)
{
    stk::mesh::EntityVector sides;
    std::vector<int> elemSideIds;

    fill_element_and_side_ids( params,
                              &ssPart,
                               parentElementBlock,
                               stk_element_topology,
                               sides,
                               elemSideIds,
                               sideOrdOffset);

    return sides.size();
}

}
}

#endif /* SIDESETTRANSLATOR_HPP_ */
