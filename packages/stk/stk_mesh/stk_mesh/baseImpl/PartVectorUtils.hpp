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
// 

#ifndef stk_mesh_base_impl_PartVectorUtils_hpp
#define stk_mesh_base_impl_PartVectorUtils_hpp

//----------------------------------------------------------------------

#include <stk_util/util/SortAndUnique.hpp>
#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

template<typename PARTVECTOR>
void fill_add_parts_and_supersets(const PARTVECTOR & add_parts,
                                  OrdinalVector &addPartsAndSupersets)
{
  const unsigned expected_min_num_supersets = 3;
  const size_t expectedSizeOfAddPartList = add_parts.size() * expected_min_num_supersets;
  addPartsAndSupersets.reserve(expectedSizeOfAddPartList);
  for(const Part* addPart : add_parts)
  {
    addPartsAndSupersets.push_back(addPart->mesh_meta_data_ordinal());
    const PartVector& supersets = addPart->supersets();
    for(const Part* superset : supersets)
    {
      addPartsAndSupersets.push_back(superset->mesh_meta_data_ordinal());
    }
  }
  stk::util::sort_and_unique(addPartsAndSupersets);
}

template<typename PARTVECTOR>
void fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(
                                  const PARTVECTOR & remove_parts,
                                  const OrdinalVector & addPartsAndSupersets,
                                  const stk::mesh::Bucket* entityBucket,
                                  OrdinalVector &removePartsAndSubsetsMinusPartsInAddPartsList)
{
  const unsigned expected_min_num_subsets = 3;
  const size_t expectedSizeOfRemovePartList = remove_parts.size() * expected_min_num_subsets;
  removePartsAndSubsetsMinusPartsInAddPartsList.reserve(expectedSizeOfRemovePartList);
  for(const Part* rmPart : remove_parts)
  {
    if(!contains_ordinal(addPartsAndSupersets.begin(),
                         addPartsAndSupersets.end(),
                         rmPart->mesh_meta_data_ordinal()))
    {
      removePartsAndSubsetsMinusPartsInAddPartsList.push_back(rmPart->mesh_meta_data_ordinal());
      if (entityBucket != nullptr) {
        const PartVector& subsets = rmPart->subsets();
        for(const Part* subset : subsets) {
          if(entityBucket->member(*subset)) {
            removePartsAndSubsetsMinusPartsInAddPartsList.push_back(subset->mesh_meta_data_ordinal());
          }
        }
      }
    }
  }
  stk::util::sort_and_unique(removePartsAndSubsetsMinusPartsInAddPartsList);
}

template<typename PARTVECTOR>
void fill_remove_parts_and_subsets_minus_parts_in_add_parts_list(
                                  const PARTVECTOR & remove_parts,
                                  const OrdinalVector & addPartsAndSupersets,
                                  const stk::mesh::BucketVector& buckets,
                                  OrdinalVector &removePartsAndSubsetsMinusPartsInAddPartsList)
{
  const unsigned expected_min_num_subsets = 3;
  const size_t expectedSizeOfRemovePartList = remove_parts.size() * expected_min_num_subsets;
  removePartsAndSubsetsMinusPartsInAddPartsList.reserve(expectedSizeOfRemovePartList);
  for(const Part* rmPart : remove_parts) {
    if(!contains_ordinal(addPartsAndSupersets.begin(),
                         addPartsAndSupersets.end(),
                         rmPart->mesh_meta_data_ordinal()))
    {
      removePartsAndSubsetsMinusPartsInAddPartsList.push_back(rmPart->mesh_meta_data_ordinal());
      const PartVector& subsets = rmPart->subsets();
      for(const Part* subset : subsets) {
        for(auto entityBucket : buckets) {
          if(entityBucket->member(*subset)) {
            removePartsAndSubsetsMinusPartsInAddPartsList.push_back(subset->mesh_meta_data_ordinal());
          }
        }
      }
    }
  }
  stk::util::sort_and_unique(removePartsAndSubsetsMinusPartsInAddPartsList);
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_base_impl_PartVectorUtils_hpp

