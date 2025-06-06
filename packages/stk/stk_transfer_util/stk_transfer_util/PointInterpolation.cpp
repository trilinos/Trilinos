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


// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_transfer_util/PointInterpolation.hpp"
#include "stk_util/util/string_case_compare.hpp"
#include <vector>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {

namespace impl {

const stk::mesh::Selector& universal_selector() {
  static stk::mesh::Selector universalSelector = stk::mesh::Selector().complement();
  return universalSelector;
}

}

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::search::SearchField& meField, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& fieldDataScratch, std::vector<double>& result)
{
  const stk::mesh::BulkData& bulkData = meField.get_mesh();
  const stk::mesh::Bucket& bucket = bulkData.bucket(entity);
  const stk::search::SearchTopology meTopo(bucket.topology(), stk::search::spmd::make_entity_key_pair(bulkData,entity), &bucket);

  unsigned numFieldComponents(0);
  unsigned numNodes(0);

  masterElement.nodal_field_data(meTopo.get_key(), meField, numFieldComponents, numNodes, fieldDataScratch);
  masterElement.evaluate_field(meTopo, paramCoords, numFieldComponents, fieldDataScratch, result);

  apply_bounds(result, lowerBound, upperBound);
}

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::mesh::FieldBase& field, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& fieldDataScratch, std::vector<double>& result)
{
  const stk::search::SearchField meField(&field);
  nodal_interpolation(entity, meField, paramCoords, masterElement, lowerBound, upperBound, fieldDataScratch, result);
}

void nodal_interpolation(stk::mesh::Entity entity,
                         const stk::mesh::FieldBase& field, const std::vector<double>& paramCoords,
                         const stk::search::MasterElementProviderInterface& masterElement, const double lowerBound,
                         const double upperBound, std::vector<double>& result)
{
  std::vector<double> fieldDataScratch;
  nodal_interpolation(entity, field, paramCoords, masterElement, lowerBound, upperBound, fieldDataScratch, result);
}

} // namespace transfer
} // namespace stk

