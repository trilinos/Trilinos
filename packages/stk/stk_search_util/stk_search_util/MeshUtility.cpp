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
#include "stk_search_util/MeshUtility.hpp"
#include "stk_mesh/base/Part.hpp"                     // for Part
#include "stk_mesh/base/Bucket.hpp"                   // for Bucket
#include "stk_mesh/base/BulkData.hpp"                 // for BulkData
#include "stk_mesh/base/CompositeRank.hpp"            // for CompositeRank
#include "stk_mesh/base/MetaData.hpp"                 // for MetaData
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/FieldBase.hpp"                // for field_data, Fie...
#include "stk_util/util/ReportHandler.hpp"            // for eval_test_condi...
#include "stk_util/parallel/ParallelReduce.hpp"       // for all_reduce_max
#include "stk_search/DistanceComparison.hpp"          // for distance_sq
#include <algorithm>                                  // for all_of, max, min
#include <array>                                      // for array
#include <cmath>                                      // for sqrt
#include <cstddef>                                    // for size_t
#include <limits>                                     // for numeric_limits
#include <memory>                                     // for __shared_ptr_ac...

#include <ctype.h>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

stk::search::Box<double> get_mesh_bounding_box(const stk::mesh::FieldBase* coords, const stk::mesh::BulkData& bulk)
{
  const auto& meta = bulk.mesh_meta_data();
  const int ndim = meta.spatial_dimension();
  auto sel = stk::mesh::selectField(*coords) & meta.locally_owned_part();
  const auto& buckets = bulk.get_buckets(coords->entity_rank(), sel);
  const double m = std::numeric_limits<double>::max();
  std::array<double, 3> min_corner = { { m, m, m } };
  std::array<double, 3> max_corner = { { -m, -m, -m } };

  for(auto&& b_ptr : buckets) {
    stk::mesh::Bucket& b = *b_ptr;
    const size_t length = b.size();
    for(size_t k = 0; k < length; ++k) {
      const double* c = static_cast<const double *>(stk::mesh::field_data(*coords, b[k]));
      for(int d = 0; d < ndim; ++d) {
        min_corner[d] = std::min(min_corner[d], c[d]);
        max_corner[d] = std::max(max_corner[d], c[d]);
      }
    }
  }

  std::array<double, 3> g_min_corner = { { 0, 0, 0 } };
  std::array<double, 3> g_max_corner = { { 0, 0, 0 } };
  stk::all_reduce_min(bulk.parallel(), min_corner.data(), g_min_corner.data(), ndim);
  stk::all_reduce_max(bulk.parallel(), max_corner.data(), g_max_corner.data(), ndim);

  return stk::search::Box<double>(g_min_corner[0], g_min_corner[1], g_min_corner[2],
                                  g_max_corner[0], g_max_corner[1], g_max_corner[2]);
}

stk::mesh::EntityRank get_objects_rank(const stk::mesh::EntityRank objectType, const stk::mesh::PartVector& parts)
{
  stk::mesh::EntityRank rank = objectType;

  if(!parts.empty()) {
    rank = stk::mesh::CompositeRank::get_rank(parts[0]);
    if(!std::all_of(parts.begin(), parts.end(),
                    [&rank](stk::mesh::Part* i) { return stk::mesh::CompositeRank::get_rank(i) == rank; })) {
      STK_ThrowRequireMsg(false, "All mesh parts must be the same entity rank");
    }
  }

  if(rank == stk::topology::INVALID_RANK) {
    if(!parts.empty() && (parts[0]->subsets().size() > 0)) {
      rank = stk::mesh::CompositeRank::get_rank(parts[0]->subsets()[0]);

      for(const stk::mesh::Part* part : parts) {
        const stk::mesh::PartVector& subset = part->subsets();

        if(!std::all_of(subset.begin(), subset.end(),
                        [&rank](stk::mesh::Part* i) { return stk::mesh::CompositeRank::get_rank(i) == rank; })) {
          STK_ThrowRequireMsg(false, "All mesh parts must be the same entity rank");
        }
      }
    }
    else {
      rank = objectType;
    }
  }

  return rank;
}

stk::mesh::Selector
get_objects_selector(const stk::mesh::MetaData& meta, const stk::mesh::PartVector& parts,
                     const stk::mesh::Selector* activeSelector, bool includeGhosts)
{
  stk::mesh::Selector ghosts = (!meta.locally_owned_part()) & (!meta.globally_shared_part());

  stk::mesh::Selector selector = meta.locally_owned_part();

  if (includeGhosts) {
    selector |= ghosts;
  }

  if(nullptr != activeSelector) {
    selector &= *activeSelector;
  }

  if(!parts.empty()) {
    selector &= stk::mesh::selectUnion(parts);
  }

  return selector;
}

unsigned get_number_of_parametric_coordinates(stk::topology topo)
{
  unsigned n = topo.dimension();
  if(topo.is_shell()) {
    --n;
  }

  return n;
}

bool is_valid_entity_rank(const stk::mesh::MetaData& meta, const stk::mesh::EntityRank rank)
{
  return (rank == stk::topology::ELEM_RANK) || (rank == stk::topology::FACE_RANK) ||
         (rank == stk::topology::EDGE_RANK && meta.spatial_dimension() == 2);
}

bool part_has_proper_entity_rank(const stk::mesh::Part* part)
{
  const stk::mesh::EntityRank rank = part->primary_entity_rank();
  const stk::mesh::MetaData& meta = part->mesh_meta_data();

  if(is_valid_entity_rank(meta, rank)) {
    return true;
  }
  else if(rank == stk::topology::INVALID_RANK) {
    bool hasAtLeastOneValidSubsetPart = false;
    for(const stk::mesh::Part* subsetPart : part->subsets()) {
      if(subsetPart->primary_entity_rank() != stk::topology::INVALID_RANK &&
         !is_valid_entity_rank(meta, subsetPart->primary_entity_rank())) {
        return false;
      }

      hasAtLeastOneValidSubsetPart = true;
    }
    return hasAtLeastOneValidSubsetPart;
  }

  return false;
}

std::vector<std::string> get_part_membership(const stk::mesh::BulkData& bulk, const stk::mesh::Entity e, const stk::mesh::PartVector& parts)
{
  std::vector<std::string> partNames;

  const stk::mesh::Bucket &bucket = bulk.bucket(e);

  for(const stk::mesh::Part* part : parts) {
    if(bucket.member(*part)) {
      std::ostringstream outputName;

      std::string partName(part->name());
      std::transform(partName.begin(), partName.end(), partName.begin(), ::toupper);
      outputName << partName;

      const stk::mesh::PartVector subsetParts = part->subsets();
      if(subsetParts.size() > 0) {
        outputName << "{ ";
        for(const stk::mesh::Part* subsetPart : subsetParts) {
          std::string subsetPartName(subsetPart->name());
          std::transform(subsetPartName.begin(), subsetPartName.end(), subsetPartName.begin(), ::toupper);

          outputName << subsetPartName;
          if(bucket.member(*subsetPart)) {
            outputName << "*";
          }
          outputName << " ";
        }
        outputName << "}";
      }

      partNames.push_back(outputName.str());
    }
  }

  return partNames;
}

std::vector<std::string> get_part_membership(const stk::mesh::BulkData& bulk, const stk::mesh::EntityKey k, const stk::mesh::PartVector& parts)
{
  stk::mesh::Entity e = bulk.get_entity(k);
  return get_part_membership(bulk, e, parts);
}

std::string get_time_stamp()
{
  const static std::string defaultStamp("------------------------------------");
  return defaultStamp;
}

} // namespace search
} // namespace stk

