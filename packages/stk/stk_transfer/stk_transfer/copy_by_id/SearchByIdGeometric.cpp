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

#include "SearchByIdGeometric.hpp"
#include <stk_util/parallel/Parallel.hpp>
#include <sstream>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>

namespace stk {
namespace transfer {

void SearchByIdGeometric::do_search(const TransferCopyByIdMeshAdapter & mesha,
                                    const TransferCopyByIdMeshAdapter & meshb,
                                    KeyToTargetProcessor & key_to_target_processor)
{
  const ParallelMachine comm = mesha.comm();
  const int p_rank = parallel_machine_rank(comm);

  key_to_target_processor.clear();
  m_remote_keys.clear();

  using MeshIDVector = TransferCopyByIdMeshAdapter::MeshIDVector;
  using Point = stk::search::Point<float>;
  using Sphere = stk::search::Sphere<float>;
  using BoundingBox = std::pair<stk::search::Sphere<float>,MeshIDProc>;

  std::vector<BoundingBox> source_bbox_vector;
  std::vector<BoundingBox> target_bbox_vector;

  {
    const MeshIDVector & source_ids = mesha.get_mesh_ids();
    for (size_t id_index = 0; id_index < source_ids.size(); ++id_index) {
      double coords[3] = {0.0,0.0,0.0};
      mesha.centroid(source_ids[id_index],coords);
      Point center;
      for (int i = 0; i < 3; ++i) {
        center[i] = static_cast<float>(coords[i]);
      }
      source_bbox_vector.emplace_back(Sphere(center,m_radius), MeshIDProc(source_ids[id_index],p_rank));
    }
  }
  {
    const MeshIDVector & target_ids = meshb.get_mesh_ids();
    for (size_t id_index = 0; id_index < target_ids.size(); ++id_index) {
      if (mesha.is_locally_owned(target_ids[id_index])) {
        key_to_target_processor.emplace_back(target_ids[id_index], p_rank);
      } else {
        m_remote_keys.insert(target_ids[id_index]);

        double coords[3] = {0.0,0.0,0.0};
        meshb.centroid(target_ids[id_index],coords);
        Point center;
        for (int i = 0; i < 3; ++i) { center[i] = coords[i]; }
        target_bbox_vector.emplace_back(Sphere(center,m_radius), MeshIDProc(target_ids[id_index],p_rank));
      }
    }
  }

  std::sort(source_bbox_vector.begin(), source_bbox_vector.end(), BoundingBoxCompare<BoundingBox>());
  std::sort(target_bbox_vector.begin(), target_bbox_vector.end(), BoundingBoxCompare<BoundingBox>());

  using EntityProcRelation = std::pair<MeshIDProc, MeshIDProc>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  EntityProcRelationVec source_to_target_vector;
  stk::search::coarse_search(source_bbox_vector,
                             target_bbox_vector,
                             stk::search::KDTREE,
                             mesha.comm(),
                             source_to_target_vector);

  // Match on Mesh_ID only and store in our special data structure.
  for (const auto & s2tEntry : source_to_target_vector) {
    Mesh_ID id = s2tEntry.first.id();
    if (s2tEntry.second.id() == id) {
      const bool idInMySourceMesh = (p_rank == s2tEntry.first.proc());
      if (idInMySourceMesh) {
        key_to_target_processor.emplace_back(id, s2tEntry.second.proc());
      }
    }
  }
}

}  } // namespace transfer stk
