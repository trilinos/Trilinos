// Copyright (c) 2015, Sandia Corporation.
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
//

#include "CopySearchGeometric.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include <sstream>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>

namespace stk {
namespace transfer {


void CopySearchGeometric::do_search(const CopyTransferMeshBase & mesha,
                                    const CopyTransferMeshBase & meshb,
                                    KeyToTargetProcessor & key_to_target_processor
                                    )
{
  const ParallelMachine comm = mesha.comm();
  const int p_rank = parallel_machine_rank(comm);

  key_to_target_processor.clear();
  m_remote_keys.clear();

  typedef CopyTransferMeshBase::Mesh_ID Mesh_ID;
  typedef CopyTransferMeshBase::MeshIDVector MeshIDVector;
  typedef stk::search::Point<float>  Point;
  typedef stk::search::Sphere<float> Sphere;
  typedef std::pair<stk::search::Sphere<float>,MeshIDProc> BoundingBox;

  std::vector<BoundingBox> source_bbox_vector;
  std::vector<BoundingBox> target_bbox_vector;

  {
    const MeshIDVector & source_ids = mesha.get_mesh_ids();
    for (size_t id_index=0 ; id_index<source_ids.size() ; ++id_index) {
      double coords[3] = {0.0,0.0,0.0};
      mesha.centroid(source_ids[id_index],coords);
      Point center;
      for (int i=0 ; i<3 ; ++i)
      {
        center[i] = static_cast<float>(coords[i]);
      }
      source_bbox_vector.push_back(std::make_pair( Sphere(center,m_radius), MeshIDProc(source_ids[id_index],p_rank)));
    }
  }
  {
    const MeshIDVector & target_ids = meshb.get_mesh_ids();
    for (size_t id_index=0 ; id_index<target_ids.size() ; ++id_index) {
      if (mesha.is_locally_owned(target_ids[id_index])) {
        key_to_target_processor[target_ids[id_index]] = p_rank;
      } else {
        m_remote_keys.insert(target_ids[id_index]);

        double coords[3] = {0.0,0.0,0.0};
        meshb.centroid(target_ids[id_index],coords);
        Point center;
        for (int i=0 ; i<3 ; ++i) { center[i] = coords[i]; }
        target_bbox_vector.push_back(std::make_pair( Sphere(center,m_radius), MeshIDProc(target_ids[id_index],p_rank)));
      }
    }
  }

  std::sort(source_bbox_vector.begin(),source_bbox_vector.end(),BoundingBoxCompare<BoundingBox>());
  std::sort(target_bbox_vector.begin(),target_bbox_vector.end(),BoundingBoxCompare<BoundingBox>());

  typedef std::pair<MeshIDProc, MeshIDProc> EntityProcRelation;
  typedef std::vector<EntityProcRelation>   EntityProcRelationVec;

  EntityProcRelationVec source_to_target_vector;
  stk::search::coarse_search(source_bbox_vector,
                             target_bbox_vector,
                             stk::search::BOOST_RTREE,
                             mesha.comm(),
                             source_to_target_vector
                             );

  // Match on Mesh_ID only and store in our special data structure.
  EntityProcRelationVec::const_iterator s2t_vec_iter = source_to_target_vector.begin();
  for ( ; s2t_vec_iter != source_to_target_vector.end() ; ++s2t_vec_iter ) {
    Mesh_ID id = s2t_vec_iter->first.id();
    if (s2t_vec_iter->second.id() == id) {
      const bool id_in_my_source_mesh = (p_rank == s2t_vec_iter->first.proc());
      if (id_in_my_source_mesh) {
        key_to_target_processor[id] = s2t_vec_iter->second.proc();
      }
    }
  }
}

}  } // namespace transfer stk
