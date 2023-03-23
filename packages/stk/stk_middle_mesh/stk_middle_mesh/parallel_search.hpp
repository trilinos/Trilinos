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

#ifndef  STK_MIDDLE_MESH_PARALLEL_SEARCH_HPP
#define  STK_MIDDLE_MESH_PARALLEL_SEARCH_HPP

#include "mesh.hpp"
#include "create_mesh.hpp"
#include "mesh_io.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"
#include "bounding_box_search.hpp"
#include "stk_search/Box.hpp"
#include "stk_transfer/GeometricTransfer.hpp"
#include "stk_util/util/SortAndUnique.hpp"

#include <limits>
#include <map>

namespace stk {
namespace middle_mesh {
namespace impl {

class MeshScatterSpec {
public:
  MeshScatterSpec() = default;

  void add_destination(mesh::MeshEntityPtr entity, int destRank)
  {
    if(nullptr != entity) {
      MapKey key(entity->get_id(), entity->get_type());

      auto iter = m_entityProcMap.find(key);
      if(iter != m_entityProcMap.end()) {
        stk::util::insert_keep_sorted_and_unique(destRank, iter->second);
      } else {
        m_entityProcMap[key].push_back(destRank);
      }
    }
  }

  void get_destinations(mesh::MeshEntityPtr entity, std::vector<int>& destRanks)
  {
    if(nullptr != entity) {
      MapKey key(entity->get_id(), entity->get_type());

      auto iter = m_entityProcMap.find(key);
      if(iter != m_entityProcMap.end()) {
        destRanks.reserve(iter->second.size());
        for(int dest : iter->second) {
          destRanks.push_back(dest);
        }
      }
    }
  }

private:
  using MapKey = std::pair<int, mesh::MeshEntityType>;
  using MapValue = std::vector<int>;

  std::map<MapKey, MapValue> m_entityProcMap;
};

class SearchMesh {
 public:
  using Entity        = mesh::MeshEntity;
  using EntityKey     = int64_t;

  using EntityProc    = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point         = stk::search::Point<double>;
  using Box           = stk::search::Box<double>;
  using BoundingBox   = std::pair<Box, EntityProc>;

  SearchMesh(std::shared_ptr<mesh::Mesh> inputMesh)
    : m_mesh(inputMesh)
    {}

  MPI_Comm get_comm() const {return m_mesh->get_comm(); }

  void fill_bounding_box(mesh::MeshEntityPtr element, stk::search::Point<double>& minCorner,
                         stk::search::Point<double>& maxCorner) const
  {
    for(unsigned j = 0; j < 3u; ++j) {
      minCorner[j] = std::numeric_limits<double>::max();
      maxCorner[j] = -std::numeric_limits<double>::max();
    }

    int numEdges = element->count_down();
    for(auto i = 0; i < numEdges; ++i) {
      auto edge = element->get_down(i);

      int numNodes = edge->count_down();
      for(auto j = 0; j < numNodes; ++j) {
        auto node = edge->get_down(j);
        auto nodeCoords = node->get_point_orig(0);

        for(unsigned k = 0; k < 3u; ++k) {
          minCorner[k] = std::min(minCorner[k], nodeCoords[k]);
          maxCorner[k] = std::max(maxCorner[k], nodeCoords[k]);
        }
      }
    }
  }

  void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const
  {
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    auto entities = m_mesh->get_elements();

    stk::search::Point<double> minCorner, maxCorner;

    for(auto entity : entities) {
      fill_bounding_box(entity, minCorner, maxCorner);
      EntityProc entityProc(entity->get_id(), proc);

      BoundingBox boundingBox(Box(minCorner, maxCorner), entityProc);
      boundingBoxes.push_back(boundingBox);
    }

    std::sort(boundingBoxes.begin(), boundingBoxes.end(),
              [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
  };

  std::shared_ptr<mesh::Mesh> get_mesh() const { return m_mesh; }

 private:
  std::shared_ptr<mesh::Mesh> m_mesh;
};

using BoundingBoxSearch =
    stk::middle_mesh::search::BoundingBoxSearch<
        stk::middle_mesh::search::BoundingBoxSearchType<SearchMesh, SearchMesh>>;

using SearchRelationVec = BoundingBoxSearch::EntityProcRelationVec;
using UnpairedRelationVec = std::vector<SearchMesh::BoundingBox>;

enum class SplitCommColor {
      RECV = 0,
      SEND,
      INVALID
};

}
}
}

#endif
