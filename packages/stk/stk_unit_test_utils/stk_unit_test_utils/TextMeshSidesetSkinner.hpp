// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshSidesetSkinner_hpp
#define TextMeshSidesetSkinner_hpp

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <ctype.h>                                   // for toupper
#include <stddef.h>                                  // for size_t
#include <algorithm>                                 // for remove, etc
#include <iterator>                                  // for insert_iterator
#include <map>
#include <set>                                       // for set
#include <sstream>                                   // for operator<<, etc
#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector
#include <unordered_map>
#include <sstream>                       // for ostringstream
#include <iostream>
#include <functional>
#include <stdexcept>
#include <numeric>
#include <strings.h>
#include <cassert>      // for assert

#include "TextMeshFuncs.hpp"
#include "TextMeshDataTypes.hpp"
#include "TextMeshEntityGroup.hpp"
#include "TextMeshAdjacencyGraph.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

template <typename EntityId, typename Topology>
class SidesetSkinner : public SideAdjacencyGraph<EntityId, Topology>
{
public:
  using BaseClass = SideAdjacencyGraph<EntityId, Topology>;

  SidesetSkinner() = default;

  ~SidesetSkinner() {}

  size_t get_num_elements() const override
  {
    assert(m_textMeshData != nullptr);
    return m_textMeshData->elementDataVec.size();
  }

  int get_element_proc(const size_t elemIndex) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];
    return elemData.proc;
  }

  bool element_has_any_node_on_proc(const size_t elemIndex, int proc) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];

    for (const EntityId &nodeId : elemData.nodeIds) {
      const std::set<int> &procsForNode = m_textMeshData->procs_for_node(nodeId);
      if (procsForNode.count(proc) > 0) {
        return true;
      }
    }

    return false;
  }

  const std::string& get_element_block_name(const size_t elemIndex) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];
    return elemData.partName;
  }

  const std::vector<EntityId>& get_element_node_ids(const size_t elemIndex) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];
    return elemData.nodeIds;
  }

  const Topology& get_element_topology(const size_t elemIndex) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];
    return elemData.topology;
  }

  EntityId get_element_id(const size_t elemIndex) const override
  {
    assert(m_textMeshData != nullptr);
    const ElementData<EntityId, Topology> &elemData = m_textMeshData->elementDataVec[elemIndex];
    return elemData.identifier;
  }

  void set_skin_blocks(const std::vector<std::string> &skinBlocks) { m_skinBlocks = skinBlocks; }

  void skin_blocks(const TextMeshData<EntityId, Topology> &textMeshData,
                   std::vector<std::pair<EntityId, int>>& elemSidePairs)
  {
    using SkinFaceConnections = typename SideAdjacencyGraph<EntityId, Topology>::FaceConnections;
    using SkinFaceConnection = typename SideAdjacencyGraph<EntityId, Topology>::FaceConnection;


    populate_skin_blocks(textMeshData.partIds);

    if (!m_skinBlocks.empty()) {
      set_text_mesh_data(textMeshData);
      BaseClass::create_graph(m_skinBlocks);
      reset_text_mesh_data();

      for (auto iter = BaseClass::begin(); iter != BaseClass::end(); iter++) {
        size_t elemIndex = iter->first;
        const SkinFaceConnections &faceConnections = iter->second;

        std::vector<bool> hasConnection(faceConnections.numSides, false);
        for (const SkinFaceConnection &connection : faceConnections.connections) {
          hasConnection[connection.thisSide - 1] = true;
        }

        for (unsigned i = 0; i < faceConnections.numSides; i++) {
          if (!hasConnection[i]) {
            EntityId elemId = textMeshData.elementDataVec[elemIndex].identifier;
            int side = i + 1;
            elemSidePairs.push_back(std::make_pair(elemId, side));
          }
        }
      }
    }
  }

private:
  void populate_skin_blocks(const PartIdMapping &partIds)
  {
    bool skinAll = false;
    for (const std::string &block : m_skinBlocks) {
      if (0 == strcasecmp("all", block.c_str())) {
        skinAll = true;
        break;
      }
    }

    if (skinAll) {
      m_skinBlocks = partIds.get_part_names();
    }
  }

  void set_text_mesh_data(const TextMeshData<EntityId, Topology> &textMeshData)
  {
    m_textMeshData = &textMeshData;
  }

  void reset_text_mesh_data()
  {
    m_textMeshData = nullptr;
  }

  const TextMeshData<EntityId, Topology> *m_textMeshData{nullptr};
  std::vector<std::string> m_skinBlocks{};
};

}  // namespace text_mesh

#endif
