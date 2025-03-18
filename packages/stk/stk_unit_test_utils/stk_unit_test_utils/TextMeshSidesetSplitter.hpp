// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshSidesetSplitter_hpp
#define TextMeshSidesetSplitter_hpp

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

#include "TextMeshFuncs.hpp"
#include "TextMeshDataTypes.hpp"
#include "TextMeshEntityGroup.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

struct SideBlockInfo {
  std::string name{};
  std::string parentName{};
  std::string sideTopology{};
  std::string elementTopology{};
  std::string touchingBlock{};
  std::vector<size_t> sideIndex{};
  unsigned numNodesPerSide{0};
};

enum SplitType { TOPOLOGY, ELEMENT_BLOCK, NO_SPLIT, INVALID_SPLIT };

inline std::ostream &operator<<(std::ostream &out, const SplitType &t)
{
  switch (t) {
    case SplitType::TOPOLOGY:
      return out << "TOPOLOGY";
      break;
    case SplitType::ELEMENT_BLOCK:
      return out << "ELEMENT_BLOCK";
      break;
    case SplitType::NO_SPLIT:
      return out << "NO_SPLIT";
      break;
    default:
      return out << "INVALID";
      break;
  }
  return out << "INVALID[" << (unsigned) t << "]";
}

template <typename EntityId>
using SidesetDataType = std::pair<EntityId, int>;

template <typename EntityId, typename Topology>
struct SidesetData;

template <typename EntityId, typename Topology>
class SidesetSplitter
{
 public:
  SidesetSplitter(SplitType splitType) : m_splitType(splitType)
  {
    ErrorHandler errorHandler = [](const std::ostringstream &errmsg) { default_error_handler(errmsg); };
    set_error_handler(errorHandler);
  }

  SidesetSplitter() : m_splitType(INVALID_SPLIT)
  {
    ErrorHandler errorHandler = [](const std::ostringstream &errmsg) { default_error_handler(errmsg); };
    set_error_handler(errorHandler);
  }

  void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  void split(
      const SidesetData<EntityId, Topology> &sideset, const std::vector<ElementData<EntityId, Topology>> &elementData)
  {
    m_splitMap.clear();
    m_sidesetName = sideset.name;

    if (get_split_type() == SplitType::TOPOLOGY) {
      split_by_topology(sideset, elementData);
    } else if (get_split_type() == SplitType::ELEMENT_BLOCK) {
      split_by_element_block(sideset, elementData);
    } else if (get_split_type() == SplitType::NO_SPLIT) {
      split_by_no_split(sideset, elementData);
    } else {
      std::ostringstream errmsg;
      errmsg << "Invalid split type: " << get_split_type();
      m_errorHandler(errmsg);
    }

    build_index_proc_map(sideset, elementData);
  }

  std::vector<SideBlockInfo> get_side_block_info() const
  {
    std::vector<SideBlockInfo> infoVec;
    infoVec.reserve(m_splitMap.size());

    for (auto iter = m_splitMap.begin(); iter != m_splitMap.end(); iter++) {
      const std::string &sideBlockName = iter->first;
      SideBlockInfo info = get_side_block_info(sideBlockName);
      infoVec.push_back(info);
    }
    return infoVec;
  }

  std::vector<size_t> get_indices_local_to_proc(const std::vector<size_t> &index, int proc) const
  {
    std::vector<size_t> indexForProc;
    indexForProc.reserve(index.size());

    for (size_t elemPairIndex : index) {
      if (is_index_local_to_proc(elemPairIndex, proc)) {
        indexForProc.push_back(elemPairIndex);
      }
    }

    indexForProc.resize(indexForProc.size());
    return indexForProc;
  }

  SideBlockInfo get_side_block_info(const std::string &name) const
  {
    SideBlockInfo info;

    auto iter = m_splitMap.find(name);
    if (iter != m_splitMap.end()) {
      const SplitData &splitData = iter->second;

      info.name = name;
      info.parentName = splitData.sidesetName;
      info.sideTopology = splitData.sideTopology;
      info.elementTopology = splitData.elemTopology;
      info.numNodesPerSide = splitData.sideNodeCount;
      info.touchingBlock = splitData.touchingBlock;
      info.sideIndex = splitData.index;
    }
    return info;
  }

  SplitType get_split_type() const { return m_splitType; }
  void set_split_type(SplitType inputSplitType) { m_splitType = inputSplitType; }

 private:
  void build_index_proc_map(
      const SidesetData<EntityId, Topology> &sideset, const std::vector<ElementData<EntityId, Topology>> &elementData)
  {
    for (size_t i = 0; i < sideset.data.size(); ++i) {
      const SidesetDataType<EntityId> &elemSidePair = sideset.data[i];
      EntityId elemId = elemSidePair.first;

      auto iter = bound_search(elementData.begin(), elementData.end(), elemId, ElementDataLess<EntityId, Topology>());
      if (iter == elementData.end()) {
        std::ostringstream errmsg;
        errmsg << "Error!  Sideset with id: " << sideset.id << " and name: " << sideset.name
               << " has reference to invalid element '" << elemId << "'.";
        m_errorHandler(errmsg);
      }

      m_indexProcMap[i] = iter->proc;
    }
  }

  bool is_index_local_to_proc(size_t elemPairIndex, int proc) const
  {
    auto iter = m_indexProcMap.find(elemPairIndex);

    if (iter == m_indexProcMap.end()) {
      std::ostringstream errmsg;
      errmsg << "Sideset with name: " << m_sidesetName << " is referencing an invalid index " << elemPairIndex;
      m_errorHandler(errmsg);
    }

    return iter->second == proc;
  }

  struct SplitData {
    bool metaDataSet{false};
    std::string sidesetName{};
    std::string touchingBlock{};
    std::string elemTopology{};
    std::string sideTopology{};
    int sideNodeCount{-1};
    std::vector<size_t> index{};

    SplitData() : metaDataSet(false), sideNodeCount(-1) {}
  };

  void fill_split_data(std::string key, size_t index, const ElementData<EntityId, Topology> &elemData, int side)
  {
    convert_to_upper_case(key);

    SplitData &splitData = m_splitMap[key];

    splitData.index.push_back(index);

    if (!splitData.metaDataSet) {
      splitData.sidesetName = m_sidesetName;
      splitData.elemTopology = elemData.topology.name();
      splitData.sideTopology = elemData.topology.side_topology_name(side);
      splitData.sideNodeCount = elemData.topology.side_topology_num_nodes(side);

      if (get_split_type() == ELEMENT_BLOCK) {
        splitData.touchingBlock = elemData.partName;
      }

      splitData.metaDataSet = true;
    }
  }

  using Criterion = std::function<std::string(
      const SidesetData<EntityId, Topology> &sideset, const ElementData<EntityId, Topology> &elemData, int side)>;

  void split_by_criterion(const SidesetData<EntityId, Topology> &sideset,
      const std::vector<ElementData<EntityId, Topology>> &elementData,
      Criterion criterion)
  {
    for (size_t index = 0; index < sideset.data.size(); ++index) {
      const SidesetDataType<EntityId> &elemSidePair = sideset.data[index];
      EntityId elemId = elemSidePair.first;
      int side = elemSidePair.second;

      auto iter = bound_search(elementData.begin(), elementData.end(), elemId, ElementDataLess<EntityId, Topology>());
      if (iter == elementData.end()) {
        std::ostringstream errmsg;
        errmsg << "Error!  Sideset with id: " << sideset.id << " and name: " << sideset.name
               << " has reference to invalid element '" << elemId << "'.";
        m_errorHandler(errmsg);
      }

      std::string key = criterion(sideset, *iter, side);
      fill_split_data(key, index, *iter, side);
    }
  }

  void split_by_topology(
      const SidesetData<EntityId, Topology> &sideset, const std::vector<ElementData<EntityId, Topology>> &elementData)
  {
    Criterion criterion = [](const SidesetData<EntityId, Topology> &sideSet,
                              const ElementData<EntityId, Topology> &elemData, int side) {
      if (sideSet.has_default_exodus_name()) {
        return "SURFACE_" + elemData.topology.name() + "_" + elemData.topology.side_topology_name(side) + "_" +
               std::to_string(sideSet.id);
      }
      return sideSet.name + "_" + elemData.topology.name() + "_" + elemData.topology.side_topology_name(side);
    };

    split_by_criterion(sideset, elementData, criterion);
  }

  void split_by_element_block(
      const SidesetData<EntityId, Topology> &sideset, const std::vector<ElementData<EntityId, Topology>> &elementData)
  {
    Criterion criterion = [](const SidesetData<EntityId, Topology> &sideSet,
                              const ElementData<EntityId, Topology> &elemData, int side) {
      if (sideSet.has_default_exodus_name()) {
        return "SURFACE_" + elemData.partName + "_" + elemData.topology.side_topology_name(side) + "_" +
               std::to_string(sideSet.id);
      }
      return sideSet.name + "_" + elemData.partName + "_" + elemData.topology.side_topology_name(side);
    };

    split_by_criterion(sideset, elementData, criterion);
  }

  void split_by_no_split(
      const SidesetData<EntityId, Topology> &sideset, const std::vector<ElementData<EntityId, Topology>> & /*elementData*/)
  {
    std::vector<size_t> splitIndex(sideset.data.size());
    std::iota(std::begin(splitIndex), std::end(splitIndex), 0);
    SplitData &splitData = m_splitMap[sideset.name];

    splitData.index = splitIndex;
    splitData.sidesetName = m_sidesetName;
    splitData.elemTopology = "unknown";
    splitData.sideTopology = "unknown";
    splitData.sideNodeCount = -1;
    splitData.metaDataSet = true;
  }

  SplitType m_splitType{INVALID_SPLIT};
  std::string m_sidesetName{};

  std::unordered_map<size_t, int> m_indexProcMap;
  std::unordered_map<std::string, SplitData> m_splitMap;
  ErrorHandler m_errorHandler;
};

}  // namespace text_mesh

#endif
