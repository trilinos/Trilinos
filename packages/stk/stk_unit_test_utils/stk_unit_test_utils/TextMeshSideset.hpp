// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshSideset_hpp
#define TextMeshSideset_hpp

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
#include "TextMeshSidesetSplitter.hpp"
#include "TextMeshSidesetSkinner.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

template <typename EntityId>
using SidesetDataType = std::pair<EntityId, int>;

template <typename EntityId, typename Topology>
struct SidesetData : public EntityGroupData<SidesetDataType<EntityId>> {
  using DataType = SidesetDataType<EntityId>;
  using BaseClass = EntityGroupData<SidesetDataType<EntityId>>;

  void set_split_type(SplitType splitType) { m_sidesetSplitter.set_split_type(splitType); }
  SplitType get_split_type() const { return m_sidesetSplitter.get_split_type(); }

  void set_skin_blocks(const std::vector<std::string> &skinBlocks)
  {
    m_sidesetSkinner.set_skin_blocks(skinBlocks);
  }

  void set_error_handler(ErrorHandler errorHandler)
  {
    m_sidesetSkinner.set_error_handler(errorHandler);
    m_sidesetSplitter.set_error_handler(errorHandler);
  }

  void split(const std::vector<ElementData<EntityId, Topology>> &elementData)
  {
    m_sidesetSplitter.split(*this, elementData);
  }

  void skin_blocks(const TextMeshData<EntityId, Topology> &textMeshData)
  {
    m_sidesetSkinner.skin_blocks(textMeshData, BaseClass::data);
  }

  std::vector<size_t> get_sideblock_indices_local_to_proc(const SideBlockInfo &info, int proc) const
  {
    return m_sidesetSplitter.get_indices_local_to_proc(info.sideIndex, proc);
  }

  SideBlockInfo get_side_block_info(const std::string &sideBlockName) const
  {
    return m_sidesetSplitter.get_side_block_info(sideBlockName);
  }

  std::vector<SideBlockInfo> get_side_block_info() const { return m_sidesetSplitter.get_side_block_info(); }

  bool has_default_exodus_name() const
  {
    if (BaseClass::has_name()) {
      std::pair<unsigned, bool> result = get_id_from_part_name(BaseClass::name, "SURFACE_");
      return result.second;
    }

    return false;
  }

  SidesetSkinner<EntityId, Topology> m_sidesetSkinner;
  SidesetSplitter<EntityId, Topology> m_sidesetSplitter;
};

template <typename EntityId, typename Topology>
class Sidesets : public EntityGroup<SidesetData<EntityId, Topology>>
{
 public:
  using BaseClass = EntityGroup<SidesetData<EntityId, Topology>>;

  Sidesets() : BaseClass("SIDESET", "SURFACE_", {"BLOCK_", "NODELIST_", "ASSEMBLY_"}) {}

  void set_error_handler(ErrorHandler errorHandler) override
  {
    BaseClass::set_error_handler(errorHandler);

    for (SidesetData<EntityId, Topology> &sidesetData : BaseClass::m_groupDataVec) {
      sidesetData.set_error_handler(errorHandler);
    }
  }

  void finalize_parse(const TextMeshData<EntityId, Topology> &data)
  {
    BaseClass::finalize_parse();

    for (SidesetData<EntityId, Topology> &sidesetData : BaseClass::m_groupDataVec) {
      sidesetData.skin_blocks(data);
      sidesetData.split(data.elementDataVec);
    }
  }
};

template <typename EntityId>
class SidesetParser
{
 public:
  SidesetParser() : m_splitType(NO_SPLIT)
  {
    ErrorHandler errorHandler = [](const std::ostringstream &errmsg) { default_error_handler(errmsg); };
    set_error_handler(errorHandler);
  }

  void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  std::string get_name() { return m_name; }

  const std::vector<std::pair<EntityId, int>> &get_sideset_data() { return m_elemSidePairs; }

  SplitType get_split_type() { return m_splitType; }

  const std::vector<std::string> &get_skin_blocks() const { return m_skinnedBlocks; }

  void parse(const std::string &parseData)
  {
    auto options = get_tokens(parseData, ";");

    for (const auto &option : options) {
      parse_option_group(option);
    }
  }

  void verify_parse() const
  {
    if (!m_skinnedBlocks.empty() && !m_elemSidePairs.empty()) {
      std::ostringstream errmsg;
      errmsg << "Error!  Sideset with name: " << m_name
             << " is attempting to set element/side pair data *AND* use skinning.";
      m_errorHandler(errmsg);
    }
  }

 private:
  void parse_option(std::string optionName, const std::string &optionValue)
  {
    convert_to_lower_case(optionName);

    if (optionName == "name") {
      parse_name(optionValue);
    } else if (optionName == "data") {
      parse_element_side_pairs(optionValue);
    } else if (optionName == "split") {
      parse_split_type(optionValue);
    } else if (optionName == "skin") {
      parse_skin_blocks(optionValue);
    } else {
      std::ostringstream errmsg;
      errmsg << "Unrecognized sideset option: " << optionName;
      m_errorHandler(errmsg);
    }
  }

  void parse_option_group(const std::string &option)
  {
    if (!option.empty()) {
      auto optionTokens = get_tokens(option, "=");

      if (optionTokens.size() != 2) {
        std::ostringstream errmsg;
        errmsg << "Unrecognized sideset option: " << option;
        m_errorHandler(errmsg);
      }

      parse_option(optionTokens[0], optionTokens[1]);
    }
  }

  void parse_name(const std::string &data) { m_name = data; }

  void parse_skin_blocks(const std::string &data) { m_skinnedBlocks = get_tokens(data, ","); }

  void parse_element_side_pairs(const std::string &data)
  {
    auto sidesetData = get_tokens(data, ",");

    if (sidesetData.size() % 2 != 0) {
      std::ostringstream errmsg;
      errmsg << "Unmatched element/ordinal pairs in sideset data: " << data;
      m_errorHandler(errmsg);
    }

    for (unsigned i = 0; i < sidesetData.size(); i += 2) {
      EntityId elem = std::stoull(sidesetData[i]);
      int side = std::stoi(sidesetData[i + 1]);

      if (side <= 0) {
        std::ostringstream errmsg;
        errmsg << "Invalid element/ordinal pair {" << sidesetData[i] << "," << sidesetData[i + 1] << "}";
        m_errorHandler(errmsg);
      }

      m_elemSidePairs.push_back(std::make_pair(elem, side));
    }
  }

  void parse_split_type(std::string splitName)
  {
    convert_to_lower_case(splitName);

    if (splitName == "none") {
      m_splitType = NO_SPLIT;
    } else if (splitName == "block") {
      m_splitType = ELEMENT_BLOCK;
    } else if (splitName == "topology") {
      m_splitType = TOPOLOGY;
    } else {
      std::ostringstream errmsg;
      errmsg << "Unrecognized sideset split type: " << splitName;
      m_errorHandler(errmsg);
    }
  }

  std::vector<std::pair<EntityId, int>> m_elemSidePairs{};
  std::string m_name{};
  std::vector<std::string> m_skinnedBlocks{};
  SplitType m_splitType{NO_SPLIT};
  ErrorHandler m_errorHandler;
};

}  // namespace text_mesh

#endif
