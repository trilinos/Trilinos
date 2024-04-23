// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshEntityGroup_hpp
#define TextMeshEntityGroup_hpp

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
#include <limits>
#include <strings.h>

#include "TextMeshFuncs.hpp"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

template <typename T>
struct EntityGroupData {
  using DataType = T;
  static constexpr unsigned INVALID_ID = std::numeric_limits<unsigned>::max();

  bool hasInputName = false;
  unsigned id = INVALID_ID;
  std::string name = "";
  std::string type = "";
  std::vector<DataType> data{};

  bool has_valid_id() const { return id != 0 && id != INVALID_ID; }
  bool has_name() const { return !name.empty(); }
};

template <typename GroupData>
class EntityGroup
{
 private:
  using DataType = typename GroupData::DataType;

 public:
  EntityGroup(const std::string &type,
      const std::string &namePrefix,
      const std::vector<std::string> &invalidNamePrefixes)
      : m_idsAssigned(false), m_type(type), m_exodusPrefix(namePrefix), m_invalidPrefixes(invalidNamePrefixes)
  {
    set_error_handler([](const std::ostringstream &errmsg) { default_error_handler(errmsg); });
  }

  virtual ~EntityGroup() {}

  virtual void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  GroupData *add_group_data(const std::string &name, const std::vector<DataType> &data)
  {
    GroupData groupData;
    groupData.data = data;
    groupData.type = m_type;

    if (!name.empty()) {
      verify_name(name);
      groupData.name = name;
      groupData.hasInputName = true;
    }

    m_groupDataVec.push_back(groupData);

    return &m_groupDataVec.back();
  }

  void finalize_parse()
  {
    assign_id_from_default_exodus_name();
    assign_id_and_name_for_empty_name();
    assign_id_for_non_default_exodus_name();

    if (m_groupDataVec.size() != m_groupDataMap.size()) {
      std::ostringstream errmsg;
      errmsg << "Error populating " << m_type << " map";
      m_errorHandler(errmsg);
    }
    m_idsAssigned = true;
  }

  size_t size() const { return m_groupDataVec.size(); }

  const std::vector<GroupData> &get_group_data() const { return m_groupDataVec; }

  const std::vector<std::string> &get_part_names() const { return m_partNames; }

  const std::string &get_group_type() const { return m_type; }

  const GroupData *get_group_data(unsigned id) const
  {
    if (is_assigned(id)) {
      auto iter = m_parts.find(id);
      return &m_groupDataVec[m_groupDataMap[iter->second]];
    }

    return nullptr;
  }

  const GroupData *get_group_data(std::string name) const
  {
    convert_to_upper_case(name);
    if (is_registered(name)) {
      return &m_groupDataVec[m_groupDataMap[name]];
    }

    return nullptr;
  }

  bool is_registered(const std::string &name) const { return m_ids.count(name) > 0; }

 protected:
  EntityGroup();

  unsigned get_unassigned_id() const
  {
    unsigned nextPartId = 1;
    while (is_assigned(nextPartId)) nextPartId++;
    return nextPartId;
  }

  void validate_group_meta_data(const GroupData &groupData)
  {
    if (!groupData.has_name()) {
      std::ostringstream errmsg;
      errmsg << m_type << " has no name";
      m_errorHandler(errmsg);
    }

    if (!groupData.has_valid_id()) {
      std::ostringstream errmsg;
      errmsg << m_type << " named " << groupData.name << " has invalid id";
      m_errorHandler(errmsg);
    }

    if (is_registered(groupData.name)) {
      std::ostringstream errmsg;
      errmsg << "Multiple declarations of " << m_type << ": " << groupData.name;
      m_errorHandler(errmsg);
    }
  }

  void assign(size_t index)
  {
    GroupData &groupData = m_groupDataVec[index];

    convert_to_upper_case(groupData.name);
    validate_group_meta_data(groupData);

    m_partNames.push_back(groupData.name);
    m_ids[groupData.name] = groupData.id;
    m_parts[groupData.id] = groupData.name;
    m_groupDataMap[groupData.name] = index;
  }

  void assign_id_from_default_exodus_name()
  {
    for (size_t i = 0; i < m_groupDataVec.size(); i++) {
      GroupData &groupData = m_groupDataVec[i];
      if (groupData.has_name()) {
        std::pair<unsigned, bool> result = get_id_from_part_name(groupData.name, m_exodusPrefix);

        if (result.second) {
          groupData.id = result.first;
          assign(i);
        }
      }
    }
  }

  void assign_id_and_name_for_empty_name()
  {
    for (size_t i = 0; i < m_groupDataVec.size(); i++) {
      GroupData &groupData = m_groupDataVec[i];
      if (!groupData.has_name()) {
        unsigned id = get_unassigned_id();

        std::ostringstream oss;
        oss << m_exodusPrefix;
        oss << id;
        std::string name = oss.str();

        groupData.id = id;
        groupData.name = name;
        assign(i);
      }
    }
  }

  void assign_id_for_non_default_exodus_name()
  {
    for (size_t i = 0; i < m_groupDataVec.size(); i++) {
      GroupData &groupData = m_groupDataVec[i];
      if (groupData.has_name()) {
        std::pair<unsigned, bool> result = get_id_from_part_name(groupData.name, m_exodusPrefix);

        if (!result.second) {
          groupData.id = get_unassigned_id();
          assign(i);
        }
      }
    }
  }

  bool is_assigned(unsigned id) const { return m_parts.count(id) > 0; }

  void verify_name(const std::string &name)
  {
    for (const std::string &invalidPrefix : m_invalidPrefixes) {
      const unsigned prefixLength = invalidPrefix.length();
      const std::string namePrefix = name.substr(0, prefixLength);

      if (strcasecmp(namePrefix.c_str(), invalidPrefix.c_str()) == 0) {
        std::ostringstream errmsg;
        errmsg << "Invalid name '" << name << "' for a " << m_type << " part";
        m_errorHandler(errmsg);
      }
    }
  }

  std::vector<std::string> m_partNames{};
  mutable std::unordered_map<std::string, unsigned> m_ids;
  mutable std::unordered_map<unsigned, std::string> m_parts;
  mutable bool m_idsAssigned{false};
  mutable std::unordered_map<std::string, size_t> m_groupDataMap;

  std::string m_type{};
  std::string m_exodusPrefix{};
  std::vector<std::string> m_invalidPrefixes{};
  std::vector<GroupData> m_groupDataVec{};

  ErrorHandler m_errorHandler;
};

}  // namespace text_mesh

#endif
