// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#ifndef TextMeshAssembly_hpp
#define TextMeshAssembly_hpp

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
#include "TextMeshEntityGroup.hpp"


// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace text_mesh {

using ErrorHandler = std::function<void(const std::ostringstream &)>;

enum AssemblyType { ASSEMBLY, BLOCK, SIDESET, NODESET, INVALID_ASSEMBLY };

inline std::ostream &operator<<(std::ostream &out, const AssemblyType &t)
{
  switch (t) {
    case AssemblyType::ASSEMBLY:
      return out << "ASSEMBLY";
      break;
    case AssemblyType::BLOCK:
      return out << "ELEMENT_BLOCK";
      break;
    case AssemblyType::SIDESET:
      return out << "SIDESET";
      break;
    case AssemblyType::NODESET:
      return out << "NODESET";
      break;
    default:
      return out << "INVALID";
      break;
  }
  return out << "INVALID[" << (unsigned) t << "]";
}

using AssemblyDataType = std::string;

struct AssemblyData : public EntityGroupData<AssemblyDataType> {
  using DataType = AssemblyDataType;

  void set_assembly_type(AssemblyType type_) { assemblyType = type_; }
  AssemblyType get_assembly_type() const {return assemblyType;}

private:
  AssemblyType assemblyType = INVALID_ASSEMBLY;
};

template <typename EntityId>
class Assemblies : public EntityGroup<AssemblyData>
{
public:
  using BaseClass = EntityGroup<AssemblyData>;

  Assemblies() : EntityGroup<AssemblyData>("ASSEMBLY", "ASSEMBLY_", {"BLOCK_", "SURFACE_", "NODELIST_"}) {}

  bool is_cyclic(const std::string& assembly) const
  {
    initialize_graph();
    return check_for_cycle(assembly);
  }

  bool is_cyclic() const
  {
    for(const std::string& assembly : BaseClass::get_part_names())
      if(is_cyclic(assembly)) {
        return true;
      }

    return false;
  }

  std::vector<std::string>&& get_forward_traversal_list(const std::string& assembly) const
  {
    initialize_graph();
    fill_traversal(assembly);

    return std::move(m_traversalList);
  }

  std::vector<std::string>&& get_reverse_traversal_list(const std::string& assembly) const
  {
    initialize_graph();
    fill_traversal(assembly);
    std::reverse(m_traversalList.begin(), m_traversalList.end());

    return std::move(m_traversalList);
  }

private:
  void fill_traversal(const std::string& assembly) const
  {
    const AssemblyData *assemblyData = BaseClass::get_group_data(assembly);
    if (nullptr != assemblyData) {
      if (m_visitedNodes[assembly] == false) {
        m_visitedNodes[assembly] = true;
        m_traversalList.push_back(assembly);

        if(assemblyData->get_assembly_type() == AssemblyType::ASSEMBLY) {
          for (const std::string& member : assemblyData->data) {
            fill_traversal(member);
          }
        }
      }
    }
  }

  bool check_for_cycle(const std::string& assembly) const
  {
    bool isCyclic = false;
    const AssemblyData *assemblyData = BaseClass::get_group_data(assembly);
    if (nullptr != assemblyData) {
      if (m_visitedNodes[assembly] == true) {
        isCyclic = true;
      } else {
        m_visitedNodes[assembly] = true;

        if(assemblyData->get_assembly_type() == AssemblyType::ASSEMBLY) {
          for (const std::string& member : assemblyData->data) {
            isCyclic |= check_for_cycle(member);
          }
        }
      }
    }
    return isCyclic;
  }

  void initialize_graph() const
  {
    m_traversalList.clear();
    m_traversalList.reserve(BaseClass::size());

    for(const std::string& name : BaseClass::get_part_names()) {
      m_visitedNodes[name] = false;
    }
  }

  mutable std::unordered_map<std::string, bool> m_visitedNodes;
  mutable std::vector<std::string> m_traversalList;
};

class AssemblyParser
{
 public:
  AssemblyParser()
  : m_assemblyType(INVALID_ASSEMBLY)
  {
    ErrorHandler errorHandler = [](const std::ostringstream &errmsg) { default_error_handler(errmsg); };
    set_error_handler(errorHandler);
  }

  void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

  std::string get_name() { return m_name; }

  AssemblyType get_assembly_type() const { return m_assemblyType; }

  const std::vector<std::string>& get_assembly_data() { return m_members; }

  void parse(const std::string &parseData)
  {
    auto options = get_tokens(parseData, ";");

    for (const auto &option : options) {
      parse_option_group(option);
    }
  }

  void verify_parse() const
  {
    if (m_assemblyType == INVALID_ASSEMBLY) {
      std::ostringstream errmsg;
      errmsg << "Error!  Assembly with name: " << m_name
             << " does not have a defined type.";
      m_errorHandler(errmsg);
    }
  }

 private:
  void parse_option(std::string optionName, const std::string &optionValue)
  {
    convert_to_lower_case(optionName);

    if (optionName == "name") {
      parse_name(optionValue);
    } else if (optionName == "type") {
      parse_assembly_type(optionValue);
    } else if (optionName == "member") {
      parse_assembly_members(optionValue);
    } else {
      std::ostringstream errmsg;
      errmsg << "Unrecognized assembly option: " << optionName;
      m_errorHandler(errmsg);
    }
  }

  void parse_option_group(const std::string &option)
  {
    if (!option.empty()) {
      auto optionTokens = get_tokens(option, "=");

      if (optionTokens.size() != 2) {
        std::ostringstream errmsg;
        errmsg << "Unrecognized assembly option: " << option;
        m_errorHandler(errmsg);
      }

      parse_option(optionTokens[0], optionTokens[1]);
    }
  }

  void parse_name(const std::string &data) { m_name = data; }

  void parse_assembly_type(std::string type)
  {
    convert_to_lower_case(type);

    if (type == "assembly") {
      m_assemblyType = ASSEMBLY;
    } else if (type == "block") {
      m_assemblyType = BLOCK;
    } else if (type == "sideset") {
      m_assemblyType = SIDESET;
    } else if (type == "nodeset") {
      m_assemblyType = NODESET;
    } else {
      std::ostringstream errmsg;
      errmsg << "Unrecognized assembly type: " << type;
      m_errorHandler(errmsg);
    }
  }

  void parse_assembly_members(const std::string &data)
  {
    std::vector<std::string> assemblyData = get_tokens(data, ",");
    for(std::string& member : assemblyData) {
      convert_to_upper_case(member);
    }

    m_members = assemblyData;
  }

  std::vector<std::string> m_members{};
  std::string m_name{};
  AssemblyType m_assemblyType{INVALID_ASSEMBLY};
  ErrorHandler m_errorHandler;
};

}  // namespace text_mesh

#endif
