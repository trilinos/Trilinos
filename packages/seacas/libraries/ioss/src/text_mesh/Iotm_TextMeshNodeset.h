// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.

#pragma once

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

#include "Iotm_TextMeshFuncs.h"
#include "Iotm_TextMeshEntityGroup.h"

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
namespace Iotm {
  namespace text_mesh {

    using ErrorHandler = std::function<void(const std::ostringstream &)>;

    template <typename EntityId> using NodesetDataType = EntityId;

    template <typename EntityId>
    struct NodesetData : public EntityGroupData<NodesetDataType<EntityId>>
    {
      using DataType = NodesetDataType<EntityId>;
    };

    template <typename EntityId> class Nodesets : public EntityGroup<NodesetData<EntityId>>
    {
    public:
      Nodesets()
          : EntityGroup<NodesetData<EntityId>>("NODESET", "NODELIST_",
                                               {"BLOCK_", "SURFACE_", "ASSEMBLY_"})
      {
      }
    };

    template <typename EntityId> class NodesetParser
    {
    public:
      NodesetParser()
      {
        ErrorHandler errorHandler = [](const std::ostringstream &errmsg) {
          default_error_handler(errmsg);
        };
        set_error_handler(errorHandler);
      }

      void set_error_handler(ErrorHandler errorHandler) { m_errorHandler = errorHandler; }

      std::string get_name() { return m_name; }

      const std::vector<EntityId> &get_nodeset_data() { return m_nodeList; }

      // Expected format for nodeset string data is:
      // "[name=<name>;] data=node_1,node_2,....,node_n"
      // Only data keyword is required
      void parse(const std::string &parseData)
      {
        auto options = get_tokens(parseData, ";");

        for (const auto &option : options) {
          parse_option_group(option);
        }
      }

    private:
      void parse_option(std::string optionName, const std::string &optionValue)
      {
        convert_to_lowercase(optionName);

        if (optionName == "name") {
          parse_name(optionValue);
        }
        else if (optionName == "data") {
          parse_node_data(optionValue);
        }
        else {
          std::ostringstream errmsg;
          errmsg << "Unrecognized nodeset option: " << optionName;
          m_errorHandler(errmsg);
        }
      }

      void parse_option_group(const std::string &option)
      {
        if (!option.empty()) {
          auto optionTokens = get_tokens(option, "=");

          if (optionTokens.size() != 2) {
            std::ostringstream errmsg;
            errmsg << "Unrecognized nodeset option: " << option;
            m_errorHandler(errmsg);
          }

          parse_option(optionTokens[0], optionTokens[1]);
        }
      }

      void parse_name(const std::string &data) { m_name = data; }

      void parse_node_data(const std::string &data)
      {
        auto nodesetData = get_tokens(data, ",");

        for (const std::string &nodeString : nodesetData) {
          if (!is_positive_number(nodeString)) {
            std::ostringstream errmsg;
            errmsg << "Unrecognized nodeset node id: " << nodeString;
            m_errorHandler(errmsg);
          }
          EntityId node = std::stoull(nodeString);
          m_nodeList.push_back(node);
        }
      }

      std::vector<EntityId> m_nodeList{};
      std::string           m_name{};
      ErrorHandler          m_errorHandler;
    };

  } // namespace text_mesh
} // namespace Iotm
