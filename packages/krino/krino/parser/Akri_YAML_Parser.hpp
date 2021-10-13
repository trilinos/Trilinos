// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_YAML_Parser_h
#define Akri_YAML_Parser_h

#include <Akri_YAML.hpp>
#include <string>
#include <vector>

namespace krino{
namespace YAML_Parser {
  
void parse();
std::string info(const YAML::Node & node);

const YAML::Node
get_if_present(const YAML::Node& node, const std::string& key);

/// these can be used to check and ensure a type of yaml node is as expected
const YAML::Node
get_type_if_present(const YAML::Node& node, const std::string& key, YAML::NodeType::value type);

inline const YAML::Node
get_null_if_present(const YAML::Node& node, const std::string& key) {return get_type_if_present(node, key, YAML::NodeType::Null); }

inline const YAML::Node
get_scalar_if_present(const YAML::Node& node, const std::string& key) {return get_type_if_present(node, key, YAML::NodeType::Scalar); }

inline const YAML::Node
get_sequence_if_present(const YAML::Node& node, const std::string& key) {return get_type_if_present(node, key, YAML::NodeType::Sequence); }

inline const YAML::Node
get_map_if_present(const YAML::Node& node, const std::string& key) {return get_type_if_present(node, key, YAML::NodeType::Map); }

/// Doesn't change @param result unless the @param key is present in the @param node
template<typename T>
bool get_if_present(const YAML::Node & node, const std::string& key, T& result)
{
  const YAML::Node value = get_if_present(node, key);
  if (value)
  {
    result = value.as<T>();
    return true;
  }
  return false;
}

template<>
bool get_if_present(const YAML::Node & node, const std::string& key, bool& result);

} // namespace parser
} // namespace krino

#endif // Akri_YAML_Parser_h
