// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PARSER_AKRI_PARSER_HPP_
#define KRINO_KRINO_PARSER_AKRI_PARSER_HPP_

#include <Akri_YAML.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <string>

namespace krino { class Simulation; }

namespace krino {
namespace Parser {

#ifdef KRINO_HAVE_YAML

class Node;

class ConstNodeIterator {
public:
  ConstNodeIterator(const YAML::Node::const_iterator & iter) : mIter(iter) {}
  Node operator*() const;
  ConstNodeIterator& operator++() { mIter++; return *this; }
  ConstNodeIterator operator++(int) { ConstNodeIterator tmp = *this; ++(*this); return tmp; }

  friend bool operator==(const ConstNodeIterator& a, const ConstNodeIterator& b) { return a.mIter == b.mIter; }
  friend bool operator!=(const ConstNodeIterator& a, const ConstNodeIterator& b) { return a.mIter != b.mIter; }
private:
  YAML::Node::const_iterator mIter;
};

class Node {
public:
  Node() : mNode(YAML::Node()) {}
  Node(const YAML::Node & node) : mNode(node) {}

  operator bool() const { return bool(mNode); }
  bool is_scalar() const;

  Node get_if_present(const std::string& key) const;
  Node get_null_if_present(const std::string& key) const;
  Node get_map_if_present(const std::string& key) const;
  Node get_sequence_if_present(const std::string& key) const;
  Node get_scalar_if_present(const std::string& key) const;

  template<typename T>
  bool get_if_present(const std::string& key, T& result) const
  {
    const Node value = get_if_present(key);
    if (value)
    {
      result = value.as<T>();
      return true;
    }
    return false;
  }

  template<typename T>
  T as() const
  {
    return mNode.as<T>();
  }

  ConstNodeIterator begin() const { return ConstNodeIterator(mNode.begin()); }
  ConstNodeIterator end() const { return ConstNodeIterator(mNode.end()); }

  std::string info() const;

private:
  std::string line_info() const;
  Node get_type_if_present(const std::string& key, YAML::NodeType::value type) const;
  YAML::NodeType::value Type() const { return mNode.Type(); }

  const YAML::Node mNode;
};

template<>
bool Node::as<bool>() const;

inline Node ConstNodeIterator::operator*() const
{
  return Node(*mIter);
}

#else

class Node {
public:
  operator bool() const { return false; }
  bool is_scalar() const { return false; }

  Node get_if_present(const std::string& key) const { return Node(); }
  Node get_null_if_present(const std::string& key) const { return Node(); }
  Node get_map_if_present(const std::string& key) const { return Node(); }
  Node get_sequence_if_present(const std::string& key) const { return Node(); }
  Node get_scalar_if_present(const std::string& key) const { return Node(); }

  template<typename T>
  bool get_if_present(const std::string& key, T& result) const { return false; }

  template<typename T>
  T as() const { return T(); }

  Node * begin() const { return nullptr; }
  Node * end() const { return nullptr; }

  std::string info() const { return std::string(); }
};

#endif

void parse(Simulation & simulation);

}
}

#endif /* KRINO_KRINO_PARSER_AKRI_PARSER_HPP_ */
