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
  Node() : mNode(YAML::Node(YAML::NodeType::Undefined)) {}
  Node(const YAML::Node & node) : mNode(node) {}

  Node(const Node &) = default;
  Node & operator=(const Node &) = default;

  operator bool() const { return bool(mNode); }
  bool is_scalar() const;
  bool is_map() const;
  bool is_sequence() const;
  bool is_empty_scalar() const;
  bool is_null() const;

  std::size_t size() const { return mNode.size(); }

  Node get_if_present(const std::string& key) const;
  Node get_null_if_present(const std::string& key) const;
  Node get_map_if_present(const std::string& key) const;
  Node get_sequence_if_present(const std::string& key) const;
  Node get_scalar_if_present(const std::string& key) const;

  // Options of a keyword, given either nested under it (the value of key is a map) or beside it in
  // this node (key has no value).  Returns the node holding the options, or an empty node if key is
  // absent; dooms if key has any other value.
  Node get_options_if_present(const std::string& key) const;

  template<typename T>
  bool get_if_present(const std::string& key, T& result) const
  {
    const Node value = get_if_present(key);
    if (value)
    {
      try
      {
        result = value.as<T>();
      }
      catch (const YAML::Exception &)
      {
        // A scalar that cannot be converted to T (e.g. "mesh_size: abc") must
        // not abort parsing.  Doom like every other input error so parse()
        // can collect and report it, leave result untouched, and report failure.
        value.report_conversion_failure(key);
        return false;
      }
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
  void report_conversion_failure(const std::string& key) const;
  Node get_type_if_present(const std::string& key, YAML::NodeType::value type) const;
  YAML::NodeType::value Type() const { return mNode.Type(); }

  YAML::Node mNode;
};

template<>
bool Node::as<bool>() const;

inline Node ConstNodeIterator::operator*() const
{
  return Node(*mIter);
}

// The message for a YAML parse error, with a hint for undefined aliases, which an unquoted value
// starting with '*' (such as a name pattern) produces.
std::string yaml_parse_error_message(const YAML::ParserException & e);

// Lines of text around a 0-based line and column (as in YAML::Mark), each prefixed with its 1-based
// line number, with '>' marking that line and a '^' under that column.  The context is clamped to
// the text.
std::string input_excerpt(const std::string & text, const int line, const int column, const int numContextLines = 3);

#else

class Node;

// No-YAML stub iterator.  Parsing is disabled in this configuration, so
// iteration is always over an empty range.  This exists only so that
// Node::begin()/end() return an iterator type matching the YAML build's
// ConstNodeIterator rather than Node*, keeping the two configurations'
// interfaces identical (e.g. for code using decltype(node.begin())).
class ConstNodeIterator {
public:
  Node operator*() const; // never invoked: begin() == end()
  ConstNodeIterator& operator++() { return *this; }
  ConstNodeIterator operator++(int) { return *this; }

  friend bool operator==(const ConstNodeIterator&, const ConstNodeIterator&) { return true; }
  friend bool operator!=(const ConstNodeIterator&, const ConstNodeIterator&) { return false; }
};

class Node {
public:
  operator bool() const { return false; }
  bool is_scalar() const { return false; }
  bool is_map() const { return false; }
  bool is_sequence() const { return false; }
  bool is_empty_scalar() const { return false; }
  bool is_null() const { return false; }
  std::size_t size() const { return 0; }

  Node get_if_present(const std::string& key) const { return Node(); }
  Node get_null_if_present(const std::string& key) const { return Node(); }
  Node get_map_if_present(const std::string& key) const { return Node(); }
  Node get_sequence_if_present(const std::string& key) const { return Node(); }
  Node get_scalar_if_present(const std::string& key) const { return Node(); }
  Node get_options_if_present(const std::string& key) const { return Node(); }

  template<typename T>
  bool get_if_present(const std::string& key, T& result) const { return false; }

  template<typename T>
  T as() const { return T(); }

  ConstNodeIterator begin() const { return ConstNodeIterator(); }
  ConstNodeIterator end() const { return ConstNodeIterator(); }

  std::string info() const { return std::string(); }
};

inline Node ConstNodeIterator::operator*() const { return Node(); }

#endif

void parse(Simulation & simulation);

}
}

#endif /* KRINO_KRINO_PARSER_AKRI_PARSER_HPP_ */
