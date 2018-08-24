// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once

#include "ROL_Ptr.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;


#include<typeinfo>
// For debug
#include <iostream>
#include <memory>

/*  Implements a unified ParameterList interface which conforms to that of
    ROL::ParameterList while using the Boost::property_tree implementation.

    This interface only implements a small subset of the available methods,
    however, at the time of creation these amount to all of the functionality
    that ROL has needed.
 */

namespace ROL {

namespace details {

using namespace std;

  // Try to get type of an object
  // FIXME: sometimes failing for std::string
  template <typename T> struct value_type
  { static std::string name() { std::string t = typeid(T).name();
      if (t == "i")
        return "int";
      else if (t == "d")
        return "double";
      else if (t == "b")
        return "bool";
      return t; }};

  template <> struct value_type<std::string>
  { static std::string name() { return "string"; }};
  template <int N> struct value_type<const char [N]>
  { static std::string name() { return "string"; }};
  template <typename T> std::string get_type(T& obj)
  { return value_type<T>::name(); }

class ParameterList {
private:

  // Keep the original ptree, but sublists, an iterator is sufficient
  pt::ptree tree_;
  pt::ptree::iterator root_;

  // For references (sublists) which must stay alive
  std::vector<std::shared_ptr<ParameterList>> refs_;

public:

  ParameterList() {
    tree_.put("ParameterList.<xmlattr>.name", "Unknown");
    root_ = tree_.begin();
  }

  ParameterList( pt::ptree tree ) : tree_(tree), root_(tree_.begin())
  {
  }

  ParameterList( pt::ptree::iterator root ) : root_(root)
  {
  }

  ParameterList( const string& name )  {
    tree_.put("ParameterList.<xmlattr>.name", name);
    root_ = tree_.begin();
  }

  virtual ~ParameterList() {
  }

  ParameterList& operator=(ParameterList& p)
  {
    // Copy tree over to this root
    root_->second = p.root_->second;

    // Copy name from "Name" Parameter in p
    std::string name = "Unknown";
    for (auto &q : p.root_->second)
      if (q.first == "Parameter" and
          q.second.get<std::string>("<xmlattr>.name") == "Name")
      {
        name = q.second.get<std::string>("<xmlattr>.value");
        break;
      }
    root_->second.put("<xmlattr>.name", name);
    return *this;
  }

  using ConstIterator = pt::ptree::const_iterator;

  ConstIterator begin() const {
    return root_->second.begin();
  }

  ConstIterator end() const {
    return root_->second.end();
  }

  std::string name(ConstIterator& it ) const {
    return it->second.get<std::string>("<xmlattr>.name");
  }

  template<class T>
  bool isType( const string& name ) const {
    T obj;
    const std::string my_type = get_type(obj);
    for (auto q : root_->second)
    {
      if (q.first == "Parameter" and
          q.second.get<std::string>("<xmlattr>.name") == name and
          q.second.get<std::string>("<xmlattr>.type") == my_type)
        return true;
    }
    return false;
  }


  template<class T>
  void set( const string& name, const T& value ) {
    // Look for existing parameter

    for (auto &q : root_->second)
    {
      if (q.first == "Parameter" and
          q.second.get<string>("<xmlattr>.name") == name)
      {
        q.second.put("<xmlattr>.value", value);
        // FIXME: check type
        return;
      }
    }

    // Make new parameter
    pt::ptree new_node;
    new_node.put("<xmlattr>.name", name);
    new_node.put("<xmlattr>.type", get_type(value));
    new_node.put("<xmlattr>.value", value);
    root_->second.add_child("Parameter", new_node);


  }

  template<class T>
  T get( const string& name ) const {
    for (auto &r : root_->second)
    {
      if (r.first == "Parameter" and
          r.second.get<std::string>("<xmlattr>.name") == name)
      {
        return r.second.get<T>("<xmlattr>.value");
      }
    }
    exit(-1);
    return T();
  }

  std::string get( const string& name, const string& default_value) {
    return get<string>(name, default_value);
  }

  template<class T>
  T get( const string& name, const T& default_value ) {
    for (auto r : root_->second)
    {
      if (r.first == "Parameter" and
          r.second.get<std::string>("<xmlattr>.name") == name)
      {
        return r.second.get<T>("<xmlattr>.value");
      }
    }
    set(name, default_value);

    return default_value;
  }

  void print()
  {
      print(tree_);
  }

  static void print(pt::ptree& r, std::string indent="")
  {
    for (auto q : r)
    {
      pt::ptree& sub = q.second;
      if (sub.size() == 0)
      {
        std::cout << indent << "[" << q.first << "] = \"";
        std::cout << r.get<std::string>(q.first) << "\"\n";
      }
      else
      {
        std::cout << indent << "[" << q.first << "]\n";
        print(sub, indent + "  ");
      }
    }
  }

  ParameterList& sublist(const string& name) {
    for (pt::ptree::iterator r = root_->second.begin(); r != root_->second.end(); ++r)
    {
      if (r->first == "ParameterList")
      {
        const std::string xml_name = r->second.get<std::string>("<xmlattr>.name");
        if (xml_name == name)
        {
          refs_.push_back(std::shared_ptr<ParameterList>(new ParameterList(r)));
          return *refs_.back();
        }
      }
    }

    // Create node and retry
    pt::ptree tr;
    tr.put("<xmlattr>.name", name);
    root_->second.add_child("ParameterList", tr);
    return sublist(name);
  }


  bool isSublist(const string& name) const
  {
    for (auto q : root_->second)
    {
      if (q.first == "ParameterList" and
          q.second.get<string>("<xmlattr>.name") == name)
        return true;
    }
    return false;
  }

  bool isParameter(const string& name) const
  {
    for (auto q : root_->second)
    {
      if (q.first == "Parameter" and
          q.second.get<string>("<xmlattr>.name") == name)
        return true;
    }
    return false;
  }

  pt::ptree& tree()
  { return tree_; }

  pt::ptree::iterator& root()
  { return root_; }

  //  friend void readParametersFromXml( const string&, ParameterList& parlist );

  void validateParametersAndSetDefaults( const ParameterList& parlist ) {
    throw Exception::NotImplemented("ParameterList validation has not been "
          "implemented for the  boost::property_tree implementation.");
  }

};

} // namespace details

  using ParameterList = details::ParameterList;

  template <class T>
  inline std::vector<T> getArrayFromStringParameter(const ParameterList& parlist,
                                                    const std::string& name)
  {
    std::string p = parlist.get<std::string>(name);

    std::vector<std::string> p_split;
    boost::split(p_split, p, boost::is_any_of("{,}"));

    std::vector<T> result;
    for (auto &q : p_split)
    {
      boost::trim(q);
      if(q.size() > 0)
        result.push_back(boost::lexical_cast<T>(q));
    }

    return result;
  }

  inline ROL::Ptr<ParameterList> getParametersFromXmlFile( const std::string& filename )
  {
    pt::ptree tr;
    boost::property_tree::read_xml(filename, tr);
    auto list = ROL::makePtr<ParameterList>(tr);
    list->root() = list->tree().begin();
    return list;
  }

  inline void readParametersFromXml( const std::string& filename,
                                     ParameterList& parlist ) {
    boost::property_tree::read_xml(filename, parlist.tree());
    parlist.root() = parlist.tree().begin();
  }

  inline void updateParametersFromXmlFile( const std::string& infile, ParameterList& inlist )
  {
    // FIXME: do something
  }

  inline void writeParameterListToXmlFile( ParameterList& parlist,
                                           const std::string& filename ) {
    boost::property_tree::write_xml(filename, parlist.tree());
  }



} // namespace ROL
