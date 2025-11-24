// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_PARAMETERLIST_HPP
#define ROL_PARAMETERLIST_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <type_traits>
#include <sstream>
#include <iostream>

#include "ROL_Stream.hpp"

namespace ROL {

#if __cplusplus >= 201703L
template<class...B>
inline constexpr bool disjunction_v = std::disjunction<B...>::value;
#else
template<class...> struct disjunction : std::false_type { };
template<class B1> struct disjunction<B1> : B1 { };
template<class B1, class... Bn>
struct disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>>  {};
template<class...B>
constexpr bool disjunction_v = disjunction<B...>::value;
#endif

namespace detail {

template<typename value_type>
std::ostream& operator << ( std::ostream& os, const std::vector<value_type>& v ) {
  os << "[ ";
  for( std::size_t k=0; k<v.size()-1; ++k ) os << v[k] << ", ";
  os << v.back() << " ]" << std::endl;
  return os;
}


template<typename ValueType>
class PList {
public:

  using key_type   = std::string;
  using value_type = ValueType;

  inline void set( key_type   key,
                   value_type value ) {
    values_[key] = value;
  }

  inline value_type get( key_type   key,
                         value_type value ) {
    if( values_.find(key) == values_.end() ) values_[key] = value;
    return values_.at(key);
  }

  inline value_type get( key_type key ) const {
    return values_.at(key);
  }

  virtual int get_level() const = 0;

  friend std::ostream& operator << ( std::ostream& os,
                                     const PList& plist ) {
    for( auto it=plist.values_.begin();
              it != plist.values_.end(); ++it )
      os << std::string(2*(plist.get_level()),' ') <<
            it->first << " : " << it->second << std::endl;
    return os;
  }

  const std::map<key_type,value_type>& get_values() const { return values_; }

private:
  std::map<key_type,value_type> values_;

}; // class detail::PList



template<typename...ValueTypes>
class ParameterList : public PList<ValueTypes>... {
public:

  using key_type = std::string;

  // True if type can be converted to any of the ParameterList value types
  template<typename value_type>
  static constexpr bool is_valid_type_v = disjunction_v<std::is_convertible<value_type,ValueTypes>...>;

  ParameterList( int level = 0 ) : level_(level) {}

  template<typename value_type>
  inline void set( key_type   key,
                   value_type value ) {
    static_assert( is_valid_type_v<value_type>, "" );
    static_cast<PList<value_type>&>(*this).set( key, value );
  }

  inline void set( key_type    key,
                   const char* value ) {
    set(key,key_type(value));
  }

  inline std::string get( key_type    key,
                          const char* value ) {
    return get(key,std::string(value));
  }

  template<typename value_type>
  inline
  std::enable_if_t<is_valid_type_v<value_type>,value_type>
  get( key_type key ) const {
    //static_assert( is_valid_type_v<value_type>, "" );
    return PList<value_type>::get(key);
  }

  template<typename value_type>
  inline value_type get( key_type key, value_type value ) {
    //static_assert( is_valid_type_v<value_type>, "" );
    return PList<value_type>::get(key,value);
  }


private:
  int level_ = 0;

}; // class detail::ParameterList

template<typename...ValueTypes>
inline void display( std::ostream& os,
                     const ParameterList<ValueTypes...>& parlist ) {
  using expander = int[];
  (void)expander{0, (void( os << static_cast<const PList<ValueTypes>&>(parlist) ),0) ... };
}

} // namespace detail


class ParameterList : public detail::ParameterList<bool,
                                                   int,
                                                   double,
                                                   std::string,
                                                   std::vector<int>,
                                                   std::vector<double>,
                                                   std::vector<std::string>> {
public:

  using key_type = std::string;

  // Iterator types for ParameterListConverters compatibility
  class ConstIterator {
  private:
    std::vector<std::string> keys_;
    std::size_t index_;

  public:
    ConstIterator(const std::vector<std::string>& keys, std::size_t index)
      : keys_(keys), index_(index) {}

    ConstIterator& operator++() { ++index_; return *this; }
    bool operator!=(const ConstIterator& other) const { return index_ != other.index_; }
    bool operator==(const ConstIterator& other) const { return index_ == other.index_; }
    std::size_t getIndex() const { return index_; }
  };

  ParameterList( int level = 0 ) : level_(level) {}

  // Constructor with name (for Teuchos compatibility) - name is ignored in simple implementation
  ParameterList( const std::string& name, int level = 0 ) : level_(level) {}

  // Deep copy constructor
  ParameterList( const ParameterList& other ) : level_(other.level_) {
    // Copy all parameter values from each base class
    copyValues<bool>(other);
    copyValues<int>(other);
    copyValues<double>(other);
    copyValues<std::string>(other);
    copyValues<std::vector<int>>(other);
    copyValues<std::vector<double>>(other);
    copyValues<std::vector<std::string>>(other);

    // Deep copy all sublists
    for (const auto& sublist_pair : other.sublists_) {
      sublists_[sublist_pair.first] = std::make_shared<ParameterList>(*sublist_pair.second);
    }
  }

  ParameterList& sublist( key_type key,
                          bool     mustAlreadyExist=false );

//  template<typename value_type>
//  inline value_type get( key_type key ) {
//    static_assert( is_valid_type_v<value_type>, "" );
//    return PList<value_type>::get( key );
//  }
//
//  template<typename value_type>
//  inline value_type get( key_type key, value_type default_value ) {
//    static_assert( is_valid_type_v<value_type>, "" );
//    return PList<value_type>::get( key,default_value );
//  }

  inline bool isSublist( key_type key ) {
    return sublists_.find(key) != sublists_.end();
  }

  // Iterator methods for ParameterListConverters compatibility
  ConstIterator begin() const {
    return ConstIterator(getAllKeys(), 0);
  }

  ConstIterator end() const {
    std::vector<std::string> keys = getAllKeys();
    return ConstIterator(keys, keys.size());
  }

  // Get parameter name from iterator
  std::string name(const ConstIterator& it) const {
    std::vector<std::string> keys = getAllKeys();
    return keys[it.getIndex()];
  }

  // Type checking methods
  template<typename T>
  bool isType(const std::string& key) const {
    return hasKey<T>(key);
  }

  // Check if parameter exists (any supported type)
  bool isParameter(const std::string& key) const {
    return hasKey<bool>(key) ||
           hasKey<int>(key) ||
           hasKey<double>(key) ||
           hasKey<std::string>(key) ||
           hasKey<std::vector<int>>(key) ||
           hasKey<std::vector<double>>(key) ||
           hasKey<std::vector<std::string>>(key);
  }

  int get_level() const override { return level_; }

  friend std::ostream& operator << ( std::ostream&  os,
                                     const ParameterList& parlist ) {
    detail::display( os , parlist );
    os << std::endl;
    for( auto it=parlist.sublists_.begin(); it != parlist.sublists_.end(); ++it )
      os << std::string(2*parlist.get_level(),' ') <<
           it->first << ":\n" << *(it->second);
    return os;
  }

private:

  // Helper method to get all parameter keys across all base classes
  std::vector<std::string> getAllKeys() const {
    std::vector<std::string> keys;
    collectKeys<bool>(keys);
    collectKeys<int>(keys);
    collectKeys<double>(keys);
    collectKeys<std::string>(keys);
    collectKeys<std::vector<int>>(keys);
    collectKeys<std::vector<double>>(keys);
    collectKeys<std::vector<std::string>>(keys);
    return keys;
  }

  // Helper template to collect keys for specific type
  template<typename T>
  void collectKeys(std::vector<std::string>& keys) const {
    const auto& values = static_cast<const detail::PList<T>&>(*this).get_values();
    for (const auto& pair : values) {
      keys.push_back(pair.first);
    }
  }

  // Helper template to check if key exists for specific type
  template<typename T>
  bool hasKey(const std::string& key) const {
    const auto& values = static_cast<const detail::PList<T>&>(*this).get_values();
    return values.find(key) != values.end();
  }

  // Helper template to copy values for specific type from another ParameterList
  template<typename T>
  void copyValues(const ParameterList& other) {
    const auto& other_values = static_cast<const detail::PList<T>&>(other).get_values();
    for (const auto& pair : other_values) {
      static_cast<detail::PList<T>&>(*this).set(pair.first, pair.second);
    }
  }

  std::map<key_type,std::shared_ptr<ParameterList>> sublists_;
  int level_ = 0;

}; // ParameterList

// Array parsing utilities similar to Teuchos::getArrayFromStringParameter
template<typename value_type>
std::vector<value_type> fromStringToArray(const std::string& arrayStr) {
  std::vector<value_type> result;

  // Remove braces and spaces
  std::string cleaned = arrayStr;
  size_t start = cleaned.find('{');
  size_t end = cleaned.find('}');
  if (start != std::string::npos && end != std::string::npos && end > start) {
    cleaned = cleaned.substr(start + 1, end - start - 1);
  }

  // Parse comma-separated values
  std::istringstream iss(cleaned);
  std::string token;

  while (std::getline(iss, token, ',')) {
    // Trim whitespace
    size_t first = token.find_first_not_of(" \t");
    if (first != std::string::npos) {
      size_t last = token.find_last_not_of(" \t");
      token = token.substr(first, last - first + 1);

      // Convert to appropriate type
      if (!token.empty()) {
        if constexpr (std::is_same_v<value_type, int>) {
          result.push_back(std::stoi(token));
        } else if constexpr (std::is_same_v<value_type, double>) {
          result.push_back(std::stod(token));
        } else if constexpr (std::is_same_v<value_type, std::string>) {
          result.push_back(token);
        }
      }
    }
  }

  return result;
}

template<typename value_type>
std::vector<value_type> getArrayFromStringParameter(const ParameterList& parlist, const std::string& key) {
  // First try to get as vector directly
  try {
    return parlist.get<std::vector<value_type>>(key);
  } catch (...) {
    // If that fails, try to get as string and parse it
    try {
      std::string arrayStr = parlist.get<std::string>(key);
      return fromStringToArray<value_type>(arrayStr);
    } catch (...) {
      throw std::runtime_error("Parameter '" + key + "' not found or not convertible to array");
    }
  }
}


std::ostream& operator << ( std::ostream& os, const ParameterList& parlist );

// Forward declaration for XML reader
class XMLParameterListReader;

/// \brief Create ParameterList from XML file using pugixml
ROL::Ptr<ParameterList> getParametersFromXmlFile(const std::string& filename);

/// \brief Write ParameterList to XML file using pugixml
void writeParameterListToXmlFile(const ParameterList& params, const std::string& filename);



} // namespace ROL

#include "ROL_ParameterList.cpp"
#include "ROL_XMLReader.hpp"

namespace ROL{
// Global function for compatibility
inline ROL::Ptr<ParameterList> getParametersFromXmlFile(const std::string& filename) {
    return XMLReader::readParametersFromFile(filename);
}

// Global function for writing parameters to XML
inline void writeParameterListToXmlFile(const ParameterList& params, const std::string& filename) {
    XMLWriter::writeParametersToFile(params, filename);
}

}
#endif // ROL_PARAMETERLIST_HPP
