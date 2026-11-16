/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <cstdint> // for int64_t
#include <string>  // for string
#include <unordered_map>
#include <variant>
#include <vector>

#include "exodusII.h" // for ex_close, etc

/*********** entity_attr.h -- Property defined in an Exodus or Database Attribute **************/

class ExodusProperty
{
public:
  enum BasicType { INVALID = -1, REAL, INTEGER, STRING, VEC_INTEGER, VEC_DOUBLE };

  ExodusProperty(std::string name, int64_t value);
  ExodusProperty(std::string name, int value);
  ExodusProperty(std::string name, double value);
  ExodusProperty(std::string name, const std::string &value);
  ExodusProperty(std::string name, const char *value);
  ExodusProperty(std::string name, const std::vector<int> &value);
  ExodusProperty(std::string name, const std::vector<double> &value);

  std::string         get_string() const;
  int64_t             get_int() const;
  double              get_real() const;
  std::vector<double> get_vec_double() const;
  std::vector<int>    get_vec_int() const;

  /** Tells whether the attribute has a valid type (currently REAL, INTEGER, or STRING)
   *
   *  \returns True if the attribute type is valid.
   */
  bool is_valid() const { return type_ != INVALID; }

  /** Tells whether the attribute has an invalid type (currently not one of REAL, INTEGER, or
   * STRING)
   *
   *  \returns True if the attribute type is invalid.
   */
  bool is_invalid() const { return type_ == INVALID; }

  /** \brief Get the attribute name.
   *
   *  \returns The attribute name.
   */
  std::string get_name() const { return name_; }

  /** \brief Get the attribute type.
   *
   *  \returns The attribute type.
   */
  BasicType get_type() const { return type_; }

  bool operator!=(const ExodusProperty &rhs) const;
  bool operator==(const ExodusProperty &rhs) const;

private:
  std::string name_{};
  BasicType   type_{INVALID};

  bool get_value(int64_t *value) const;
  bool get_value(double *value) const;
  bool get_value(std::string *value) const;
  bool get_value(void *&value) const;
  bool get_value(std::vector<double> *value) const;
  bool get_value(std::vector<int> *value) const;

  /// The actual value of the attribute.  Use 'type_' to
  /// discriminate the actual type of the attribute.
  std::variant<std::string, double, int64_t, std::vector<double>, std::vector<int>> data_;
};

using PropMapType = std::unordered_map<std::string, ExodusProperty>;
using ValuePair   = PropMapType::value_type;

/** \brief A collection of Property objects  */
using NameList = std::vector<std::string>;

class ExodusPropertyManager
{
public:
  ExodusPropertyManager() = default;
  ExodusPropertyManager(const ExodusPropertyManager &from)
      : m_properties(from.m_properties)
  { /* Do not make this `=default` since that breaks the thread-safe build */ }
  ExodusPropertyManager &operator=(const ExodusPropertyManager &from) = delete;

  // Add the specified attribute to the list.
  void add(const ExodusProperty &new_prop);

  // Assumes: Attribute 'name' must exist.
  void erase(const std::string &attribute_name);

  // Checks if a attribute with 'attribute_name' exists in the database.
  bool exists(const std::string &attribute_name) const;

  ExodusProperty get(const std::string &attribute_name) const;
  double         get_optional(const std::string &attribute_name, double optional_value) const;
  int64_t        get_optional(const std::string &attribute_name, int64_t optional_value) const;
  int            get_optional(const std::string &attribute_name, int optional_value) const;
  std::string    get_optional(const std::string &attribute_name,
                              const std::string &optional_value) const;

  // Returns the names of all attributes
  int      describe(NameList *names) const;
  NameList describe() const;

  size_t count() const;

  void copy(const ExodusPropertyManager &from);

private:
  PropMapType m_properties{};
};

bool read_exodus_entity_properties(int mesh_exoid, ex_entity_type type, int id,
                                   ExodusPropertyManager &attributes);
bool write_exodus_entity_properties(int exoid, ex_entity_type type, int id,
                                    const ExodusPropertyManager &attributes);
