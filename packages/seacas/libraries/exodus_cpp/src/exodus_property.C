/*
 * Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "exodus_property.h"
#include "vector_data.h"
#include <assert.h>
#include <cstdio>  // for stderr
#include <cstdlib> // for exit
#include <fmt/ostream.h>
#include <ostream>
#include <sstream>

namespace {
  std::string type_string(ExodusProperty::BasicType type)
  {
    switch (type) {
    case ExodusProperty::INVALID: return {"invalid"};
    case ExodusProperty::REAL: return {"real"};
    case ExodusProperty::INTEGER: return {"integer"};
    case ExodusProperty::STRING: return {"string"};
    case ExodusProperty::VEC_INTEGER: return {"vector<int>"};
    case ExodusProperty::VEC_DOUBLE: return {"vector<double>"};
    }
    return {"internal error"};
  }

  void exit_with_error(const std::ostringstream &errmsg) { throw std::runtime_error(errmsg.str()); }

  void error_message(const ExodusProperty &property, const std::string &requested_type)
  {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: For property named '{}', code requested value of type '{}', but property "
               "type is '{}'. Types must match\n",
               property.get_name(), requested_type, type_string(property.get_type()));
    exit_with_error(errmsg);
  }
} // namespace

/** \brief Create an INTEGER type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, int value)
    : name_(std::move(name)), type_(INTEGER), data_(static_cast<int64_t>(value))
{
}

/** \brief Create an INTEGER type property using an int64_t variable.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, int64_t value)
    : name_(std::move(name)), type_(INTEGER), data_(value)
{
}

/** \brief Create a REAL type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, double value)
    : name_(std::move(name)), type_(REAL), data_(value)
{
}

/** \brief Create a STRING type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, const std::string &value)
    : name_(std::move(name)), type_(STRING), data_(value)
{
}

/** \brief Create a VEC_INTEGER type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, const std::vector<int> &value)
    : name_(std::move(name)), type_(VEC_INTEGER), data_(value)
{
}

/** \brief Create a VEC_DOUBLE type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, const std::vector<double> &value)
    : name_(std::move(name)), type_(VEC_DOUBLE), data_(value)
{
}

/** \brief Create a STRING type property from const char* argument.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 */
ExodusProperty::ExodusProperty(std::string name, const char *value)
    : name_(std::move(name)), type_(STRING), data_(std::string(value))
{
}

bool ExodusProperty::operator==(const ExodusProperty &rhs) const
{
  if (this->name_ != rhs.name_) {
    return false;
  }

  if (this->type_ != rhs.type_) {
    return false;
  }

  switch (this->type_) {
  case INVALID: break;
  case REAL:
    double r_lhs, r_rhs;
    this->get_value(&r_lhs);
    rhs.get_value(&r_rhs);
    if (r_lhs != r_rhs) {
      return false;
    }
    break;
  case INTEGER:
    int64_t i_lhs, i_rhs;
    this->get_value(&i_lhs);
    rhs.get_value(&i_rhs);
    if (i_lhs != i_rhs) {
      return false;
    }
    break;
  case VEC_DOUBLE: {
    auto rh = rhs.get_vec_double();
    auto lh = get_vec_double();
    if (lh != rh) {
      return false;
    }
  } break;
  case VEC_INTEGER: {
    auto rh = rhs.get_vec_int();
    auto lh = get_vec_int();
    if (lh != rh) {
      return false;
    }
  } break;
  case STRING:
    std::string s_lhs, s_rhs;
    this->get_value(&s_lhs);
    rhs.get_value(&s_rhs);
    if (s_lhs != s_rhs) {
      return false;
    }
    break;
  }
  return true;
}

bool ExodusProperty::operator!=(const ExodusProperty &rhs) const { return !(*this == rhs); }

/** \brief Get the property value if it is of type STRING.
 *
 *  \returns The STRING-type property value
 */
std::string ExodusProperty::get_string() const
{
  std::string value;
  bool        valid = get_value(&value);
  if (!valid) {
    error_message(*this, "string");
  }
  return value;
}

/** \brief Get the property value if it is of type VEC_DOUBLE.
 *
 *  \returns The VEC_DOUBLE-type property value
 */
std::vector<double> ExodusProperty::get_vec_double() const
{
  std::vector<double> value;
  bool                valid = get_value(&value);
  if (!valid) {
    error_message(*this, "vector<double>");
  }
  return value;
}

/** \brief Get the property value if it is of type VEC_INT.
 *
 *  \returns The VEC_INT-type property value
 */
std::vector<int> ExodusProperty::get_vec_int() const
{
  std::vector<int> value;
  bool             valid = get_value(&value);
  if (!valid) {
    error_message(*this, "vector<int>");
  }
  return value;
}

/** \brief Get the property value if it is of type INTEGER.
 *
 *  \returns The INTEGER-type property value
 */
int64_t ExodusProperty::get_int() const
{
  int64_t value;
  bool    valid = get_value(&value);
  if (!valid) {
    error_message(*this, "int");
  }
  return value;
}

/** \brief Get the property value if it is of type REAL.
 *
 *  \returns The REAL-type property value.
 */
double ExodusProperty::get_real() const
{
  double value;
  bool   valid = get_value(&value);
  if (!valid) {
    error_message(*this, "real");
  }
  return value;
}

bool ExodusProperty::get_value(int64_t *value) const
{
  bool valid_request = type_ == INTEGER;

  assert(std::holds_alternative<int64_t>(data_));
  *value = std::get<int64_t>(data_);

  return valid_request;
}

bool ExodusProperty::get_value(double *value) const
{
  bool valid_request = type_ == REAL;

  assert(std::holds_alternative<double>(data_));
  *value = std::get<double>(data_);

  return valid_request;
}

bool ExodusProperty::get_value(std::string *value) const
{
  bool valid_request = type_ == STRING;

  assert(std::holds_alternative<std::string>(data_));
  *value = std::get<std::string>(data_);

  return valid_request;
}

bool ExodusProperty::get_value(std::vector<int> *value) const
{
  bool valid_request = type_ == VEC_INTEGER;

  assert(std::holds_alternative<std::vector<int>>(data_));
  auto ivec = std::get<std::vector<int>>(data_);
  std::copy(ivec.begin(), ivec.end(), std::back_inserter(*value));

  return valid_request;
}

bool ExodusProperty::get_value(std::vector<double> *value) const
{
  bool valid_request = type_ == VEC_DOUBLE;

  assert(std::holds_alternative<std::vector<double>>(data_));
  auto dvec = std::get<std::vector<double>>(data_);
  std::copy(dvec.begin(), dvec.end(), std::back_inserter(*value));

  return valid_request;
}

/** \brief Add a property to the property manager.
 *
 *  \param[in] new_attr The property to add.
 */
void ExodusPropertyManager::add(const ExodusProperty &new_property)
{
  auto iter = m_properties.find(new_property.get_name());
  if (iter != m_properties.end()) {
    m_properties.erase(iter);
  }
  m_properties.emplace(new_property.get_name(), new_property);
}

/** \brief Checks if a property exists in the database.
 *
 *  \param[in] property_name The property to check
 *  \returns True if the property exists, false otherwise.
 */
bool ExodusPropertyManager::exists(const std::string &property_name) const
{
  return m_properties.find(property_name) != m_properties.end();
}

/** \brief Get a property object from the property manager.
 *
 *  \param[in] property_name The name of the property to get.
 *  \returns The property object.
 */
ExodusProperty ExodusPropertyManager::get(const std::string &property_name) const
{
  auto iter = m_properties.find(property_name);
  if (iter == m_properties.end()) {
    std::ostringstream errmsg;
    fmt::print(errmsg, "ERROR: Could not find property '{}'\n", property_name);
    exit_with_error(errmsg);
  }
  return (*iter).second;
}

/** \brief Get an optional property object from the property manager.
 *
 *  \param[in] property_name The name of the property to get.
 *  \param[in] optional_value The value to return if the property does not exist.
 *  \returns The property object.
 */
int ExodusPropertyManager::get_optional(const std::string &property_name, int optional_value) const
{
  auto iter = m_properties.find(property_name);
  if (iter == m_properties.end()) {
    return optional_value;
  }
  return (*iter).second.get_int();
}

int64_t ExodusPropertyManager::get_optional(const std::string &property_name,
                                            int64_t            optional_value) const
{
  auto iter = m_properties.find(property_name);
  if (iter == m_properties.end()) {
    return optional_value;
  }
  return (*iter).second.get_int();
}

double ExodusPropertyManager::get_optional(const std::string &property_name,
                                           double             optional_value) const
{
  auto iter = m_properties.find(property_name);
  if (iter == m_properties.end()) {
    return optional_value;
  }
  return (*iter).second.get_real();
}

std::string ExodusPropertyManager::get_optional(const std::string &property_name,
                                                const std::string &optional_value) const
{
  auto iter = m_properties.find(property_name);
  if (iter == m_properties.end()) {
    return optional_value;
  }
  return (*iter).second.get_string();
}

/** \brief Remove a property from the property manager.
 *
 *  Assumes that the property with the given name already exists in the property manager.
 *
 *  \param[in] property_name The name of the property to remove.
 *
 */
void ExodusPropertyManager::erase(const std::string &property_name)
{
  auto iter = m_properties.find(property_name);
  if (iter != m_properties.end()) {
    m_properties.erase(iter);
  }
}

/** \brief Get the names of all properties in the property manager
 *
 *  \returns All the property names in the property manager.
 */
NameList ExodusPropertyManager::describe() const
{
  NameList names;
  describe(&names);
  return names;
}

/** \brief Get the names of all properties in the property manager
 *
 *  \param[out] names All the property names in the property manager.
 *  \returns The number of properties extracted from the property manager.
 */
int ExodusPropertyManager::describe(NameList *names) const
{
  int                         the_count = 0;
  PropMapType::const_iterator I;
  for (I = m_properties.begin(); I != m_properties.end(); ++I) {
    names->push_back((*I).first);
    the_count++;
  }
  return the_count;
}

/** Get the number of properties in the property manager
 *
 *  \returns The number of properties in the property manager.
 */
size_t ExodusPropertyManager::count() const { return m_properties.size(); }

/** \brief copies all the properties from an property manager
 *
 *  \param[in] from The property manager to duplicate.
 */
void ExodusPropertyManager::copy(const ExodusPropertyManager &from)
{
  NameList properties = from.describe();

  for (const auto &property_name : properties) {
    auto prop = from.get(property_name);
    add(prop);
  }
}

namespace {
  void exit_on_exodus_error(int exoid, int error, const char *function_name)
  {
    if (error < 0) {
      fmt::print(stderr, "ERROR returned from {}!\n", function_name);
      ex_close(exoid);
      exit(1);
    }
  }

  // Memory allocated in `ex_get_attributes`, this makes deletion cleaner...
  class EX_attribute : public ex_attribute
  {
  public:
    EX_attribute() { values = nullptr; }
    ~EX_attribute() { free(values); }
  };
} // namespace

bool read_exodus_entity_properties(int exoid, ex_entity_type type, int id,
                                   ExodusPropertyManager &attributes)
{
  int att_count = ex_get_attribute_count(exoid, type, id);

  if (att_count > 0) {
    std::vector<EX_attribute> attr(att_count);
    exit_on_exodus_error(exoid, ex_get_attribute_param(exoid, type, id, suplib_cpp::Data(attr)),
                         "ex_get_attribute_param");
    exit_on_exodus_error(exoid, ex_get_attributes(exoid, att_count, suplib_cpp::Data(attr)),
                         "ex_get_attributes");

    // Create a property on `entity` for each `attribute`
    for (const auto &att : attr) {
      if (att.value_count == 0) {
        // Just an attribute name.  Give it an empty value...
        attributes.add(ExodusProperty(att.name, ""));
        continue;
      }
      assert(att.values != nullptr);

      switch (att.type) {
      case EX_INTEGER: {
        const auto *idata = static_cast<int *>(att.values);
        if (att.value_count == 1) {
          attributes.add(ExodusProperty(att.name, *idata));
        }
        else {
          std::vector<int> tmp(att.value_count);
          std::copy(idata, idata + att.value_count, tmp.begin());
          attributes.add(ExodusProperty(att.name, tmp));
        }
      } break;
      case EX_DOUBLE: {
        const auto *ddata = static_cast<double *>(att.values);
        if (att.value_count == 1) {
          attributes.add(ExodusProperty(att.name, *ddata));
        }
        else {
          std::vector<double> tmp(att.value_count);
          std::copy(ddata, ddata + att.value_count, tmp.begin());
          attributes.add(ExodusProperty(att.name, tmp));
        }
      } break;
      case EX_CHAR: {
        const auto *cdata = static_cast<char *>(att.values);
        attributes.add(ExodusProperty(att.name, cdata));
      } break;
      }
    }
  }
  return true;
}

bool write_exodus_entity_properties(int exoid, ex_entity_type type, int id,
                                    const ExodusPropertyManager &attributes)
{
  NameList properties = attributes.describe();

  double  rval = 0.0;
  int64_t ival = 0;
  for (const auto &property_name : properties) {
    auto prop = attributes.get(property_name);

    switch (prop.get_type()) {
    case ExodusProperty::BasicType::REAL:
      rval = prop.get_real();
      exit_on_exodus_error(
          exoid, ex_put_double_attribute(exoid, type, id, property_name.c_str(), 1, &rval),
          "ex_put_double_attribute");
      break;
    case ExodusProperty::BasicType::INTEGER:
      ival = prop.get_int();
      exit_on_exodus_error(
          exoid, ex_put_integer_attribute(exoid, type, id, property_name.c_str(), 1, &ival),
          "ex_put_integer_attribute");
      break;
    case ExodusProperty::BasicType::STRING:
      exit_on_exodus_error(
          exoid,
          ex_put_text_attribute(exoid, type, id, property_name.c_str(), prop.get_string().c_str()),
          "ex_put_text_attribute");
      break;
    case ExodusProperty::BasicType::VEC_INTEGER:
      exit_on_exodus_error(exoid,
                           ex_put_integer_attribute(exoid, type, id, property_name.c_str(),
                                                    prop.get_vec_int().size(),
                                                    suplib_cpp::Data(prop.get_vec_int())),
                           "ex_put_integer_attribute");
      break;
    case ExodusProperty::BasicType::VEC_DOUBLE:
      exit_on_exodus_error(exoid,
                           ex_put_double_attribute(exoid, type, id, property_name.c_str(),
                                                   prop.get_vec_double().size(),
                                                   suplib_cpp::Data(prop.get_vec_double())),
                           "ex_put_double_attribute");
      break;
    default:; // Do nothing
    }
  }

  return true;
}
