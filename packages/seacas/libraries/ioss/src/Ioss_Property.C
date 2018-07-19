// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include <Ioss_GroupingEntity.h>
#include <Ioss_Property.h>
#include <Ioss_Utils.h>
#include <cstddef>
#include <ostream>
#include <string>

namespace {
  std::string type_string(Ioss::Property::BasicType type)
  {
    switch (type) {
    case Ioss::Property::INVALID: return std::string("invalid");
    case Ioss::Property::REAL: return std::string("real");
    case Ioss::Property::INTEGER: return std::string("integer");
    case Ioss::Property::POINTER: return std::string("pointer");
    case Ioss::Property::STRING: return std::string("string");
    default: return std::string("internal error");
    }
  }

  void error_message(const Ioss::Property &property, const std::string &requested_type)
  {
    std::ostringstream errmsg;
    errmsg << "ERROR: For property named '" << property.get_name()
           << "', code requested value of type '" << requested_type << "', but property type is '"
           << type_string(property.get_type()) << "'. Types must match\n";
    IOSS_ERROR(errmsg);
  }
} // namespace

/** \brief Create an explicit, empty property having no name, INVALID type
 */
Ioss::Property::Property() : name_(""), type_(INVALID), isImplicit_(false) { data_.pval = nullptr; }

/** \brief Create a property with a specified type.
 *
 *  \param[in] name The property name.
 *  \param[in] type The property type.
 *  \param[in] data The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, const BasicType type, void *data, bool is_it_implicit)
    : name_(std::move(name)), type_(type), isImplicit_(is_it_implicit)
{
  data_.pval = data;
}

/** \brief Create an INTEGER type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, int value, bool is_it_implicit)
    : name_(std::move(name)), type_(INTEGER), isImplicit_(is_it_implicit)
{
  data_.ival = value;
}

/** \brief Create an INTEGER type property using an int64_t variable.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, int64_t value, bool is_it_implicit)
    : name_(std::move(name)), type_(INTEGER), isImplicit_(is_it_implicit)
{
  data_.ival = value;
}

/** \brief Create a REAL type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, double value, bool is_it_implicit)
    : name_(std::move(name)), type_(REAL), isImplicit_(is_it_implicit)
{
  data_.rval = value;
}

/** \brief Create a STRING type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, const std::string &value, bool is_it_implicit)
    : name_(std::move(name)), type_(STRING), isImplicit_(is_it_implicit)
{
  data_.sval = new std::string(value);
}

/** \brief Create a POINTER type property.
 *
 *  \param[in] name The property name.
 *  \param[in] value The property value.
 *  \param[in] is_it_implicit True if the property is calculated, rather than stored directly
 */
Ioss::Property::Property(std::string name, void *value, bool is_it_implicit)
    : name_(std::move(name)), type_(POINTER), isImplicit_(is_it_implicit)
{
  data_.pval = value;
}

/** \brief Set implicit property with a specified type.
 *
 *  \param[in] ge The property value.
 *  \param[in] name The property name.
 *  \param[in] type The property type.
 */
Ioss::Property::Property(const Ioss::GroupingEntity *ge, std::string name, const BasicType type)
    : name_(std::move(name)), type_(type), isImplicit_(true)
{
  data_.ge = ge;
}

/** \brief Copy constructor.
 *
 *  \param[in] from The Ioss::Property to copy
 */
Ioss::Property::Property(const Ioss::Property &from)
{
  name_       = from.name_;
  type_       = from.type_;
  isImplicit_ = from.isImplicit_;
  if (!isImplicit_ && type_ == STRING) {
    data_.sval = new std::string(*(from.data_.sval));
  }
  else {
    data_ = from.data_;
  }
}

Ioss::Property::~Property()
{
  if (!isImplicit_ && type_ == STRING) {
    delete data_.sval;
  }
}

Ioss::Property &Ioss::Property::operator=(Ioss::Property rhs)
{
  std::swap(*this, rhs);
  return *this;
}

/** \brief Get the property value if it is of type STRING.
 *
 *  \returns The STRING-type property value
 */
std::string Ioss::Property::get_string() const
{
  std::string value;
  bool        valid = get_value(&value);
  if (!valid) {
    error_message(*this, "string");
  }
  return value;
}

/** \brief Get the property value if it is of type INTEGER.
 *
 *  \returns The INTEGER-type property value
 */
int64_t Ioss::Property::get_int() const
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
double Ioss::Property::get_real() const
{
  double value;
  bool   valid = get_value(&value);
  if (!valid) {
    error_message(*this, "real");
  }
  return value;
}

/** \brief Get the property value if it is of type POINTER.
 *
 *  \returns The POINTER-type property value.
 */
void *Ioss::Property::get_pointer() const
{
  void *value = nullptr;
  bool  valid = get_value(value);
  if (!valid) {
    error_message(*this, "pointer");
  }
  return value;
}

bool Ioss::Property::get_value(int64_t *value) const
{
  bool valid_request = false;
  if (type_ == INTEGER) {
    valid_request = true;
  }
  if (is_explicit()) {
    *value = data_.ival;
  }
  else {
    const Ioss::GroupingEntity *ge       = data_.ge;
    const Ioss::Property        implicit = ge->get_implicit_property(name_);
    valid_request                        = implicit.get_value(value);
  }
  return valid_request;
}

bool Ioss::Property::get_value(double *value) const
{
  bool valid_request = false;
  if (type_ == REAL) {
    valid_request = true;
  }
  if (is_explicit()) {
    *value = data_.rval;
  }
  else {
    const Ioss::GroupingEntity *ge       = data_.ge;
    const Ioss::Property        implicit = ge->get_implicit_property(name_);
    valid_request                        = implicit.get_value(value);
  }
  return valid_request;
}

bool Ioss::Property::get_value(std::string *value) const
{
  bool valid_request = false;
  if (type_ == STRING) {
    valid_request = true;
  }
  if (is_explicit()) {
    *value = *(data_.sval);
  }
  else {
    const Ioss::GroupingEntity *ge       = data_.ge;
    const Ioss::Property        implicit = ge->get_implicit_property(name_);
    valid_request                        = implicit.get_value(value);
  }
  return valid_request;
}

bool Ioss::Property::get_value(void *&value) const
{
  bool valid_request = false;
  if (type_ == POINTER) {
    valid_request = true;
  }
  if (is_explicit()) {
    value = data_.pval;
  }
  else {
    const Ioss::GroupingEntity *ge       = data_.ge;
    const Ioss::Property        implicit = ge->get_implicit_property(name_);
    valid_request                        = implicit.get_value(value);
  }
  return valid_request;
}
