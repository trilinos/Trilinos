/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <Ioss_Property.h>

#include <Ioss_GroupingEntity.h>
#include <Ioss_Utils.h>

#include <string>

namespace {
  std::string type_string(Ioss::Property::BasicType type)
  {
    switch (type) {
    case Ioss::Property::INVALID:
      return std::string("invalid");
    case Ioss::Property::REAL:
      return std::string("real");
    case Ioss::Property::INTEGER:
      return std::string("integer");
    case Ioss::Property::POINTER:
      return std::string("pointer");
    case Ioss::Property::STRING:
      return std::string("string");
    default:
      return std::string("internal error");
    }
  }
  
  void error_message(const Ioss::Property &property, const std::string &requested_type)
  {
    std::ostringstream errmsg;
    errmsg << "ERROR: For property named '" << property.get_name()
	   << "', code requested value of type '" << requested_type
	   << "', but property type is '" << type_string(property.get_type())
	   << "'. Types must match\n";
    IOSS_ERROR(errmsg);
  }
}

Ioss::Property::Property() :
  name_(""), type_(INVALID), storage_(UNKNOWN_VAR_TYPE), isImplicit_(false)
{ data_.pval = NULL; }

Ioss::Property::Property(const std::string &name,
			 const BasicType type,
			 const VariableType storage,
			 void *data, bool is_it_implicit) :
  name_(name), type_(type), storage_(storage), isImplicit_(is_it_implicit)
{ data_.pval = data; }

Ioss::Property::Property(const std::string &name, const int    value,
			     bool is_it_implicit) :
  name_(name), type_(INTEGER), storage_(SCALAR), isImplicit_(is_it_implicit)
{ data_.ival = value; }

Ioss::Property::Property(const std::string &name, const double   value,
			     bool is_it_implicit) :
  name_(name), type_(REAL), storage_(SCALAR), isImplicit_(is_it_implicit)
{ data_.rval = value; }

Ioss::Property::Property(const std::string &name, const std::string &value,
			     bool is_it_implicit) :
  name_(name), type_(STRING), storage_(SCALAR), isImplicit_(is_it_implicit)
{ data_.sval = new std::string(value); }

// To set implicit property
Ioss::Property::Property(const Ioss::GroupingEntity* ge,
			     const std::string &name, const BasicType type) :
  name_(name), type_(type), storage_(SCALAR), isImplicit_(true)
{ data_.ge = ge; }

Ioss::Property::Property(const Ioss::Property& from)
{
  name_ = from.name_;
  type_ = from.type_;
  storage_ = from.storage_;
  isImplicit_ = from.isImplicit_;
  if (!isImplicit_ && type_ == STRING)
    data_.sval = new std::string(*(from.data_.sval));
  else
    data_ = from.data_;
}

Ioss::Property::~Property()
{ if (!isImplicit_ && type_ == STRING) delete data_.sval; }


std::string Ioss::Property::get_string()  const
{
  std::string value;
  bool valid = get_value(&value);
  if (!valid) {
    error_message(*this, "string");
  }
  return value;
}

int    Ioss::Property::get_int()     const
{
  int value;
  bool valid = get_value(&value);
  if (!valid) {
    error_message(*this, "int");
  }
  return value;
}

double   Ioss::Property::get_real()    const
{
  double value;
  bool valid = get_value(&value);
  if (!valid) {
    error_message(*this, "real");
  }
  return value;
}

void*  Ioss::Property::get_pointer() const
{
  void* value = NULL;
  bool valid = get_value(value);
  if (!valid) {
    error_message(*this, "pointer");
  }
  return value;
}

bool Ioss::Property::get_value(int *value) const
{
  bool valid_request = false;
  if (type_ == INTEGER) {
    valid_request = true;
  }
  if (is_explicit())
    *value = data_.ival;
  else {
    const Ioss::GroupingEntity* ge = data_.ge;
    const Ioss::Property implicit = ge->get_implicit_property(name_);
    valid_request = implicit.get_value(value);
  }
  return valid_request;
}

bool Ioss::Property::get_value(double *value) const
{
  bool valid_request = false;
  if (type_ == REAL) {
    valid_request = true;
  }
  if (is_explicit())
    *value = data_.rval;
  else {
    const Ioss::GroupingEntity* ge = data_.ge;
    const Ioss::Property implicit = ge->get_implicit_property(name_);
    valid_request = implicit.get_value(value);
  }
  return valid_request;
}


bool Ioss::Property::get_value(std::string *value) const
{
  bool valid_request = false;
  if (type_ == STRING) {
    valid_request = true;
  }
  if (is_explicit())
    *value = *(data_.sval);
  else {
    const Ioss::GroupingEntity* ge = data_.ge;
    const Ioss::Property implicit = ge->get_implicit_property(name_);
    valid_request = implicit.get_value(value);
  }
  return valid_request;
}

bool Ioss::Property::get_value(void*   value) const
{
  bool valid_request = false;
  if (type_ == POINTER) {
    valid_request = true;
  }
  if (is_explicit())
    value = data_.pval;
  else {
    const Ioss::GroupingEntity* ge = data_.ge;
    const Ioss::Property implicit = ge->get_implicit_property(name_);
    valid_request = implicit.get_value(value);
  }
  return valid_request;
}
