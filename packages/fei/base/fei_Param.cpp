/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <fei_Param.hpp>

fei::Param::Param(const char* name,
		  const char* value)
  : type_(STRING),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
  if (value != 0) string_value_ = value;
}

fei::Param::Param(const char* name,
		  double value)
  : type_(DOUBLE),
    name_(),
    string_value_(),
    double_value_(value),
    int_value_(0),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  int value)
  : type_(INT),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(value),
    bool_value_(false),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  bool value)
  : type_(BOOL),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(value),
    void_value_(NULL)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const char* name,
		  const void* value)
  : type_(VOID),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    bool_value_(false),
    void_value_(value)
{
  if (name != 0) name_ = name;
}

fei::Param::Param(const fei::Param& src)
  : type_(src.type_),
    name_(),
    string_value_(),
    double_value_(0.0),
    int_value_(0),
    void_value_(NULL)
{
  *this = src;
}

fei::Param::~Param()
{
}

fei::Param& fei::Param::operator=(const fei::Param& src)
{
  name_ = src.name_;
  switch(type_) {
  case STRING:
    string_value_ = src.string_value_; break;
  case DOUBLE:
    double_value_ = src.double_value_; break;
  case INT:
    int_value_ = src.int_value_; break;
  case VOID:
    void_value_ = src.void_value_; break;
  case BOOL:
    bool_value_ = src.bool_value_; break;
  case BAD_TYPE:
    break;
  }

  return(*this);
}
