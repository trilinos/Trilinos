/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_COUPLING_IMPL_VALUE_HPP
#define STK_COUPLING_IMPL_VALUE_HPP

#include <stk_util/util/ReportHandler.hpp>
#include <Teuchos_any.hpp>
#include <cstdint>
#include <string>
#include <vector>


namespace stk {
namespace coupling {
namespace impl {

enum ValueTypes : uint8_t
{
  BOOL = 0,
  INT,
  FLOAT,
  DOUBLE,
  STRING,
  VECTOR_INT,
  VECTOR_FLOAT,
  VECTOR_DOUBLE,
  VECTOR_STRING,
  VECTOR_PAIR_STRING_INT,
  INVALID_TYPE
};

template<typename T> ValueTypes to_value_type();
template<>    inline ValueTypes to_value_type<bool>() { return BOOL; }
template<>    inline ValueTypes to_value_type<int>() { return INT; }
template<>    inline ValueTypes to_value_type<float>() { return FLOAT; }
template<>    inline ValueTypes to_value_type<double>() { return DOUBLE; }
template<>    inline ValueTypes to_value_type<std::string>() { return STRING; }
template<>    inline ValueTypes to_value_type<std::vector<int>>() { return VECTOR_INT; }
template<>    inline ValueTypes to_value_type<std::vector<float>>() { return VECTOR_FLOAT; }
template<>    inline ValueTypes to_value_type<std::vector<double>>() { return VECTOR_DOUBLE; }
template<>    inline ValueTypes to_value_type<std::vector<std::string>>() { return VECTOR_STRING; }
template<>    inline ValueTypes to_value_type<std::vector<std::pair<std::string,int>>>() { return VECTOR_PAIR_STRING_INT; }
template<>    inline ValueTypes to_value_type<char*>() { return INVALID_TYPE; }

struct Value
{
  Value() : value(), type(INVALID_TYPE) {}

  template<typename ValueType>
  Value(const ValueType& val) : value(val), type(to_value_type<ValueType>()) {}

  bool operator==(const Value& rhs) const
  {
    return type == rhs.type && value == rhs.value;
  }

  Teuchos::any value;
  ValueTypes type;
};

template <typename CommandType>
void execute(CommandType& command, Value& value)
{
  switch(value.type) {
    case BOOL: command.template execute<bool>(value); break;
    case INT: command.template execute<int>(value); break;
    case FLOAT: command.template execute<float>(value); break;
    case DOUBLE: command.template execute<double>(value); break;
    case STRING: command.template execute<std::string>(value); break;
    case VECTOR_INT: command.template execute<std::vector<int>>(value); break;
    case VECTOR_FLOAT: command.template execute<std::vector<float>>(value); break;
    case VECTOR_DOUBLE: command.template execute<std::vector<double>>(value); break;
    case VECTOR_STRING: command.template execute<std::vector<std::string>>(value); break;
    case VECTOR_PAIR_STRING_INT: command.template execute<std::vector<std::pair<std::string,int>>>(value); break;
    default: ThrowErrorMsg("Found unsupported type: " << value.type << " while executing command '" << command.name() << "'"); break;
  }
}

} // namespace impl
} // namespace coupling
} // namespace stk

#endif
