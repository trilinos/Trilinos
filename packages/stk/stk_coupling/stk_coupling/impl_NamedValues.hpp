/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_COUPLING_IMPL_NAMED_VALUES_HPP
#define STK_COUPLING_IMPL_NAMED_VALUES_HPP

#include <stk_coupling/impl_Value.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <map>

namespace stk
{
namespace coupling
{
namespace impl
{

class NamedValues
{
public:

template<typename ValueType>
void set_value(const std::string & name, const ValueType & value)
{
  m_values[name] = Value(value);
}

template<typename ValueType>
bool has_value(const std::string & name) const
{
  MapType::const_iterator iter = m_values.find(name);
  if (iter == m_values.end()) {
    return false;
  }

  return (to_value_type<ValueType>() != INVALID_TYPE)
      && (to_value_type<ValueType>() == (iter->second).type);
}

template<typename ValueType>
ValueType get_value(const std::string & name) const
{
  MapType::const_iterator iter = m_values.find(name);
  STK_ThrowRequireMsg(iter != m_values.end(), "get_value: name='" << name << "' not found.");
  return std::any_cast<ValueType>((iter->second).value);
}

bool operator==(const NamedValues& rhs) const
{
  return m_values == rhs.m_values;
}

void pack(stk::CommBuffer& buf) const;
void unpack(stk::CommBuffer& buf);

using MapType = std::map<std::string, Value>;

MapType::const_iterator begin() const { return m_values.begin(); }
MapType::const_iterator end() const { return m_values.end(); }

private:
  MapType m_values;
};

} // namespace impl
} // namespace coupling
} // namespace stk
#endif /* STK_COUPLING_IMPL_NAMED_VALUES_HPP */
