/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_COUPLING_SYNC_INFO_HPP
#define STK_COUPLING_SYNC_INFO_HPP

#include <map>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp> // for CommBuffer
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/impl_NamedValues.hpp>

namespace stk
{
namespace coupling
{

class SyncInfo
{
public:
  using ColorToSyncInfoMap = std::map<int, SyncInfo>;

  SyncInfo(const std::string& name = "")
    : m_name(name)
  { }

  template <typename ValueType>
  void set_value(const std::string & parameterName, const ValueType& value)
  {
    m_vals.set_value<ValueType>(parameterName, value);
  }

  template <typename ValueType>
  ValueType get_value(const std::string & parameterName) const
  {
    return m_vals.get_value<ValueType>(parameterName);
  }

  template <typename ValueType>
  ValueType get_value(const std::string & parameterName, const ValueType& defaultValue) const
  {
    if (has_value<ValueType>(parameterName)) return get_value<ValueType>(parameterName);
    else return defaultValue;
  }

  template <typename ValueType>
  bool has_value(const std::string & parameterName) const
  {
    return m_vals.has_value<ValueType>(parameterName);
  }

  SyncInfo exchange(const SplitComms & splitComms, int otherColor) const;

  ColorToSyncInfoMap exchange(const SplitComms & splitComms) const;

  const std::string& get_name() const { return m_name; }

  const impl::NamedValues& get_values() const { return m_vals; }

protected:
  void pack(stk::CommBuffer & b) const;
  void unpack(stk::CommBuffer & b);

  impl::NamedValues m_vals;
  std::string m_name;
};

template <>
inline void SyncInfo::set_value<SyncMode>(const std::string & parameterName, const SyncMode & value)
{
  set_value<int>(parameterName, static_cast<int>(value));
}

template <>
inline SyncMode SyncInfo::get_value<SyncMode>(const std::string & parameterName) const
{
  return static_cast<SyncMode>(get_value<int>(parameterName));
}

template <>
inline bool SyncInfo::has_value<SyncMode>(const std::string & parameterName) const
{
  return has_value<int>(parameterName);
}

} // namespace coupling
} // namespace stk

#endif /* STK_COUPLING_SYNC_INFO_HPP */

