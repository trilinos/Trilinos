/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_COUPLING_OLDSYNCINFO_HPP
#define STK_COUPLING_OLDSYNCINFO_HPP

#include <map>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp> // for CommBuffer
#include <stk_coupling/Constants.hpp>
#include <stk_util/stk_config.h>

#ifndef STK_HIDE_DEPRECATED_CODE
namespace stk
{
namespace coupling
{


class STK_DEPRECATED OldSyncInfo
{
public:
  OldSyncInfo()
  {}

  template <typename ValueType>
  void set_value(const std::string & parameterName, const ValueType& value);

  template <typename ValueType>
  ValueType get_value(const std::string & parameterName) const;

  template <typename ValueType>
  ValueType get_value(const std::string & parameterName, const ValueType& defaultValue) const
  {
    if (has_value<ValueType>(parameterName)) return get_value<ValueType>(parameterName);
    else return defaultValue;
  }

  template <typename ValueType>
  bool has_value(const std::string & parameterName) const;

  OldSyncInfo exchange(stk::ParallelMachine global, stk::ParallelMachine local);

private:
  template<typename MapType>
  typename MapType::mapped_type get_value_from_map(const MapType& vals, const std::string& parameterName) const
  {
    typename MapType::const_iterator iter = vals.find(parameterName);
    ThrowRequireMsg(iter != vals.end(), "OldSyncInfo::get_value didn't find parameterName " << parameterName);
    return iter->second;
  }

  std::map<std::string, int> ivals;
  std::map<std::string, double> dvals;
  std::map<std::string, std::string> svals;

  void pack(stk::CommBuffer & b) const;
  void unpack(stk::CommBuffer & b);

};

template <>
inline void OldSyncInfo::set_value<int>(const std::string & parameterName, const int & value)
{
  ivals[parameterName] = value;
}

template <>
inline void OldSyncInfo::set_value<SyncMode>(const std::string & parameterName, const SyncMode & value)
{
  set_value<int>(parameterName, static_cast<int>(value));
}

template <>
inline void OldSyncInfo::set_value<double>(const std::string & parameterName, const double & value)
{
  dvals[parameterName] = value;
}

template <>
inline void OldSyncInfo::set_value<std::string>(const std::string & parameterName, const std::string & value)
{
  svals[parameterName] = value;
}

template <>
inline int OldSyncInfo::get_value<int>(const std::string & parameterName) const
{
  return get_value_from_map(ivals, parameterName);
}

template <>
inline SyncMode OldSyncInfo::get_value<SyncMode>(const std::string & parameterName) const
{
  return static_cast<SyncMode>(get_value<int>(parameterName));
}

template <>
inline double OldSyncInfo::get_value<double>(const std::string & parameterName) const
{
  return get_value_from_map(dvals, parameterName);
}

template <>
inline std::string OldSyncInfo::get_value<std::string>(const std::string & parameterName) const
{
  return get_value_from_map(svals, parameterName);
}

template <>
inline bool OldSyncInfo::has_value<int>(const std::string & parameterName) const
{
  return ivals.count(parameterName) != 0;
}

template <>
inline bool OldSyncInfo::has_value<SyncMode>(const std::string & parameterName) const
{
  return has_value<int>(parameterName);
}

template <>
inline bool OldSyncInfo::has_value<double>(const std::string & parameterName) const
{
  return dvals.count(parameterName) != 0;
}

template <>
inline bool OldSyncInfo::has_value<std::string>(const std::string & parameterName) const
{
  return svals.count(parameterName) != 0;
}

} // namespace coupling
} // namespace stk

#endif  /* STK_HIDE_DEPRECATED_CODE */

#endif /* STK_COUPLING_OLDSYNCINFO_HPP */
