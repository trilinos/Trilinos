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

namespace stk
{
namespace coupling
{

class SyncInfo
{
public:
  template <typename ValueType>
  ValueType get_value(const std::string & parameterName) const;

  void set_value(const std::string & parameterName, const bool &value) { bvals[parameterName] = value; }
  void set_value(const std::string & parameterName, const int &value) { ivals[parameterName] = value; }
  void set_value(const std::string & parameterName, const double &value) { dvals[parameterName] = value; }
  void set_value(const std::string & parameterName, const std::string &value) { svals[parameterName] = value; }

  template <typename ValueType>
  bool has_value(const std::string & parameterName) const;

  SyncInfo exchange(stk::ParallelMachine global, stk::ParallelMachine local);

private:
  std::map<std::string, bool> bvals;
  std::map<std::string, int> ivals;
  std::map<std::string, double> dvals;
  std::map<std::string, std::string> svals;

  void pack(stk::CommBuffer & b) const;
  void unpack(stk::CommBuffer & b);
};

namespace impl {
  template<typename MapType>
  typename MapType::mapped_type get_value_from_map(const MapType& vals, const std::string& parameterName)
  {
    typename MapType::const_iterator iter = vals.find(parameterName);
    ThrowRequireMsg(iter != vals.end(), "get_value didn't find parameterName " << parameterName);
    return iter->second;
  }
}

template <>
inline bool SyncInfo::get_value<bool>(const std::string & parameterName) const
{
  return impl::get_value_from_map(bvals, parameterName);
}

template <>
inline int SyncInfo::get_value<int>(const std::string & parameterName) const
{
  return impl::get_value_from_map(ivals, parameterName);
}

template <>
inline double SyncInfo::get_value<double>(const std::string & parameterName) const
{
  return impl::get_value_from_map(dvals, parameterName);
}

template <>
inline std::string SyncInfo::get_value<std::string>(const std::string & parameterName) const
{
  return impl::get_value_from_map(svals, parameterName);
}

template <>
inline bool SyncInfo::has_value<bool>(const std::string & parameterName) const
{
  return bvals.count(parameterName) != 0;
}

template <>
inline bool SyncInfo::has_value<double>(const std::string & parameterName) const
{
  return dvals.count(parameterName) != 0;
}

template <>
inline bool SyncInfo::has_value<int>(const std::string & parameterName) const
{
  return ivals.count(parameterName) != 0;
}

template <>
inline bool SyncInfo::has_value<std::string>(const std::string & parameterName) const
{
  return svals.count(parameterName) != 0;
}

template<typename ValueType>
bool check_consistency(const SyncInfo & localSyncInfo, 
                       const SyncInfo & remoteSyncInfo, 
                       const std::string & parameterName)
{
  return (localSyncInfo.has_value<ValueType>(parameterName) && remoteSyncInfo.has_value<ValueType>(parameterName)) 
      && (localSyncInfo.get_value<ValueType>(parameterName) == remoteSyncInfo.get_value<ValueType>(parameterName));
}

template<typename ValueType>
bool copy_value(const SyncInfo & source, 
                const std::string & sourceParameterName,
                SyncInfo & destination, 
                const std::string & destinationParameterName)
{
  ThrowAssertMsg(source.has_value<ValueType>(sourceParameterName), "Parameter " << sourceParameterName << " is missing from source SyncInfo");
  if(source.has_value<ValueType>(sourceParameterName)) {
    std::cout<<source.get_value<std::string>(stk::coupling::AppName)<<": setting " << destinationParameterName << " from " << sourceParameterName << std::endl;
    destination.set_value(destinationParameterName, source.get_value<ValueType>(sourceParameterName));
    return true;
  }

  return false;
}

void check_sync_mode_consistency(const SyncInfo & myInfo,
                                 const SyncInfo & otherInfo);

double choose_value(const SyncInfo & myInfo,
                    const SyncInfo & otherInfo,
                    const std::string & parameterName,
                    stk::coupling::SyncMode syncMode);

} // namespace coupling
} // namespace stk

#endif /* STK_COUPLING_CONFIGURATION_INFO_HPP */
