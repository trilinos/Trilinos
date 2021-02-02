/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_COUPLING_CONFIGURATION_INFO_HPP
#define STK_COUPLING_CONFIGURATION_INFO_HPP

#include <map>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp> // for CommBuffer
#include <stk_coupling/Constants.hpp>

namespace stk
{
namespace coupling
{

class ConfigurationInfo
{
public:
  template <typename ValueType>
  ValueType get_value(const std::string & parameterName) const;

  void set_value(const std::string & parameterName, const double &value) { dvals[parameterName] = value; }
  void set_value(const std::string & parameterName, const int &value) { ivals[parameterName] = value; }
  void set_value(const std::string & parameterName, const std::string &value) { svals[parameterName] = value; }

  template <typename ValueType>
  bool has_value(const std::string & paramterName) const;

  ConfigurationInfo exchange(stk::ParallelMachine global, stk::ParallelMachine local);

private:
  std::map<std::string, double> dvals;
  std::map<std::string, int> ivals;
  std::map<std::string, std::string> svals;

  void pack(stk::CommBuffer & b) const;
  void unpack(stk::CommBuffer & b);
};

template <>
inline std::string ConfigurationInfo::get_value<std::string>(const std::string & parameterName) const
{
  ThrowAssertMsg(svals.count(parameterName) != 0, "get_value<std::string> didn't find parameterName "<<parameterName);
  return svals.at(parameterName);
}

template <>
inline double ConfigurationInfo::get_value<double>(const std::string & parameterName) const
{
  ThrowAssertMsg(dvals.count(parameterName) != 0, "get_value<double> didn't find parameterName "<<parameterName);
  return dvals.at(parameterName);
}

template <>
inline int ConfigurationInfo::get_value<int>(const std::string & parameterName) const
{
  ThrowAssertMsg(ivals.count(parameterName) != 0, "get_value<int> didn't find parameterName "<<parameterName);
  return ivals.at(parameterName);
}

template <>
inline bool ConfigurationInfo::has_value<double>(const std::string & parameterName) const
{
  return dvals.count(parameterName) != 0;
}

template <>
inline bool ConfigurationInfo::has_value<int>(const std::string & parameterName) const
{
  return ivals.count(parameterName) != 0;
}

template <>
inline bool ConfigurationInfo::has_value<std::string>(const std::string & parameterName) const
{
  return svals.count(parameterName) != 0;
}

template<typename ValueType>
bool check_consistency(const ConfigurationInfo & localConfigurationInfo, 
                       const ConfigurationInfo & remoteConfigurationInfo, 
                       const std::string & parameterName)
{
  return (localConfigurationInfo.has_value<ValueType>(parameterName) && remoteConfigurationInfo.has_value<ValueType>(parameterName)) 
      && (localConfigurationInfo.get_value<ValueType>(parameterName) == remoteConfigurationInfo.get_value<ValueType>(parameterName));
}

template<typename ValueType>
bool copy_value(const ConfigurationInfo & source, 
                const std::string & sourceParameterName,
                ConfigurationInfo & destination, 
                const std::string & destinationParameterName)
{
  ThrowAssertMsg(source.has_value<ValueType>(sourceParameterName), "Parameter " << sourceParameterName << " is missing from source ConfigurationInfo");
  if(source.has_value<ValueType>(sourceParameterName)) {
    std::cout<<source.get_value<std::string>(stk::coupling::AppName)<<": setting " << destinationParameterName << " from " << sourceParameterName << std::endl;
    destination.set_value(destinationParameterName, source.get_value<ValueType>(sourceParameterName));
    return true;
  }

  return false;
}

void check_sync_mode_consistency(const ConfigurationInfo & myInfo,
                                 const ConfigurationInfo & otherInfo);

double choose_value(const ConfigurationInfo & myInfo,
                    const ConfigurationInfo & otherInfo,
                    const std::string & parameterName,
                    stk::coupling::SyncMode syncMode);

} // namespace coupling
} // namespace stk

#endif /* STK_COUPLING_CONFIGURATION_INFO_HPP */
