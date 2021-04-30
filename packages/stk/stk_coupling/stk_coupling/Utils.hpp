/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_UTILS_HPP
#define STK_COUPLING_UTILS_HPP

#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/Constants.hpp>
#include <string>
#include <iostream>
#include <utility>

namespace stk
{
namespace coupling
{

template<typename ValueType>
bool check_consistency(const SyncInfo & localSyncInfo,
                       const SyncInfo & remoteSyncInfo,
                       const std::string & parameterName,
                       stk::coupling::SyncMode syncMode)
{
  bool hasLocal = localSyncInfo.has_value<ValueType>(parameterName);
  bool hasRemote = remoteSyncInfo.has_value<ValueType>(parameterName);

  bool consistent = false;
  consistent |= (syncMode == Receive && hasRemote);
  consistent |= (syncMode == Send && hasLocal);
  consistent |= (hasRemote && hasLocal &&
      (localSyncInfo.get_value<ValueType>(parameterName) == remoteSyncInfo.get_value<ValueType>(parameterName)));

  return consistent;
}

void check_sync_mode_consistency(const SyncInfo & myInfo,
                                 const SyncInfo & otherInfo);

double choose_value(const SyncInfo & myInfo,
                    const SyncInfo & otherInfo,
                    const std::string & parameterName,
                    stk::coupling::SyncMode syncMode);

SyncMode get_time_sync_mode(int argc, char** argv, const std::string& argName);

int string_to_color(const std::string& appString);

SyncMode string_to_sync_mode(const std::string& syncModeString);

std::ostream& operator<<(std::ostream& os, const SyncMode & mode);

}
}

#endif /* STK_COUPLING_UTILS_HPP */
