/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/Utils.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/string_case_compare.hpp>
#include <stk_util/command_line/CommandLineParser.hpp>

#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include <limits>
#include <cctype>

namespace stk
{
namespace coupling
{

void check_sync_mode_consistency(const SyncInfo & myInfo,
                                 const SyncInfo & otherInfo)
{
  bool myInfoHasMode = myInfo.has_value<int>(TimeSyncMode);
  bool otherInfoHasMode = otherInfo.has_value<int>(TimeSyncMode);
  ThrowRequireMsg(myInfoHasMode || otherInfoHasMode,
                  "neither myInfo nor otherInfo contains value for '" << TimeSyncMode << "'");
  
  if (myInfoHasMode) {
    SyncMode myMode = myInfo.get_value<SyncMode>(TimeSyncMode);
    if (myMode == Minimum || myMode == Receive) {
      ThrowRequireMsg(otherInfo.has_value<int>(TimeSyncMode),
          "otherInfo doesn't contain value for '" << TimeSyncMode << "' and my mode is '" << myMode << "'");
    }
  }

  if (otherInfoHasMode) {
    SyncMode otherMode = otherInfo.get_value<SyncMode>(TimeSyncMode);
    if (otherMode == Minimum || otherMode == Receive) {
      ThrowRequireMsg(myInfo.has_value<int>(TimeSyncMode),
          "myInfo doesn't contain value for '" << TimeSyncMode << "' and other mode is '" << otherMode << "'");
    }
  }

  if (myInfoHasMode && otherInfoHasMode) {
    SyncMode myMode = myInfo.get_value<SyncMode>(TimeSyncMode);
    SyncMode otherMode = otherInfo.get_value<SyncMode>(TimeSyncMode);
    bool consistent = ((myMode==Any || otherMode==Any) && (myMode != otherMode)) ||
                    (myMode==Minimum && otherMode==Minimum) ||
                    (myMode==Send && otherMode==Receive) ||
                    (myMode==Receive && otherMode==Send);
    ThrowRequireMsg(consistent, "Inconsistent TimeSyncMode (my mode="<<myMode<<", other mode="<<otherMode
                               <<"). Required to both be Minimum, or one Send and one Receive or either can be 'Any'.");
  }
}

double choose_value(const SyncInfo & myInfo,
                    const SyncInfo & otherInfo,
                    const std::string & parameterName,
                    stk::coupling::SyncMode syncMode)
{
  double myValue = 0.0;
  double otherValue = 0.0;

  if (syncMode == Any) {
    if (otherInfo.has_value<stk::coupling::SyncMode>(stk::coupling::TimeSyncMode)) {
      SyncMode otherSyncMode = otherInfo.get_value<stk::coupling::SyncMode>(stk::coupling::TimeSyncMode);
      switch(otherSyncMode) {
        case Send: syncMode = Receive; break;
        case Receive: syncMode = Send; break;
        case Minimum: syncMode = Minimum; break;
        case Any:
        default: ThrowErrorMsg("choose_value: if syncMode is Any, then other syncMode must be Send, Receive or Minimum.");
          break;
      }
    }
    else {
      ThrowRequireMsg(myInfo.has_value<double>(parameterName),"choose_value: syncMode=Any, but no value for '"<<parameterName<<"' in myInfo.");
      return myInfo.get_value<double>(parameterName);
    }
  }

  if (syncMode != Receive) {
    ThrowRequireMsg(myInfo.has_value<double>(parameterName),
        "choose_value: myInfo "<<myInfo.get_name()<<" doesn't contain " << parameterName
         << " and sync mode is " << syncMode);
    myValue = myInfo.get_value<double>(parameterName);
  }
  if (syncMode != Send) {
    ThrowRequireMsg(otherInfo.has_value<double>(parameterName),
        "choose_value: otherInfo "<<otherInfo.get_name()<<" doesn't contain " << parameterName
         << " and sync mode is " << syncMode);
    otherValue = otherInfo.get_value<double>(parameterName);
  }

  if (syncMode == Minimum) {
    return std::min(myValue, otherValue);
  }

  return (syncMode == Receive) ? otherValue : myValue;
}

int string_to_color(const std::string& appString)
{
  ThrowRequireMsg(!appString.empty(), "string_to_color: empty input string, probable programming error.");

  std::hash<std::string> hasher;
  size_t hashedValue = hasher(appString);
  int intValue = static_cast<int>(hashedValue % std::numeric_limits<int>::max());
  ThrowRequireMsg(intValue >= 0, "string_to_color is required to produce non-negative int, but produced '"<<intValue<<"'");
  return intValue;
}

SyncMode string_to_sync_mode(const std::string& syncModeString)
{
  if (stk::equal_case(syncModeString, "Minimum")) return Minimum;
  if (stk::equal_case(syncModeString, "Receive")) return Receive;
  if (stk::equal_case(syncModeString, "Send")) return Send;
  ThrowErrorMsg("string_to_sync_mode: invalid sync mode: " << syncModeString);
  return Send;
}

std::ostream& operator<<(std::ostream& os, const SyncMode & mode)
{
  switch(mode) {
  case Minimum: os << "Minimum"; break;
  case Receive: os << "Receive"; break;
  case Send: os << "Send"; break;
  case Any: os << "Any"; break;
  default: break;
  }
  return os;
}

}
}
