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
  default: break;
  }
  return os;
}

}
}
