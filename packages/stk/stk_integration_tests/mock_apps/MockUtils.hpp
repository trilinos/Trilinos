/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_UTILS_HPP
#define MOCK_UTILS_HPP

#include <stk_coupling/Constants.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace mock_utils {

inline
void exchange_and_print_info(const stk::coupling::SplitComms & splitComms,
                             const std::string & appName, int color)
{
  stk::coupling::SyncInfo myInfo;
  myInfo.set_value(stk::coupling::AppName, appName);
  myInfo.set_value("App Color", color);
  stk::coupling::SyncInfo::ColorToSyncInfoMap colorToSyncInfo = myInfo.exchange(splitComms);
  if (0 == stk::parallel_machine_rank(splitComms.get_split_comm())) {
    std::ostringstream os;
    const std::vector<int>& otherColors = splitComms.get_other_colors();
    os << appName << " detected " << otherColors.size()+1 << "-way coupling with other apps: ";
    for(int otherColor : otherColors) {
      os << "{" << colorToSyncInfo[otherColor].get_value<std::string>(stk::coupling::AppName)
         << ",color=" << colorToSyncInfo[otherColor].get_value<int>("App Color") << "} ";
    }
    os << std::endl;
    std::cout << os.str();
  }
}


}

#endif // MOCK_UTILS_HPP
