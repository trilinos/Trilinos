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

#include <stk_coupling/Constants.hpp>
#include <string>
#include <iostream>
#include <utility>

namespace stk
{
namespace coupling
{

SyncMode get_time_sync_mode(int argc, char** argv, const std::string& argName);

int string_to_color(const std::string& appString);

SyncMode string_to_sync_mode(const std::string& syncModeString);

std::ostream& operator<<(std::ostream& os, const SyncMode & mode);

}
}

#endif /* STK_COUPLING_UTILS_HPP */
