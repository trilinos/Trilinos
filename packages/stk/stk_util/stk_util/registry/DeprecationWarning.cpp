/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "stk_util/registry/DeprecationWarning.hpp"
#include "stk_util/environment/Env.hpp"             // for parallel_rank
#include "stk_util/environment/RuntimeDoomed.hpp"   // for RuntimeDoomed
#include "stk_util/environment/RuntimeWarning.hpp"  // for RuntimeWarning
#include <iostream>                                 // for cerr
#include <string>                                   // for operator<<

namespace stk
{
namespace util
{

DeprecationWarning::DeprecationWarning(const VersionNumber & removal_version)
    : my_removal_version(removal_version)
{
  if( my_removal_version.is_release() ) {
    throw std::runtime_error("Cannot set deprecations on release versions");
  }

  message << "Deprecated feature removed in " << removal_version << " detected.\n";
}

DeprecationWarning::~DeprecationWarning()
{
  if (VersionNumber::current_version() < my_removal_version)
  {
    if (sierra::Env::parallel_rank() == 0)
    {
      std::cerr << "WARNING: " << message.str() << std::endl;
    }
    stk::RuntimeWarningP0() << message.str();
  }
  else
  {
    stk::RuntimeDoomed() << message.str();
  }
}
}
}
