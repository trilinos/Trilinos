/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/SplitCommsSingleton.hpp>

namespace stk
{
namespace coupling
{

namespace {
  static SplitComms splitCommsSingleton;
}

SplitComms& get_split_comms_singleton()
{
  return splitCommsSingleton;
}

void set_split_comms_singleton(const SplitComms& splitComms)
{
  splitCommsSingleton = splitComms;
}

}
}
