/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_SPLITCOMMS_SINGLETON_HPP
#define STK_COUPLING_SPLITCOMMS_SINGLETON_HPP

#include <stk_coupling/SplitComms.hpp>

namespace stk
{
namespace coupling
{

SplitComms& get_split_comms_singleton();

void set_split_comms_singleton(const SplitComms& splitComms);

}
}

#endif /* STK_COUPLING_SPLITCOMMS_SINGLETON_HPP */
