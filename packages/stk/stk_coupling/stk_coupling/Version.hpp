/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_VERSION_HPP
#define STK_COUPLING_VERSION_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/stk_config.h>
#include <string>

namespace stk
{
namespace coupling
{
 #ifndef STK_HIDE_DEPRECATED_CODE // delete after June 2022
STK_DEPRECATED int version();
#endif

}
}


#endif /* STK_COUPLING_VERSION_HPP */
