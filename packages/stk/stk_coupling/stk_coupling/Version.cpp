/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/Version.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_util/parallel/CouplingVersions.hpp>

namespace stk
{
namespace coupling
{

int version()
{
  return stk::util::get_common_coupling_version();
}

} // namespace coupling
} // namespace stk
