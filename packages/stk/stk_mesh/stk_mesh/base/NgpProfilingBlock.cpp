/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Kokkos_Macros.hpp>
#include <stk_mesh/base/NgpProfilingBlock.hpp>
#include <impl/Kokkos_Profiling.hpp>

namespace stk {
namespace mesh {

ProfilingBlock::ProfilingBlock(const std::string & name) { Kokkos::Profiling::pushRegion(name); }

ProfilingBlock::~ProfilingBlock()
{
  Kokkos::Profiling::popRegion();
}

}
}
