/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Kokkos_Macros.hpp>
#include <stk_ngp/NgpProfilingBlock.hpp>
// mfh 27 Sep 2017: Note that it is possible to disable Kokkos'
// profiling interface at configure time.  There's no good reason for
// this as far as I know, but it's possible.
#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#endif // defined(KOKKOS_ENABLE_PROFILING)

namespace ngp
{

#if defined(KOKKOS_ENABLE_PROFILING)
ProfilingBlock::ProfilingBlock(const std::string & name) { Kokkos::Profiling::pushRegion(name); }
#else
ProfilingBlock::ProfilingBlock(const std::string & /* name */) {}
#endif // defined(KOKKOS_ENABLE_PROFILING)

ProfilingBlock::~ProfilingBlock()
{
#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::popRegion();
#endif // defined(KOKKOS_ENABLE_PROFILING)
}
}
