/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_diag_PrintTimer_hpp
#define stk_util_diag_PrintTimer_hpp

#include <iosfwd>

#include <stk_util/diag/Timer.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace stk {
namespace diag {

std::ostream &printXML(std::ostream& os, MetricsMask metrics_mask, bool checkpoint);

std::ostream &printTimersTable(std::ostream& os, Timer root_timer, MetricsMask metrics_mask, bool timer_checkpoint);

std::ostream &printTimersTable(std::ostream& os, Timer root_timer, MetricsMask metrics_mask, bool timer_checkpoint, ParallelMachine parallel_machine);

} // namespace diag
} // namespace stk

#endif // stk_util_diag_PrintTimer_hpp
