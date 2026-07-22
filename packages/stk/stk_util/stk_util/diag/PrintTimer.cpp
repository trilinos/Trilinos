// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "stk_util/diag/PrintTimer.hpp"
#include "stk_util/diag/PrintTable.hpp"           // for operator<<, PrintTable, end_col, PrintT...
#include "stk_util/diag/Timer.hpp"                // for getEnabledTimerMetricsMask, Timer, Time...
#include "stk_util/diag/WriterExt.hpp"            // for operator<<
#include "stk_util/diag/ParallelTimerImpl.hpp"
#include "stk_util/environment/WallTime.hpp"      // for wall_time
#include "stk_util/parallel/Parallel.hpp"         // for parallel_machine_rank, MPI_Gather, para...
#include "stk_util/stk_config.h"                  // for STK_HAS_MPI
#include "stk_util/util/Marshal.hpp"              // for operator>>, Marshal, operator<<
#include "stk_util/util/Writer.hpp"               // for operator<<, Writer, dendl, pop, push
#include "stk_util/util/WriterManip.hpp"          // for hex
#include <cstddef>                                // for size_t
#include <algorithm>                              // for find_if, max, min
#include <functional>                             // for unary_function
#include <iomanip>                                // for setw, operator<<, _Setw, setprecision
#include <limits>                                 // for numeric_limits
#include <list>                                   // for list, _List_const_iterator, list<>::con...
#include <ostream>                                // for operator<<, basic_ostream, basic_ostrea...
#include <stdexcept>                              // for runtime_error
#include <string>                                 // for string, char_traits, operator<<
#include <vector>                                 // for vector


namespace stk {
namespace diag {

namespace {

/**
 * Class <b>Percent</b> is a functor which display the percentage of the numerator
 * to the denominator.  The value is displayed at (xx.xx%) for values in the range 0.01
 * to 99.99, (0.00%) if zero, (<0.01%) is less than 0.01, and (100.0%) for 100 percent.
 *
 */
struct Percent
{
  Percent(double numerator, double denominator)
    : m_numerator(numerator),
      m_denominator(denominator)
  {}

  /**
   * Member function <b>operator()</b> writes the percentage as a string to the
   * output stream.
   *
   * @param os    a <b>std::ostream</b> reference to the output stream
   *        to write to.
   *
   * @return      a <b>std::ostream</b> reference to the output stream
   *        written to.
   */
  std::ostream &operator()(std::ostream &os) const;

private:
  double    m_numerator;
  double    m_denominator;
};


std::ostream &
Percent::operator()(
  std::ostream &  os) const
{
  std::ostringstream strout;

  if (m_numerator == 0.0)
    strout << "(0.00%)";
  else if (m_denominator == 0.0)
    strout << "( NaN)";
  else {
    double ratio = m_numerator/m_denominator*100.0;
    if (ratio < 0.01)
      strout << "(<0.01%)";
    else if (ratio >= 100.0)
      strout << "(" << std::setw(5) << std::setprecision(1) << std::fixed << ratio << "%)";
    else if (ratio >= 10.0)
      strout << "(" << std::setw(5) << std::setprecision(2) << std::fixed << ratio << "%)";
    else
      strout << "(" << std::setw(5) << std::setprecision(3) << std::fixed << ratio << "%)";
  }

  return os << strout.str();
}

inline std::ostream &operator<<(std::ostream &os, const Percent &p) {
  return p(os);
}

// PrintTable &printTable(PrintTable &table, MPI_Comm mpi_comm, MetricsMask metrics_mask) const;

PrintTable &
printSubtable(
  PrintTable &      table,
  const Timer &                 root_timer,
  const Timer &                 timer,
  MetricsMask      metrics_mask,
  int        depth,
  bool        timer_checkpoint)
{
  if (timer.getSubtimerLapCount() != 0.0) {
    if (timer.shouldRecord()) {
      if (timer.getTimerMask() == 0 || timer.getMetric<LapCount>().getAccumulatedLap(timer_checkpoint) > 0) {
        table << justify(PrintTable::Cell::LEFT) << indent(depth) << timer.getName() << end_col
              << justify(PrintTable::Cell::RIGHT) << timer.getMetric<LapCount>().getAccumulatedLap(timer_checkpoint) << end_col;

        if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<CPUTime>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<CPUTime>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<WallTime>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<WallTime>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPICount>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<MPICount>().getAccumulatedLap(timer_checkpoint)) << end_col;
        if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
          table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint))
                << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint), root_timer.getMetric<MPIByteCount>().getAccumulatedLap(timer_checkpoint)) << end_col;
      }
      else
        table << justify(PrintTable::Cell::LEFT) << indent(depth) << span << timer.getName() << end_col;

      table << end_row;
      depth++;
    }

    for (TimerList::const_iterator it = timer.begin(); it != timer.end(); ++it)
      printSubtable(table, root_timer, *it, metrics_mask, depth, timer_checkpoint);
  }

  return table;
}


PrintTable &
printSubtable(
  PrintTable &      table,
  const impl::ParallelTimer &         root_timer,
  const impl::ParallelTimer &         timer,
  MetricsMask      metrics_mask,
  int        depth,
  bool        timer_checkpoint)
{
  if (timer.m_subtimerLapCount != 0.0) {
    if (timer.m_timerMask == 0 || timer.getMetric<LapCount>().m_sum > 0) {
      table << justify(PrintTable::Cell::LEFT) << indent(depth) << timer.m_name << end_col
            << justify(PrintTable::Cell::RIGHT) << timer.getMetric<LapCount>().m_sum << end_col;

      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_sum, root_timer.getMetric<CPUTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_min, root_timer.getMetric<CPUTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<CPUTime>::format(timer.getMetric<CPUTime>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<CPUTime>().m_max, root_timer.getMetric<CPUTime>().m_sum) << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_sum, root_timer.getMetric<WallTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_min, root_timer.getMetric<WallTime>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<WallTime>::format(timer.getMetric<WallTime>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<WallTime>().m_max, root_timer.getMetric<WallTime>().m_sum) << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_sum, root_timer.getMetric<MPICount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_min, root_timer.getMetric<MPICount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPICount>::format(timer.getMetric<MPICount>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<MPICount>().m_max, root_timer.getMetric<MPICount>().m_sum) << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_sum)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_sum, root_timer.getMetric<MPIByteCount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_min)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_min, root_timer.getMetric<MPIByteCount>().m_sum) << end_col
              << justify(PrintTable::Cell::RIGHT) << std::setw(12) << MetricTraits<MPIByteCount>::format(timer.getMetric<MPIByteCount>().m_max)
              << " " << std::setw(8) << Percent(timer.getMetric<MPIByteCount>().m_max, root_timer.getMetric<MPIByteCount>().m_sum) << end_col;
    }
    else
      table << justify(PrintTable::Cell::LEFT) << indent(depth) << span << timer.m_name << end_col;

    table << end_row;
    depth++;
  }

  for (std::list<impl::ParallelTimer>::const_iterator it = timer.m_subtimerList.begin(); it != timer.m_subtimerList.end(); ++it)
    printSubtable(table, root_timer, *it, metrics_mask, depth, timer_checkpoint);

  return table;
}


PrintTable &
printTable(
  PrintTable &          table,
  Timer &               root_timer,
  MetricsMask           metrics_mask,
  size_t                name_width,
  bool                  timer_checkpoint)
{
  updateRootTimer(root_timer);

  root_timer.accumulateSubtimerLapCounts();

  if (metrics_mask & getEnabledTimerMetricsMask()) {
    table.setAutoEndCol(false);

    table << cell_width(name_width) << justify(PrintTable::Cell::CENTER) << "Timer" << (timer_checkpoint ? " (delta time)" : "") << end_col
          << justify(PrintTable::Cell::CENTER) << "Count"  << end_col;

    if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<CPUTime>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col;
    if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<WallTime>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col;
    if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPICount>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col;
    if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
      table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col;

    table << end_header;

    printSubtable(table, root_timer, root_timer, metrics_mask, 0, timer_checkpoint);

    if (timer_checkpoint)
      root_timer.checkpoint();
  }

  return table;
}


PrintTable &
printTable(
  PrintTable &          table,
  Timer &               root_timer,
  MetricsMask           metrics_mask,
  size_t                name_width,
  bool                  timer_checkpoint,
  ParallelMachine       parallel_machine)
{
  updateRootTimer(root_timer);

  root_timer.accumulateSubtimerLapCounts();

  impl::ParallelTimer parallel_timer = stk::diag::impl::collect_timers(root_timer, timer_checkpoint, parallel_machine);

  int parallel_rank = parallel_machine_rank(parallel_machine);
  if (parallel_rank == 0) {
    if (metrics_mask & getEnabledTimerMetricsMask()) {
      table.setAutoEndCol(false);

      table << end_col << end_col;

      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<CPUTime>::table_header() << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<WallTime>::table_header() << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPICount>::table_header() << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col
              << justify(PrintTable::Cell::CENTER) << MetricTraits<MPIByteCount>::table_header() << end_col;

      table << end_header;
      table << cell_width(name_width) << justify(PrintTable::Cell::CENTER) << "Timer" << (timer_checkpoint ? " (delta time)" : "") << end_col
            << justify(PrintTable::Cell::CENTER) << "Count"  << end_col;

      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<CPUTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<WallTime>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPICount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;
      if (metrics_mask & getEnabledTimerMetricsMask() & MetricTraits<MPIByteCount>::METRIC)
        table << justify(PrintTable::Cell::CENTER) << "Sum (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Min (% of System)" << end_col
              << justify(PrintTable::Cell::CENTER) << "Max (% of System)" << end_col;

      table << end_header;

      printSubtable(table, parallel_timer, parallel_timer, metrics_mask, 0, timer_checkpoint);
    }

    if (timer_checkpoint)
      root_timer.checkpoint();
  }

  return table;
}

} // namespace <empty>

void printTimeToPrintTable(std::ostream& os, double durationToPrintTable)
{
    os << "Took " << durationToPrintTable << " seconds to generate the table above." << std::endl;
}

std::ostream &printTimersTable(std::ostream& os, Timer root_timer, MetricsMask metrics_mask, bool timer_checkpoint)
{
  double startTimeToPrintTable = stk::wall_time();
  stk::PrintTable print_table;

  printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint);

  os << print_table;

  double durationToPrintTable = stk::wall_time() - startTimeToPrintTable;
  printTimeToPrintTable(os, durationToPrintTable);
  return os;
}


std::ostream &printTimersTable(std::ostream& os, Timer root_timer, MetricsMask metrics_mask, bool timer_checkpoint, ParallelMachine parallel_machine)
{
  double startTimeToPrintTable = stk::wall_time();
  stk::PrintTable print_table;

  int parallel_size = parallel_machine_size(parallel_machine);
  if (parallel_size == 1)
    printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint);
  else
    printTable(print_table, root_timer, metrics_mask, 40, timer_checkpoint, parallel_machine);

  os << print_table;

  double durationToPrintTable = stk::wall_time() - startTimeToPrintTable;
  if (parallel_machine_rank(parallel_machine) == 0)
    printTimeToPrintTable(os, durationToPrintTable);
  return os;
}

} // namespace diag


} // namespace stk
