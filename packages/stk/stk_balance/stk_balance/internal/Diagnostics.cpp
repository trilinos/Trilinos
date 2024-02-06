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
#include "Diagnostics.hpp"
#include "DiagnosticsContainer.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_balance/balanceUtils.hpp"
#include <iomanip>

namespace stk {
namespace balance {

void
UnsignedDiagnostic::collect_data(stk::ParallelMachine pm, int numRanks)
{
  stk::CommSparse comm(pm);
  stk::pack_and_communicate(comm, [this, &comm]()
  {
    for (const auto & localValue : m_localValues) {
      comm.send_buffer(0).pack<int>(localValue.first);
      comm.send_buffer(0).pack<unsigned>(localValue.second);
    }
  });

  if (stk::parallel_machine_rank(pm) == 0) {
    m_values.resize(numRanks, 0);

    for (int proc = 0; proc < stk::parallel_machine_size(pm); ++proc) {
      while (comm.recv_buffer(proc).remaining()) {
        int rank;
        unsigned value;
        comm.recv_buffer(proc).unpack<int>(rank);
        comm.recv_buffer(proc).unpack<unsigned>(value);
        STK_ThrowRequire(rank < numRanks);
        m_values[rank] = value;
      }
    }
  }
}

void
UnsignedDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    m_min = compute_min(m_values);
    m_max = compute_max(m_values);
    m_avg = compute_avg(m_values);
  }
}

unsigned
UnsignedDiagnostic::compute_min(const std::vector<unsigned> & data)
{
  return *std::min_element(data.begin(), data.end());
}

unsigned
UnsignedDiagnostic::compute_max(const std::vector<unsigned> & data)
{
  return *std::max_element(data.begin(), data.end());
}

unsigned
UnsignedDiagnostic::compute_avg(const std::vector<unsigned> & data)
{
  if (data.empty()) {
    return 0;
  }
  return std::accumulate(data.begin(), data.end(), 0ull) / data.size();
}


void
UnsignedWithPercentDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    m_min = compute_min(m_values);
    m_max = compute_max(m_values);
    m_avg = compute_avg(m_values);

    const double maxPercent = (m_avg > 0) ? 100. * (static_cast<int>(m_max) - static_cast<int>(m_avg)) / static_cast<double>(m_avg)
                                          : 0.0;
    const double minPercent = (m_avg > 0) ? 100. * (static_cast<int>(m_min) - static_cast<int>(m_avg)) / static_cast<double>(m_avg)
                                          : 0.0;

    std::ostringstream osMax;
    osMax << std::fixed << std::setprecision(1) << std::showpos << maxPercent;

    std::ostringstream osMin;
    osMin << std::fixed << std::setprecision(1) << std::showpos << minPercent;

    m_percentSize = std::max(osMax.str().size(), osMin.str().size());
  }
}

std::string
UnsignedWithPercentDiagnostic::print_value_with_percent(unsigned value)
{
  const double percent = (m_avg > 0) ? 100. * (static_cast<int>(value) - static_cast<int>(m_avg)) / static_cast<double>(m_avg)
                                     : 0.0;
  std::ostringstream os;
  os << value << " (" << std::setw(m_percentSize) << std::fixed << std::setprecision(1) << std::showpos
     << percent << "%)";

  return os.str();
}

std::string
UnsignedWithPercentDiagnostic::print_rank_value(unsigned , int rank)
{
  return print_value_with_percent(m_values[rank]);
}

std::string
UnsignedWithPercentDiagnostic::print_min(unsigned )
{
  return print_value_with_percent(m_min);
}

std::string
UnsignedWithPercentDiagnostic::print_max(unsigned )
{
  return print_value_with_percent(m_max);
}

std::string
UnsignedWithPercentDiagnostic::print_avg(unsigned )
{
  std::ostringstream os;
  os << m_avg << std::string(m_percentSize+4, ' ');

  return os.str();
}


MultiUnsignedDiagnostic::MultiUnsignedDiagnostic(unsigned numColumns)
  : m_localValues(numColumns),
    m_values(numColumns),
    m_min(numColumns),
    m_max(numColumns),
    m_avg(numColumns),
    m_numColumns(numColumns)
{
}

void
MultiUnsignedDiagnostic::collect_data(stk::ParallelMachine pm, int numRanks)
{
  stk::CommSparse comm(pm);
  stk::pack_and_communicate(comm, [this, &comm]()
  {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      for (const auto & localValue : m_localValues[column]) {
        comm.send_buffer(0).pack<unsigned>(column);
        comm.send_buffer(0).pack<int>(localValue.first);
        comm.send_buffer(0).pack<unsigned>(localValue.second);
      }
    }
  });

  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_values[column].resize(numRanks, 0);
    }

    for (int proc = 0; proc < stk::parallel_machine_size(pm); ++proc) {
      while (comm.recv_buffer(proc).remaining()) {
        unsigned column;
        int rank;
        unsigned value;
        comm.recv_buffer(proc).unpack<unsigned>(column);
        comm.recv_buffer(proc).unpack<int>(rank);
        comm.recv_buffer(proc).unpack<unsigned>(value);
        STK_ThrowRequire(rank < numRanks);
        m_values[column][rank] = value;
      }
    }
  }
}

void
MultiUnsignedDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_min[column] = compute_min(m_values[column]);
      m_max[column] = compute_max(m_values[column]);
      m_avg[column] = compute_avg(m_values[column]);
    }
  }
}

unsigned
MultiUnsignedDiagnostic::compute_min(const std::vector<unsigned> & data)
{
  return *std::min_element(data.begin(), data.end());
}

unsigned
MultiUnsignedDiagnostic::compute_max(const std::vector<unsigned> & data)
{
  return *std::max_element(data.begin(), data.end());
}

unsigned
MultiUnsignedDiagnostic::compute_avg(const std::vector<unsigned> & data)
{
  if (data.empty()) {
    return 0;
  }
  return std::accumulate(data.begin(), data.end(), 0ull) / data.size();
}


void
MultiUnsignedWithPercentDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_min[column] = compute_min(m_values[column]);
      m_max[column] = compute_max(m_values[column]);
      m_avg[column] = compute_avg(m_values[column]);

      const double maxPercent = (m_avg[column] > 0) ? 100. * (static_cast<int>(m_max[column]) - static_cast<int>(m_avg[column])) /
                                                      static_cast<double>(m_avg[column])
                                                    : 0.0;
      const double minPercent = (m_avg[column] > 0) ? 100. * (static_cast<int>(m_min[column]) - static_cast<int>(m_avg[column])) /
                                                      static_cast<double>(m_avg[column])
                                                    : 0.0;

      std::ostringstream osMax;
      osMax << std::fixed << std::setprecision(1) << std::showpos << maxPercent;

      std::ostringstream osMin;
      osMin << std::fixed << std::setprecision(1) << std::showpos << minPercent;

      m_percentSize[column] = std::max(osMax.str().size(), osMin.str().size());
    }
  }
}

std::string
MultiUnsignedWithPercentDiagnostic::print_value_with_percent(unsigned column, unsigned value)
{
  const double percent = (m_avg[column] > 0) ? 100. * (static_cast<int>(value) - static_cast<int>(m_avg[column])) /
                                               static_cast<double>(m_avg[column])
                                             : 0.0;
  std::ostringstream os;
  os << value << " (" << std::setw(m_percentSize[column]) << std::fixed << std::setprecision(1) << std::showpos
     << percent << "%)";

  return os.str();
}

std::string
MultiUnsignedWithPercentDiagnostic::print_rank_value(unsigned column, int rank)
{
  return print_value_with_percent(column, m_values[column][rank]);
}

std::string
MultiUnsignedWithPercentDiagnostic::print_min(unsigned column)
{
  return print_value_with_percent(column, m_min[column]);
}

std::string
MultiUnsignedWithPercentDiagnostic::print_max(unsigned column)
{
  return print_value_with_percent(column, m_max[column]);
}

std::string
MultiUnsignedWithPercentDiagnostic::print_avg(unsigned column)
{
  std::ostringstream os;
  os << m_avg[column] << std::string(m_percentSize[column]+4, ' ');

  return os.str();
}


void
DoubleWithPercentDiagnostic::collect_data(stk::ParallelMachine pm, int numRanks)
{
  stk::CommSparse comm(pm);
  stk::pack_and_communicate(comm, [this, &comm]()
  {
    for (const auto & localValue : m_localValues) {
      comm.send_buffer(0).pack<int>(localValue.first);
      comm.send_buffer(0).pack<double>(localValue.second);
    }
  });

  if (stk::parallel_machine_rank(pm) == 0) {
    m_values.resize(numRanks, 0);

    for (int proc = 0; proc < stk::parallel_machine_size(pm); ++proc) {
      while (comm.recv_buffer(proc).remaining()) {
        int rank;
        double value;
        comm.recv_buffer(proc).unpack<int>(rank);
        comm.recv_buffer(proc).unpack<double>(value);
        STK_ThrowRequire(rank < numRanks);
        m_values[rank] = value;
      }
    }
  }
}

void
DoubleWithPercentDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    m_min = compute_min(m_values);
    m_max = compute_max(m_values);
    m_avg = compute_avg(m_values);

    const double maxPercent = (m_avg > 0.0) ? 100. * (m_max - m_avg) / m_avg
                                            : 0.0;
    const double minPercent = (m_avg > 0.0) ? 100. * (m_min - m_avg) / m_avg
                                            : 0.0;

    std::ostringstream osMax;
    osMax << std::fixed << std::setprecision(1) << std::showpos << maxPercent;

    std::ostringstream osMin;
    osMin << std::fixed << std::setprecision(1) << std::showpos << minPercent;

    m_percentSize = std::max(osMax.str().size(), osMin.str().size());
  }
}

double
DoubleWithPercentDiagnostic::compute_min(const std::vector<double> & data)
{
  return *std::min_element(data.begin(), data.end());
}

double
DoubleWithPercentDiagnostic::compute_max(const std::vector<double> & data)
{
  return *std::max_element(data.begin(), data.end());
}

double
DoubleWithPercentDiagnostic::compute_avg(const std::vector<double> & data)
{
  if (data.empty()) {
    return 0.0;
  }
  return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

std::string
DoubleWithPercentDiagnostic::print_value_with_percent(double value)
{
  const double percent = (m_avg > 0.0) ? 100. * (value - m_avg) / m_avg
                                       : 0.0;
  std::ostringstream os;
  os << std::fixed << std::setprecision(m_decimalOutput) << value
     << " (" << std::setw(m_percentSize) << std::fixed << std::setprecision(1) << std::showpos << percent << "%)";

  return os.str();
}

std::string
DoubleWithPercentDiagnostic::print_rank_value(unsigned , int rank)
{
  return print_value_with_percent(m_values[rank]);
}

std::string
DoubleWithPercentDiagnostic::print_min(unsigned )
{
  return print_value_with_percent(m_min);
}

std::string
DoubleWithPercentDiagnostic::print_max(unsigned )
{
  return print_value_with_percent(m_max);
}

std::string
DoubleWithPercentDiagnostic::print_avg(unsigned )
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(m_decimalOutput) << m_avg << std::string(m_percentSize+4, ' ');

  return os.str();
}


MultiDoubleDiagnostic::MultiDoubleDiagnostic(unsigned numColumns, unsigned decimalOutput)
  : m_localValues(numColumns),
    m_values(numColumns),
    m_min(numColumns),
    m_max(numColumns),
    m_avg(numColumns),
    m_numColumns(numColumns),
    m_decimalOutput(decimalOutput)
{
}

void
MultiDoubleDiagnostic::collect_data(stk::ParallelMachine pm, int numRanks)
{
  stk::CommSparse comm(pm);
  stk::pack_and_communicate(comm, [this, &comm]()
  {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      for (const auto & localValue : m_localValues[column]) {
        comm.send_buffer(0).pack<unsigned>(column);
        comm.send_buffer(0).pack<int>(localValue.first);
        comm.send_buffer(0).pack<double>(localValue.second);
      }
    }
  });

  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_values[column].resize(numRanks, 0);
    }

    for (int proc = 0; proc < stk::parallel_machine_size(pm); ++proc) {
      while (comm.recv_buffer(proc).remaining()) {
        unsigned column;
        int rank;
        double value;
        comm.recv_buffer(proc).unpack<unsigned>(column);
        comm.recv_buffer(proc).unpack<int>(rank);
        comm.recv_buffer(proc).unpack<double>(value);
        STK_ThrowRequire(rank < numRanks);
        m_values[column][rank] = value;
      }
    }
  }
}

void
MultiDoubleDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_min[column] = compute_min(m_values[column]);
      m_max[column] = compute_max(m_values[column]);
      m_avg[column] = compute_avg(m_values[column]);
    }
  }
}

std::string
MultiDoubleDiagnostic::print_value(double value)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(m_decimalOutput) << value;

  return os.str();
}

std::string
MultiDoubleDiagnostic::print_rank_value(unsigned column, int rank)
{
  return print_value(m_values[column][rank]);
}

std::string
MultiDoubleDiagnostic::print_min(unsigned column)
{
  return print_value(m_min[column]);
}

std::string
MultiDoubleDiagnostic::print_max(unsigned column)
{
  return print_value(m_max[column]);
}

std::string
MultiDoubleDiagnostic::print_avg(unsigned column)
{
  return print_value(m_avg[column]);
}

double
MultiDoubleDiagnostic::compute_min(const std::vector<double> & data)
{
  return *std::min_element(data.begin(), data.end());
}

double
MultiDoubleDiagnostic::compute_max(const std::vector<double> & data)
{
  return *std::max_element(data.begin(), data.end());
}

double
MultiDoubleDiagnostic::compute_avg(const std::vector<double> & data)
{
  if (data.empty()) {
    return 0.0;
  }
  return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}


void
MultiDoubleWithPercentDiagnostic::process_data(stk::ParallelMachine pm)
{
  if (stk::parallel_machine_rank(pm) == 0) {
    for (unsigned column = 0; column < m_numColumns; ++column) {
      m_min[column] = compute_min(m_values[column]);
      m_max[column] = compute_max(m_values[column]);
      m_avg[column] = compute_avg(m_values[column]);

      const double maxPercent = (m_avg[column] > 0.0) ? 100.0 * (m_max[column] - m_avg[column]) / m_avg[column]
                                                      : 0.0;
      const double minPercent = (m_avg[column] > 0.0) ? 100.0 * (m_min[column] - m_avg[column]) / m_avg[column]
                                                      : 0.0;

      std::ostringstream osMax;
      osMax << std::fixed << std::setprecision(1) << std::showpos << maxPercent;

      std::ostringstream osMin;
      osMin << std::fixed << std::setprecision(1) << std::showpos << minPercent;

      m_percentSize[column] = std::max(osMax.str().size(), osMin.str().size());
    }
  }
}

std::string
MultiDoubleWithPercentDiagnostic::print_value_with_percent(unsigned column, double value)
{
  const double percent = (m_avg[column] > 0.0) ? 100.0 * (value - m_avg[column]) / m_avg[column]
                                               : 0.0;
  std::ostringstream os;
  os << std::fixed << std::setprecision(m_decimalOutput) << value
     << " (" << std::setw(m_percentSize[column]) << std::fixed << std::setprecision(1)
     << std::showpos << percent << "%)";

  return os.str();
}

std::string
MultiDoubleWithPercentDiagnostic::print_rank_value(unsigned column, int rank)
{
  return print_value_with_percent(column, m_values[column][rank]);
}

std::string
MultiDoubleWithPercentDiagnostic::print_min(unsigned column)
{
  return print_value_with_percent(column, m_min[column]);
}

std::string
MultiDoubleWithPercentDiagnostic::print_max(unsigned column)
{
  return print_value_with_percent(column, m_max[column]);
}

std::string
MultiDoubleWithPercentDiagnostic::print_avg(unsigned column)
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(m_decimalOutput) << m_avg[column]
     << std::string(m_percentSize[column]+4, ' ');

  return os.str();
}


void set_up_diagnostics(const stk::balance::BalanceSettings & balanceSettings)
{
  if ((impl::g_diagnosticsContainer.size() >= 2) &&
      (get_diagnostic<TotalElementWeightDiagnostic>()->num_columns() !=
       static_cast<unsigned>(balanceSettings.getNumCriteria()))) {
    impl::g_diagnosticsContainer.clear();
  }

  if (impl::g_diagnosticsContainer.size() == 0) {
    register_diagnostic<ElementCountDiagnostic>();
    register_diagnostic<TotalElementWeightDiagnostic>(static_cast<unsigned>(balanceSettings.getNumCriteria()));
    if (not balanceSettings.get_is_rebalancing()) {
      register_diagnostic<RelativeNodeInterfaceSizeDiagnostic>();
      register_diagnostic<ConnectivityWeightDiagnostic>();
    }
  }
}




} // namespace balance
} // namespace stk
