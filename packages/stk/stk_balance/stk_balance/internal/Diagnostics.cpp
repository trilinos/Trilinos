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
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/CommSparse.hpp"
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
        ThrowRequire(rank < numRanks);
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

    const double maxPercent = 100. * (static_cast<int>(m_max) - static_cast<int>(m_avg)) / static_cast<double>(m_avg);
    const double minPercent = 100. * (static_cast<int>(m_min) - static_cast<int>(m_avg)) / static_cast<double>(m_avg);

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
  const double percent = 100. * (static_cast<int>(value) - static_cast<int>(m_avg)) / static_cast<double>(m_avg);
  std::ostringstream os;
  os << value << " (" << std::setw(m_percentSize) << std::fixed << std::setprecision(1) << std::showpos
     << percent << "%)";

  return os.str();
}

std::string
UnsignedWithPercentDiagnostic::print_rank_value(int rank)
{
  return print_value_with_percent(m_values[rank]);
}

std::string
UnsignedWithPercentDiagnostic::print_min()
{
  return print_value_with_percent(m_min);
}

std::string
UnsignedWithPercentDiagnostic::print_max()
{
  return print_value_with_percent(m_max);
}

std::string
UnsignedWithPercentDiagnostic::print_avg()
{
  std::ostringstream os;
  os << m_avg << std::string(m_percentSize+4, ' ');

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
        ThrowRequire(rank < numRanks);
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

    const double maxPercent = 100. * (m_max - m_avg) / m_avg;
    const double minPercent = 100. * (m_min - m_avg) / m_avg;

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
  const double percent = 100. * (value - m_avg) / m_avg;
  std::ostringstream os;
  os << std::fixed << std::setprecision(3) << value
     << " (" << std::setw(m_percentSize) << std::fixed << std::setprecision(1) << std::showpos << percent << "%)";

  return os.str();
}

std::string
DoubleWithPercentDiagnostic::print_rank_value(int rank)
{
  return print_value_with_percent(m_values[rank]);
}

std::string
DoubleWithPercentDiagnostic::print_min()
{
  return print_value_with_percent(m_min);
}

std::string
DoubleWithPercentDiagnostic::print_max()
{
  return print_value_with_percent(m_max);
}

std::string
DoubleWithPercentDiagnostic::print_avg()
{
  std::ostringstream os;
  os << std::fixed << std::setprecision(3) << m_avg << std::string(m_percentSize+4, ' ');

  return os.str();
}


} // namespace balance
} // namespace stk
