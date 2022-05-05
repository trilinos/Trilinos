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
#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP

#include "stk_util/parallel/Parallel.hpp"
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <type_traits>
#include <numeric>
#include <memory>
#include <typeinfo>  // for type_info

namespace stk {
namespace balance {

class Diagnostic
{
public:
  Diagnostic() = default;
  virtual ~Diagnostic() = default;

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) = 0;
  virtual void process_data(stk::ParallelMachine comm) = 0;

  virtual std::string print_header1() = 0;
  virtual std::string print_header2() = 0;

  virtual std::string print_rank_value(int rank) = 0;

  virtual std::string print_min() = 0;
  virtual std::string print_max() = 0;
  virtual std::string print_avg() = 0;

protected:
  template <typename T> T compute_min(const std::vector<T> & data);
  template <typename T> T compute_max(const std::vector<T> & data);

  template <typename T>
  std::enable_if_t<std::is_integral<T>::value, T> compute_avg(const std::vector<T> & data);

  template <typename T>
  std::enable_if_t<std::is_floating_point<T>::value, T> compute_avg(const std::vector<T> & data);
};

class UnsignedDiagnostic : public Diagnostic
{
public:
  UnsignedDiagnostic() = default;
  virtual ~UnsignedDiagnostic() override = default;

  void store_value(int rank, unsigned value) { m_localValues[rank] = value; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(int rank) override { return std::to_string(m_values[rank]); }
  virtual std::string print_min() override { return std::to_string(m_min); }
  virtual std::string print_max() override { return std::to_string(m_max); }
  virtual std::string print_avg() override { return std::to_string(m_avg); }

protected:
  unsigned compute_min(const std::vector<unsigned> & data);
  unsigned compute_max(const std::vector<unsigned> & data);
  unsigned compute_avg(const std::vector<unsigned> & data);

  std::map<int, unsigned> m_localValues;
  std::vector<unsigned> m_values;
  unsigned m_min = 0;
  unsigned m_max = 0;
  unsigned m_avg = 0;
};

class UnsignedWithPercentDiagnostic : public UnsignedDiagnostic
{
public:
  UnsignedWithPercentDiagnostic() = default;
  virtual ~UnsignedWithPercentDiagnostic() override = default;

  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(int rank) override;
  virtual std::string print_min() override;
  virtual std::string print_max() override;
  virtual std::string print_avg() override;

protected:
  std::string print_value_with_percent(unsigned value);

  unsigned m_percentSize = 0;
};

class DoubleWithPercentDiagnostic : public Diagnostic
{
public:
  DoubleWithPercentDiagnostic() = default;
  virtual ~DoubleWithPercentDiagnostic() override = default;

  void store_value(int rank, double value) { m_localValues[rank] = value; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(int rank) override;
  virtual std::string print_min() override;
  virtual std::string print_max() override;
  virtual std::string print_avg() override;

protected:
  double compute_min(const std::vector<double> & data);
  double compute_max(const std::vector<double> & data);
  double compute_avg(const std::vector<double> & data);
  std::string print_value_with_percent(double value);

  std::map<int, double> m_localValues;
  std::vector<double> m_values;
  double m_min = 0;
  double m_max = 0;
  double m_avg = 0;
  double m_percentSize = 0;
};

} // namespace balance
} // namespace stk

#endif // DIAGNOSTICS_HPP
