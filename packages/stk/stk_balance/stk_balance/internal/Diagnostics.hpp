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

class BalanceSettings;

class Diagnostic
{
public:
  Diagnostic() = default;
  virtual ~Diagnostic() = default;

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) = 0;
  virtual void process_data(stk::ParallelMachine comm) = 0;

  virtual unsigned num_columns() { return 1; }

  virtual std::string print_header1(unsigned column) = 0;
  virtual std::string print_header2(unsigned column) = 0;

  virtual std::string print_rank_value(unsigned column, int rank) = 0;

  virtual std::string print_min(unsigned column) = 0;
  virtual std::string print_max(unsigned column) = 0;
  virtual std::string print_avg(unsigned column) = 0;
};

class UnsignedDiagnostic : public Diagnostic
{
public:
  UnsignedDiagnostic() = default;
  virtual ~UnsignedDiagnostic() override = default;

  void store_value(int rank, unsigned value) { m_localValues[rank] = value; }
  unsigned get_rank_value(int rank) { return m_values[rank]; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned , int rank) override { return std::to_string(m_values[rank]); }
  virtual std::string print_min(unsigned ) override { return std::to_string(m_min); }
  virtual std::string print_max(unsigned ) override { return std::to_string(m_max); }
  virtual std::string print_avg(unsigned ) override { return std::to_string(m_avg); }

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

  virtual std::string print_rank_value(unsigned , int rank) override;
  virtual std::string print_min(unsigned ) override;
  virtual std::string print_max(unsigned ) override;
  virtual std::string print_avg(unsigned ) override;

protected:
  std::string print_value_with_percent(unsigned value);

  unsigned m_percentSize = 0;
};

class MultiUnsignedDiagnostic : public Diagnostic
{
public:
  explicit MultiUnsignedDiagnostic(unsigned numColumns);
  virtual ~MultiUnsignedDiagnostic() override = default;

  void store_value(unsigned column, int rank, unsigned value) { m_localValues[column][rank] = value; }
  double get_rank_value(unsigned column, int rank) { return m_values[column][rank]; }

  virtual unsigned num_columns() override { return m_numColumns; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned column, int rank) override { return std::to_string(m_values[column][rank]); }
  virtual std::string print_min(unsigned column) override { return std::to_string(m_min[column]); }
  virtual std::string print_max(unsigned column) override { return std::to_string(m_max[column]); }
  virtual std::string print_avg(unsigned column) override { return std::to_string(m_avg[column]); }

protected:
  unsigned compute_min(const std::vector<unsigned> & data);
  unsigned compute_max(const std::vector<unsigned> & data);
  unsigned compute_avg(const std::vector<unsigned> & data);

  std::vector<std::map<int, unsigned>> m_localValues;
  std::vector<std::vector<unsigned>> m_values;
  std::vector<unsigned> m_min;
  std::vector<unsigned> m_max;
  std::vector<unsigned> m_avg;
  unsigned m_numColumns;
};

class MultiUnsignedWithPercentDiagnostic : public MultiUnsignedDiagnostic
{
public:
  MultiUnsignedWithPercentDiagnostic(unsigned numColumns)
    : MultiUnsignedDiagnostic(numColumns),
      m_percentSize(numColumns)
  {}
  virtual ~MultiUnsignedWithPercentDiagnostic() override = default;

  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned column, int rank) override;
  virtual std::string print_min(unsigned column) override;
  virtual std::string print_max(unsigned column) override;
  virtual std::string print_avg(unsigned column) override;

protected:
  std::string print_value_with_percent(unsigned column, unsigned value);

  std::vector<unsigned> m_percentSize;
};


class DoubleWithPercentDiagnostic : public Diagnostic
{
public:
  explicit DoubleWithPercentDiagnostic(unsigned decimalOutput = 3)
    : m_min(0),
      m_max(0),
      m_avg(0),
      m_percentSize(0),
      m_decimalOutput(decimalOutput)
  {}
  virtual ~DoubleWithPercentDiagnostic() override = default;

  void store_value(int rank, double value) { m_localValues[rank] = value; }
  double get_rank_value(int rank) { return m_values[rank]; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned , int rank) override;
  virtual std::string print_min(unsigned ) override;
  virtual std::string print_max(unsigned ) override;
  virtual std::string print_avg(unsigned ) override;

protected:
  double compute_min(const std::vector<double> & data);
  double compute_max(const std::vector<double> & data);
  double compute_avg(const std::vector<double> & data);
  std::string print_value_with_percent(double value);

  std::map<int, double> m_localValues;
  std::vector<double> m_values;
  double m_min;
  double m_max;
  double m_avg;
  double m_percentSize;
  unsigned m_decimalOutput;
};

class MultiDoubleDiagnostic : public Diagnostic
{
public:
  explicit MultiDoubleDiagnostic(unsigned numColumns, unsigned decimalOutput = 3);
  virtual ~MultiDoubleDiagnostic() override = default;

  void store_value(unsigned column, int rank, double value) { m_localValues[column][rank] = value; }
  double get_rank_value(unsigned column, int rank) { return m_values[column][rank]; }

  virtual unsigned num_columns() override { return m_numColumns; }

  virtual void collect_data(stk::ParallelMachine comm, int numRanks) override;
  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned column, int rank) override;
  virtual std::string print_min(unsigned column) override;
  virtual std::string print_max(unsigned column) override;
  virtual std::string print_avg(unsigned column) override;

protected:
  double compute_min(const std::vector<double> & data);
  double compute_max(const std::vector<double> & data);
  double compute_avg(const std::vector<double> & data);
  std::string print_value(double value);

  std::vector<std::map<int, double>> m_localValues;
  std::vector<std::vector<double>> m_values;
  std::vector<double> m_min;
  std::vector<double> m_max;
  std::vector<double> m_avg;
  unsigned m_numColumns;
  unsigned m_decimalOutput;
};


class MultiDoubleWithPercentDiagnostic : public MultiDoubleDiagnostic
{
public:
  explicit MultiDoubleWithPercentDiagnostic(unsigned numColumns, unsigned decimalOutput = 3)
    : MultiDoubleDiagnostic(numColumns, decimalOutput),
      m_percentSize(numColumns)
  {}
  virtual ~MultiDoubleWithPercentDiagnostic() override = default;

  virtual void process_data(stk::ParallelMachine comm) override;

  virtual std::string print_rank_value(unsigned column, int rank) override;
  virtual std::string print_min(unsigned column) override;
  virtual std::string print_max(unsigned column) override;
  virtual std::string print_avg(unsigned column) override;

protected:
  std::string print_value_with_percent(unsigned column, double value);

  std::vector<unsigned> m_percentSize;
};



class ElementCountDiagnostic : public UnsignedWithPercentDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Number of"; }
  virtual std::string print_header2(unsigned ) override { return "Elements"; }
};


class TotalElementWeightDiagnostic : public MultiUnsignedWithPercentDiagnostic
{
public:
  TotalElementWeightDiagnostic(unsigned numColumns)
    : MultiUnsignedWithPercentDiagnostic(numColumns)
  {}
  virtual std::string print_header1(unsigned ) override { return "Total Element"; }
  virtual std::string print_header2(unsigned column) override {
    if (m_numColumns == 1) {
      return "Weight";
    }
    else {
      return "Weight " + std::to_string(column);
    }
  }
};


class RelativeNodeInterfaceSizeDiagnostic : public DoubleWithPercentDiagnostic
{
public:
  virtual std::string print_header1(unsigned ) override { return "Relative Node"; }
  virtual std::string print_header2(unsigned ) override { return "Interface Size"; }
};


class ConnectivityWeightDiagnostic : public DoubleWithPercentDiagnostic
{
public:
  ConnectivityWeightDiagnostic()
    : DoubleWithPercentDiagnostic(0)
  {}

  virtual std::string print_header1(unsigned ) override { return "Connectivity"; }
  virtual std::string print_header2(unsigned ) override { return "Weight"; }
};


void set_up_diagnostics(const stk::balance::BalanceSettings & balanceSettings);

} // namespace balance
} // namespace stk

#endif // DIAGNOSTICS_HPP
