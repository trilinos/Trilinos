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
#ifndef DIAGNOSTICSPRINTER_HPP
#define DIAGNOSTICSPRINTER_HPP

#include "stk_util/parallel/Parallel.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cstddef>

namespace stk {
namespace balance {

class DiagnosticsPrinter {
public:
  DiagnosticsPrinter(stk::ParallelMachine comm, int numLogicalRanks);
  ~DiagnosticsPrinter() = default;

  void print(std::ostream & stream);

private:
  std::string print_rank_header1() { return "MPI"; }
  std::string print_rank_header2() { return "Rank"; }
  std::string print_rank(int rank) { return std::to_string(rank); }
  std::string print_summary_header1() { return ""; }
  std::string print_summary_header2() { return "Qty"; }
  std::string print_summary_min() { return "Min"; }
  std::string print_summary_max() { return "Max"; }
  std::string print_summary_avg() { return "Avg"; }
  std::vector<size_t> required_column_header_widths();
  std::vector<size_t> required_column_data_widths();
  std::vector<size_t> required_column_total_widths(const std::vector<size_t> & columnWidths,
                                                   const std::vector<size_t> & dataWidths);

  void print_horizontal_divider(std::ostream & stream, size_t totalWidth);

  void print_data_header1_row(std::ostream & stream, const std::vector<size_t> & columnWidths);
  void print_data_header2_row(std::ostream & stream, const std::vector<size_t> & columnWidths);
  void print_data_rank_rows(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                   const std::vector<size_t> & columnWidths);

  void print_summary_header1_row(std::ostream & stream, const std::vector<size_t> & columnWidths);
  void print_summary_header2_row(std::ostream & stream, const std::vector<size_t> & columnWidths);
  void print_summary_min_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                    const std::vector<size_t> & columnWidths);
  void print_summary_max_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                    const std::vector<size_t> & columnWidths);
  void print_summary_avg_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                    const std::vector<size_t> & columnWidths);

  stk::ParallelMachine m_comm;
  int m_numLogicalRanks;
};

} // namespace balance
} // namespace stk

#endif // DIAGNOSTICSPRINTER_HPP
