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
#include "DiagnosticsPrinter.hpp"
#include "DiagnosticsContainer.hpp"
#include "Diagnostics.hpp"

namespace stk {
namespace balance {

const std::string ColumnDivider = " | ";

DiagnosticsPrinter::DiagnosticsPrinter(stk::ParallelMachine comm, int numLogicalRanks)
  : m_comm(comm),
    m_numLogicalRanks(numLogicalRanks)
{
}

std::vector<size_t>
DiagnosticsPrinter::required_column_header_widths()
{
  std::vector<size_t> headerWidths;

  headerWidths.push_back(std::max({print_rank_header1().size(),
                                   print_rank_header2().size()}));

  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      headerWidths.push_back(std::max({diag->print_header1(column).size(),
                                       diag->print_header2(column).size()}));
    }
  }

  return headerWidths;
}

std::vector<size_t>
DiagnosticsPrinter::required_column_data_widths()
{
  std::vector<size_t> dataWidths;

  dataWidths.push_back(print_rank(m_numLogicalRanks-1).size());

  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      dataWidths.push_back(std::max({diag->print_min(column).size(),
                                     diag->print_max(column).size()}));
    }
  }

  return dataWidths;
}

std::vector<size_t>
DiagnosticsPrinter::required_column_total_widths(const std::vector<size_t> & headerWidths,
                                                 const std::vector<size_t> & dataWidths)
{
  std::vector<size_t> totalWidths(headerWidths.size(), 0);

  for (size_t i = 0; i < headerWidths.size(); ++i) {
    totalWidths[i] = std::max(headerWidths[i], dataWidths[i]);
  }

  return totalWidths;
}

std::string centered_in_width(const std::string & str, size_t width)
{
  size_t strLen = str.length();
  if (width <= strLen) { return str; }

  size_t diff = width - strLen;
  size_t padRight = diff/2;
  size_t padLeft = diff - padRight;
  return std::string(padLeft, ' ') + str + std::string(padRight, ' ');
}

std::string centered_data_in_width(const std::string & str, size_t dataWidth, size_t totalWidth)
{
  size_t strLen = str.length();
  if (totalWidth <= strLen) { return str; }
  size_t dataPad = (dataWidth > strLen) ? dataWidth - strLen : 0;
  std::string rightJustifiedData = std::string(dataPad, ' ') + str;

  if (totalWidth <= dataWidth) { return rightJustifiedData; }
  size_t totalPad = totalWidth - dataWidth;
  size_t padRight = totalPad / 2;
  size_t padLeft = totalPad - padRight;
  return std::string(padLeft, ' ') + rightJustifiedData + std::string(padRight, ' ');
}

size_t total_output_width(const std::vector<size_t> & columnWidths)
{
  const size_t dividerWidth = ColumnDivider.size();
  size_t totalWidth = columnWidths[0] + 2;
  for (size_t i = 1; i < columnWidths.size(); ++i) {
    totalWidth += dividerWidth + columnWidths[i];
  }
  return totalWidth;
}

void DiagnosticsPrinter::print_data_header1_row(std::ostream & stream, const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_rank_header1(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      stream << ColumnDivider << centered_in_width(diag->print_header1(column), columnWidths[++index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_data_header2_row(std::ostream & stream, const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_rank_header2(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      stream << ColumnDivider << centered_in_width(diag->print_header2(column), columnWidths[++index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_summary_header1_row(std::ostream & stream, const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_summary_header1(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      stream << ColumnDivider << centered_in_width(diag->print_header1(column), columnWidths[++index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_summary_header2_row(std::ostream & stream, const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_summary_header2(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      stream << ColumnDivider << centered_in_width(diag->print_header2(column), columnWidths[++index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_summary_min_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                                      const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_summary_min(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      ++index;
      stream << ColumnDivider << centered_data_in_width(diag->print_min(column),
                                                        dataWidths[index], columnWidths[index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_summary_max_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                                      const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_summary_max(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      ++index;
      stream << ColumnDivider << centered_data_in_width(diag->print_max(column),
                                                        dataWidths[index], columnWidths[index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_summary_avg_row(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                                      const std::vector<size_t> & columnWidths)
{
  stream << " " << centered_in_width(print_summary_avg(), columnWidths[0]);
  int index = 0;
  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    for (unsigned column = 0; column < diag->num_columns(); ++column) {
      ++index;
      stream << ColumnDivider << centered_data_in_width(diag->print_avg(column),
                                                        dataWidths[index], columnWidths[index]);
    }
  }
  stream << " " << std::endl;
}

void DiagnosticsPrinter::print_data_rank_rows(std::ostream & stream, const std::vector<size_t> & dataWidths,
                                                                     const std::vector<size_t> & columnWidths)
{
  for (int rank = 0; rank < m_numLogicalRanks; ++rank) {
    stream << " " << centered_data_in_width(print_rank(rank), dataWidths[0], columnWidths[0]);
    int index = 0;
    for (Diagnostic * diag : impl::g_diagnosticsContainer) {
      for (unsigned column = 0; column < diag->num_columns(); ++column) {
        ++index;
        stream << ColumnDivider << centered_data_in_width(diag->print_rank_value(column, rank),
                                                          dataWidths[index], columnWidths[index]);
      }
    }
    stream << " " << std::endl;
  }
}

void DiagnosticsPrinter::print_horizontal_divider(std::ostream & stream, size_t totalWidth)
{
  stream << std::string(totalWidth, '=') << std::endl;
}

void
DiagnosticsPrinter::print(std::ostream & stream)
{
  if (impl::g_diagnosticsContainer.size() == 0) {
    return;
  }

  for (Diagnostic * diag : impl::g_diagnosticsContainer) {
    diag->collect_data(m_comm, m_numLogicalRanks);
    diag->process_data(m_comm);
  }

  if (stk::parallel_machine_rank(m_comm) == 0) {
    std::vector<size_t> headerWidths = required_column_header_widths();
    std::vector<size_t> dataWidths = required_column_data_widths();
    std::vector<size_t> columnWidths = required_column_total_widths(headerWidths, dataWidths);
    const size_t totalWidth = total_output_width(columnWidths);

    stream << std::endl;
    print_horizontal_divider(stream, totalWidth);
    print_data_header1_row(stream, columnWidths);
    print_data_header2_row(stream, columnWidths);
    print_horizontal_divider(stream, totalWidth);
    print_data_rank_rows(stream, dataWidths, columnWidths);
    print_horizontal_divider(stream, totalWidth);

    stream << std::endl;
    print_horizontal_divider(stream, totalWidth);
    print_summary_header1_row(stream, columnWidths);
    print_summary_header2_row(stream, columnWidths);
    print_horizontal_divider(stream, totalWidth);
    print_summary_min_row(stream, dataWidths, columnWidths);
    print_summary_max_row(stream, dataWidths, columnWidths);
    print_summary_avg_row(stream, dataWidths, columnWidths);
    print_horizontal_divider(stream, totalWidth);
  }
}

} // namespace balance
} // namespace stk
