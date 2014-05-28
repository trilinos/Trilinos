/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/diag/Writer.hpp>
#include <iostream>

namespace stk {

void
PrintTable::transpose_table() const
{}


void
PrintTable::calculate_column_widths() const
{
  ColumnWidthVector	column_width_set;

  // loop over the headers and find the longest field for each column by size and from m_columnWidth
  for (Table::const_iterator row_it = m_header.begin(); row_it != m_header.end(); ++row_it) {
    if ((*row_it).size() > m_columnWidth.size())
      m_columnWidth.resize((*row_it).size(), 0);
    if ((*row_it).size() > column_width_set.size())
      column_width_set.resize((*row_it).size(), 0);

    int i = 0;
    for (Row::const_iterator cell_it = (*row_it).begin(); cell_it != (*row_it).end(); ++cell_it, ++i) {
      m_columnWidth[i] = std::max(m_columnWidth[i], (*cell_it).m_string.size());
      column_width_set[i] = std::max(column_width_set[i], (*cell_it).m_width);
    }
  }

  // loop over the table and find the longest field for each column by size and from m_columnWidth
  for (Table::const_iterator row_it = m_table.begin(); row_it != m_table.end(); ++row_it) {
    if ((*row_it).size() > m_columnWidth.size())
      m_columnWidth.resize((*row_it).size(), 0);
    if ((*row_it).size() > column_width_set.size())
      column_width_set.resize((*row_it).size(), 0);

    int i = 0;
    for (Row::const_iterator cell_it = (*row_it).begin(); cell_it != (*row_it).end(); ++cell_it, ++i) {
      m_columnWidth[i] = std::max(m_columnWidth[i], (*cell_it).m_string.size());
      column_width_set[i] = std::max(column_width_set[i], (*cell_it).m_width);
    }
  }

  // choose m_width width over size() width for each column
  m_tableWidth = 0;
  for (ColumnWidthVector::size_type i = 0; i < m_columnWidth.size(); ++i) {
    if (column_width_set[i] != 0)
      m_columnWidth[i] = column_width_set[i];
    m_tableWidth += m_columnWidth[i] + 1;
  }
}


PrintTable &
PrintTable::at(
  size_t        row,
  size_t        col)
{
  for (Table::size_type i = m_table.size(); i <= row; ++i)
    m_table.push_back(Row());
  for (Row::size_type i = m_table[row].size(); i <= col; ++i) 
    m_table[row].push_back(Cell());
  
  m_currentCell.m_string = std::string(m_currentCell.m_indent*2, ' ') + m_currentString.str();  
  m_table[row][col] = m_currentCell;
  
  m_currentCell = Cell();
  m_currentString.str("");

  return *this;
}


PrintTable &
PrintTable::end_col()
{
  m_currentCell.m_string = std::string(m_currentCell.m_indent*2, ' ') + m_currentString.str();
  m_table.back().push_back(m_currentCell);
  if (m_table.size() > 1 && m_table[0].size() <= m_table.back().size()) {
    m_currentCell.m_string = "";
    m_currentCell.m_flags = 0;
    m_currentCell.m_justification = m_table[0][m_table[0].size() - 1].m_justification;
    m_currentCell.m_width = m_table[0][m_table[0].size() - 1].m_width;
    m_currentCell.m_indent = m_table[0][m_table[0].size() - 1].m_indent;
 }
  else {
    m_currentCell = Cell();
  }
  m_currentString.str("");

  return *this;
}


PrintTable &
PrintTable::end_row()
{
  if (!m_currentString.str().empty())
    end_col();
  m_table.push_back(Row());
  return *this;
}


std::ostream &
PrintTable::print(
  std::ostream &	os) const
{
  if (m_flags & COMMA_SEPARATED_VALUES)
    csvPrint(os);

  else {
    if (m_flags & PRINT_TRANSPOSED)
      transpose_table();

    calculate_column_widths();

    if (!m_title.empty()) {
      int prespaces = 0;

      if(m_title.length() < m_tableWidth)
	prespaces = (m_tableWidth - m_title.length())/2;;

      os << m_commentPrefix;
      os << std::left << std::setw(prespaces) << "" << m_title << '\n';
    }

    for (Table::const_iterator row_it = m_header.begin(); row_it != m_header.end(); ++row_it) {
      os << m_commentPrefix;
      printRow(os, *row_it);
      os << '\n';
    }

    if (m_header.size() > 0) {
      os << m_commentPrefix;
      printHeaderBar(os);
      os << '\n';
    }

    for (Table::const_iterator row_it = m_table.begin(); row_it != m_table.end(); ++row_it) {
      os << std::left << std::setw(m_commentPrefix.size()) << "";
      printRow(os, *row_it);
      os << '\n';
    }
  }

  return os;
}


std::ostream &
PrintTable::printRow(
  std::ostream &	os,
  const Row &		row) const
{
  int i = 0;
  int postspaces = 0;
  for (Row::const_iterator cell_it = row.begin(); cell_it != row.end(); ++cell_it, ++i) {
    os // << postspaces << ", "
       << std::left << std::setw(postspaces) << "";
    postspaces = 0;

    if (cell_it != row.begin())
      os << " ";

    if ((*cell_it).m_flags & Cell::SPAN)
      os << (*cell_it).m_string;
    else if ((*cell_it).m_string.length() > m_columnWidth[i]) {
      if ((*cell_it).m_justification & Cell::ENDS) {
	int front_end = m_columnWidth[i]/4;
	int back_begin = (*cell_it).m_string.size() - (m_columnWidth[i] - front_end);
	os << (*cell_it).m_string.substr(0, front_end - 3) + "..." + (*cell_it).m_string.substr(back_begin, (*cell_it).m_string.size());
      }
      else { // if ((*cell_it).m_justification & Cell::TRUNC) {
	os << (*cell_it).m_string.substr(0, m_columnWidth[i]);
      }
    }
    else {
      if ((*cell_it).m_string.length() == 0)
	postspaces = m_columnWidth[i];
      else if (((*cell_it).m_justification & Cell::JUSTIFY_MASK) == Cell::LEFT) {
	postspaces = m_columnWidth[i] - (*cell_it).m_string.length();
	os // << m_columnWidth[i] << ", " << postspaces << ", "
	  << std::left << (*cell_it).m_string;
      }
      else if (((*cell_it).m_justification & Cell::JUSTIFY_MASK) == Cell::CENTER) {
	int prespaces = (m_columnWidth[i] - (*cell_it).m_string.length())/2;
	postspaces = m_columnWidth[i] - (*cell_it).m_string.length() - prespaces;
	os // << prespaces << " " << postspaces << ", "
	  << std::left << std::setw(prespaces) << "" << (*cell_it).m_string;
      }
      else // if (((*cell_it).m_justification & Cell::JUSTIFY_MASK) == Cell::RIGHT)
	os // << m_columnWidth[i] << ", "
	   << std::right << std::setw(m_columnWidth[i]) << (*cell_it).m_string;
    }
  }

  return os;
}


std::ostream &
PrintTable::printHeaderBar(
  std::ostream &	os) const
{
  os << std::setfill('-');

  for (ColumnWidthVector::size_type i = 0; i < m_columnWidth.size(); ++i) {
    if (i != 0)
      os << " ";
    os << std::setw(m_columnWidth[i]) << "";
  }
  os << std::setfill(' ');

  return os;
}


std::ostream &
PrintTable::csvPrint(
  std::ostream &	os) const
{
  if (!m_title.empty())
    os << m_title << '\n';

  for (Table::const_iterator row_it = m_header.begin(); row_it != m_header.end(); ++row_it) {
    const Row &row = (*row_it);
    for (Row::const_iterator cell_it = row.begin(); cell_it != row.end(); ++cell_it) {
      if (cell_it != row.begin())
	os << ",";
      os << (*cell_it).m_string;
    }
    os << '\n';
  }

  for (Table::const_iterator row_it = m_table.begin(); row_it != m_table.end(); ++row_it) {
    const Row &row = (*row_it);
    for (Row::const_iterator cell_it = row.begin(); cell_it != row.end(); ++cell_it) {
      if (cell_it != row.begin())
	os << ",";
      os << (*cell_it).m_string;
    }
    os << '\n';
  }

  return os;
}


diag::Writer &
PrintTable::verbose_print(
  diag::Writer &	dout) const
{
//   const ColumnWidthVector &column_width = calculate_column_widths();

//   for (Table::const_iterator row_it = m_header.begin(); row_it != m_header.end(); ++row_it) {
//     printRow(os, *row_it);
//     os << '\n';
//   }

//   if (m_header.size() > 0)
//     printHeaderBar(os);

//   for (Table::const_iterator row_it = m_table.begin(); row_it != m_table.end(); ++row_it) {
//     int i = 0;
//     for (Row::const_iterator cell_it = (*row_it).begin(); cell_it != (*row_it).end(); ++cell_it, ++i)
//       if ((*cell_it).m_flags & Cell::SPAN)
//	dout << (*cell_it).m_string;
//       else
//	dout << std::setw(column_width[i]) << (*cell_it).m_string;
//     dout << dendl;
//   }

  calculate_column_widths();

  dout << m_title << std::endl;

  for (Table::const_iterator row_it = m_header.begin(); row_it != m_header.end(); ++row_it) {
    dout << "";
    printRow(dout.getStream(), *row_it);
    dout << diag::dendl;
  }

  if (m_header.size() > 0) {
    dout << "";
    printHeaderBar(dout.getStream());
    dout << diag::dendl;
  }

  for (Table::const_iterator row_it = m_table.begin(); row_it != m_table.end(); ++row_it) {
    dout << "";
    printRow(dout.getStream(), *row_it);
    dout << diag::dendl;
  }

  return dout;
}

} // namespace stk
