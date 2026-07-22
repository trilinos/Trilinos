// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include "Teuchos_TableFormat.hpp"
#include "Teuchos_Assert.hpp"


namespace Teuchos {


std::string TableFormat::thinline() const
{
  std::ostringstream toss;
  for (int i=0; i<pageWidth_; i++)
  {
    toss << "-";
  }
  return toss.str();
}



std::string TableFormat::thickline() const
{
  std::ostringstream toss;
  for (int i=0; i<pageWidth_; i++)
  {
    toss << "=";
  }
  return toss.str();
}


std::string TableFormat::blanks(int size) const
{
  std::ostringstream toss;
  for (int i=0; i<size; i++)
  {
    toss << " ";
  }
  return toss.str();
}


int TableFormat::computeRequiredColumnWidth(
  const std::string& name,
  const TableColumn& column
  ) const
{
  int rtn = name.length();

  for (int i=0; i<column.numRows(); i++)
  {
    int x = column.entry(i)->toString().length();
    rtn = std::max(rtn, x);
  }

  return rtn + columnSpacing_;
}


void TableFormat::writeRow(
  std::ostream& out,
  const Array<RCP<TableEntry> >& entries
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(entries.size() != columnWidths_.size()
    && columnWidths_.size() != 0);

  std::ios::fmtflags f( out.flags() );
  for (Array<RCP<TableEntry> >::size_type i=0; i<entries.size(); i++)
  {
    int cw = defaultColumnWidth();
    if (columnWidths_.size() != 0) cw = columnWidths_[i];

    out << std::left << std::setw(cw) << entries[i]->toString();
  }
  out << std::endl;
  out.flags(f);
}


void TableFormat::writeRow(
  std::ostream& out,
  int rowIndex,
  const Array<TableColumn>& columns
  ) const
{
  Array<RCP<TableEntry> > entries(columns.size());
  for (Array<TableColumn>::size_type i=0; i<columns.size(); i++)
  {
    entries[i] = columns[i].entry(rowIndex);
  }

  writeRow(out, entries);
}


void TableFormat::writeWholeTable(
  std::ostream& out,
  const std::string& header,
  const Array<std::string>& columnNames,
  const Array<TableColumn>& columns
  ) const
{
  std::ios::fmtflags f(out.flags());

  /* compute the total width */
  int pgWidth = 0;
  for (Array<TableColumn>::size_type i=0; i<columnNames.size(); i++)
  {
    int cw = defaultColumnWidth();
    if (columnWidths_.size() != 0) cw = columnWidths_[i];
    pgWidth += cw;
  }
  setPageWidth(std::max(pageWidth_, pgWidth));

  /* write the header */
  out << thickline() << std::endl;
  out << std::endl;
  int numBlanks = (pageWidth_ - header.length())/2;
  out << blanks(numBlanks) << header << std::endl;
  out << std::endl;

  /* write the column titles */
  for (Array<std::string>::size_type i=0; i<columnNames.size(); i++)
  {
    int cw = defaultColumnWidth();
    if (columnWidths_.size() != 0) cw = columnWidths_[i];

    out << std::left << std::setw(cw) << columnNames[i];
  }
  out << std::endl;

  /* ensure that all columns have the same number of rows */
  int numRows = columns[0].numRows();
  for (Array<TableColumn>::size_type i=1; i<columns.size(); i++)
  {
    TEUCHOS_ASSERT_EQUALITY(columns[i].numRows(), numRows);
  }

  /* write the table data */
  for (int i=0; i<numRows; i++)
  {
    if (i % lineInterval_ == 0)
      out << std::left << thinline() << std::endl;
    writeRow(out, i, columns);
  }

  /* write the footer */
  out << thickline() << std::endl;

  // Restore flags
  out.flags(f);
}


} // namespace Teuchos
