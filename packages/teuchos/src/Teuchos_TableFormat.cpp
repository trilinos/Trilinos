// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


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
  TEST_FOR_EXCEPT(entries.size() != columnWidths_.size() 
    && columnWidths_.size() != 0);

  for (Array<RCP<TableEntry> >::size_type i=0; i<entries.size(); i++)
  {
    int cw = defaultColumnWidth();
    if (columnWidths_.size() != 0) cw = columnWidths_[i];

    out << std::left << std::setw(cw) << entries[i]->toString();
  }
  out << std::endl;
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

}


} // namespace Teuchos
