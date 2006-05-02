// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_TableFormat.hpp"

using namespace Teuchos;


string TableFormat::thinline() const
{
  TeuchosOStringStream toss;
  for (int i=0; i<pageWidth_; i++)
    {
      toss << "-";
    }
  return toss.str();
}



string TableFormat::thickline() const
{
  TeuchosOStringStream toss;
  for (int i=0; i<pageWidth_; i++)
    {
      toss << "=";
    }
  return toss.str();
}


string TableFormat::blanks(int size) const
{
  TeuchosOStringStream toss;
  for (int i=0; i<size; i++)
    {
      toss << " ";
    }
  return toss.str();
}



int TableFormat
::computeRequiredColumnWidth(const string& name,
                             const TableColumn& column) const
{
  int rtn = name.length();
  
  for (int i=0; i<column.numRows(); i++)
    {
      int x = column.entry(i)->toString().length();
      rtn = max(rtn, x);
    }
  
  return rtn + columnSpacing_;
}




void TableFormat::writeRow(RefCountPtr<std::ostream>& out,
                           const Array<RefCountPtr<TableEntry> >& entries) const
{
  TEST_FOR_EXCEPT(entries.size() != columnWidths_.size() 
                  && columnWidths_.size() != 0);

  for (unsigned int i=0; i<entries.size(); i++)
    {
      int cw = defaultColumnWidth();
      if (columnWidths_.size() != 0) cw = columnWidths_[i];

      *out << std::left << std::setw(cw) << entries[i]->toString();
    }
  *out << std::endl;
}


void TableFormat::writeRow(RefCountPtr<std::ostream>& out,
                           int rowIndex,
                           const Array<TableColumn>& columns) const
{
  Array<RefCountPtr<TableEntry> > entries(columns.size());
  for (unsigned int i=0; i<columns.size(); i++)
    {
      entries[i] = columns[i].entry(rowIndex);
    }
  
  writeRow(out, entries);
}


void TableFormat::writeWholeTable(RefCountPtr<std::ostream>& out,
                                  const string& header,
                                  const Array<string>& columnNames,
                                  const Array<TableColumn>& columns) const
{
  /* compute the total width */
  int pageWidth = 0;
  for (unsigned int i=0; i<columnNames.size(); i++)
    {
      int cw = defaultColumnWidth();
      if (columnWidths_.size() != 0) cw = columnWidths_[i];
      pageWidth += cw;
    }
  setPageWidth(max(pageWidth_, pageWidth));
  
  /* write the header */
  *out << thickline() << endl;
  *out << endl;
  int numBlanks = (pageWidth_ - header.length())/2;
  *out << blanks(numBlanks) << header << endl;
  *out << endl;

  /* write the column titles */
  for (unsigned int i=0; i<columnNames.size(); i++)
    {
      int cw = defaultColumnWidth();
      if (columnWidths_.size() != 0) cw = columnWidths_[i];

      *out << std::left << std::setw(cw) << columnNames[i];
    }
  *out << endl;

  /* ensure that all columns have the same number of rows */
  int numRows = columns[0].numRows();
  for (unsigned int i=1; i<columns.size(); i++)
    {
      TEST_FOR_EXCEPTION(columns[i].numRows() != numRows, runtime_error,
                         "inconsistent column sizes");
    }
  
  /* write the table data */
  for (int i=0; i<numRows; i++)
    {
      if (i % lineInterval_ == 0) *out << std::left << thinline() << endl;   
      writeRow(out, i, columns);
    }
  
  /* write the footer */
  *out << thickline() << endl;
}



