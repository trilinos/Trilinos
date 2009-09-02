//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#include <ispatest_read_distribution.hpp>
#include <Isorropia_Exception.hpp>

/** ispatest is the namespace that contains isorropia's test-utilities
*/
namespace ispatest {

void read_distribution(const char* filename,
                       std::vector<int>& rows,
                       std::vector<int>& cols,
                       std::vector<int>& partitions)
{
  std::ifstream infile(filename);
  if (!infile) {
    throw Isorropia::Exception("ispatest::read_distribution ERROR opening file.");
  }

  rows.clear();
  cols.clear();
  partitions.clear();

  int row, col, partition;

  infile >> row;

  while(!infile.eof()) {
    infile >> col;
    infile >> partition;

    rows.push_back(row);
    cols.push_back(col);
    partitions.push_back(partition);

    infile >> row;
  }
}

}//namespace ispatest

