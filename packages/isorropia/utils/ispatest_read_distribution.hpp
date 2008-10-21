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

#ifndef _ispatest_read_distribution_hpp_
#define _ispatest_read_distribution_hpp_

#include <Isorropia_ConfigDefs.hpp>

namespace ispatest {

/** Read a matrix distribution from a file. The file is expected to
  follow a very simple format where each line corresponds to a single
  nonzero of the matrix, and each line contains 3 numbers which specify
  a row number, a column-number, and a partition number. When this
  function returns, each of the three vectors should have the same length,
  which is the number of nonzeros.
*/
void read_distribution(const char* filename,
                       std::vector<int>& rows,
                       std::vector<int>& cols,
                       std::vector<int>& partitions);

}//namespace ispatest

#endif

