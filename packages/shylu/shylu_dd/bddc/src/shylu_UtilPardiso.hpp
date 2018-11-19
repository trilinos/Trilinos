
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_UTILPARDISO_H
#define BDDC_UTILPARDISO_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <vector>

namespace bddc {
  
  template <class LO, class SX> 
  class UtilPardiso
{
public:
  UtilPardiso()
  {
  }

  static void constructPardisoMatrix(const LO numRows, 
				     const LO* rowBegin, 
				     const LO* columns, 
				     const SX* values,
				     const bool matrixIsSymmetric,
				     LO* & rowBeginP,
				     LO* & columnsP,
				     SX* & valuesP)
  {
    //
    // will switch to Fortran-based numbering for Pardiso at the very end
    // (i.e., start with 1 instead of 0)
    //
    rowBeginP = new LO[numRows+1]; rowBeginP[0] = 0;
    LO maxNumCols(0);
    for (LO i=0; i<numRows; i++) {
      LO numCols(0), dflag(0);
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	LO col = columns[j];
	if ((col >= i) || (!matrixIsSymmetric)) numCols++;
	if (col == i) dflag = 1;
      }
      if ((dflag == 0) && (matrixIsSymmetric)) numCols++;
      rowBeginP[i+1] = rowBeginP[i] + numCols;
      if (numCols > maxNumCols) maxNumCols = numCols;
    }
    LO nnz = rowBeginP[numRows];
    columnsP = new LO[nnz];
    valuesP = new SX[nnz];
    for (LO i=0; i<numRows; i++) {
      LO numCols(0), dflag(0);
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	LO col = columns[j];
	if ((col >= i) || (!matrixIsSymmetric)) {
	  LO index = rowBeginP[i] + numCols;
	  columnsP[index] = col;
	  valuesP[index] = values[j];
	  numCols++;
	}
	if (col == i) dflag = 1;
      }
      if ((dflag == 0) && (matrixIsSymmetric)) {
	LO index = rowBeginP[i] + numCols;
	columnsP[index] = i;
	valuesP[index] = 0;
      }
    }
    //
    // make sure columns in each row are in ascending order
    //
    std::vector<SX> rowValues(maxNumCols);
    for (LO i=0; i<numRows; i++) {
      LO numCols = rowBeginP[i+1] - rowBeginP[i];
      std::vector< std::pair<LO, LO> > pairColumns(numCols);
      for (LO j=rowBeginP[i]; j<rowBeginP[i+1]; j++) {
	LO index = j - rowBeginP[i];
	pairColumns[index] = std::make_pair(columnsP[j], index);
	rowValues[index] = valuesP[j];
      }
      std::sort(pairColumns.begin(), pairColumns.end());
      typename std::vector< std::pair<LO, LO> >::const_iterator iter;
      for (LO j=0; j<numCols; j++) {
	iter = pairColumns.begin() + j;
	LO index = j + rowBeginP[i];
	columnsP[index] = iter->first;
	valuesP[index] = rowValues[iter->second];
      }
    }
    //
    // switch to Fortran numbering
    //
    for (int i=0; i<rowBeginP[numRows]; i++) columnsP[i]++;
    for (int i=0; i<=numRows; i++) rowBeginP[i]++;
  }

private:

};

} // namespace bddc

#endif // BDDC_UTILPARDISO_H
  
