
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

#ifndef BDDC_BANDEDSOLVER_H
#define BDDC_BANDEDSOLVER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "Teuchos_LAPACK.hpp"
#include "shylu_errorBDDC.hpp"

namespace bddc {

template <class SX> class BandedSolver 
{
public:

  BandedSolver(const int numRows,
	       const int* rowBegin, 
	       const int* columns, 
	       const SX* values) :
    m_numRows(numRows),
    m_rowBegin(rowBegin),
    m_columns(columns),
    m_values(values)
  {
  }

  ~BandedSolver()
  {
  }

  int initialize()
  {
    if (m_numRows == 0) return 0;
    determineRKMOrdering(m_numRows, m_rowBegin, m_columns, m_perm);
    int INFO;
    factorMatrix(m_numRows, m_rowBegin, m_columns, m_values, INFO);
    return INFO;
  }

  int solve(int NRHS,
	    SX* Rhs, 
	    SX* Sol)
  {
    if (m_numRows == 0) return 0;
    Teuchos::LAPACK<int, SX> LAPACK;
    int INFO(0);
    m_work.resize(m_numRows*NRHS);
    m_work.shrink_to_fit();
    if (m_usingDense) {
      memcpy(m_work.data(), Rhs, m_numRows*NRHS*sizeof(SX));
      LAPACK.GETRS('N', m_numRows, NRHS, m_A.data(), m_numRows, m_IPIV.data(), 
		   m_work.data(), m_numRows, &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRS error");
      memcpy(Sol, m_work.data(), m_numRows*NRHS*sizeof(SX));
    }
    else {
      // permute rows of rhs
      for (int j=0; j<NRHS; j++) {
	SX* rhsSource = &Rhs[m_numRows*j];
	SX* rhs = &m_work[m_numRows*j];
	for (int i=0; i<m_numRows; i++) {
	  const int row = m_perm[i];
	  rhs[i] = rhsSource[row];
	}
      }
      LAPACK.GBTRS('N', m_numRows, m_KL, m_KU, NRHS, m_A.data(), m_LDA, 
		   m_IPIV.data(), m_work.data(), m_numRows, &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GBTRS error");
      // permute rows of sol
      for (int j=0; j<NRHS; j++) {
	SX* sol = &m_work[m_numRows*j];
	SX* solSource = &Sol[m_numRows*j];
	for (int i=0; i<m_numRows; i++) {
	  const int row = m_perm[i];
	  solSource[row] = sol[i];
	}
      }
    }
    return INFO;
  }

private:
  int m_numRows{0}, m_LDA{0}, m_KL{0}, m_KU{0};
  std::vector<SX> m_A, m_work;
  std::vector<int> m_perm, m_IPIV;
  const int *m_rowBegin{nullptr}, *m_columns{nullptr};
  const SX *m_values{nullptr};
  bool m_usingDense{false};

  void factorMatrix(const int numRows, 
		    const int* rowBegin, 
		    const int* columns, 
		    const SX* values,
		    int & INFO)
  {
    std::vector<int> invPerm(numRows);
    for (int i=0; i<numRows; i++) invPerm[m_perm[i]] = i;
    int KL(0), KU(0);
    for (int row=0; row<numRows; row++) {
      const int rowSource = m_perm[row];
      for (int j=rowBegin[rowSource]; j<rowBegin[rowSource+1]; j++) {
	const int colSource = columns[j];
	const int col = invPerm[colSource];
	const int belowDiag = row - col;
	if (belowDiag > KL) KL = belowDiag;
	const int aboveDiag = col - row;
	if (aboveDiag > KU) KU = aboveDiag;
      }
    }
    m_IPIV.resize(numRows);
    Teuchos::LAPACK<int, SX> LAPACK;
    const int LDA = 2*KL + KU + 1;
    if (LDA > numRows) { // simply use dense LAPACK solver to save memory
      m_usingDense = true;
      m_A.assign(numRows*numRows, 0);
      for (int i=0; i<numRows; i++) {
	for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	  const int col = columns[j];
	  m_A[i+col*numRows] = values[j];
	}
      }
      LAPACK.GETRF(numRows, numRows, m_A.data(), numRows, m_IPIV.data(), 
		   &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GETRF error");
    }
    else { // use banded LAPACK solver
      m_usingDense = false;
      const int numTerms = LDA*numRows;
      m_A.assign(numTerms, 0);
      for (int rowSource=0; rowSource<numRows; rowSource++) {
	const int row = invPerm[rowSource];
	for (int j=rowBegin[rowSource]; j<rowBegin[rowSource+1]; j++) {
	  const int colSource = columns[j];
	  const int col = invPerm[colSource];
	  const int localRow = KL + KU + row - col;
	  BDDC_TEST_FOR_EXCEPTION(localRow < KL, std::runtime_error, 
				  "localRow too small");
	  BDDC_TEST_FOR_EXCEPTION(localRow >= LDA, std::runtime_error, 
				  "localRow too large");
	  m_A[localRow+LDA*col] = values[j];
	}
      }
      LAPACK.GBTRF(numRows, numRows, KL, KU, m_A.data(), LDA, m_IPIV.data(), 
		   &INFO);
      BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GBTRF error");
      m_KL = KL;
      m_KU = KU;
      m_LDA = LDA;
    }
  }

  void getRowsToAddToQueue(const int parent,
			   const int* rowBegin,
			   const int* columns,
			   const std::vector<int> & degree,
			   const std::vector<int> & perm,
			   const std::vector<bool> & inQueue,
			   std::vector< std::pair<int,int> > & rowsToAddToQueue,
			   int & numAdjacent)
  {
    numAdjacent = 0;
    for (int i=rowBegin[parent]; i<rowBegin[parent+1]; i++) {
      const int col = columns[i];
      if (inQueue[col] == false) {
	if ((col != parent) && (perm[col] == -1)) {
	  rowsToAddToQueue[numAdjacent].first = degree[col];
	  rowsToAddToQueue[numAdjacent++].second = col;
	}
      }
    }
    std::sort(rowsToAddToQueue.begin(), rowsToAddToQueue.begin()+numAdjacent);
  }

  void determineDegreeAndSort(const int numRows,
			      const int* rowBegin,
			      int & maxDegree,
			      std::vector<int> & degree, 
			      std::vector<int> & sortedDegree)
  {
    degree.resize(numRows);
    std::vector< std::pair<int,int> > degreePair(numRows);
    maxDegree = -1;
    for (int i=0; i<numRows; i++) {
      degree[i] = rowBegin[i+1] - rowBegin[i];
      degreePair[i].first = degree[i];
      degreePair[i].second = i;
      if (degree[i] > maxDegree) maxDegree = degree[i];
    }
    std::sort(degreePair.begin(), degreePair.end());
    sortedDegree.resize(numRows);
    for (int i=0; i<numRows; i++) {
      sortedDegree[i] = degreePair[i].second;
    }
  }

  int getRowWithMinDegree(const std::vector<int> & sortedDegree, 
			  const std::vector<int> & perm, 
			  int & startLooking)
  {
    bool foundRow(false);
    int row(-1);
    while (foundRow == false) {
      const int sortRow = sortedDegree[startLooking];
      if (perm[sortRow] == -1) {
	row = sortRow;
	foundRow = true;
      }
      startLooking++;
    }
    BDDC_TEST_FOR_EXCEPTION(row == -1, std::runtime_error, "invalid row");
    return row;
  }

  void determineRKMOrdering(const int numRows, 
			    const int* rowBegin, 
			    const int* columns,
			    std::vector<int> & perm)
  {
    int maxDegree;
    std::vector<int> degree, sortedDegree;
    determineDegreeAndSort(numRows, rowBegin, maxDegree, degree, sortedDegree);
    std::vector< std::pair<int, int> > rowsToAddToQueue(maxDegree);
    perm.assign(numRows, -1);
    std::vector<int> queue(numRows), results(numRows);
    std::vector<bool> inQueue(numRows, false);
    int startLooking(0), rEnd(0);
    while (rEnd < numRows) {
      int parent = getRowWithMinDegree(sortedDegree, perm, startLooking);
      perm[parent] = rEnd;
      results[rEnd++] = parent;
      int numAdjacent;
      getRowsToAddToQueue(parent, rowBegin, columns, degree, perm,
			  inQueue, rowsToAddToQueue, numAdjacent);
      for (int i=0; i<numAdjacent; i++) {
	queue[i] = rowsToAddToQueue[i].second;
	inQueue[queue[i]] = true;
      }
      int qStart(0), qEnd(numAdjacent);
      while (qEnd > qStart) {
	const int child = queue[qStart];
	BDDC_TEST_FOR_EXCEPTION(perm[child] != -1, std::runtime_error, 
				"perm[child] should be -1");
	perm[child] = rEnd;
	results[rEnd++] = child;
	qStart++;
	getRowsToAddToQueue(child, rowBegin, columns, degree, perm,
			    inQueue, rowsToAddToQueue, numAdjacent);
	for (int i=0; i<numAdjacent; i++) {
	  queue[qEnd+i] = rowsToAddToQueue[i].second;
	  inQueue[rowsToAddToQueue[i].second] = true;
	}
	qEnd += numAdjacent;
      }
    }
    // reverse order
    for (int i=0; i<numRows; i++) perm[i] = results[numRows-1-i];
  }

};
  
} // namespace bddc

#endif // BDDC_BANDEDSOLVER_H
  
