
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

#ifndef BDDC_ORTHOGGCR_H
#define BDDC_ORTHOGGCR_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

#include "shylu_PreconditionerBase.hpp"
#include "shylu_errorBDDC.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
  template <class SX, class SM, class LO, class GO> 
  class OrthogGCR
{
 public:
  OrthogGCR()
  {
  }

  OrthogGCR
    (RCP< PreconditionerBase<SX,SM,LO,GO> > Preconditioner,
     RCP<Teuchos::ParameterList> Parameters,
     LO vectorLength) :
  m_Preconditioner(Preconditioner),
    m_Parameters(Parameters),
    m_vectorLength(vectorLength),
    m_numOrthogSteps
    (Parameters->get("Number Orthogonalization Steps GCR", 2)),
    m_maxNumVectors(Parameters->get("Maximum Stored Directions", 0)),
    m_numVectors(0),
    m_numVectorsStart(0),
    m_projectionTime(0),
    m_orthogonalizationTime(0)
  {
    reserveMemory();
  }
   
  ~OrthogGCR()
  {
  }

  void ApplyProjection(SX* r,
		       SX* u)
  {
    if (m_numVectors > 0) {
      double startTime = m_Preconditioner->GetTime();
      LO nRow = m_vectorLength;
      int nCol = m_numVectors;
      SX ALPHA(1), BETA(0), dotProd(0), dotProdSum(0);
      int INCX(1), INCY(1);
      LO start = nRow*(nCol-1);
      if (nRow > 0) {
	m_BLAS.GEMV(Teuchos::CONJ_TRANS, nRow, 1, ALPHA, &m_QAPhi[start], 
		    nRow, r, INCX, BETA, &dotProd, INCY);
      }
      m_Preconditioner->ReduceSum(&dotProd, 1, &dotProdSum);
      SX scaleFacu =  dotProdSum;
      SX scaleFacr = -dotProdSum;
      if (nRow > 0) {
	m_BLAS.AXPY(nRow, scaleFacu, &m_Phi[start],  INCX, u, INCY);
	m_BLAS.AXPY(nRow, scaleFacr, &m_QAPhi[start], INCX, r, INCY);
      }
      m_projectionTime += m_Preconditioner->GetTime() - startTime;
    }
  }

  int ProjectionCorrection(SX* residualVector,
			   SX* solutionVector,
			   SX ALPHAr,
			   SX BETAr,
			   SX ALPHA,
			   SX BETA)
  {
    if (m_numVectorsStart > 0) {
      double startTime = m_Preconditioner->GetTime();
      LO nRow = m_vectorLength;
      int nCol = m_numVectorsStart;
      int INCX(1), INCY(1);
      SX ZERO(0);
      m_dotProd.assign(nCol, 0);
      if (nRow > 0) {
	m_BLAS.GEMV(Teuchos::CONJ_TRANS, nRow, nCol, ALPHA, &m_QAPhi[0], nRow, 
		    residualVector, INCX, ZERO, &m_dotProd[0], INCY);
      }
      m_Preconditioner->ReduceSum(&m_dotProd[0], nCol, &m_dotProdSum[0]);
      // update residual to account for correction
      if (nRow > 0) {
	m_BLAS.GEMV(Teuchos::NO_TRANS, nRow, nCol, ALPHAr, &m_QAPhi[0], nRow,
		    &m_dotProdSum[0], INCX, BETAr, residualVector, INCY);
      }
      // update solution
      if (nRow > 0) {
	m_BLAS.GEMV(Teuchos::NO_TRANS, nRow, nCol, ALPHA, &m_Phi[0], nRow,
		    &m_dotProdSum[0], INCX, BETA, solutionVector, INCY);
      }
      m_projectionTime += m_Preconditioner->GetTime() - startTime;
      return 1;
    }
    else return 0;
  }

  void StoreSearchDirection(SX* x,
			    SX* Ax)
  {
    if (m_maxNumVectors <= 0) return;
    if (m_numVectors+1 <= m_maxNumVectors) {
      double startTime = m_Preconditioner->GetTime();
      LO nRow = m_vectorLength;
      int nCol = m_numVectors;
      LO istart = nRow*nCol;
      for (LO i=0; i<nRow; i++) {
	m_Phi[istart+i] = x[i];
	m_QAPhi[istart+i] = Ax[i];
      }
      int INCX(1), INCY(1);
      if (nCol > 0) {
	m_dotProd.assign(nCol, 0);
	for (int step=0; step<m_numOrthogSteps; step++) {
	  SX ALPHA(1), BETA(0);
	  if (nRow > 0) {
	    m_BLAS.GEMV(Teuchos::CONJ_TRANS, nRow, nCol, ALPHA, 
			&m_QAPhi[0], nRow, &m_QAPhi[istart], INCX, 
			BETA, &m_dotProd[0], INCY);
	  }
	  m_Preconditioner->ReduceSum(&m_dotProd[0], nCol, &m_dotProdSum[0]);
	  ALPHA = -1; BETA = 1;
	  if (nRow > 0) {
	    m_BLAS.GEMV(Teuchos::NO_TRANS, nRow, nCol, ALPHA, 
			&m_QAPhi[0], nRow, &m_dotProdSum[0], INCX, 
			BETA, &m_QAPhi[istart], INCY);
	    m_BLAS.GEMV(Teuchos::NO_TRANS, nRow, nCol, ALPHA, 
			&m_Phi[0], nRow, &m_dotProdSum[0], INCX, 
			BETA, &m_Phi[istart], INCY);
	  }
	}
      }
      SM normQ = m_Preconditioner->Norm2(&m_QAPhi[istart], nRow);
      for (LO i=0; i<nRow; i++) {
	m_Phi[istart+i] /= normQ;
	m_QAPhi[istart+i] /= normQ;
      }
      m_numVectors++;
      m_orthogonalizationTime += m_Preconditioner->GetTime() - startTime;
    }
  }

  bool maxSearchDirectionsUsed() const
  {
    if (m_maxNumVectors == 0) return false;
    if (m_numVectors == m_maxNumVectors) {
      return true;
    }
    else {
      return false;
    }
  }

  void UpdateSearchSpace(int avgNumIter)
  {
    if (m_maxNumVectors == 0) return;
    BDDC_TEST_FOR_EXCEPTION(m_numVectors != m_maxNumVectors, 
			    std::runtime_error, "invalid m_numVectors");
    int numKeep = std::max(m_maxNumVectors/4,
			   m_maxNumVectors - 2*avgNumIter);
    m_numVectors = numKeep;
    m_numVectorsStart = numKeep;
  }

  int GetNumVectors() {
    return m_numVectors;
  }

  void SetNumVectorsStart(int numVectorsStart) {
    m_numVectorsStart = numVectorsStart;
  }

  double getProjectionTime() const {
    return m_projectionTime;
  }

  double getOrthogonalizationTime() const {
    return m_orthogonalizationTime;
  }

 private: // data
  RCP< PreconditionerBase<SX,SM,LO,GO> > m_Preconditioner;
  RCP<Teuchos::ParameterList> m_Parameters;
  LO m_vectorLength;
  int m_numOrthogSteps, m_maxNumVectors, m_numVectors, m_numVectorsStart;
  double m_projectionTime, m_orthogonalizationTime;
  std::vector<SX> m_Phi, m_QAPhi, m_dotProd, m_dotProdSum;
  Teuchos::BLAS<int, SX> m_BLAS;

 private: // methods

  void reserveMemory()
  {
    LO numRow = m_vectorLength;
    int numCol = m_maxNumVectors;
    m_Phi.resize(numRow*numCol);
    m_QAPhi.resize(numRow*numCol);
    m_dotProd.resize(numCol+1);
    m_dotProdSum.resize(numCol+1);
  }

};

} // namespace bddc

#endif // BDDC_ORTHOGGCR_H
  
