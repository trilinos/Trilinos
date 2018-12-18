
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

#ifndef BDDC_SOLVERGCR_H
#define BDDC_SOLVERGCR_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

#include "shylu_PreconditionerBDDC.hpp"
#include "shylu_OrthogGCR.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class SolverGCR
{
public:
  SolverGCR()
  {
  }

  SolverGCR
    (RCP< PreconditionerBDDC<SX,SM,LO,GO> > Preconditioner,
     RCP<Teuchos::ParameterList> Parameters) :
  m_Preconditioner(Preconditioner),
    m_Parameters(Parameters),
    m_interfacePreconditioner
    (Parameters->get("Interface Preconditioner", false)),
    m_calculateRelativeResidualActual(false),
    m_vectorLength(Preconditioner->NumMyRows()),
    m_krylovLength(Preconditioner->NumMyRowsKrylov()),
    m_maxIter
    (Parameters->get("Maximum Iterations", 200)),
    m_myPID(Preconditioner->MyPID()),
    m_printFlag(m_Parameters->get("Print Summary", 0)),
    m_numIterTotal(0),
    m_numSolvesTotal(0),
    m_convergenceCheckInterval
    (Parameters->get("Convergence Check Interval", 10)),
    m_solverTol
    (Parameters->get("Convergence Tolerance", 1e-6)),
    m_totalTime(0),
    m_relativeResidual(0),
    m_relativeResidualActual(0),
    m_originalNorm(0),
    m_currentNorm(0)
  {
    if (m_printFlag > 1) m_calculateRelativeResidualActual = true;
    if (m_myPID != 0) m_printFlag = 0;
    reserveMemory();
    m_Orthog = Teuchos::rcp( new OrthogGCR<SX,SM,LO,GO>
			     (m_Preconditioner,
			      m_Parameters,
			      m_krylovLength) );
  }
   
  ~SolverGCR()
  {
  }

  double getProjectionTime() const {
    return m_Orthog->getProjectionTime();
  }

  double getOrthogonalizationTime() const {
    return m_Orthog->getOrthogonalizationTime();
  }

  int Solve(SX* rightHandSide,
	    SX* solution)
  {
    if (m_printFlag > 0) {
      m_outputFileKrylov.open("krylov_solver.dat",
			      std::ios::out | std::ios::app);
    }
    double startTime = GetTime();
    LO numIter(0);
    int solverStatus = SolveEquations(rightHandSide,
				      solution,
				      numIter);
    double endTime = GetTime();
    double deltaTime = endTime - startTime;
    m_totalTime += deltaTime;
    m_numIterTotal += numIter;
    m_numSolvesTotal++;
    if (m_printFlag > 0) {
      m_outputFileKrylov << "elapsed time for this solve (seconds) = "
			 << deltaTime << std::endl;
      m_outputFileKrylov.close();
      OutputSolverData(deltaTime, solverStatus, numIter);
    }
    return solverStatus;
  }

 private: // data
  RCP< PreconditionerBDDC<SX,SM,LO,GO> > m_Preconditioner;
  RCP<Teuchos::ParameterList> m_Parameters;
  bool m_interfacePreconditioner, m_calculateRelativeResidualActual;
  LO m_vectorLength, m_krylovLength, m_maxIter, 
    m_myPID, m_printFlag, m_numIterTotal, m_numSolvesTotal,
    m_convergenceCheckInterval;
  double m_solverTol, m_totalTime, m_relativeResidual,
    m_relativeResidualActual;
  SM m_originalNorm, m_currentNorm;
  std::vector<SX> m_residualVector, m_deltaSol, m_Pr, m_APr;
  std::ofstream m_outputFileKrylov, m_outputFileDD;
  Teuchos::RCP< OrthogGCR<SX,SM,LO,GO> > m_Orthog;
  Teuchos::BLAS<int, SX> m_BLAS;

 private: // methods
  void reserveMemory() 
  {
    m_residualVector.resize(m_vectorLength, 0);
    m_Pr.resize(m_krylovLength, 0); // preconditioned residual
    m_APr.resize(m_krylovLength); // opertor time preconditioned residual
    m_deltaSol.resize(m_vectorLength);
  }
  
  double GetTime()
  {
    return m_Preconditioner->GetTime();
  }

  int SolveEquations(SX* rightHandSide,
		     SX* solution,
		     LO & numIter)
  {
    int solverStatus(-1);
    m_originalNorm = m_Preconditioner->Norm2(rightHandSide, m_vectorLength);
    if (m_printFlag > 1) {
      m_outputFileKrylov << "initial residual = "
			 << m_originalNorm  << std::endl;
    }
    for (LO i=0; i<m_vectorLength; i++) {
      m_residualVector[i] = rightHandSide[i];
      solution[i] = 0;
    }
    for (LO i=0; i<m_krylovLength; i++) m_deltaSol[i] = 0;
    m_currentNorm = m_originalNorm;
    if (m_originalNorm == 0) {
      m_relativeResidual = 0;
      m_relativeResidualActual = 0;
      return 0;
    }
    //
    // initial corrections
    //
    int valuesChanged  =
      m_Preconditioner->InitialStaticCondensation(&m_residualVector[0],
						  solution);
    valuesChanged += m_Orthog->ProjectionCorrection(&m_residualVector[0],
						    &m_deltaSol[0],
						    -1, 1, 1, 1);
    if (valuesChanged) {
      m_currentNorm = 
	m_Preconditioner->Norm2(&m_residualVector[0], m_krylovLength);
      if (m_printFlag > 1) {
	/*
	  m_outputFileKrylov << "number of search directions used         = "
	  << m_Orthog->NumSearchDirectionsUsed() << std::endl;
	*/
	m_outputFileKrylov << "residual after using initial corrections = " 
			   << m_currentNorm << std::endl;
      }
    }
    m_relativeResidual = m_currentNorm/m_originalNorm;
    if (m_relativeResidual <= m_solverTol) {   
      m_Preconditioner->StaticExpansion(&m_deltaSol[0], solution); 
      if (m_calculateRelativeResidualActual) {
	m_currentNorm = CalculateResidualNorm(rightHandSide, solution);
	m_relativeResidualActual = m_currentNorm/m_originalNorm;
      }
      return 0;
    }
    solverStatus = IterateGCR(numIter);
    m_Preconditioner->StaticExpansion(&m_deltaSol[0], solution);
    if (m_calculateRelativeResidualActual) {
      m_currentNorm = CalculateResidualNorm(rightHandSide, solution);
      m_relativeResidualActual = m_currentNorm/m_originalNorm;
      if (m_printFlag > 1) {
	m_outputFileKrylov << "recursive final residual              = "
			   << m_relativeResidual*m_originalNorm << std::endl;
	m_outputFileKrylov << "actual final residual                 = "
			   << m_currentNorm << std::endl;
      }
    }
    return solverStatus;
  }

  SM CalculateResidualNorm(SX* rightHandSide,
			   SX* solution)
  {
    m_Preconditioner->ApplyFullOperator(solution, &m_residualVector[0]);
    for (LO i=0; i<m_vectorLength; i++) {
      m_residualVector[i] = rightHandSide[i] - m_residualVector[i];
    }
    return(m_Preconditioner->Norm2(&m_residualVector[0], m_vectorLength));
  }

  int IterateGCR(LO & numIter)
  {
    numIter = 0;
    //
    // GCR iterations
    //
    double previousRelativeResidual(0);
    for (int iter=0; iter<m_maxIter; iter++) {
      SX* r   = &m_residualVector[0];
      SX* Pr  = &m_Pr[0];
      SX* APr = &m_APr[0];
      m_Preconditioner->Apply(r, Pr, APr);
      m_Orthog->StoreSearchDirection(Pr, APr);
      m_Orthog->ApplyProjection(r, &m_deltaSol[0]);
      m_currentNorm = m_Preconditioner->Norm2(r, m_krylovLength);
      m_relativeResidual = m_currentNorm/m_originalNorm;
      cullSearchDirections(iter);
      numIter++;
      if (m_printFlag >= 3) {
	m_outputFileKrylov << "GCR iteration " << (iter%m_maxIter)+1 
			   << " of maxIter = " 
			   << m_maxIter << ": relative residual = "
			   << m_relativeResidual << std::endl;
      }
      if (m_relativeResidual <= m_solverTol) break;
      if (isConverging(m_relativeResidual, 
		       previousRelativeResidual, iter) == false) {
	return -1;
      }
    }
    m_Orthog->SetNumVectorsStart(m_Orthog->GetNumVectors());
    if (m_relativeResidual <= m_solverTol) return 0;
    else return -1;
  }

  bool isConverging(double relativeResidual,
		    double & previousRelativeResidual,
		    int iteration)
  {
    double epsilon(0.01);
    if (iteration == 0) {
      previousRelativeResidual = relativeResidual;
      return true;
    }
    else {
      if (iteration%m_convergenceCheckInterval == 0) {
	double ratio = relativeResidual/previousRelativeResidual;
	if ((ratio >= 1-epsilon) && (ratio <= 1+epsilon)) {
	  return false;
	}
	previousRelativeResidual = relativeResidual;
      }
    }
    return true;
  }

  void cullSearchDirections(int iter)
  {
    if (m_Orthog->maxSearchDirectionsUsed() == true) {
      int iterPlus1 = iter + 1;
      int avgNumIter = (m_numIterTotal+iterPlus1)/(m_numSolvesTotal+1);
      avgNumIter = std::max(avgNumIter, 2);
      m_Orthog->UpdateSearchSpace(avgNumIter);
    }
  }

  void OutputSolverData(double deltaTime,
			int solverStatus,
			int numIter)
  {
    m_outputFileDD.open("dd_solver.dat", std::ios::out | std::ios::app);
    if (m_numSolvesTotal == 1) {
      m_outputFileDD.width(36); m_outputFileDD << "Recursive";
      if (m_printFlag > 1) {
	m_outputFileDD.width(12); m_outputFileDD << "Actual";
      }
      m_outputFileDD << std::endl;
      m_outputFileDD.width(36); m_outputFileDD << "Relative";
      if (m_printFlag > 1) {
	m_outputFileDD.width(12); m_outputFileDD << "Relative";
      }
      m_outputFileDD << std::endl;
      m_outputFileDD.width(5); m_outputFileDD << "Solve";
      m_outputFileDD.width(6); m_outputFileDD << "Iter";
      m_outputFileDD.width(7); m_outputFileDD << "Total";
      m_outputFileDD.width(5); m_outputFileDD << "Avg";
      m_outputFileDD.width(13); m_outputFileDD << "Residual";
      if (m_printFlag > 1) {
	m_outputFileDD.width(12); m_outputFileDD << "Residual";
      }
      m_outputFileDD.width(12); m_outputFileDD << "CPU (s)";
      m_outputFileDD.width(12); m_outputFileDD << "Total (s)";
      m_outputFileDD.width(13); m_outputFileDD << "Avg (s)";
      m_outputFileDD << std::endl;
    }
    m_outputFileDD.width(5); m_outputFileDD << m_numSolvesTotal;
    m_outputFileDD.width(6); m_outputFileDD << numIter;
    m_outputFileDD.width(7); m_outputFileDD << m_numIterTotal;
    m_outputFileDD.width(5); m_outputFileDD << m_numIterTotal/m_numSolvesTotal;
    m_outputFileDD.width(13); m_outputFileDD << m_relativeResidual;
    if (m_printFlag > 1) {
      m_outputFileDD.width(12); m_outputFileDD << m_relativeResidualActual;
    }
    m_outputFileDD.width(12); m_outputFileDD << deltaTime;
    m_outputFileDD.width(12); m_outputFileDD << m_totalTime;
    m_outputFileDD.width(13); m_outputFileDD << m_totalTime/m_numSolvesTotal;
    m_outputFileDD << std::endl;
    if (solverStatus != 0) {
      m_outputFileDD << "SolverGCR: iterative solver failed to converge" 
		     << std::endl;
    }
    m_outputFileDD.close();
  }

};

} // namespace bddc

#endif // BDDC_SOLVERGCR_H
  
