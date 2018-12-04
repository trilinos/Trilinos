
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

#ifndef BDDC_KRYLOVSOLVER_H
#define BDDC_KRYLOVSOLVER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

#include "shylu_PreconditionerBase.hpp"
#include "shylu_OrthogGCR.hpp"
#include "shylu_UtilBDDC.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
enum KrylovMethod {GCR, PCG};

template <class SX, class SM, class LO, class GO> 
  class KrylovSolver
{
public:
  KrylovSolver()
  {
  }

  KrylovSolver
    (PreconditionerBase<SX,SM,LO,GO>* Preconditioner,
     Teuchos::ParameterList Parameters)
  {
    m_Preconditioner = rcp(Preconditioner);
    m_Parameters = rcp( new Teuchos::ParameterList(Parameters) );
    initialize();
  }

  KrylovSolver
    (RCP< PreconditionerBase<SX,SM,LO,GO> > Preconditioner,
     RCP<Teuchos::ParameterList> Parameters) :
  m_Preconditioner(Preconditioner),
    m_Parameters(Parameters)
  {
    initialize();
  }
   
  ~KrylovSolver() = default;

  double getProjectionTime() const {
    if (m_Orthog == Teuchos::null) return 0;
    return m_Orthog->getProjectionTime();
  }

  double getOrthogonalizationTime() const {
    if (m_Orthog == Teuchos::null) return 0;
    return m_Orthog->getOrthogonalizationTime();
  }

  int getNumIterations()
  {
    return m_numIter;
  }

  SM getCurrentResidual()
  {
    return m_currentNorm;
  }

  int Solve(SX* rightHandSide,
	    SX* solution)
  {
    if (m_printFlag > 0) {
      m_outputFileKrylov.open("krylov_solver.dat",
			      std::ios::out | std::ios::app);
    }
    double startTime = GetTime();
    m_numIter = 0;
    int solverStatus = SolveEquations(rightHandSide,
				      solution,
				      m_numIter);
    double endTime = GetTime();
    double deltaTime = endTime - startTime;
    m_totalTime += deltaTime;
    m_numIterTotal += m_numIter;
    m_numSolvesTotal++;
    if (m_printFlag > 0) {
      m_outputFileKrylov << "elapsed time for this solve (seconds) = "
			 << deltaTime << std::endl;
      m_outputFileKrylov.close();
      OutputSolverData(deltaTime, solverStatus, m_numIter);
    }
    return solverStatus;
  }

  SX getMaxEigenvalue() const
  {
    return m_maxEigenvalue;
  }

 private: // data
  RCP< PreconditionerBase<SX,SM,LO,GO> > m_Preconditioner;
  RCP<Teuchos::ParameterList> m_Parameters;
  bool m_interfacePreconditioner, m_calculateRelativeResidualActual,
    m_avoidInitialStaticCondensation;
  LO m_vectorLength, m_krylovLength, m_maxIter,
    m_myPID, m_printFlag, m_numIter, m_numIterTotal, m_numSolvesTotal,
    m_convergenceCheckInterval, m_krylovMethod;
  bool m_estimateConditionNumber{false};
  double m_solverTol{0}, m_totalTime{0}, m_relativeResidual{0};
  double m_relativeResidualActual{0}, m_maxEigenvalue{0};
  SM m_originalNorm{0}, m_currentNorm{0};
  std::vector<SX> 
    m_residualVector, m_deltaSol, m_Pr, m_APr, m_p, m_Ap, 
    m_rhoArray, m_betaArray, m_pApArray;
  std::ofstream m_outputFileKrylov, m_outputFileDD;
  Teuchos::RCP< OrthogGCR<SX,SM,LO,GO> > m_Orthog;
  Teuchos::BLAS<int, SX> m_BLAS;

 private: // methods
  void initialize()
  {
    initializeParameters();
    if ((m_printFlag > 1) || 
	m_Parameters->get("Calculate Actual Residual", false)) {
      m_calculateRelativeResidualActual = true;
    }
    if (m_myPID != 0) m_printFlag = 0;
    reserveMemory();
    if (m_krylovMethod == GCR) {
      m_Orthog = Teuchos::rcp( new OrthogGCR<SX,SM,LO,GO>
			       (m_Preconditioner,
				m_Parameters,
				m_krylovLength) );
    }
  }

  void initializeParameters()
  {
    m_interfacePreconditioner = 
      m_Parameters->get("Interface Preconditioner", false);
    m_calculateRelativeResidualActual = false;
    m_avoidInitialStaticCondensation = 
      m_Parameters->get("Avoid Initial Static Condensation", false);
    m_vectorLength = m_Preconditioner->NumMyRows();
    m_krylovLength = m_Preconditioner->NumMyRowsKrylov();
    m_maxIter = m_Parameters->get("Maximum Iterations", 200);
    m_myPID = m_Preconditioner->MyPID();
    m_printFlag = m_Parameters->get("Print Summary", 0);
    m_numIter = m_numIterTotal = m_numSolvesTotal = 0;
    m_convergenceCheckInterval = 
      m_Parameters->get("Convergence Check Interval", 50);
    m_krylovMethod = m_Parameters->get("Krylov Method", 0);
    m_estimateConditionNumber = 
      m_Parameters->get("Estimate Condition Number", false);
    m_solverTol = m_Parameters->get("Convergence Tolerance", 1e-6);
    m_Orthog = Teuchos::null;
  }

  void reserveMemory() 
  {
    m_residualVector.resize(m_vectorLength, 0);
    m_Pr.resize(m_krylovLength, 0); // preconditioned residual
    m_APr.resize(m_krylovLength); // opertor time preconditioned residual
    m_deltaSol.resize(m_vectorLength);
    if (m_krylovMethod == PCG) {
      m_p.resize(m_krylovLength);
      m_Ap.resize(m_krylovLength);
      if (m_estimateConditionNumber) {
	const int length = m_maxIter + 1;
	m_rhoArray.resize(length);
	m_betaArray.resize(length);
	m_pApArray.resize(length);
      }
    }
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
    int valuesChanged(0);
    if (m_avoidInitialStaticCondensation == false) {
      valuesChanged  =
	m_Preconditioner->InitialStaticCondensation(&m_residualVector[0],
						    solution);
    }
    if (m_krylovMethod == GCR) {
      valuesChanged += m_Orthog->ProjectionCorrection(&m_residualVector[0],
						      &m_deltaSol[0],
						      -1, 1, 1, 1);
    }
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
      //      m_Preconditioner->StaticExpansion(&m_deltaSol[0], solution); 
      expandSolution(&m_deltaSol[0], solution); 
      if (m_calculateRelativeResidualActual) {
	m_currentNorm = CalculateResidualNorm(rightHandSide, solution);
	m_relativeResidualActual = m_currentNorm/m_originalNorm;
      }
      return 0;
    }
    switch(m_krylovMethod) {
    case PCG:
      solverStatus = IteratePCG(numIter);
      break;
    case GCR:
      solverStatus = IterateGCR(numIter);
      break;
    default:
      solverStatus = IterateGCR(numIter);
      if (m_myPID == 0) {
	std::cout << "Krylov method defaulted to GCR\n";
      }
      break;
    }
    //    m_Preconditioner->StaticExpansion(&m_deltaSol[0], solution);
    expandSolution(&m_deltaSol[0], solution);
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

  void expandSolution(SX* deltaSolution, 
		      SX* solution)
  {
    if (m_interfacePreconditioner) {
      const bool interfaceValuesOnlyInput(true);
      const bool useCurrentInteriorValues(false);
      m_Preconditioner->StaticExpansion
	(&m_deltaSol[0], interfaceValuesOnlyInput, useCurrentInteriorValues,
	 solution);
    }
    else {
      for (LO i=0; i<m_vectorLength; i++) {
	solution[i] += deltaSolution[i];
      }
    }
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
    SX* r   = m_residualVector.data();
    SX* Pr  = m_Pr.data();
    SX* APr = m_APr.data();
    numIter = 0;
    // matrix output for diagnostic purposes
    outputMatrices(Pr, APr);
    //
    // GCR iterations
    //
    double previousRelativeResidual(0);
    for (int iter=0; iter<m_maxIter; iter++) {
      m_Preconditioner->Apply(r, Pr, APr, m_interfacePreconditioner);
      m_Orthog->StoreSearchDirection(Pr, APr);
      m_Orthog->ApplyProjection(r, m_deltaSol.data());
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

  int IteratePCG(LO & numIter)
  {
    SX* u = m_deltaSol.data();
    SX* r = m_residualVector.data();
    SX* Pr = m_Pr.data();
    SX* APr = m_APr.data();
    SX* p = m_p.data();
    SX* Ap = m_Ap.data();
    SX rPrevDotzPrev(0), beta(0);
    numIter = 0;
    // matrix output for diagnostic purposes
    outputMatrices(Pr, APr);
    //
    // PCG iterations
    //
    double previousRelativeResidual(0);
    for (int iter=0; iter<m_maxIter; iter++) {
      m_Preconditioner->Apply(r, Pr, APr, m_interfacePreconditioner);
      SX rDotz = m_Preconditioner->DotProd(r, Pr, m_krylovLength);
      if (iter == 0) beta = 0;
      else beta = rDotz/rPrevDotzPrev;
      for (int i=0; i<m_krylovLength; i++) {
	p[i ] = Pr[i]  + beta*p[i];
	Ap[i] = APr[i] + beta*Ap[i];
      }
      rPrevDotzPrev = rDotz;
      SX pAp = m_Preconditioner->DotProd(p, Ap, m_krylovLength);
      SX alpha = rDotz/pAp;
      for (int i=0; i<m_krylovLength; i++) {
	u[i] += alpha*p[i];
	r[i] -= alpha*Ap[i];
      }
      if (m_estimateConditionNumber) {
	m_rhoArray[iter] = std::sqrt(std::abs(rDotz));
	m_betaArray[iter] = beta;
	m_pApArray[iter] = pAp;
      }
      m_currentNorm = m_Preconditioner->Norm2(r, m_krylovLength);
      m_relativeResidual = m_currentNorm/m_originalNorm;
      numIter = iter + 1;
      if (m_printFlag >= 3) {
	m_outputFileKrylov << "PCG iteration " << (iter%m_maxIter)+1 
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
    estimateConditionNumber(numIter);
    if (m_relativeResidual <= m_solverTol) return 0;
    else return -1;
  }

  void outputMatrices(SX* Pr,
		      SX* APr)
  {
    const bool printOperator = 
      m_Parameters->get("Print Operator Matrix", false);
    const bool printPrecond = 
      m_Parameters->get("Print Preconditioner Matrix", false);
    const int numProc = m_Preconditioner->numProc();
    if ((printOperator || printPrecond) && (numProc == 1)) {
      const int N = m_krylovLength;
      std::vector<SX> A(N*N), Minv(N*N), AMinv(N*N), rr(N);
      for (int i=0; i<N; i++) {
	rr.assign(N,0);
	rr[i] = 1;
	if (printOperator) {
	  m_Preconditioner->ApplyFullOperator(rr.data(), Pr);
	  memcpy(&A[i*N], Pr, N*sizeof(SX));
	}
	if (printPrecond) {
	  m_Preconditioner->Apply(rr.data(), Pr, APr, m_interfacePreconditioner);
	  memcpy(&Minv[i*N], Pr, N*sizeof(SX));
	  memcpy(&AMinv[i*N], APr, N*sizeof(SX));
	}
      }
      UtilBDDC<SX,SM>::printDenseMatrix(N, N, A.data(), "A.dat");
      UtilBDDC<SX,SM>::printDenseMatrix(N, N, Minv.data(), "Minv.dat");
      UtilBDDC<SX,SM>::printDenseMatrix(N, N, AMinv.data(), "AMinv.dat");
    }
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
      m_outputFileDD << "KrylovSolver: iterative solver failed to converge" 
		     << std::endl;
    }
    m_outputFileDD.close();
  }

  void calculateEigenvalues(int & numRows, 
			    SX* D, 
			    SX* E, 
			    int & INFO)
  {
    char COMPZ = 'N'; // calculate eigenvalues only
    SX *Z(0), *WORK(0);
    int LDZ(1);
    Teuchos::LAPACK<int, SX> LAPACK;
    LAPACK.STEQR(COMPZ, numRows, D, E, Z, LDZ, WORK, &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "STEQR error");
  }

  void estimateConditionNumber(const int numIter)
  {
    if ((m_estimateConditionNumber == false) || (m_myPID != 0)) return;
    const int length = numIter + 1;
    std::vector<SX> DtriArray(length), EtriArray(length), econArray(length);
    if (numIter == 1) {
      econArray[0] = 1;
    }
    else {
      for (int i=0; i<numIter; i++) {
	int numRows = i + 1;
	DtriArray[0] = m_pApArray[0]/m_rhoArray[0]/m_rhoArray[0];
	for (int j=1; j<numRows; j++) {
	  DtriArray[j] = (m_pApArray[j-1]*m_betaArray[j]*m_betaArray[j] +
			  m_pApArray[j])/m_rhoArray[j]/m_rhoArray[j];
	  EtriArray[j-1] =
	    -m_pApArray[j-1]*m_betaArray[j]/m_rhoArray[j-1]/m_rhoArray[j];
	}
	int INFO(0);
	calculateEigenvalues(numRows, DtriArray.data(), EtriArray.data(), INFO);
	if (INFO != 0) {
	  std::cout << "error in KrylovSolver::estimateConditionNumber\n"; 
	  std::cout << "i, numRows, INFO = " << i << " " << numRows << " "
		    << INFO << std::endl;
	  continue;
	}
	econArray[i] = DtriArray[i]/DtriArray[0];
      }
    }
    std::cout <<"eigenvalue estimates = " << std::endl;
    for (int i=0; i<numIter; i++) {
      std::cout << econArray[i] << std::endl;
    }
  }

};

} // namespace bddc

#endif // BDDC_KRYLOVSOLVER_H
  
