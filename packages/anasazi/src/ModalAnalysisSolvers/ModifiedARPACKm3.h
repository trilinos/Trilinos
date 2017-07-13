// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#ifndef MODIFIED_ARPACK_MODE3_H
#define MODIFIED_ARPACK_MODE3_H

#include "Epetra_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#endif

#include "CheckingTools.h"
#include "FortranRoutines.h"
#include "ModalAnalysisSolver.h"
#include "MyMemory.h"
#include "ModalTools.h"
#include "SortingTools.h"

class ModifiedARPACKm3 : public ModalAnalysisSolver {

  private:

    const CheckingTools myVerify;
    const Epetra_BLAS callBLAS;
    const FortranRoutines callFortran;
    ModalTools modalTool;
    const SortingTools mySort;

    const Epetra_Comm &MyComm;
    const Epetra_Operator *K;
    const Epetra_Operator *M;
    const Epetra_Time MyWatch;

    double tolEigenSolve;
    int maxIterEigenSolve;

    double *normWeight;

    int verbose;

    int historyCount;
    double *resHistory;

    double memRequested;
    double highMem;

    int massOp;
    int numResidual;
    int orthoOp;
    int outerIter;
    int stifOp;

    double timeMassOp;
    double timeOuterLoop;
    double timePostProce;
    double timeResidual;
    double timeStifOp;

    // Don't define these functions
    ModifiedARPACKm3(const ModifiedARPACKm3 &ref);
    ModifiedARPACKm3& operator=(const ModifiedARPACKm3 &ref);

  public:

    ModifiedARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK, 
                        double _tol = 1.0e-08, int _maxIter = 100, int _verb = 0);

    ModifiedARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK, 
                        const Epetra_Operator *MM,
                        double _tol = 1.0e-08, int _maxIter = 100, int _verb = 0,
                        double *_weight = 0);

    ~ModifiedARPACKm3();

    int solve(int numEigen, Epetra_MultiVector &Q, double *lambda);

    int reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV = 0);

    int reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV,
                const Epetra_MultiVector *orthoVec);

    int minimumSpaceDimension(int nev) const         { return nev+1; }

    void initializeCounters();

    void algorithmInfo() const;
    void historyInfo() const;
    void memoryInfo() const;
    void operationInfo() const;
    void timeInfo() const;

};

#endif
