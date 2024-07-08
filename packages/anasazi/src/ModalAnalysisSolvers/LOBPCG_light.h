// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#ifndef LOBPCG_LIGHT_H
#define LOBPCG_LIGHT_H

#include "Epetra_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"

#include "CheckingTools.h"
#include "FortranRoutines.h"
#include "ModalAnalysisSolver.h"
#include "MyMemory.h"
#include "ModalTools.h"
#include "SortingTools.h"

class LOBPCG_light : public ModalAnalysisSolver {

  private:

    const CheckingTools myVerify;
    const Epetra_BLAS callBLAS;
    const FortranRoutines callFortran;
    ModalTools modalTool;
    const SortingTools mySort;

    const Epetra_Comm &MyComm;
    const Epetra_Operator *K;
    const Epetra_Operator *M;
    const Epetra_Operator *Prec;
    const Epetra_Time MyWatch;

    double tolEigenSolve;
    int maxIterEigenSolve;

    int blockSize;

    double *normWeight;

    int verbose;

    int historyCount;
    double *resHistory;

    double memRequested;
    double highMem;

    int massOp;
    int numRestart;
    int outerIter;
    int precOp;
    int residual;
    int stifOp;

    double timeLocalProj;
    double timeLocalSolve;
    double timeLocalUpdate;
    double timeMassOp;
    double timeNorm;
    double timeOrtho;
    double timeOuterLoop;
    double timePostProce;
    double timePrecOp;
    double timeResidual;
    double timeRestart;
    double timeStifOp;

    // Private functions
    void accuracyCheck(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
                       const Epetra_MultiVector *R, const Epetra_MultiVector *Q,
                       const Epetra_MultiVector *H, const Epetra_MultiVector *P) const;

    // Don't define these functions
    LOBPCG_light(const LOBPCG_light &ref);
    LOBPCG_light& operator=(const LOBPCG_light &ref);

  public:

    LOBPCG_light(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
              const Epetra_Operator *PP, int _blk,
              double _tol = 1.0e-08, int _maxIter = 100, int _verb = 0);

    LOBPCG_light(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
              const Epetra_Operator *MM, const Epetra_Operator *PP, int _blk,
              double _tol = 1.0e-08, int _maxIter = 100, int _verb = 0, double *_weight = 0);

    ~LOBPCG_light();

    int solve(int numEigen, Epetra_MultiVector &Q, double *lambda)
        { return LOBPCG_light::reSolve(numEigen, Q, lambda); }

    int reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV = 0);

    int minimumSpaceDimension(int nev) const         { return nev+blockSize; }

    void initializeCounters();

    void algorithmInfo() const;
    void historyInfo() const;
    void memoryInfo() const;
    void operationInfo() const;
    void timeInfo() const;

};

#endif
