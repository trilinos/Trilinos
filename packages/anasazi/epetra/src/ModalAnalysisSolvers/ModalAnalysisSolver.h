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

#ifndef MODAL_ANALYSIS_SOLVER_H
#define MODAL_ANALYSIS_SOLVER_H

class ModalAnalysisSolver {

  public:

    virtual ~ModalAnalysisSolver()                { }

    virtual int solve(int numEigen, Epetra_MultiVector &Q, double *lambda) = 0;

    virtual int reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV=0)=0;

    virtual int minimumSpaceDimension(int nev) const  = 0;

    virtual void initializeCounters() { }

    virtual void algorithmInfo() const = 0;
    virtual void historyInfo() const { }
    virtual void memoryInfo() const { }
    virtual void operationInfo() const  { }
    virtual void timeInfo() const    { }

};


#endif
