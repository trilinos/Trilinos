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

#ifndef MODAL_PROBLEM_H
#define MODAL_PROBLEM_H

class Epetra_MultiVector;

class ModalProblem {

  public:

    virtual ~ModalProblem() { }

    virtual const Epetra_Operator* getStiffness() const = 0;
    virtual const Epetra_Operator* getMass()      const = 0;

    virtual int eigenCheck(const Epetra_MultiVector &Q, double *lambda,
                           double *normWeight) const { return 0; };

    virtual void memoryInfo() const { };
    virtual void problemInfo() const { };

};

#endif
