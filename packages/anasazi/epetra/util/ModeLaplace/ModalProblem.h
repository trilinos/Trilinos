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

#ifndef ANASAZI_MODAL_PROBLEM_H
#define ANASAZI_MODAL_PROBLEM_H

#include "Anasaziepetra_ModeLaplace_DLLExportMacro.h"

class Epetra_MultiVector;

class ANASAZIEPETRA_MODELAPLACE_LIB_DLL_EXPORT ModalProblem {

  public:

    virtual ~ModalProblem() { }

    virtual const Epetra_CrsMatrix* getStiffness() const = 0;
    virtual const Epetra_CrsMatrix* getMass()      const = 0;

    virtual int eigenCheck(const Epetra_MultiVector &Q, double *lambda,
                           double *normWeight, bool smallest = 0) const { 
      (void)Q;
      (void)lambda;
      (void)normWeight;
      (void)smallest;
      return 0; 
    };

    virtual void memoryInfo() const { };
    virtual void problemInfo() const { };

};

#endif
