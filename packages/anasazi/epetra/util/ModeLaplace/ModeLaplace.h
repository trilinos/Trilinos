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

#ifndef ANASAZI_MODE_LAPLACE_H
#define ANASAZI_MODE_LAPLACE_H

#include "ModalProblem.h"
#include "Anasaziepetra_ModeLaplace_DLLExportMacro.h"

class ANASAZIEPETRA_MODELAPLACE_LIB_DLL_EXPORT ModeLaplace : public ModalProblem {

  public:

    virtual ~ModeLaplace() { }

    virtual double getFirstMassEigenValue() const  = 0;

};

#endif
