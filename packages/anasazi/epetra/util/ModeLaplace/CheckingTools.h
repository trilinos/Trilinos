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

#ifndef ANASAZI_CHECKING_TOOLS_H
#define ANASAZI_CHECKING_TOOLS_H

#include "Epetra_ConfigDefs.h"
#include "Anasaziepetra_ModeLaplace_DLLExportMacro.h"

#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"

class ANASAZIEPETRA_MODELAPLACE_LIB_DLL_EXPORT CheckingTools {

  private:

    const Epetra_Comm &MyComm;

  public:

    CheckingTools(const Epetra_Comm &_Comm);

    double errorOrthogonality(const Epetra_MultiVector *X, const Epetra_MultiVector *R,
                              const Epetra_Operator *M = 0) const;

    double errorOrthonormality(const Epetra_MultiVector *X, const Epetra_Operator *M = 0) const;

    double errorEquality(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
                         const Epetra_Operator *M = 0) const;

    int errorSubspaces(const Epetra_MultiVector &Q, const Epetra_MultiVector &Qex,
                       const Epetra_Operator *M) const;

    void errorEigenResiduals(const Epetra_MultiVector &Q, double *lambda,
                             const Epetra_Operator *K, const Epetra_Operator *M,
                             double *normWeight = 0) const;

    void errorEigenResiduals(const Epetra_MultiVector &Q, double *lambda,
                             const Epetra_Operator *K, const Epetra_Operator *M,
                             const Epetra_Operator *Msolver) const;

    int errorLambda(double *continuous, double *discrete, int numDiscrete, double *lambda,
                    int nev, bool ascend=true) const;

    int inputArguments(const int &numEigen, const Epetra_Operator *K,
                       const Epetra_Operator *M, const Epetra_Operator *P,
                       const Epetra_MultiVector &Q, const int &minSize) const;

};


#endif

