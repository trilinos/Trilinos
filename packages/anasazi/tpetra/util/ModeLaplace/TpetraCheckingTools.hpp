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

#ifndef ANASAZI_TPETRA_CHECKING_TOOLS_H
#define ANASAZI_TPETRA_CHECKING_TOOLS_H

#include "Anasazitpetra_ModeLaplace_DLLExportMacro.h"

#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

template <class Scalar, class LO, class GO, class Node>
class ANASAZITPETRA_MODELAPLACE_LIB_DLL_EXPORT TpetraCheckingTools {

  private:

    typedef Tpetra::MultiVector<Scalar,LO,GO,Node>   MV;
    typedef Tpetra::Operator<Scalar,LO,GO,Node>      OP;
    Teuchos::RCP<const Teuchos::Comm<int> >&     MyComm;

  public:

    TpetraCheckingTools(Teuchos::RCP<const Teuchos::Comm<int> > &_Comm)
    : MyComm( _Comm )
    {}

/*
    // These methods need to be filled out.
    double errorOrthogonality(const MV *X, const MV *R,
                              const OP *M = 0) const { return 0.0; }

    double errorOrthonormality(const MV *X, const OP *M = 0) const { return 0.0; }

    double errorEquality(const MV *X, const MV *MX,
                         const OP *M = 0) const { return 0.0; }

    int errorSubspaces(const MV &Q, const MV &Qex,
                       const OP *M) const { return 0; }

    void errorEigenResiduals(const MV &Q, double *lambda,
                             const OP *K, const OP *M,
                             double *normWeight = 0) const {}

    void errorEigenResiduals(const MV &Q, double *lambda,
                             const OP *K, const OP *M,
                             const OP *Msolver) const {}

    int errorLambda(double *continuous, double *discrete, int numDiscrete, double *lambda,
                    int nev, bool ascend=true) const { return 0; }

    int inputArguments(const int &numEigen, const OP *K,
                       const OP *M, const OP *P,
                       const MV &Q, const int &minSize) const { return 0; }
*/
};


#endif
