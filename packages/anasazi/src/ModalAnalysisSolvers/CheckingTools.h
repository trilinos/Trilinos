//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on any reproduction of this software, in whole or
// in part.
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#ifndef CHECKING_TOOLS_H
#define CHECKING_TOOLS_H

#include "Epetra_ConfigDefs.h"

#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"

class CheckingTools {

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

    int errorLambda(double *continuous, double *discrete, int numDiscrete, double *lambda,
                    int nev) const;

    int inputArguments(const int &numEigen, const Epetra_Operator *K,
                       const Epetra_Operator *M, const Epetra_Operator *P,
                       const Epetra_MultiVector &Q, const int &minSize) const;

};


#endif

